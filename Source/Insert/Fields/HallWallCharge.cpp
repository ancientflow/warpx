#include "Insert/Fields/HallWallCharge.H"

#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <ablastr/utils/Communication.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

#include <fstream>
#include <iomanip>
#include <numeric>

namespace Insert {

using namespace amrex::literals;

WallChargeGrid
MakeWallChargeGrid (amrex::Geometry const& geom)
{
#if defined(WARPX_DIM_3D)
    amrex::Box domain = geom.Domain();
    domain.surroundingNodes();
    auto const prob_lo = geom.ProbLoArray();
    auto const cell_size = geom.CellSizeArray();

    return WallChargeGrid{
        amrex::lbound(domain), amrex::ubound(domain),
        {AMREX_D_DECL(prob_lo[0], prob_lo[1], prob_lo[2])},
        {AMREX_D_DECL(1.0_rt / cell_size[0], 1.0_rt / cell_size[1],
                      1.0_rt / cell_size[2])},
        1.0_rt / (cell_size[0] * cell_size[1] * cell_size[2])};
#else
    amrex::ignore_unused(geom);
    return {};
#endif
}

HallWallCharge&
HallWallCharge::GetInstance ()
{
    static HallWallCharge wall_charge;
    return wall_charge;
}

void
HallWallCharge::prepare (
    ablastr::fields::MultiLevelScalarField const& rho, int const max_level)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        static_cast<int>(rho.size()) == max_level + 1,
        "Persistent analytic-wall charge requires rho on every grid level.");
    for (auto const* field : rho) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            field && field->nComp() == 1 && field->ixType().nodeCentered(),
            "Persistent analytic-wall charge requires single-component nodal rho.");
    }
    if (max_level > 0) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!WarpX::use_filter,
            "Multilevel wall charge requires warpx.use_filter = 0.");
        bool ecdi_enabled = false;
        amrex::ParmParse("insert.ecdi_control").query("enabled", ecdi_enabled);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!ecdi_enabled,
            "Multilevel wall charge cannot be combined with ECDI filtering.");
    }

    if (m_defined) {
        assertLayout(rho, max_level);
        return;
    }

    m_wall_charge.resize(max_level + 1);
    for (int lev = 0; lev <= max_level; ++lev) {
        m_wall_charge[lev] = std::make_unique<amrex::MultiFab>(
            rho[lev]->boxArray(), rho[lev]->DistributionMap(), 1,
            rho[lev]->nGrowVect());
        m_wall_charge[lev]->setVal(0.0_rt);
    }
    amrex::ParmParse("insert").query("wall_charge_diagnostic", m_diagnostic);
    if (m_diagnostic && amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream output("wall_charge.dat", std::ios::trunc);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(output.is_open(), "Cannot create wall_charge.dat.");
        output << "# step time deposited_C stored_C level_charge_C... "
                  "composite_C composite_minus_deposited_C\n";
    }
    m_defined = true;
}

void
HallWallCharge::addTo (
    ablastr::fields::MultiLevelScalarField const& rho, int const max_level)
{
    prepare(rho, max_level);

    if (max_level == 0) {
        auto const& raw = *m_wall_charge[0];
        amrex::MultiFab assembled(raw.boxArray(), raw.DistributionMap(), 1, raw.nGrowVect());
        amrex::MultiFab::Copy(assembled, raw, 0, 0, 1, raw.nGrowVect());
        ablastr::utils::communication::SumBoundary(
            assembled, 0, 1, assembled.nGrowVect(), assembled.nGrowVect(),
            WarpX::do_single_precision_comms, WarpX::GetInstance().Geom(0).periodicity());
        assembled.setBndry(0.0_rt);
        amrex::MultiFab::Add(*rho[0], assembled, 0, 0, 1, rho[0]->nGrowVect());
        return;
    }

    for (int lev = 0; lev <= max_level; ++lev) {
        amrex::MultiFab::Add(
            *rho[lev], *m_wall_charge[lev], 0, 0, 1, rho[lev]->nGrowVect());
    }
}

amrex::MultiFab&
HallWallCharge::get (int const lev)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_defined && lev >= 0 && lev < static_cast<int>(m_wall_charge.size()),
        "Persistent analytic-wall charge was accessed before it was initialized.");
    return *m_wall_charge[lev];
}

void
HallWallCharge::depositImpacts (ImpactContainer& impacts)
{
#if defined(WARPX_DIM_3D)
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_defined, "Wall charge must be prepared before deposition.");
    auto const& warpx = WarpX::GetInstance();
    if (m_diagnostic) {
        // Read the staging tiles directly: they have not yet acquired MPI ownership.
        for (auto const& entry : impacts.GetParticles(0)) {
            auto const& tile = entry.second;
            auto const* particles = tile.GetArrayOfStructs().data();
            m_deposited_charge += amrex::Reduce::Sum<amrex::Real>(tile.numParticles(),
                [=] AMREX_GPU_DEVICE (int i) noexcept -> amrex::Real {
                    return particles[i].id() > 0 ? particles[i].rdata(0) : 0.0_rt;
                });
        }
    }
    // All ranks participate, including ranks with no local impacts. Each event
    // is assigned once by its intersection, independently of the incident tile.
    impacts.Redistribute();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev) {
        auto const grid = MakeWallChargeGrid(warpx.Geom(lev));
        for (ImpactContainer::ParIterType pti(impacts, lev); pti.isValid(); ++pti) {
            auto const* particles = pti.GetArrayOfStructs().data();
            auto const charge = m_wall_charge[lev]->array(pti);
            // Scatter-add must not use ParallelFor's CPU SIMD promise.
            amrex::For(pti.numParticles(), [=] AMREX_GPU_DEVICE (int i) noexcept {
                auto const& p = particles[i];
                DepositWallChargeToNodes(charge, grid,
                    p.rdata(1), p.rdata(2), p.rdata(3), p.rdata(0));
            });
        }
    }
    amrex::Gpu::streamSynchronize();
    if (m_diagnostic) {
        auto expected = m_deposited_charge;
        amrex::ParallelDescriptor::ReduceRealSum(expected);
        auto const levels = chargeByLevel();
        auto const actual = std::accumulate(levels.begin(), levels.end(), amrex::Real(0.0));
        auto const composite = compositeCharge();
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream output("wall_charge.dat", std::ios::app);
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(output.is_open(), "Cannot append wall_charge.dat.");
            output << std::setprecision(17) << warpx.getistep(0) << ' '
                   << warpx.gett_new(0) << ' ' << expected << ' ' << actual;
            for (auto const value : levels) { output << ' ' << value; }
            output << ' ' << composite << ' ' << composite-expected << '\n';
        }
    }
#else
    amrex::ignore_unused(impacts);
#endif
}

amrex::Vector<amrex::Real>
HallWallCharge::chargeByLevel () const
{
    amrex::Vector<amrex::Real> result(m_wall_charge.size(), 0.0);
    auto const& warpx = WarpX::GetInstance();
    for (int lev = 0; lev < static_cast<int>(m_wall_charge.size()); ++lev) {
        auto const& field = *m_wall_charge[lev];
        // Raw FAB entries are distinct contributions, including shared nodes and
        // ghosts. Summing them before assembly counts each hat weight once and
        // includes no fine-to-coarse solver copies. Their normalization is q/dV.
        auto const volume = 1.0 / MakeWallChargeGrid(warpx.Geom(lev)).inv_cell_volume;
        amrex::ReduceOps<amrex::ReduceOpSum> ops;
        amrex::ReduceData<amrex::Real> data(ops);
        using Tuple = decltype(data)::Type;
        for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
            auto const source = field.const_array(mfi);
            ops.eval(mfi.fabbox(), data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> Tuple {
                    return {source(i,j,k)};
                });
        }
        result[lev] = amrex::get<0>(data.value()) * volume;
    }
    amrex::ParallelDescriptor::ReduceRealSum(result.data(), static_cast<int>(result.size()));
    return result;
}

amrex::Real
HallWallCharge::compositeCharge () const
{
    auto& warpx = WarpX::GetInstance();
    int const nlev = static_cast<int>(m_wall_charge.size());
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> fine(nlev), coarse(nlev);
    ablastr::fields::MultiLevelScalarField buffer(nlev, nullptr);
    for (int lev = 0; lev < nlev; ++lev) {
        auto const& raw = *m_wall_charge[lev];
        fine[lev] = std::make_unique<amrex::MultiFab>(
            raw.boxArray(), raw.DistributionMap(), 1, raw.nGrowVect());
        amrex::MultiFab::Copy(*fine[lev], raw, 0, 0, 1, raw.nGrowVect());
        if (lev > 0) {
            auto ba = raw.boxArray();
            ba.coarsen(warpx.refRatio(lev-1));
            coarse[lev] = std::make_unique<amrex::MultiFab>(
                ba, raw.DistributionMap(), 1, raw.nGrowVect());
        }
    }
    // Never write assembly/restriction results back into the persistent source.
    // Report before PEC image-charge processing and conductor rho removal.
    warpx.SyncRho(amrex::GetVecOfPtrs(fine), amrex::GetVecOfPtrs(coarse), buffer);
    amrex::Real result = 0.0;
    for (int lev = 0; lev < nlev; ++lev) {
        amrex::iMultiFab uncovered(warpx.boxArray(lev), warpx.DistributionMap(lev), 1, 0);
        if (lev+1 < nlev) {
            auto mask = amrex::makeFineMask(warpx.boxArray(lev), warpx.DistributionMap(lev),
                warpx.boxArray(lev+1), warpx.refRatio(lev), 1, 0);
            amrex::iMultiFab::Copy(uncovered, mask, 0, 0, 1, 0);
        } else {
            uncovered.setVal(1);
        }
        amrex::ReduceOps<amrex::ReduceOpSum> ops;
        amrex::ReduceData<amrex::Real> data(ops);
        using Tuple = decltype(data)::Type;
        for (amrex::MFIter mfi(uncovered, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto const source = fine[lev]->const_array(mfi);
            auto const mask = uncovered.const_array(mfi);
            // Cell-wise quadrature gives each corner dV/8 in 3D. This is the
            // nodal dual-volume weight, including patch and physical edges;
            // each cell has one MPI owner, so shared nodes are not double counted.
            ops.eval(mfi.tilebox(), data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> Tuple {
                    if (mask(i,j,k) == 0) { return {0.0}; }
                    amrex::Real value = 0.0;
                    for (int corner = 0; corner < (1 << AMREX_SPACEDIM); ++corner) {
                        amrex::IntVect p(AMREX_D_DECL(i,j,k));
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            p[d] += (corner >> d) & 1;
                        }
                        value += source(p);
                    }
                    return {value / (1 << AMREX_SPACEDIM)};
                });
        }
        result += amrex::get<0>(data.value()) /
            MakeWallChargeGrid(warpx.Geom(lev)).inv_cell_volume;
    }
    amrex::ParallelDescriptor::ReduceRealSum(result);
    return result;
}

void
HallWallCharge::assertLayout (
    ablastr::fields::MultiLevelScalarField const& rho, int const max_level) const
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        static_cast<int>(m_wall_charge.size()) == max_level + 1,
        "Persistent analytic-wall charge grid-level count changed after initialization.");

    for (int lev = 0; lev <= max_level; ++lev) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_wall_charge[lev]->boxArray() == rho[lev]->boxArray() &&
                m_wall_charge[lev]->DistributionMap() == rho[lev]->DistributionMap() &&
                m_wall_charge[lev]->nGrowVect() == rho[lev]->nGrowVect() &&
                rho[lev]->nComp() == 1,
            "Persistent analytic-wall charge layout changed after initialization.");
    }
}

} // namespace Insert
