#include "Insert/Fields/HallWallCharge.H"

#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <ablastr/utils/Communication.H>

#include <AMReX_Geometry.H>

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
        max_level == 0,
        "Persistent analytic-wall charge currently supports a single grid level only.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        static_cast<int>(rho.size()) == max_level + 1 && rho[0] != nullptr &&
            rho[0]->nComp() == 1,
        "Persistent analytic-wall charge requires valid rho on every grid level.");

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
    m_defined = true;
}

void
HallWallCharge::addTo (
    ablastr::fields::MultiLevelScalarField const& rho, int const max_level)
{
    prepare(rho, max_level);
    synchronize();

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
    m_dirty = true;
    return *m_wall_charge[lev];
}

void
HallWallCharge::synchronize ()
{
    if (!m_dirty) {
        return;
    }

    WarpX const& warpx = WarpX::GetInstance();
    for (int lev = 0; lev < static_cast<int>(m_wall_charge.size()); ++lev) {
        amrex::MultiFab& wall_charge = *m_wall_charge[lev];
        ablastr::utils::communication::SumBoundary(
            wall_charge, 0, wall_charge.nComp(), wall_charge.nGrowVect(),
            wall_charge.nGrowVect(), WarpX::do_single_precision_comms,
            warpx.Geom(lev).periodicity());
        // SumBoundary leaves the source values in ghost cells intact. Clear
        // them so that a later persistent deposition cannot add them again.
        wall_charge.setBndry(0.0_rt);
    }
    m_dirty = false;
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
            m_wall_charge[lev]->boxArray().CellEqual(rho[lev]->boxArray()) &&
                m_wall_charge[lev]->DistributionMap() == rho[lev]->DistributionMap() &&
                m_wall_charge[lev]->nGrowVect() == rho[lev]->nGrowVect() &&
                rho[lev]->nComp() == 1,
            "Persistent analytic-wall charge layout changed after initialization.");
    }
}

} // namespace Insert
