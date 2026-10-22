#include "BackgroundCoupledDensity.h"

#include "BoundaryConditions/WarpX_PEC.H"
#include "Insert/Utils/InsertUtils.h"
#include "Parallelization/WarpXSumGuardCells.H"
#include "WarpX.H"

#include <AMReX_ParmParse.H>

#if defined(MCC_DENSITY_AVERAGE_CALC) || defined(MCC_DENSITY_AVERAGE_USE)
#include <AMReX_VisMF.H>
#endif

#if defined(MCC_DENSITY_AVERAGE_CALC) || defined(MCC_DENSITY_AVERAGE_USE)
#include <cctype>
#endif
#include <cmath>
#ifdef MCC_DENSITY_AVERAGE_CALC
#include <iomanip>
#endif
#include <limits>
#if defined(MCC_DENSITY_AVERAGE_CALC) || defined(MCC_DENSITY_AVERAGE_USE)
#include <sstream>
#endif

#if defined(MCC_DENSITY_AVERAGE_CALC) || defined(MCC_DENSITY_AVERAGE_USE)
namespace {

std::string
SafeName (std::string name) {
    for (char& c : name) {
        if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' &&
            c != '-') {
            c = '_';
        }
    }
    return name;
}

std::string
FabOutputParent (std::string const& base_dir, std::string const& species,
                 std::string const& quantity, int lev) {
    std::ostringstream os;
    os << base_dir << "/" << SafeName(species) << "/" << quantity << "/lev"
       << lev;
    return os.str();
}

std::string
FabOutputPath (std::string const& base_dir, std::string const& species,
               std::string const& quantity, int lev) {
    std::ostringstream os;
    os << FabOutputParent(base_dir, species, quantity, lev) << "/" << quantity;
    return os.str();
}

} // namespace
#endif

#ifdef MCC_DENSITY_AVERAGE_CALC
namespace {

std::string
FabStepOutputPath (std::string const& base_dir, std::string const& species,
                   std::string const& quantity, int lev, int step) {
    std::ostringstream os;
    os << FabOutputParent(base_dir, species, quantity, lev) << "/step"
       << std::setw(8) << std::setfill('0') << step;
    return os.str();
}

void
WriteSingleFabMultiFab (amrex::MultiFab const& mf, std::string const& path,
                        std::string const& species, int lev) {
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        mf.boxArray().size() == 1,
        "Background density FAB output requires exactly one FAB for species " +
            species + " at level " + std::to_string(lev) + ".");
    Insert::CreateDirectoryTree(Insert::ParentPath(path));
    amrex::VisMF::Write(mf, path);
}

} // namespace
#endif

/**
 * @brief deposit the atom density by WarpXParticleContainer
 */
void
AtomDepositAPI (WarpXParticleContainer& pc, amrex::MultiFab& rho,
                const int lev) {
    for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {
        auto np = pti.numParticles();
        amrex::AllPrint() << "rank " << amrex::ParallelDescriptor::MyProc()
                          << ": deposit " << np << " particles" << "\n";
        // Extract particle data
        auto& attribs = pti.GetAttribs();
        auto& wp = attribs[PIdx::w];
        pc.DepositCharge(pti, wp, nullptr, &rho, 0, 0, np, 0, lev, lev);
    }

    auto& warpx_instance = WarpX::GetInstance();
    WarpXSumGuardCells(rho, warpx_instance.Geom(lev).periodicity(),
                       rho.nGrowVect(), 0, 1);

    // 处理边界密度
    // 尝试使用内置函数进行修改，便于进行高阶插值
    // 考虑到更多的层级，必须采用这种内置函数
    // 考虑到PEC中的修改，这会导致边界密度被大大低估，但是一般来说这不会明细影响放电
    PEC::ApplyReflectiveBoundarytoRhofield(
        &rho, WarpX::field_boundary_lo, WarpX::field_boundary_hi,
        WarpX::particle_boundary_lo, WarpX::particle_boundary_hi,
        warpx_instance.Geom(lev), lev, PatchType::fine,
        warpx_instance.refRatio());
}

/**
 * @brief init the vector of background density species
 */
void
BackgroundCoupledDensity::backgroundDensityInit () {
    WarpX& warpx_instance = WarpX::GetInstance();
    if (warpx_instance.maxLevel() > 0) {
        int regrid_interval = -1;
        int subcycling = 0;
        amrex::ParmParse("warpx").query("regrid_int", regrid_interval);
        amrex::ParmParse("warpx").query("do_subcycling", subcycling);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            regrid_interval <= 0 && subcycling == 0 && !WarpX::do_moving_window &&
                !warpx_instance.get_load_balance_intervals().isActivated(),
            "Coupled MCC with mesh refinement requires a static mesh: disable regrid, "
            "moving window, AMR time subcycling and load balancing.");
    }

#ifndef MCC_DENSITY_AVERAGE_USE
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();
    auto& background_species =
        mypc.GetParticleContainerFromName(m_ground_species);
    // resize the vector
    auto const flvl = background_species.finestLevel();
#else
    // In USE mode no neutral particle container is required;
    // the averaged density is loaded directly from file.
    int const flvl = 0;
#endif
    m_background_density_fabs.resize(flvl + 1);
    m_ionization_consumption_fabs.resize(flvl + 1);
#ifdef MCC_DENSITY_AVERAGE_CALC
    m_background_density_sum_fabs.resize(flvl + 1);
#endif
    m_tiles.resize(flvl + 1);

    amrex::Print() << "Background density initialized with size: \n";

    for (int lev = 0; lev <= flvl; lev++) {
        // Private nodal density: do not borrow a Poisson charge FAB or its contents.
        auto const ba = amrex::convert(warpx_instance.boxArray(lev), amrex::IntVect(1));
        auto const& dm = warpx_instance.DistributionMap(lev);
        auto ng = warpx_instance.get_ng_depos_rho();
        ng.max(amrex::IntVect(WarpX::nox + 1));
        m_background_density_fabs[lev] = amrex::MultiFab(ba, dm, WarpX::ncomps, ng);
        m_background_density_fabs[lev].setVal(0.0);
        m_ionization_consumption_fabs[lev] = amrex::MultiFab(ba, dm, 1, ng);
        m_ionization_consumption_fabs[lev].setVal(0.0);
#ifdef MCC_DENSITY_AVERAGE_CALC
        m_background_density_sum_fabs[lev] =
            amrex::MultiFab(ba, dm, WarpX::ncomps, ng);
        m_background_density_sum_fabs[lev].setVal(0.0);
#endif
    }
#ifdef MCC_DENSITY_AVERAGE_CALC
    amrex::ParmParse const pp_bg("background_density");
    pp_bg.query("average_steps_per_period", m_average_steps_per_period);
    pp_bg.query("average_periods", m_average_periods);
    pp_bg.query("raw_output_interval", m_raw_output_interval);
    pp_bg.query("output_dir", m_output_dir);
    pp_bg.query("output_fab", m_output_fab);

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_average_steps_per_period > 0,
        "background_density.average_steps_per_period must be positive.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_average_periods > 0,
        "background_density.average_periods must be positive.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_raw_output_interval >= 0,
        "background_density.raw_output_interval must be non-negative.");

    m_average_window_steps = m_average_steps_per_period * m_average_periods;
#endif
#ifdef MCC_DENSITY_AVERAGE_USE
    amrex::ParmParse const pp_bg("background_density");
    pp_bg.query("input_fab", m_input_fab);
    if (m_input_fab.empty()) {
        m_input_fab = FabOutputPath("background_density_fab", m_ground_species,
                                    "average", 0);
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_background_density_fabs[0].boxArray().size() == 1,
        "MCC_DENSITY_AVERAGE_USE requires the current background density "
        "MultiFab to contain exactly one FAB.");
    amrex::Print() << "Reading averaged background density FAB from "
                   << m_input_fab << "\n";
    amrex::VisMF::Read(m_background_density_fabs[0], m_input_fab);
    //m_background_density_fabs[0].mult(0.5);
#endif
}

/**
 * @brief update the background density data
 */
void
BackgroundCoupledDensity::backgroundDensityUpdate (
    MultiParticleContainer& mypc, amrex::ParticleReal elec_weight, int step) {
    using namespace amrex::literals;
#ifndef MCC_DENSITY_AVERAGE_CALC
    static_cast<void>(step);
#endif

    WarpX& warpx_instance = WarpX::GetInstance();
    auto& background_species =
        mypc.GetParticleContainerFromName(m_ground_species);
    auto const flvl = background_species.finestLevel();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(elec_weight) && elec_weight > 0,
        "Coupled MCC requires a finite positive my_constants.elec_weight.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        static_cast<int>(m_background_density_fabs.size()) == flvl + 1,
        "Coupled MCC background level count changed after initialization.");
    amrex::AllPrint() << "rank " << amrex::ParallelDescriptor::MyProc()
                      << ":Start updating background density\n";
    for (int lev = 0; lev <= flvl; lev++) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_background_density_fabs[lev].boxArray() ==
                amrex::convert(warpx_instance.boxArray(lev), amrex::IntVect(1)) &&
            m_background_density_fabs[lev].DistributionMap() ==
                warpx_instance.DistributionMap(lev),
            "Coupled MCC background layout changed after initialization.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            background_species.ParticleBoxArray(lev) == warpx_instance.boxArray(lev) &&
            background_species.ParticleDistributionMap(lev) ==
                warpx_instance.DistributionMap(lev),
            "Coupled MCC background particles must use the field mesh layout.");
        m_background_density_fabs[lev].setVal(0.0_prt);
        AtomDepositAPI(background_species, m_background_density_fabs[lev], lev);
#ifdef MCC_DENSITY_AVERAGE_CALC
        if (m_raw_output_interval > 0 && step % m_raw_output_interval == 0) {
            WriteSingleFabMultiFab(m_background_density_fabs[lev],
                                   FabStepOutputPath(m_output_dir,
                                                     m_ground_species, "raw",
                                                     lev, step),
                                   m_ground_species, lev);
        }
        if (m_average_sample_count < m_average_window_steps) {
            amrex::MultiFab::Add(m_background_density_sum_fabs[lev],
                                 m_background_density_fabs[lev], 0, 0,
                                 m_background_density_fabs[lev].nComp(),
                                 m_background_density_fabs[lev].nGrowVect());
        }
#endif
        auto const& geo = warpx_instance.Geom(lev);
        auto& tiles = m_tiles[lev];
        amrex::Gpu::streamSynchronize();
        tiles.clear();
        // Include empty tiles: a projectile can occupy a tile with no neutrals.
        for (auto mfi = background_species.MakeMFIter(lev); mfi.isValid(); ++mfi) {
            auto const key = std::make_pair(mfi.index(), mfi.LocalTileIndex());
            auto& ptile = background_species.DefineAndReturnParticleTile(
                lev, key.first, key.second);
            auto& tile = tiles[key];
            tile.bins = ParticleUtils::findParticlesInEachCell(geo, mfi, ptile);
            auto const* offsets = tile.bins.offsetsPtr();
            auto const* indices = tile.bins.permutationPtr();
            auto const numbins = tile.bins.numBins();
            tile.available.resize(numbins);
            auto* available = tile.available.dataPtr();
            auto const& soa = ptile.GetStructOfArrays();
            auto const* idcpu = soa.GetIdCPUData().data();
            auto const* pw = soa.GetRealData()[PIdx::w].dataPtr();
            amrex::ParallelFor(numbins, [=] AMREX_GPU_DEVICE(long ibin) {
                amrex::Long count = 0;
                for (int i = offsets[ibin]; i < offsets[ibin + 1]; ++i) {
                    int const ip = indices[i];
                    if (amrex::ConstParticleIDWrapper{idcpu[ip]}.is_valid() && pw[ip] > 0) {
                        // Count consumable weight units per particle, not fractions
                        // pooled across particles that cannot individually be consumed.
                        count += static_cast<amrex::Long>(pw[ip] / elec_weight + 1.e-6);
                    }
                }
                available[ibin] = static_cast<int>(count);
            });
        }
    }
#ifdef MCC_DENSITY_AVERAGE_CALC
    if (m_average_sample_count < m_average_window_steps) {
        ++m_average_sample_count;
    }
#endif
    amrex::Print() << "rank " << amrex::ParallelDescriptor::MyProc()
                   << ": Updated background species: " << m_ground_species
                   << " end\n";
}

void
BackgroundCoupledDensity::resetIonizationConsumption (int lev) {
    m_ionization_consumption_fabs[lev].setVal(0.0);
}

void
BackgroundCoupledDensity::applyIonizationConsumption (int lev) {
    auto& warpx_instance = WarpX::GetInstance();
    auto& consumption = m_ionization_consumption_fabs[lev];
    auto& background_density = m_background_density_fabs[lev];

    // Deposits in guard cells and overlapping nodal regions are additive.
    // WarpXSumGuardCells delegates the MPI exchange to AMReX SumBoundary.
    WarpXSumGuardCells(consumption, warpx_instance.Geom(lev).periodicity(),
                       consumption.nGrowVect(), 0, 1);

    // Apply the same linear physical-boundary operation as the full background
    // density deposition before subtracting this incremental consumption.
    PEC::ApplyReflectiveBoundarytoRhofield(
        &consumption, WarpX::field_boundary_lo, WarpX::field_boundary_hi,
        WarpX::particle_boundary_lo, WarpX::particle_boundary_hi,
        warpx_instance.Geom(lev), lev, PatchType::fine,
        warpx_instance.refRatio());

    amrex::MultiFab::Subtract(background_density, consumption, 0, 0, 1,
                              background_density.nGrowVect());
}

void
BackgroundCoupledDensity::backgroundDensityFinalize () {
#ifdef MCC_DENSITY_AVERAGE_CALC
    using namespace amrex::literals;

    if (m_average_sample_count <= 0) {
        amrex::Print() << "No background density samples collected for species "
                       << m_ground_species << "; skip averaged FAB output.\n";
        return;
    }

    auto const sample_count = static_cast<amrex::Real>(m_average_sample_count);
    amrex::Vector<amrex::MultiFab> average_fabs(
        m_background_density_sum_fabs.size());
    for (int lev = 0;
         lev < static_cast<int>(m_background_density_sum_fabs.size()); ++lev) {
        average_fabs[lev] = amrex::MultiFab(
            m_background_density_sum_fabs[lev].boxArray(),
            m_background_density_sum_fabs[lev].DistributionMap(),
            m_background_density_sum_fabs[lev].nComp(),
            m_background_density_sum_fabs[lev].nGrowVect());
        amrex::MultiFab::Copy(average_fabs[lev],
                              m_background_density_sum_fabs[lev], 0, 0,
                              m_background_density_sum_fabs[lev].nComp(),
                              m_background_density_sum_fabs[lev].nGrowVect());
        average_fabs[lev].mult(1.0_rt / sample_count);

        std::string const output_path =
            (lev == 0 && !m_output_fab.empty())
                ? m_output_fab
                : FabOutputPath(m_output_dir, m_ground_species, "average", lev);
        WriteSingleFabMultiFab(average_fabs[lev], output_path, m_ground_species,
                               lev);
    }

    amrex::Print() << "Wrote averaged background density FAB for species "
                   << m_ground_species << " from " << m_average_sample_count
                   << " samples.\n";
#endif
}

/**
 * @brief delete the particles with zero weight
 */
void
BackgroundCoupledDensity::backgroundSpeciesClean (
    MultiParticleContainer& mypc) const {
    using namespace amrex::literals;

    auto& background_species =
        mypc.GetParticleContainerFromName(m_ground_species);
    auto const flvl = background_species.finestLevel();
    for (int lev = 0; lev <= flvl; lev++) {
        for (WarpXParIter pti(background_species, lev); pti.isValid(); ++pti) {
            amrex::AllPrint()
                << "rank " << amrex::ParallelDescriptor::MyProc()
                << ": clean lev " << lev << ": box " << pti.index() << "\n";
            auto& ptile = background_species.ParticlesAt(lev, pti);
            long const np = ptile.numParticles();
            auto& soa = ptile.GetStructOfArrays();
            uint64_t* const AMREX_RESTRICT idcpu = soa.GetIdCPUData().data();
            auto& soa_arr = soa.GetRealData();
            amrex::ParticleReal const* pw = soa_arr[PIdx::w].dataPtr();
            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
                if (std::abs(pw[ip]) < 10.0_prt) {
                    auto pidw = amrex::ParticleIDWrapper{idcpu[ip]};
                    if (pidw.is_valid()) {
                        pidw.make_invalid();
                    }
                }
            });
        }
    }
}
