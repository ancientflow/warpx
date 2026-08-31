/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "CoupledBackgroundMCCCollision.H"

#include "Insert/Collisions/FilterCopyTransformCoupled.H"
#include "Insert/Collisions/IonizationSourceTable.h"
#include "Insert/Config/WarpXFunctionConfig.h"

#include "Particles/Algorithms/KineticEnergy.H"
#include "Particles/Collision/BackgroundMCC/ImpactIonization.H"
#include "Particles/Collision/BinaryCollision/BinaryCollisionUtils.H"
#include "Particles/Collision/BinaryCollision/TwoProductUtil.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/ParticleUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "WarpX.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <ablastr/profiler/ProfilerWrapper.H>
#include <string>

#include "Insert/Background/InsertBackgroundDensity.h"
#include "Insert/Core/WarpXFunc.h"

namespace {

class CoupledImpactIonizationFilterFunc {
  public:
    CoupledImpactIonizationFilterFunc (
        ScatteringProcess const& mcc_process, double const mass,
        amrex::ParticleReal const total_collision_prob,
        amrex::ParticleReal const nu_max)
        : m_mcc_process(mcc_process.executor()), m_mass(mass),
          m_total_collision_prob(total_collision_prob), m_nu_max(nu_max) {}

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
    operator()(const PData& ptd, int const i, amrex::RandomEngine const& engine,
               amrex::Real local_density) const noexcept {
        using namespace amrex;
        using std::sqrt;

        if (Random(engine) > m_total_collision_prob) {
            return false;
        }

        const ParticleReal ux = ptd.m_rdata[PIdx::ux][i];
        const ParticleReal uy = ptd.m_rdata[PIdx::uy][i];
        const ParticleReal uz = ptd.m_rdata[PIdx::uz][i];

        constexpr auto eV = PhysConst::q_e;
        const double E_coll =
            Algorithms::KineticEnergy<double>(ux, uy, uz, m_mass) / eV;

        const ParticleReal sigma_E =
            m_mcc_process.getCrossSection(static_cast<ParticleReal>(E_coll));
        const ParticleReal u_coll2 = ux * ux + uy * uy + uz * uz;
        const ParticleReal nu_i =
            local_density * sigma_E * sqrt(u_coll2) / m_nu_max;

        return (Random(engine) <= nu_i);
    }

  private:
    ScatteringProcess::Executor m_mcc_process;
    double m_mass;
    amrex::ParticleReal m_total_collision_prob = 0;
    amrex::ParticleReal m_nu_max;
};

} // namespace

CoupledBackgroundMCCCollision::CoupledBackgroundMCCCollision (
    std::string const& collision_name)
    : CollisionBase(collision_name) {
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_species_names.size() == 1,
        "Background MCC must have exactly one species.");

    const amrex::ParmParse pp_collision_name(collision_name);

    amrex::ParticleReal background_temperature;
    if (utils::parser::queryWithParser(pp_collision_name,
                                       "background_temperature",
                                       background_temperature)) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            (background_temperature >= 0),
            "The background temperature must be positive.");
        m_background_temperature_parser = utils::parser::makeParser(
            std::to_string(background_temperature), {"x", "y", "z", "t"});
    } else {
        std::string background_temperature_str;
        utils::parser::Store_parserString(pp_collision_name,
                                          "background_temperature(x,y,z,t)",
                                          background_temperature_str);
        m_background_temperature_parser = utils::parser::makeParser(
            background_temperature_str, {"x", "y", "z", "t"});
    }

    // The coupled model reads density from BackgroundCoupledDensity fields.
    m_background_temperature_func =
        m_background_temperature_parser.compile<4>();

    // if the neutral mass is specified use it, but if ionization is
    // included the mass of the secondary species of that interaction
    // will be used. If no neutral mass is specified and ionization is not
    // included the mass of the colliding species will be used
    m_background_mass = -1;
    utils::parser::queryWithParser(pp_collision_name, "background_mass",
                                   m_background_mass);

    // Parse the list of scattering processes (these could be elastic,
    // excitation, charge_exchange, etc.) and create a vector of
    // ScatteringProcess objects from each scattering process name. The
    // angular behavior of each process is set with the per-process
    // `<process>_scattering_angle_model` input (isotropic | forward |
    // backward).
    amrex::Vector<ScatteringProcess> scattering_processes =
        BinaryCollisionUtils::parse_scattering_processes(collision_name);

    for (auto& process : scattering_processes) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            process.type() != ScatteringProcessType::INVALID,
            "Cannot add an unknown scattering process type");

        // if the scattering process is ionization get the secondary species
        // only one ionization process is supported, the vector
        // m_ionization_processes is only used to make it simple to calculate
        // the maximum collision frequency with the same function used for
        // particle conserving processes
        if (process.type() == ScatteringProcessType::IONIZATION) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                !ionization_flag,
                "Background MCC only supports a single ionization process");
            ionization_flag = true;

            std::string secondary_species;
            pp_collision_name.get("ionization_species", secondary_species);
            m_species_names.push_back(secondary_species);

            m_ionization_processes.push_back(std::move(process));
        } else {
            m_scattering_processes.push_back(std::move(process));
        }
    }

#ifdef AMREX_USE_GPU
    amrex::Gpu::HostVector<ScatteringProcess::Executor>
        h_scattering_processes_exe;
    amrex::Gpu::HostVector<ScatteringProcess::Executor>
        h_ionization_processes_exe;
    for (auto const& p : m_scattering_processes) {
        h_scattering_processes_exe.push_back(p.executor());
    }
    for (auto const& p : m_ionization_processes) {
        h_ionization_processes_exe.push_back(p.executor());
    }
    m_scattering_processes_exe.resize(h_scattering_processes_exe.size());
    m_ionization_processes_exe.resize(h_ionization_processes_exe.size());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, h_scattering_processes_exe.begin(),
        h_scattering_processes_exe.end(), m_scattering_processes_exe.begin());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, h_ionization_processes_exe.begin(),
        h_ionization_processes_exe.end(), m_ionization_processes_exe.begin());
    amrex::Gpu::streamSynchronize();
#else
    for (auto const& p : m_scattering_processes) {
        m_scattering_processes_exe.push_back(p.executor());
    }
    for (auto const& p : m_ionization_processes) {
        m_ionization_processes_exe.push_back(p.executor());
    }
#endif

    pp_collision_name.get("ground_rho_index", m_ground_rho_index);
#ifdef MCC_EXCITATION
    m_have_excitation = false;
    m_excitation_product = m_ground_species;
    m_excitation_rho_index = m_ground_rho_index;

    pp_collision_name.query("have_excitation", m_have_excitation);
    pp_collision_name.query("excitation_product", m_excitation_product);
    pp_collision_name.query("excitation_rho_index", m_excitation_rho_index);
#endif
}

/** Calculate the maximum collision frequency using a fixed energy grid that
 *  ranges from 1e-4 to 5000 eV in 0.2 eV increments
 */
amrex::ParticleReal
CoupledBackgroundMCCCollision::get_nu_max (
    amrex::Vector<ScatteringProcess> const& mcc_processes) const {
    using namespace amrex::literals;
    amrex::ParticleReal nu, nu_max = 0.0;
    amrex::ParticleReal E_start = 1e-4_prt;
    amrex::ParticleReal E_end = 5000._prt;
    amrex::ParticleReal E_step = 0.2_prt;

    // set the energy limits and step size for calculating nu_max based
    // on the given cross-section inputs
    for (const auto& process : mcc_processes) {
        auto energy_lo = process.getMinEnergyInput();
        E_start = (energy_lo < E_start) ? energy_lo : E_start;
        auto energy_hi = process.getMaxEnergyInput();
        E_end = (energy_hi > E_end) ? energy_hi : E_end;
        auto energy_step = process.getEnergyInputStep();
        E_step = (energy_step < E_step) ? energy_step : E_step;
    }

    amrex::ParticleReal E = E_start;
    while (E < E_end) {
        amrex::ParticleReal sigma_E = 0.0;

        // loop through all collision pathways
        for (const auto& scattering_process : mcc_processes) {
            // get collision cross-section
            sigma_E += scattering_process.getCrossSection(E);
        }

        // calculate collision frequency
        nu = (m_max_background_density *
              std::sqrt(2.0_prt / m_mass1 * PhysConst::q_e) * sigma_E *
              std::sqrt(E));
        nu_max = std::max(nu_max, nu);
        E += E_step;
    }
    return nu_max;
}

void
CoupledBackgroundMCCCollision::doCollisions (amrex::Real cur_time,
                                             amrex::Real dt,
                                             MultiParticleContainer* mypc) {
    ABLASTR_PROFILE("CoupledBackgroundMCCCollision::doCollisions()");
    using namespace amrex::literals;

    auto& species1 = mypc->GetParticleContainerFromName(m_species_names[0]);
    // this is a very ugly hack to have species2 be a reference and be
    // defined in the scope of doCollisions
    auto& species2 =
        ((m_species_names.size() == 2)
             ? mypc->GetParticleContainerFromName(m_species_names[1])
             : mypc->GetParticleContainerFromName(m_species_names[0]));

    WarpX& warpx_instance = WarpX::GetInstance();
    constexpr int base_lev = 0;
    auto const& geom = warpx_instance.Geom(base_lev);
    const int ncell = geom.Domain().length(0);
    const amrex::ParticleReal inv_gap = static_cast<amrex::ParticleReal>(
        static_cast<amrex::Real>(ncell) / geom.ProbLength(0));

    amrex::ParticleReal elec_weight;
    amrex::ParmParse pp_mc("my_constants");
    pp_mc.getWithParser("elec_weight", elec_weight);
    BackgroundCoupledDensity& m_background_density =
        global_background_density[m_ground_rho_index];

    if (!init_flag) {
        m_mass1 = species1.getMass();

        // calculate maximum collision frequency without ionization
        m_max_background_density = 1.0_prt;
        m_nu_max = get_nu_max(m_scattering_processes);
        m_sigma_max = m_nu_max;

        // calculate total collision probability
        auto coll_n = m_nu_max * dt;
        m_total_collision_prob = 1.0_prt - std::exp(-coll_n);

        // dt has to be small enough that a linear expansion of the collision
        // probability is sufficiently accurately, otherwise the MCC results
        // will be very heavily affected by small changes in the timestep
        if (coll_n > 0.1_prt) {
            ablastr::warn_manager::WMRecordWarning(
                "BackgroundMCC Collisions",
                "dt is too large to ensure accurate MCC results , coll_n: " +
                    std::to_string(coll_n) +
                    " is > 0.1 and collision probability is = " +
                    std::to_string(m_total_collision_prob) + "\n");
        }

        if (ionization_flag) {
            // calculate maximum collision frequency for ionization
            m_nu_max_ioniz = get_nu_max(m_ionization_processes);
            m_ioni_sigma_max = m_nu_max_ioniz;
            // calculate total ionization probability
            auto coll_n_ioniz = m_nu_max_ioniz * dt;
            m_total_collision_prob_ioniz = 1.0_prt - std::exp(-coll_n_ioniz);

            if (coll_n_ioniz > 0.1_prt) {
                ablastr::warn_manager::WMRecordWarning(
                    "BackgroundMCC Collisions",
                    "dt is too large to ensure accurate MCC ionization , "
                    "coll_n_ionization: " +
                        std::to_string(coll_n_ioniz) +
                        " is > 0.1 and ionization probability is = " +
                        std::to_string(m_total_collision_prob_ioniz) + "\n");
            }

            // if an ionization process is included the secondary species mass
            // is taken as the background mass
            m_background_mass = species2.getMass();
        }
        // if no neutral species mass was specified and ionization is not
        // included assume that the collisions will be with neutrals of the
        // same mass as the colliding species (as in ion-neutral collisions)
        else if (m_background_mass == -1) {
            m_background_mass = species1.getMass();
        }

        amrex::Print() << Utils::TextMsg::Info(
            "Setting up Monte-Carlo collisions for " + m_species_names[0] +
            " with:\n" + "     total non-ionization collision probability: " +
            std::to_string(m_total_collision_prob) +
            "\n     total ionization collision probability: " +
            std::to_string(m_total_collision_prob_ioniz));

        init_flag = true;
    }

#if !defined(MCC_DENSITY_AVERAGE_USE) && !defined(IONIZATION_SOURCE_INJECT)
    auto& ground_pc = mypc->GetParticleContainer(mypc->getSpeciesID(
        global_background_density[m_ground_rho_index].m_ground_species));
#endif

#if defined(MCC_EXCITATION) && !defined(MCC_DENSITY_AVERAGE_USE) &&            \
    !defined(IONIZATION_SOURCE_INJECT)
    auto& excitation_pc =
        mypc->GetParticleContainer(mypc->getSpeciesID(m_excitation_product));
    MultiFab& m_excitation_rho = global_rho[m_excitation_rho_index];

    const SmartCopyFactory copy_factory_exc(ground_pc, excitation_pc);
    const auto CopyExc = copy_factory_exc.getSmartCopy();
#endif

    amrex::MultiFab& max_density =
        m_background_density.m_background_density_fabs[0];
    m_max_background_density = max_density.max(0);

    //  calculate maximum collision frequency without ionization
    m_nu_max = m_sigma_max * m_max_background_density;

    // calculate total collision probability
    auto coll_n = m_nu_max * dt;
    m_total_collision_prob = 1.0_prt - std::exp(-coll_n);

    amrex::Print() << "max local atom density: " << m_max_background_density
                   << " m^3; total collision probability: "
                   << m_total_collision_prob << "\n";

    if (ionization_flag) {
        // calculate maximum collision frequency for ionization
        m_nu_max_ioniz = m_ioni_sigma_max * m_max_background_density;

        // calculate total ionization probability
        auto coll_n_ioniz = m_nu_max_ioniz * dt;
        m_total_collision_prob_ioniz = 1.0_prt - std::exp(-coll_n_ioniz);

        if (coll_n_ioniz > 0.1_prt) {
            ablastr::warn_manager::WMRecordWarning(
                "BackgroundMCC Collisions",
                "dt is too large to ensure accurate MCC ionization , "
                "coll_n_ionization: " +
                    std::to_string(coll_n_ioniz) +
                    " is > 0.1 and ionization probability is = " +
                    std::to_string(m_total_collision_prob_ioniz) + "\n");
        }
    }

    // Loop over refinement levels
    const int depos_order = WarpX::nox;
#ifdef IONIZATION_SOURCE_RECORD
    if (ionization_flag) {
        Insert::IonizationSourceRecordCollisionSample(
            geom, dt, elec_weight, warpx_instance.maxLevel());
    }
#endif
#if defined(MCC_EXCITATION) && !defined(MCC_DENSITY_AVERAGE_USE) &&            \
    !defined(IONIZATION_SOURCE_INJECT)
    amrex::Vector<amrex::Vector<amrex::ParticleReal>> pdata(7);
#endif
    auto const flvl = species1.finestLevel();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        flvl == 0,
        "CoupledBackgroundMCCCollision only supports one mesh level.");
    for (int lev = 0; lev <= flvl; ++lev) {
        auto* cost = WarpX::getCosts(lev);

        // firstly loop over particles box by box and do all particle conserving
        // scattering
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

        /*
        auto binIter = m_background_density.m_background_bins[lev].begin();
        auto npIter =
            m_background_density.m_n_particle_in_each_cell[lev].begin();*/
#ifndef MCC_DENSITY_AVERAGE_USE
        auto& background_bin = m_background_density.m_background_bins[lev];
        auto& background_np =
            m_background_density.m_n_particle_in_each_cell[lev];
#endif
        auto cell_size = warpx_instance.CellSize(lev);
        amrex::XDim3 const inv_cell_size{1.0_rt / cell_size[0],
                                         1.0_rt / cell_size[1],
                                         1.0_rt / cell_size[2]};

        amrex::MultiFab& ground_density =
            m_background_density.m_background_density_fabs[lev];

        for (WarpXParIter pti(species1, lev); pti.isValid(); ++pti) {
            if (cost && WarpX::load_balance_costs_update_algo ==
                            LoadBalanceCostsUpdateAlgo::Timers) {
                amrex::Gpu::synchronize();
            }
            auto wt = static_cast<amrex::Real>(amrex::second());

            // 调用耦合函数
            /**********************准备需要用到的数组***************************/
#ifndef MCC_DENSITY_AVERAGE_USE
            const int box_index = pti.index();
            auto& bin = background_bin[box_index];
            auto& ptile = ground_pc.ParticlesAt(lev, pti);
            const int* offsets = bin.offsetsPtr();
            int* indices = bin.permutationPtr();
            int np = ptile.numParticles();
            long const numbins = bin.numBins();
#else
            long const numbins = 1;
#endif
            amrex::Gpu::DeviceVector<int> num_delete(numbins, 0);
            int* p_delete = num_delete.dataPtr();
            // 记录每个cell中的原子数
#ifndef MCC_DENSITY_AVERAGE_USE
            int* p_particle_num = background_np[box_index].dataPtr();
#else
            int* p_particle_num = p_delete;
#endif
            // npIter++;
            // binIter++;
#if defined(MCC_EXCITATION) && !defined(MCC_DENSITY_AVERAGE_USE) &&            \
    !defined(IONIZATION_SOURCE_INJECT)
            amrex::Gpu::DeviceVector<int> mask(np, 0);
            int* p_mask = mask.dataPtr();
#endif
            /*****************************************************************/
            if (depos_order == 1) {
                doBackgroundCollisionsWithinTileCouple<1>(
                    pti, cur_time, ground_density, p_delete, p_particle_num,
                    ncell, inv_cell_size);
#if defined(MCC_EXCITATION) && !defined(MCC_DENSITY_AVERAGE_USE) &&            \
    !defined(IONIZATION_SOURCE_INJECT)
                if (m_have_excitation) {
                    ReplaceParticlesEachCell<1>(p_delete, offsets, indices,
                                                numbins, ground_pc, pti,
                                                m_ground_rho, m_excitation_rho,
                                                excitation_pc, inv_gap, p_mask);
                }
#endif
            } else if (depos_order == 2) {
                doBackgroundCollisionsWithinTileCouple<2>(
                    pti, cur_time, ground_density, p_delete, p_particle_num,
                    ncell, inv_cell_size);
#if defined(MCC_EXCITATION) && !defined(MCC_DENSITY_AVERAGE_USE) &&            \
    !defined(IONIZATION_SOURCE_INJECT)
                if (m_have_excitation) {
                    ReplaceParticlesEachCell<2>(p_delete, offsets, indices,
                                                numbins, ground_pc, pti,
                                                m_ground_rho, m_excitation_rho,
                                                excitation_pc, inv_gap, p_mask);
                }
#endif
            } else if (depos_order == 3) {
                doBackgroundCollisionsWithinTileCouple<3>(
                    pti, cur_time, ground_density, p_delete, p_particle_num,
                    ncell, inv_cell_size);
#if defined(MCC_EXCITATION) && !defined(MCC_DENSITY_AVERAGE_USE) &&            \
    !defined(IONIZATION_SOURCE_INJECT)
                if (m_have_excitation) {
                    ReplaceParticlesEachCell<3>(p_delete, offsets, indices,
                                                numbins, ground_pc, pti,
                                                m_ground_rho, m_excitation_rho,
                                                excitation_pc, inv_gap, p_mask);
                }
#endif
            } else if (depos_order == 4) {
                doBackgroundCollisionsWithinTileCouple<4>(
                    pti, cur_time, ground_density, p_delete, p_particle_num,
                    ncell, inv_cell_size);
#if defined(MCC_EXCITATION) && !defined(MCC_DENSITY_AVERAGE_USE) &&            \
    !defined(IONIZATION_SOURCE_INJECT)
                if (m_have_excitation) {
                    ReplaceParticlesEachCell<4>(p_delete, offsets, indices,
                                                numbins, ground_pc, pti,
                                                m_ground_rho, m_excitation_rho,
                                                excitation_pc, inv_gap, p_mask);
                }
#endif
            } else {
                doBackgroundCollisionsWithinTileCouple<1>(
                    pti, cur_time, ground_density, p_delete, p_particle_num,
                    ncell, inv_cell_size);
#if defined(MCC_EXCITATION) && !defined(MCC_DENSITY_AVERAGE_USE) &&            \
    !defined(IONIZATION_SOURCE_INJECT)
                if (m_have_excitation) {
                    ReplaceParticlesEachCell<1>(p_delete, offsets, indices,
                                                numbins, ground_pc, pti,
                                                m_ground_rho, m_excitation_rho,
                                                excitation_pc, inv_gap, p_mask);
                }
#endif
            }
#if defined(MCC_EXCITATION) && !defined(MCC_DENSITY_AVERAGE_USE) &&            \
    !defined(IONIZATION_SOURCE_INJECT)
            if (m_have_excitation) {
                auto& exc_ptile = excitation_pc.ParticlesAt(lev, pti);
                int np_exc = exc_ptile.numParticles();

                const auto num_added = filterCopyTransformParticles<1>(
                    excitation_pc, exc_ptile, ptile, p_mask, np_exc, CopyExc);
                setNewParticleIDs(exc_ptile, np_exc, num_added);
            }
#endif
            if (cost && WarpX::load_balance_costs_update_algo ==
                            LoadBalanceCostsUpdateAlgo::Timers) {
                amrex::Gpu::synchronize();
                wt = static_cast<amrex::Real>(amrex::second()) - wt;
                amrex::HostDevice::Atomic::Add(&(*cost)[pti.index()], wt);
            }
        }
        // secondly perform ionization through the SmartCopyFactory if needed
#if defined(MCC_EXCITATION) && !defined(MCC_DENSITY_AVERAGE_USE) &&            \
    !defined(IONIZATION_SOURCE_INJECT)
        ground_pc.deleteInvalidParticles();
#endif
        if (ionization_flag) {
            if (depos_order == 1) {
                doBackgroundIonizationCouple<1>(lev, cost, species1, species2,
                                                cur_time, ground_density,
                                                inv_cell_size, elec_weight);
            } else if (depos_order == 2) {
                doBackgroundIonizationCouple<2>(lev, cost, species1, species2,
                                                cur_time, ground_density,
                                                inv_cell_size, elec_weight);
            } else if (depos_order == 3) {
                doBackgroundIonizationCouple<3>(lev, cost, species1, species2,
                                                cur_time, ground_density,
                                                inv_cell_size, elec_weight);
            } else if (depos_order == 4) {
                doBackgroundIonizationCouple<4>(lev, cost, species1, species2,
                                                cur_time, ground_density,
                                                inv_cell_size, elec_weight);
            } else {
                doBackgroundIonizationCouple<1>(lev, cost, species1, species2,
                                                cur_time, ground_density,
                                                inv_cell_size, elec_weight);
            }
        }
    }
}

// 耦合到本地密度
template <int depos_order>
void
CoupledBackgroundMCCCollision::doBackgroundCollisionsWithinTileCouple (
    WarpXParIter& pti, amrex::Real t, amrex::MultiFab& ground_rho,
    [[maybe_unused]] int* p_delete, [[maybe_unused]] int* p_particle_num,
    [[maybe_unused]] int ncell, [[maybe_unused]] amrex::XDim3 inv_cell_size) {
    using namespace amrex::literals;

    // So that CUDA code gets its intrinsic, not the host-only C++ library
    // version
    using std::sqrt;

    // get particle count
    const long np = pti.numParticles();

    // get parsers for the background density and temperature
    auto T_a_func = m_background_temperature_func;

    // get collision parameters
    auto* scattering_processes = m_scattering_processes_exe.data();
    auto const process_count =
        static_cast<int>(m_scattering_processes_exe.size());

    auto const total_collision_prob = m_total_collision_prob;
    auto const nu_max = m_nu_max;

    // store projectile and target masses
    auto const m = m_mass1;
    auto const M = m_background_mass;

    // we need particle positions in order to calculate the local density
    // and temperature
    auto GetPosition = GetParticlePosition<PIdx>(pti);

    // get Struct-Of-Array particle data, also called attribs
    auto& attribs = pti.GetAttribs();
    amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

    // 统计原子数密度 仅单个box
    auto& fab = ground_rho[pti.index()];
    auto const& ground_rho_arr = fab.array();

    // this box is after grow
    amrex::Box const box = fab.box();
    const amrex::XDim3 xyzmin = WarpX::LowerCorner(box, 0, 0._rt); // 待修改
    const amrex::Dim3 lo = lbound(box);

    amrex::ParallelForRNG(np, [=] AMREX_GPU_HOST_DEVICE(
                                  long ip, amrex::RandomEngine const& engine) {
        // determine if this particle should collide
        if (amrex::Random(engine) > total_collision_prob) {
            return;
        }

        // The background temperature parser takes Cartesian coordinates as
        // arguments, in all geometries. The coupled density interpolation
        // below only supports Cartesian geometries, where these are also the
        // stored coordinates used for the grid indexing.
        // 1D 物理坐标z方向，实际用到x
        // 2D 物理xz，实际xy
        amrex::ParticleReal x, y, z;
        GetPosition(ip, x, y, z);
        const amrex::Real rpx = (x - xyzmin.x) * inv_cell_size.x,
                          rpy = (y - xyzmin.y) * inv_cell_size.y,
                          rpz = (z - xyzmin.z) * inv_cell_size.z;
        amrex::ParticleReal n_a = 0;
        Compute_shape_factor<depos_order> const compute_shape_factor;
        amrex::Real sx[depos_order + 1] = {0._rt},
                                     sy[depos_order + 1] = {0._rt},
                                     sz[depos_order + 1] = {0._rt};

        int px = compute_shape_factor(sx, rpx),
            py = compute_shape_factor(sy, rpy),
            pz = compute_shape_factor(sz, rpz);

#if defined(WARPX_DIM_1D_Z)
        for (int ix = 0; ix <= depos_order; ix++) {
            n_a += sx[ix] * ground_rho_arr(lo.x + px + ix, 0, 0);
        }
#elif defined(WARPX_DIM_XZ)
            for (int iy = 0; iy <= depos_order; iy++) {
                for (int ix = 0; ix <= depos_order; ix++) {
                    n_a += sx[ix] * sy[iy] *
                           ground_rho_arr(lo.x + px + ix, lo.y + py + iy, 0);
                }
            }
#elif defined(WARPX_DIM_3D)
            for (int iz = 0; iz <= depos_order; iz++) {
                for (int iy = 0; iy <= depos_order; iy++) {
                    for (int ix = 0; ix <= depos_order; ix++) {
                        n_a += sx[ix] * sy[iy] * sz[iz] *
                               ground_rho_arr(lo.x + px + ix, lo.y + py + iy,
                                              lo.z + pz + iz);
                    }
                }
            }
#endif

        const amrex::ParticleReal T_a = T_a_func(x, y, z, t);

        amrex::ParticleReal v_coll, v_coll2, sigma_E, nu_i = 0;
        double gamma, E_coll;
        amrex::ParticleReal ua_x, ua_y, ua_z, vx, vy, vz;
        const amrex::ParticleReal col_select = amrex::Random(engine);

        // get velocities of gas particles from a Maxwellian distribution
        auto const vel_std = sqrt(PhysConst::kb * T_a / M);
        ua_x = vel_std * amrex::RandomNormal(0_prt, 1.0_prt, engine);
        ua_y = vel_std * amrex::RandomNormal(0_prt, 1.0_prt, engine);
        ua_z = vel_std * amrex::RandomNormal(0_prt, 1.0_prt, engine);

        // we assume the target particle is not relativistic (in
        // the lab frame) and therefore we can transform the projectile
        // velocity to a frame in which the target is stationary with
        // a simple Galilean boost
        // not doing the full Lorentz boost here saves us computation
        // since most particles will not actually collide
        vx = ux[ip] - ua_x;
        vy = uy[ip] - ua_y;
        vz = uz[ip] - ua_z;
        v_coll2 = (vx * vx + vy * vy + vz * vz);
        v_coll = std::sqrt(v_coll2);

        // calculate the collision energy in eV
        ParticleUtils::getCollisionEnergy(v_coll2, m, M, gamma, E_coll);

        // loop through all collision pathways
        for (int i = 0; i < process_count; i++) {
            auto const& scattering_process = *(scattering_processes + i);

            // get collision cross-section
            sigma_E = scattering_process.getCrossSection(
                static_cast<amrex::ParticleReal>(E_coll));

            // calculate normalized collision frequency
            nu_i += n_a * sigma_E * v_coll / nu_max;

            // check if this collision should be performed
            if (col_select > nu_i) {
                continue;
            }

            // For the coupled excitation model (MCC_EXCITATION without
            // MCC_DENSITY_AVERAGE_USE), the projectile velocity is only
            // scattered (and the cell bookkeeping updated) if a background
            // particle is available in the current cell.
            bool update_velocity = true;
#if defined(MCC_EXCITATION) && !defined(MCC_DENSITY_AVERAGE_USE)
            if (scattering_process.m_type ==
                ScatteringProcessType::EXCITATION) {
                int pos = Insert::CollisionDetail::DenseBinIndex(
                    amrex::Dim3{px, py, pz}, amrex::Dim3{ncell, ncell, ncell});
                update_velocity = (p_particle_num[pos] > 0);
                if (update_velocity) {
                    amrex::Gpu::Atomic::Add(&p_particle_num[pos], -1);
                    amrex::Gpu::Atomic::Add(&p_delete[pos], 1);
                    // 需要CAS来解决冲突，需要处理网格索引偏差
                }
            }
#endif

            if (update_velocity) {
                // At this point the given particle has been chosen for a
                // collision with a background-gas particle of velocity
                // (ua_x, ua_y, ua_z). Compute the post-collision momentum of
                // the projectile using conservation of energy and momentum.
                // The angular distribution in the center-of-mass frame is set
                // by the process's scattering angle model, and any inelastic
                // energy loss is passed as the (released) reaction energy.
                // The background particle is treated as a reservoir: its recoil
                // is computed as the second product but discarded.
                amrex::ParticleReal u1x_out, u1y_out, u1z_out;
                amrex::ParticleReal u2x_out, u2y_out, u2z_out;
                TwoProductComputeProductMomenta(
                    ux[ip], uy[ip], uz[ip], m, ua_x, ua_y, ua_z, M, u1x_out,
                    u1y_out, u1z_out, m, u2x_out, u2y_out, u2z_out, M,
                    -scattering_process.m_energy_penalty * PhysConst::q_e,
                    // TwoProductComputeProductMomenta expects the *released*
                    // energy here, hence the negative sign; the energy penalty
                    // is also converted from eV to Joules.
                    scattering_process.m_scattering_angle_model, engine);

                // update projectile velocity with new components in labframe
                // (the background-gas recoil u2*_out is discarded)
                ux[ip] = u1x_out;
                uy[ip] = u1y_out;
                uz[ip] = u1z_out;
            }
            break;
        }
    });
}

template <int depos_order>
void
CoupledBackgroundMCCCollision::doBackgroundIonizationCouple (
    int lev, amrex::LayoutData<amrex::Real>* cost,
    WarpXParticleContainer& species1, WarpXParticleContainer& species2,
    [[maybe_unused]] amrex::Real t, amrex::MultiFab& ground_rho,
    const amrex::XDim3& inv_cell_size, amrex::ParticleReal elec_weight) {
    using namespace amrex::literals;
#ifdef MCC_DENSITY_AVERAGE_USE
    amrex::ignore_unused(elec_weight);
#endif

    ABLASTR_PROFILE(
        "CoupledBackgroundMCCCollision::doBackgroundIonizationCouple()");

    const SmartCopyFactory copy_factory_elec(species1, species1);
    const SmartCopyFactory copy_factory_ion(species1, species2);
    const auto CopyElec = copy_factory_elec.getSmartCopy();
    const auto CopyIon = copy_factory_ion.getSmartCopy();

    const auto Filter = CoupledImpactIonizationFilterFunc(
        m_ionization_processes[0], m_mass1, m_total_collision_prob_ioniz,
        m_nu_max_ioniz);

#ifndef IONIZATION_SOURCE_INJECT
    const amrex::ParticleReal sqrt_kb_m =
        std::sqrt(PhysConst::kb / m_background_mass);
#endif

#if !defined(MCC_DENSITY_AVERAGE_USE) && !defined(IONIZATION_SOURCE_INJECT)
    BackgroundCoupledDensity& m_background_density =
        global_background_density[m_ground_rho_index];
    m_background_density.resetIonizationConsumption(lev);
    WarpX& warpx_instance = WarpX::GetInstance();
    auto& pc = warpx_instance.GetPartContainer().GetParticleContainerFromName(
        m_background_density.m_ground_species);

    auto& background_bin = m_background_density.m_background_bins[lev];
    auto& background_np = m_background_density.m_n_particle_in_each_cell[lev];
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
        for (WarpXParIter pti(species1, lev); pti.isValid(); ++pti) {

            if (cost && WarpX::load_balance_costs_update_algo ==
                            LoadBalanceCostsUpdateAlgo::Timers) {
                amrex::Gpu::synchronize();
            }
            auto wt = static_cast<amrex::Real>(amrex::second());

            auto& elec_tile = species1.ParticlesAt(lev, pti);
            auto& ion_tile = species2.ParticlesAt(lev, pti);

            const auto np_elec = elec_tile.numParticles();
            const auto np_ion = ion_tile.numParticles();

#ifndef IONIZATION_SOURCE_INJECT
            auto Transform = ImpactIonizationTransformFunc(
                m_ionization_processes[0].getEnergyPenalty(), m_mass1,
                sqrt_kb_m, m_background_temperature_func, t);
#else
            auto Transform =
                Insert::CollisionDetail::IonizationSourceTransformFunc(
                    m_ionization_processes[0].getEnergyPenalty(), m_mass1);
#endif

            // 获取粒子网格信息
            const int box_index = pti.index();
            long numbins = 1;
#if !defined(MCC_DENSITY_AVERAGE_USE) && !defined(IONIZATION_SOURCE_INJECT)
            auto& ptile = pc.ParticlesAt(lev, pti);
            auto& bin = background_bin[box_index];
            const int* offsets = bin.offsetsPtr();
            int* indices = bin.permutationPtr();
            numbins = bin.numBins();

            auto& soa = ptile.GetStructOfArrays();
            auto get_atom_postion = GetParticlePosition<PIdx>(ptile);
            uint64_t* const AMREX_RESTRICT idcpu = soa.GetIdCPUData().data();
            auto& soa_arr = soa.GetRealData();
            amrex::Real* pw = soa_arr[PIdx::w].dataPtr();
#endif

            amrex::Gpu::DeviceVector<int> num_delete(numbins, 0);
            int* p_delete = num_delete.dataPtr();

#if !defined(MCC_DENSITY_AVERAGE_USE) && !defined(IONIZATION_SOURCE_INJECT)
            // 记录每个cell中的原子数，存储原始数据
            amrex::Gpu::DeviceVector<int> particle_num_in_cell_origin(numbins,
                                                                      0);
            int *p_particle_num = background_np[box_index].dataPtr(),
                *p_particle_num_origin = particle_num_in_cell_origin.dataPtr();
            amrex::Gpu::copy(amrex::Gpu::deviceToDevice,
                             background_np[box_index].begin(),
                             background_np[box_index].end(),
                             particle_num_in_cell_origin.begin());
#else
            int* p_particle_num = p_delete;
#endif

            auto& fab = ground_rho[box_index];
            auto const& rho_arr = fab.array();

#if !defined(MCC_DENSITY_AVERAGE_USE) && !defined(IONIZATION_SOURCE_INJECT)
            auto& consumption_fab =
                m_background_density
                    .m_ionization_consumption_fabs[lev][box_index];
            auto const& consumption_arr = consumption_fab.array();
#endif

#ifdef MCC_ION_CACHED_RHO
            // Cached ion charge density target for this tile. The cache is
            // synchronized after the neutral-consumption loop below by
            // depositing the newly created ions at their own positions and
            // weights; the consumption loop itself only updates the
            // neutral-consumption MultiFab.
            auto& ion_fab = species2.m_cached_rho[lev][box_index];
            auto const& ion_rho_arr = ion_fab.array();
            amrex::ParticleReal const ion_charge = species2.getCharge();
#endif

            /*
             * The box of ground[pti] is correct after verifying.
             * Maybe it is gotten by the number index.
             * Note that ground[pti].box() have ghost cell.
             */

            amrex::Box box = fab.box();
            const amrex::XDim3 xyzmin = WarpX::LowerCorner(box, lev, 0._rt);
#if !defined(MCC_DENSITY_AVERAGE_USE) && !defined(IONIZATION_SOURCE_INJECT)
            const amrex::Dim3 lo = lbound(box);
#endif

            const auto num_added = filterCopyTransformParticles<1, depos_order>(
                species1, species2, elec_tile, ion_tile, elec_tile, np_elec,
                np_ion, Filter, CopyElec, CopyIon, Transform, p_delete, rho_arr,
                p_particle_num, xyzmin, box, ground_rho.nGrowVect(),
                inv_cell_size
#if defined(MCC_DENSITY_AVERAGE_USE) || defined(IONIZATION_SOURCE_INJECT)
                ,
                false
#endif
            );

#ifdef IONIZATION_SOURCE_RECORD
            Insert::IonizationSourceDepositNewElectrons(elec_tile, np_elec,
                                                        num_added, elec_weight);
#endif

#if !defined(MCC_DENSITY_AVERAGE_USE) && !defined(IONIZATION_SOURCE_INJECT)
            amrex::Gpu::DeviceScalar<int> all_deleted(0);
            int* p_num = all_deleted.dataPtr();

            amrex::Real const inv_vol =
                inv_cell_size.x * inv_cell_size.y * inv_cell_size.z;
            if (num_added > 0) {
                amrex::ParallelForRNG(
                    numbins, [=] AMREX_GPU_DEVICE(
                                 long ibin, const amrex::RandomEngine& engine) {
                        const int offset_start = offsets[ibin],
                                  offset_end = offsets[ibin + 1],
                                  np_in_cell = offset_end - offset_start,
                                  np_by_weight = p_particle_num_origin[ibin];
                        int rest = amrex::min(p_delete[ibin], np_by_weight);
                        amrex::Gpu::Atomic::Add(p_num, rest);

                        while (rest > 0) {
                            int const indices_pos =
                                offset_start +
                                amrex::Random_int(np_in_cell, engine);
                            int const pos = indices[indices_pos];
                            auto pidw = amrex::ParticleIDWrapper{idcpu[pos]};
                            if (pidw.is_valid() && pw[pos] > 10.0_prt) {
                                // 减去权重，原则上在细胞内进行操作，每个线程处理完全不同的集合
                                pw[pos] -= elec_weight;
                                rest--;

                                // Deposit positive neutral consumption. It is
                                // summed across boxes after all tile kernels
                                // have completed.
                                amrex::ParticleReal x, y, z;
                                amrex::ParticleReal const consumed_density =
                                    elec_weight * inv_vol;
                                get_atom_postion(pos, x, y, z);

                                const amrex::ParticleReal rpx = (x - xyzmin.x) *
                                                                inv_cell_size.x,
                                                          rpy = (y - xyzmin.y) *
                                                                inv_cell_size.y,
                                                          rpz = (z - xyzmin.z) *
                                                                inv_cell_size.z;
                                Compute_shape_factor<depos_order> const
                                    compute_shape_factor;

                                amrex::Real sx[depos_order +
                                               1] = {0._rt},
                                               sy[depos_order + 1] = {0._rt},
                                               sz[depos_order + 1] = {0._rt};
                                int px = compute_shape_factor(sx, rpx),
                                    py = compute_shape_factor(sy, rpy),
                                    pz = compute_shape_factor(sz, rpz);

#if defined(WARPX_DIM_1D_Z)
                                for (int ix = 0; ix <= depos_order; ix++) {
                                    amrex::Gpu::Atomic::AddNoRet(
                                        &consumption_arr(lo.x + px + ix, 0, 0),
                                        sx[ix] * consumed_density);
                                }
#elif defined(WARPX_DIM_XZ)
                        for (int iy = 0; iy <= depos_order; iy++) {
                            for (int ix = 0; ix <= depos_order; ix++) {
                                amrex::Gpu::Atomic::AddNoRet(
                                    &consumption_arr(lo.x + px + ix,
                                                     lo.y + py + iy, 0),
                                    sx[ix] * sy[iy] * consumed_density);
                            }
                        }
#elif defined(WARPX_DIM_3D)
                        for (int iz = 0; iz <= depos_order; iz++) {
                            for (int iy = 0; iy <= depos_order; iy++) {
                                for (int ix = 0; ix <= depos_order; ix++) {
                                    amrex::Gpu::Atomic::AddNoRet(
                                        &consumption_arr(
                                            lo.x + px + ix, lo.y + py + iy,
                                            lo.z + pz + iz),
                                        sx[ix] * sy[iy] * sz[iz] *
                                            consumed_density);
                                }
                            }
                        }
#endif
                            }
                        }
                    });
            }

#ifdef MCC_ION_CACHED_RHO
            // Deposit the charge of the newly created ions into the cached ion
            // charge density, using their own positions and weights, so that
            // the cache stays consistent with a fresh deposit between refresh
            // steps. The new ions occupy the contiguous slot range
            // [np_ion, np_ion + num_added) of ion_tile; the position/weight
            // accessors must be constructed after filterCopyTransformParticles
            // because the tile resize may reallocate the underlying arrays.
            if (num_added > 0) {
                auto& ion_soa = ion_tile.GetStructOfArrays();
                amrex::ParticleReal const* const AMREX_RESTRICT pw_ion =
                    ion_soa.GetRealData()[PIdx::w].dataPtr() + np_ion;
                auto const get_ion_position =
                    GetParticlePosition<PIdx>(ion_tile, np_ion);
                amrex::ParticleReal const q_inv_vol = ion_charge * inv_vol;

                amrex::ParallelFor(num_added, [=] AMREX_GPU_DEVICE(long i) {
                    amrex::ParticleReal x, y, z;
                    get_ion_position(i, x, y, z);

                    const amrex::ParticleReal rpx = (x - xyzmin.x) *
                                                    inv_cell_size.x,
                                              rpy = (y - xyzmin.y) *
                                                    inv_cell_size.y,
                                              rpz = (z - xyzmin.z) *
                                                    inv_cell_size.z;
                    Compute_shape_factor<depos_order> const
                        compute_shape_factor;
                    amrex::Real sx[depos_order + 1] = {0._rt},
                                                 sy[depos_order + 1] = {0._rt},
                                                 sz[depos_order + 1] = {0._rt};
                    int const px = compute_shape_factor(sx, rpx),
                              py = compute_shape_factor(sy, rpy),
                              pz = compute_shape_factor(sz, rpz);

                    amrex::ParticleReal const ion_charge_density =
                        pw_ion[i] * q_inv_vol;

#if defined(WARPX_DIM_1D_Z)
                    for (int ix = 0; ix <= depos_order; ix++) {
                        amrex::Gpu::Atomic::AddNoRet(
                            &ion_rho_arr(lo.x + px + ix, 0, 0),
                            sx[ix] * ion_charge_density);
                    }
#elif defined(WARPX_DIM_XZ)
                        for (int iy = 0; iy <= depos_order; iy++) {
                            for (int ix = 0; ix <= depos_order; ix++) {
                                amrex::Gpu::Atomic::AddNoRet(
                                    &ion_rho_arr(lo.x + px + ix,
                                                 lo.y + py + iy, 0),
                                    sx[ix] * sy[iy] * ion_charge_density);
                            }
                        }
#elif defined(WARPX_DIM_3D)
                        for (int iz = 0; iz <= depos_order; iz++) {
                            for (int iy = 0; iy <= depos_order; iy++) {
                                for (int ix = 0; ix <= depos_order; ix++) {
                                    amrex::Gpu::Atomic::AddNoRet(
                                        &ion_rho_arr(
                                            lo.x + px + ix, lo.y + py + iy,
                                            lo.z + pz + iz),
                                        sx[ix] * sy[iy] * sz[iz] *
                                            ion_charge_density);
                                }
                            }
                        }
#endif
                });
            }
#endif

            amrex::AllPrint()
                << "rank " << amrex::ParallelDescriptor::MyProc() << ": lev "
                << lev << " box " << box_index << " ionization: delete "
                << all_deleted.dataValue() << " particle"
                << pc.getSpeciesId() + 1 << "\n";
#endif

#ifndef IONIZATION_SOURCE_INJECT
            setNewParticleIDs(elec_tile, np_elec, num_added);
            setNewParticleIDs(ion_tile, np_ion, num_added);
#endif

            if (cost && WarpX::load_balance_costs_update_algo ==
                            LoadBalanceCostsUpdateAlgo::Timers) {
                amrex::Gpu::synchronize();
                wt = static_cast<amrex::Real>(amrex::second()) - wt;
                amrex::HostDevice::Atomic::Add(&(*cost)[pti.index()], wt);
            }
        }
    }

#if !defined(MCC_DENSITY_AVERAGE_USE) && !defined(IONIZATION_SOURCE_INJECT)
    m_background_density.applyIonizationConsumption(lev);
#endif
}
