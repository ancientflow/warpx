/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_INSERT_COUPLEDBACKGROUNDMCCCOLLISION_H_
#define WARPX_INSERT_COUPLEDBACKGROUNDMCCCOLLISION_H_

#include "Insert/Config/WarpXFunctionConfig.h"
#include "Particles/Collision/CollisionBase.H"
#include "Particles/Collision/ScatteringProcess.H"
#include "Particles/MultiParticleContainer.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>

class CoupledBackgroundMCCCollision final : public CollisionBase {
  public:
    explicit CoupledBackgroundMCCCollision (std::string const& collision_name);

    ~CoupledBackgroundMCCCollision () override = default;

    CoupledBackgroundMCCCollision (CoupledBackgroundMCCCollision const&) =
        delete;
    CoupledBackgroundMCCCollision&
    operator=(CoupledBackgroundMCCCollision const&) = delete;
    CoupledBackgroundMCCCollision (CoupledBackgroundMCCCollision&&) = delete;
    CoupledBackgroundMCCCollision&
    operator=(CoupledBackgroundMCCCollision&&) = delete;

    [[nodiscard]] amrex::ParticleReal
    get_nu_max (amrex::Vector<ScatteringProcess> const& mcc_processes) const;

    /** Perform the collisions
     *
     * @param cur_time Current time
     * @param dt Time step size
     * @param mypc Container of species involved
     *
     */
    void doCollisions (amrex::Real cur_time, amrex::Real dt,
                       MultiParticleContainer* mypc) override;

    /** Perform particle conserving MCC collisions within a tile
     *
     * @param pti particle iterator
     * @param t current time
     * @param ground_rho local background density
     * @param p_delete pointer of particles to be deleted
     * @param p_particle_num pointer of particles in each cell
     * @param ncell number of cell at an axis
     */
    template <int depos_order>
    void doBackgroundCollisionsWithinTileCouple (WarpXParIter& pti,
                                                 amrex::Real t,
                                                 amrex::MultiFab& ground_rho,
                                                 int* p_delete,
                                                 int* p_particle_num, int ncell,
                                                 amrex::XDim3 inv_cell_size);

    /** Perform MCC ionization interactions with local density
     *
     * @param[in] lev the mesh-refinement level
     * @param[in,out] cost pointer to (load balancing) cost corresponding to box
     * where present particles are ionized.
     * @param[in,out] species1,species2 reference to species container used to
     * inject new particles
     * @param t current time
     * @param ground_rho local background density
     * @param inv_cell_size the reciprocal of the cell size
     * @param elec_weight the weight of electron
     */
    template <int depos_order>
    void doBackgroundIonizationCouple (
        int lev, amrex::LayoutData<amrex::Real>* cost,
        WarpXParticleContainer& species1, WarpXParticleContainer& species2,
        amrex::Real t, amrex::MultiFab& ground_rho,
        const amrex::XDim3& inv_cell_size, amrex::ParticleReal elec_weight);

  private:
    amrex::Vector<ScatteringProcess> m_scattering_processes;
    amrex::Vector<ScatteringProcess> m_ionization_processes;
    amrex::Gpu::DeviceVector<ScatteringProcess::Executor>
        m_scattering_processes_exe;
    amrex::Gpu::DeviceVector<ScatteringProcess::Executor>
        m_ionization_processes_exe;

    bool init_flag = false;
    bool ionization_flag = false;

    amrex::ParticleReal m_mass1;

    amrex::ParticleReal m_max_background_density = 0;
    amrex::ParticleReal m_background_mass;
    amrex::ParticleReal m_total_collision_prob;
    amrex::ParticleReal m_total_collision_prob_ioniz = 0;
    amrex::ParticleReal m_nu_max;
    amrex::ParticleReal m_nu_max_ioniz;

    amrex::Parser m_background_temperature_parser;

    amrex::ParserExecutor<4> m_background_temperature_func;

    // Coupled background-density state.
    int m_ground_rho_index;
    amrex::ParticleReal m_sigma_max;
    amrex::ParticleReal m_ioni_sigma_max;

#ifdef MCC_EXCITATION
    std::string m_excitation_product;
    bool m_have_excitation;
    int m_excitation_rho_index;
#endif
};

#endif // WARPX_INSERT_COUPLEDBACKGROUNDMCCCOLLISION_H_
