#pragma once

#include "Insert/Config/WarpXFunctionConfig.h"
#include "Particles/MultiParticleContainer.H"
#include "Utils/ParticleUtils.H"

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <map>
#include <string>
#include <utility>

class BackgroundCoupledDensity {
  public:
    // Collision-owned neutral-density state. Rebuilt at actual direct-MCC invocations
    // so migration, injection and particle reordering cannot leave stale bins.
    std::string m_ground_species;
    amrex::Vector<amrex::MultiFab> m_background_density_fabs;
    amrex::Vector<amrex::MultiFab> m_ionization_consumption_fabs;
#ifdef MCC_DENSITY_AVERAGE_CALC
    amrex::Vector<amrex::MultiFab> m_background_density_sum_fabs;
    int m_average_sample_count = 0;
    int m_average_steps_per_period = 1;
    int m_average_periods = 1;
    int m_average_window_steps = 1;
    int m_raw_output_interval = 0;
    std::string m_output_dir = "background_density_fab";
    std::string m_output_fab;
#endif
#ifdef MCC_DENSITY_AVERAGE_USE
    std::string m_input_fab;
#endif
    struct TileData {
        amrex::DenseBins<ParticleUtils::ParticleTileDataType> bins;
        amrex::Gpu::DeviceVector<int> available;
    };
    // Particle tiles can be empty or share a FAB on CPU. Key by both indices.
    amrex::Vector<std::map<std::pair<int, int>, TileData>> m_tiles;

    /**
     * @brief init the vector of background density species
     */
    void backgroundDensityInit ();

    /**
     * @brief update the background density data
     */
    void backgroundDensityUpdate (MultiParticleContainer& mypc,
                                  amrex::ParticleReal elec_weight, int step);

    /**
     * @brief reset the ionization-consumption density before deposition
     */
    void resetIonizationConsumption (int lev);

    /**
     * @brief sum ionization consumption across boxes and subtract it from the
     * background density
     *
     * This updates only the neutral background density. Cached ion charge
     * density is intentionally kept on its existing update path.
     */
    void applyIonizationConsumption (int lev);

    /**
     * @brief finalize background density data
     */
    void backgroundDensityFinalize ();

    /**
     * @brief delete the particles with zero weight
     */
    void backgroundSpeciesClean (MultiParticleContainer& mypc) const;
};
