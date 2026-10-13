#include "InsertBackgroundDensity.h"

#include "WarpX.H"

#include "Particles/MultiParticleContainer.H"
#include "Insert/Config/WarpXFunctionConfig.h"
#include "Insert/Core/WarpXInsert.h"

#include <AMReX_ParmParse.H>

#ifdef MCC_DENSITY
amrex::Vector<BackgroundCoupledDensity> global_background_density;
#endif

namespace Insert {

void
GlobalBackgroundDensityInit ()
{
#ifdef MCC_DENSITY
    amrex::ParmParse const pp_coll("collisions");
    amrex::Vector<std::string> species_names;
    pp_coll.queryarr("background_species", species_names);

    global_background_density.resize(species_names.size());
    for (int i = 0; i < species_names.size(); i++) {
        global_background_density[i].m_ground_species = species_names[i];
        global_background_density[i].backgroundDensityInit();
    }
#endif
}

void
GlobalBackgroundDensityUpdate (const int step)
{
#ifdef MCC_DENSITY
    if (global_background_density.empty()) {
        return;
    }
#ifdef MCC_DENSITY_AVERAGE_USE
    amrex::ignore_unused(step);
    return;
#else
    amrex::ParmParse pp_mc("my_constants");
    amrex::ParticleReal elec_weight;
    pp_mc.getWithParser("elec_weight", elec_weight);

    WarpX& warpx_instance = WarpX::GetInstance();
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();
    const bool if_update_sort = warpx_instance.sort_intervals.contains(step);
    for (auto& density : global_background_density) {
        int const ndt = ParticleSubcyclingNdt(density.m_ground_species);
        const bool if_update_push = (step % ndt == 1);
        if (if_update_push || if_update_sort) {
            density.backgroundDensityUpdate(mypc, elec_weight, step);
        }
    }
#endif
#else
    amrex::ignore_unused(step);
#endif
}

void
GlobalBackgroundDensityClean (const int step)
{
#ifdef MCC_DENSITY
    if (global_background_density.empty()) {
        return;
    }
#ifdef MCC_DENSITY_AVERAGE_USE
    amrex::ignore_unused(step);
    return;
#else
    WarpX& warpx_instance = WarpX::GetInstance();
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();
    for (auto& density : global_background_density) {
        int const ndt = ParticleSubcyclingNdt(density.m_ground_species);
        const bool if_clean = (step % ndt == 0);
        if (if_clean) {
            density.backgroundSpeciesClean(mypc);
        }
    }
#endif
#else
    amrex::ignore_unused(step);
#endif
}

void
GlobalBackgroundDensityFinalize ()
{
#ifdef MCC_DENSITY
    for (auto& density : global_background_density) {
        density.backgroundDensityFinalize();
    }
#endif
}

} // namespace Insert
