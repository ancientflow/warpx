#pragma once

#include "Insert/Injection/HallCoordinateDistribution.h"
#include "Insert/Injection/HallPositionSampler.h"
#include "Insert/Injection/HallRateModel.h"

#include <AMReX_REAL.H>

#include <memory>
#include <string>
#include <vector>

class WarpX;

namespace Insert {

struct HallSpeciesVelocityConfig
{
    std::string species_name;
    amrex::ParticleReal macro_weight = amrex::ParticleReal(1.0);
    std::unique_ptr<HallCoordinateDistribution> velocity_space;
};

struct HallHoleArrayPlaneConfig
{
    int hole_count = 48;
    amrex::ParticleReal ring_radius = amrex::ParticleReal(0.0);
    amrex::ParticleReal hole_radius = amrex::ParticleReal(0.0);
    amrex::ParticleReal z = amrex::ParticleReal(0.0);
};

class HallInjectionSource
{
public:
    HallInjectionSource (
        std::string source_name,
        std::unique_ptr<HallRateModel> rate_model,
        std::unique_ptr<HallPositionSampler> position_space,
        std::vector<HallSpeciesVelocityConfig> species,
        amrex::ParticleReal x_offset = amrex::ParticleReal(0.0),
        amrex::ParticleReal y_offset = amrex::ParticleReal(0.0));

    void setHoleArrayPlane (HallHoleArrayPlaneConfig config);
    void setTimeStepOverride (amrex::Real dt);
    void setBatchMultiplier (amrex::Real batch_multiplier);

    void Inject (WarpX& warpx, amrex::Real dt);

    [[nodiscard]] std::string const& name () const noexcept;

private:
    [[nodiscard]] int consumeParticleCount (amrex::Real expected_macro_particles);

    [[nodiscard]] EmissionSample samplePosition (
        int particle_index, amrex::RandomEngine& engine);

    void addSpeciesParticles (
        WarpX& warpx, HallSpeciesVelocityConfig const& species,
        std::vector<EmissionSample> const& positions,
        amrex::RandomEngine& normal_engine) const;

    std::string m_source_name;
    std::unique_ptr<HallRateModel> m_rate_model;
    std::unique_ptr<HallPositionSampler> m_position_space;
    std::vector<HallSpeciesVelocityConfig> m_species;
    HallFractionalParticleAccumulator m_accumulator;
    amrex::ParticleReal m_x_offset = amrex::ParticleReal(0.0);
    amrex::ParticleReal m_y_offset = amrex::ParticleReal(0.0);
    bool m_use_hole_array_plane = false;
    HallHoleArrayPlaneConfig m_hole_array_plane;
    int m_hole_start = 0;
    bool m_has_dt_override = false;
    amrex::Real m_dt_override = 0.0;
    amrex::Real m_batch_multiplier = 0.0;
    amrex::Real m_batch_remainder = 0.0;
};

} // namespace Insert
