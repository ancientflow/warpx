#pragma once

#include "Insert/Injection/CoordinateDistribution.h"
#include "Insert/Injection/PositionSampler.h"
#include "Insert/Injection/InjectionRateModel.h"

#include <AMReX_REAL.H>

#include <memory>
#include <string>
#include <vector>

class WarpX;

namespace Insert {

struct SpeciesVelocityConfig
{
    std::string species_name;
    amrex::ParticleReal macro_weight = amrex::ParticleReal(1.0);
    std::unique_ptr<CoordinateDistribution> velocity_space;
};

/** Higher-layer hole assignment: implements no distribution of its own and
 *  must be paired with a regular position distribution that defines the
 *  within-hole local coordinates. Each particle is dispatched to one hole
 *  on the ring and its sampled position is offset to that hole's center. */
struct HoleArrayPlaneConfig
{
    int hole_count = 48;
    amrex::ParticleReal ring_radius = amrex::ParticleReal(0.0);
};

class ParticleInjectionSource
{
public:
    ParticleInjectionSource (
        std::string source_name,
        std::unique_ptr<InjectionRateModel> rate_model,
        std::unique_ptr<PositionSampler> position_space,
        std::vector<SpeciesVelocityConfig> species,
        amrex::ParticleReal x_offset = amrex::ParticleReal(0.0),
        amrex::ParticleReal y_offset = amrex::ParticleReal(0.0));

    void setHoleArrayPlane (HoleArrayPlaneConfig config);
    void setTimeStepOverride (amrex::Real dt);
    void setBatchMultiplier (amrex::Real batch_multiplier);

    void Inject (WarpX& warpx, amrex::Real dt);

    [[nodiscard]] std::string const& name () const noexcept;

private:
    [[nodiscard]] int consumeParticleCount (amrex::Real expected_macro_particles);

    [[nodiscard]] EmissionSample samplePosition (
        int particle_index, amrex::RandomEngine& engine);

    void addSpeciesParticles (
        WarpX& warpx, SpeciesVelocityConfig const& species,
        std::vector<EmissionSample> const& positions,
        amrex::RandomEngine& normal_engine) const;

    std::string m_source_name;
    std::unique_ptr<InjectionRateModel> m_rate_model;
    std::unique_ptr<PositionSampler> m_position_space;
    std::vector<SpeciesVelocityConfig> m_species;
    FractionalParticleAccumulator m_accumulator;
    amrex::ParticleReal m_x_offset = amrex::ParticleReal(0.0);
    amrex::ParticleReal m_y_offset = amrex::ParticleReal(0.0);
    bool m_use_hole_array_plane = false;
    HoleArrayPlaneConfig m_hole_array_plane;
    int m_hole_start = 0;
    bool m_has_dt_override = false;
    amrex::Real m_dt_override = 0.0;
    amrex::Real m_batch_multiplier = 0.0;
    amrex::Real m_batch_remainder = 0.0;
};

} // namespace Insert
