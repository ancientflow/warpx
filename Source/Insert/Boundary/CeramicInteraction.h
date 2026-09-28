#pragma once

#include "Insert/Boundary/WallBehavior.h"
#include "Insert/Boundary/WallInteractionOperators.h"
#include "Insert/Math/ParticleVector.h"
#include "Utils/WarpXConst.H"

#include <AMReX_Array.H>
#include <AMReX_Random.H>

#include <cmath>
#include <type_traits>

namespace Insert {

/** Electron/ceramic collision model from Script/3d_hall_dielectric_slot.
 * Defines event probabilities and outgoing velocity distributions only.
 * Particle creation, charge deposition, deletion and diagnostics are handled
 * by the common boundary driver. Input velocities are physical m/s.
 */
class CeramicInteraction
{
public:
    /** Thermal parameter kT in eV, not mean emitted kinetic energy. */
    explicit CeramicInteraction (amrex::ParticleReal secondary_temperature_eV = 3.0);

    /** Select the event once; its multiplicity is given by WallEmissionCount. */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    WallBehavior SelectEvent (
        ParticleVector const& hit_velocity, amrex::RandomEngine const& engine) const noexcept;

    /** Sample the whole emission event, retaining model-specific spectra and
     * correlations. Only WallEmissionCount(event) entries are used. This ceramic
     * preset draws each secondary independently from a Maxwellian flux at kT,
     * with mean kinetic energy 2 kT and no event-wise energy truncation.
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuArray<ParticleVector, max_wall_emissions> SampleEmission (
        WallBehavior event, ParticleVector const& hit_velocity,
        ParticleVector const& normal, amrex::RandomEngine const& engine) const noexcept;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleVector SampleReflection (
        ParticleVector const& hit_velocity, ParticleVector const& normal) const noexcept;

private:
    amrex::ParticleReal m_thermal_velocity;
};

static_assert(std::is_trivially_copyable_v<CeramicInteraction>);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
WallBehavior CeramicInteraction::SelectEvent (
    ParticleVector const& hit_velocity, amrex::RandomEngine const& engine) const noexcept
{
    amrex::ParticleReal const speed_squared = Math::NormSquared(hit_velocity);

    amrex::ParticleReal const energy_eV = PhysConst::m_e * speed_squared / (2.0 * PhysConst::q_e);
    amrex::ParticleReal const a = energy_eV / 43.4592;
    amrex::ParticleReal const r = energy_eV / 30.0;
    amrex::ParticleReal const s = energy_eV / 127.8958;
    amrex::ParticleReal const p_absorb = 0.5 * std::exp(-a*a);
    amrex::ParticleReal const p_specular = 0.5 * std::exp(-r*r);
    amrex::ParticleReal const p_see2 = -std::expm1(-s*s);
    amrex::ParticleReal const draw = amrex::Random(engine);
    if (draw < p_absorb) {
        return WallBehavior::absorb;
    }
    if (draw < p_absorb + p_specular) {
        return WallBehavior::specular;
    }
    return draw < p_absorb + p_specular + p_see2 ?
        WallBehavior::secondary_electron_2 : WallBehavior::secondary_electron_1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<ParticleVector, max_wall_emissions> CeramicInteraction::SampleEmission (
    WallBehavior const event, ParticleVector const& hit_velocity,
    ParticleVector const& normal_to_domain, amrex::RandomEngine const& engine) const noexcept
{
    ParticleVector normal = normal_to_domain;
    Math::Normalize(normal_to_domain, normal);
    amrex::GpuArray<ParticleVector, max_wall_emissions> velocities{};
    int const count = WallEmissionCount(event);
    for (int j = 0; j < count; ++j) {
        DiffuseReemissionOperator{m_thermal_velocity}(
            normal, hit_velocity, velocities[j], engine);
    }
    return velocities;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ParticleVector CeramicInteraction::SampleReflection (
    ParticleVector const& hit_velocity, ParticleVector const& normal_to_domain) const noexcept
{
    ParticleVector normal = normal_to_domain;
    Math::Normalize(normal_to_domain, normal);
    return Math::Reflect(hit_velocity, normal);
}

} // namespace Insert
