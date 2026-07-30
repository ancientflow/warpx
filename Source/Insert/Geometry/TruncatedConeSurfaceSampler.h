#pragma once

#include <AMReX_Dim3.H>
#include <AMReX_REAL.H>
#include <AMReX_RandomEngine.H>

namespace Insert {

/**
 * Pure geometric sampler for an axisymmetric truncated-cone surface generated
 * from a one-dimensional radial interval.
 *
 * The generating line is
 *
 *     z(r) = z_reference + slope * (r - r_reference).
 *
 * The radial coordinate is sampled with density proportional to r so that
 * equal-weight samples are uniform per unit cone surface area. The azimuthal
 * coordinate only rotates the sampled generating-line point about the z axis.
 */
class TruncatedConeSurfaceSampler
{
public:
    TruncatedConeSurfaceSampler (
        amrex::ParticleReal slope,
        amrex::ParticleReal r_min,
        amrex::ParticleReal r_max,
        amrex::ParticleReal r_reference,
        amrex::ParticleReal z_reference,
        amrex::ParticleReal theta_min,
        amrex::ParticleReal theta_max);

    [[nodiscard]] amrex::ParticleReal
    sampleRadius (amrex::RandomEngine const& engine) const;

    [[nodiscard]] amrex::ParticleReal
    sampleTheta (amrex::RandomEngine const& engine) const;

    [[nodiscard]] amrex::XDim3
    coordinates (
        amrex::ParticleReal radius,
        amrex::ParticleReal theta) const noexcept;

    [[nodiscard]] amrex::XDim3
    sampleCoordinates (amrex::RandomEngine const& engine) const;

private:
    amrex::ParticleReal m_slope;
    amrex::ParticleReal m_r2_min;
    amrex::ParticleReal m_r2_max;
    amrex::ParticleReal m_r_reference;
    amrex::ParticleReal m_z_reference;
    amrex::ParticleReal m_theta_min;
    amrex::ParticleReal m_theta_max;
};

} // namespace Insert
