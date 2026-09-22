#include "TruncatedConeSurfaceSampler.h"

#include "Insert/Math/Sampling.h"
#include "Utils/TextMsg.H"

#include <cmath>

namespace Insert {

TruncatedConeSurfaceSampler::TruncatedConeSurfaceSampler (
    amrex::ParticleReal slope,
    amrex::ParticleReal r_min,
    amrex::ParticleReal r_max,
    amrex::ParticleReal r_reference,
    amrex::ParticleReal z_reference,
    amrex::ParticleReal theta_min,
    amrex::ParticleReal theta_max)
    : m_slope(slope),
      m_r_min(r_min),
      m_r_max(r_max),
      m_r_reference(r_reference),
      m_z_reference(z_reference),
      m_theta_min(theta_min),
      m_theta_max(theta_max)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(m_slope),
        "TruncatedConeSurfaceSampler requires a finite slope.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        r_min >= amrex::ParticleReal(0.0) && r_max > r_min,
        "TruncatedConeSurfaceSampler requires 0 <= r_min < r_max.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(m_r_reference) && std::isfinite(m_z_reference),
        "TruncatedConeSurfaceSampler requires a finite reference point.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(m_theta_min) && std::isfinite(m_theta_max) &&
            m_theta_max > m_theta_min,
        "TruncatedConeSurfaceSampler requires theta_max > theta_min.");
}

amrex::XDim3
TruncatedConeSurfaceSampler::coordinates (
    amrex::ParticleReal radius,
    amrex::ParticleReal theta) const noexcept
{
    const auto planar = Math::PolarToCartesian(radius, theta);
    return amrex::XDim3{
        planar.x,
        planar.y,
        m_z_reference + m_slope * (radius - m_r_reference)};
}

amrex::ParticleReal
TruncatedConeSurfaceSampler::sampleRadius (
    amrex::RandomEngine const& engine) const
{
    return Math::SampleAreaUniformRadius(m_r_min, m_r_max, engine);
}

amrex::ParticleReal
TruncatedConeSurfaceSampler::sampleTheta (
    amrex::RandomEngine const& engine) const
{
    return Math::SampleUniform(m_theta_min, m_theta_max, engine);
}

amrex::XDim3
TruncatedConeSurfaceSampler::sampleCoordinates (
    amrex::RandomEngine const& engine) const
{
    const auto radius = sampleRadius(engine);
    const auto theta = sampleTheta(engine);
    return coordinates(radius, theta);
}

} // namespace Insert
