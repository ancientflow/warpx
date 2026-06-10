#pragma once

#include <AMReX_REAL.H>
#include <AMReX_RandomEngine.H>
#include <AMReX_Vector.H>

#include <functional>

namespace Insert {

/**
 * One-dimensional sampler for a non-negative distribution.
 *
 * The constructor tabulates the cumulative distribution function from the
 * supplied distribution function. The sample routine then uses inverse-CDF
 * sampling. The original integral is kept so callers can decide whether the
 * distribution should only redistribute particles or also scale the total
 * particle count.
 */
class DistributionSampler1D
{
public:
    DistributionSampler1D (
        amrex::ParticleReal xmin, amrex::ParticleReal xmax, int num_bins,
        std::function<amrex::ParticleReal(amrex::ParticleReal)> distribution);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const noexcept;

    [[nodiscard]] amrex::ParticleReal integral () const noexcept;
    [[nodiscard]] amrex::ParticleReal xmin () const noexcept;
    [[nodiscard]] amrex::ParticleReal xmax () const noexcept;

private:
    amrex::ParticleReal m_xmin = 0.0;
    amrex::ParticleReal m_xmax = 0.0;
    amrex::ParticleReal m_dx = 0.0;
    amrex::ParticleReal m_integral = 0.0;
    amrex::Vector<amrex::ParticleReal> m_cdf;
};

} // namespace Insert
