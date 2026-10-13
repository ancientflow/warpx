#pragma once

#include <AMReX_REAL.H>
#include <AMReX_RandomEngine.H>
#include <AMReX_Vector.H>

#include <functional>

namespace Insert {

/**
 * Numerical inverse-CDF sampler for a one-dimensional non-negative PDF.
 *
 * Use this class for distributions whose cumulative distribution function, or
 * inverse cumulative distribution function, is not available in a convenient
 * closed form. Simple distributions with elementary inverse CDFs should sample
 * directly instead of building this table.
 *
 * The constructor tabulates the cumulative distribution function from the
 * supplied PDF. The sample routine then uses inverse-CDF sampling. The original
 * integral is kept so callers can decide whether the distribution should only
 * redistribute particles or also scale the total particle count.
 */
class NumericalInverseCDFSampler1D
{
public:
    NumericalInverseCDFSampler1D (
        amrex::ParticleReal xmin, amrex::ParticleReal xmax, int num_bins,
        std::function<amrex::ParticleReal(amrex::ParticleReal)> pdf);

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

using DistributionSampler1D = NumericalInverseCDFSampler1D;

} // namespace Insert
