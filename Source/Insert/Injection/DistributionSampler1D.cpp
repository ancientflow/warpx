#include "DistributionSampler1D.H"

#include "Utils/TextMsg.H"

#include <AMReX_Random.H>

#include <algorithm>
#include <cmath>

namespace Insert {

DistributionSampler1D::DistributionSampler1D (
    amrex::ParticleReal xmin, amrex::ParticleReal xmax, int num_bins,
    std::function<amrex::ParticleReal(amrex::ParticleReal)> distribution)
    : m_xmin(xmin), m_xmax(xmax)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        xmax > xmin, "DistributionSampler1D requires xmax > xmin.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        num_bins > 0, "DistributionSampler1D requires num_bins > 0.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        static_cast<bool>(distribution),
        "DistributionSampler1D requires a valid distribution function.");

    m_dx = (m_xmax - m_xmin) / static_cast<amrex::ParticleReal>(num_bins);
    m_cdf.assign(num_bins + 1, static_cast<amrex::ParticleReal>(0.0));

    // Use the midpoint rule to build a tabulated, unnormalized CDF.
    for (int i = 0; i < num_bins; ++i) {
        const amrex::ParticleReal x =
            m_xmin +
            (static_cast<amrex::ParticleReal>(i) + static_cast<amrex::ParticleReal>(0.5)) *
                m_dx;
        const amrex::ParticleReal density = distribution(x);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(density >= amrex::ParticleReal(0.0) &&
                                             std::isfinite(density),
                                         "DistributionSampler1D distribution "
                                         "must be finite and non-negative.");

        m_cdf[i + 1] = m_cdf[i] + density * m_dx;
    }

    m_integral = m_cdf.back();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_integral > amrex::ParticleReal(0.0),
        "DistributionSampler1D distribution integral must be positive.");

    // Normalize the CDF to [0, 1], while preserving the original integral.
    for (auto& cdf_value : m_cdf) {
        cdf_value /= m_integral;
    }
    m_cdf.back() = amrex::ParticleReal(1.0);
}

amrex::ParticleReal
DistributionSampler1D::sample (
    amrex::RandomEngine const& engine) const noexcept
{
    const amrex::ParticleReal selector = amrex::Random(engine);
    // Locate the tabulated CDF interval and linearly interpolate within it.
    auto const upper = std::upper_bound(m_cdf.begin(), m_cdf.end(), selector);

    const auto upper_index = static_cast<int>(upper - m_cdf.begin());
    const int bin = upper_index - 1;

    const amrex::ParticleReal cdf_lo = m_cdf[bin];
    const amrex::ParticleReal cdf_hi = m_cdf[upper_index];
    const amrex::ParticleReal x_lo =
        m_xmin + static_cast<amrex::ParticleReal>(bin) * m_dx;

    const amrex::ParticleReal fraction =
        (selector - cdf_lo) / (cdf_hi - cdf_lo);
    return x_lo + fraction * m_dx;
}

amrex::ParticleReal
DistributionSampler1D::integral () const noexcept
{
    return m_integral;
}

amrex::ParticleReal
DistributionSampler1D::xmin () const noexcept
{
    return m_xmin;
}

amrex::ParticleReal
DistributionSampler1D::xmax () const noexcept
{
    return m_xmax;
}

} // namespace Insert
