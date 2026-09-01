#ifndef WARPX_INSERT_MATH_CDFUTILS_H_
#define WARPX_INSERT_MATH_CDFUTILS_H_

#include "Utils/TextMsg.H"

#include <AMReX_REAL.H>

#include <algorithm>
#include <cstddef>
#include <vector>

/** \brief Helpers to build tabulated cumulative distribution functions and
 *         interpolate tabulated functions for the Insert injection code
 *         paths. Host-only: the tabulated CDFs are built once at setup
 *         time. */
namespace Insert::Math {

/** \brief Build a normalized cumulative distribution from non-negative
 *         discrete weights.
 *
 *  On return, cdf[i] = sum_{j <= i} weights[j] / total_weight, so the CDF
 *  starts at cdf[0] = weights[0] / total_weight and ends exactly at
 *  cdf.back() = 1. The total weight (the unnormalized integral) is
 *  returned so callers can keep track of the distribution's overall
 *  normalization.
 */
template <typename T>
T
BuildNormalizedCumulativeWeights (std::vector<T> const& weights,
                                  std::vector<T>& cdf) {
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !weights.empty(),
        "BuildNormalizedCumulativeWeights requires at least one weight.");
    cdf.assign(weights.size(), T(0.0));

    T cumulative = T(0.0);
    for (std::size_t i = 0; i < weights.size(); ++i) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            weights[i] >= T(0.0),
            "BuildNormalizedCumulativeWeights weights must be non-negative.");
        cumulative += weights[i];
        cdf[i] = cumulative;
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        cumulative > T(0.0),
        "BuildNormalizedCumulativeWeights requires a positive total weight.");
    for (auto& value : cdf) {
        value /= cumulative;
    }
    cdf.back() = T(1.0);
    return cumulative;
}

/** \brief Linear interpolation in a tabulated function with endpoint
 *         clamping.
 *
 *  values must be strictly increasing and have the same size as y. Queries
 *  at or below values.front() return y.front(); queries at or above
 *  values.back() return y.back(). In between, the enclosing interval is
 *  located with std::upper_bound and y is linearly interpolated.
 */
template <typename T>
T
InterpolateTableClamped (std::vector<T> const& values,
                         std::vector<T> const& y, T x_query) noexcept {
    if (x_query <= values.front()) {
        return y.front();
    }
    if (x_query >= values.back()) {
        return y.back();
    }

    auto const upper = std::upper_bound(values.begin(), values.end(), x_query);
    const auto hi = static_cast<std::size_t>(upper - values.begin());
    const auto lo = hi - 1;
    const T fraction = (x_query - values[lo]) / (values[hi] - values[lo]);
    return y[lo] + fraction * (y[hi] - y[lo]);
}

} // namespace Insert::Math

#endif // WARPX_INSERT_MATH_CDFUTILS_H_
