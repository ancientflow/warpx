#ifndef WARPX_INSERT_MATH_SAMPLING_H_
#define WARPX_INSERT_MATH_SAMPLING_H_

#include <AMReX_Dim3.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_Random.H>

#include <cmath>

/** \brief Shared sampling primitives for the Insert injection and geometry
 *         code paths. All random draws go through the supplied
 *         amrex::RandomEngine so call sites keep control over their random
 *         stream. */
namespace Insert::Math {

/** \brief Sample uniformly on the interval [lo, hi):
 *         lo + (hi - lo) * U with U ~ U[0, 1). */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T
SampleUniform (T lo, T hi, amrex::RandomEngine const& engine) noexcept {
    return static_cast<T>(lo + (hi - lo) * amrex::Random(engine));
}

/** \brief Sample a radius on [rmin, rmax] with probability density
 *         proportional to r, i.e. uniformly per unit area of the
 *         corresponding annulus: sqrt(rmin^2 + (rmax^2 - rmin^2) * U). */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T
SampleAreaUniformRadius (T rmin, T rmax,
                         amrex::RandomEngine const& engine) noexcept {
    using std::sqrt;
    return static_cast<T>(sqrt(rmin * rmin + (rmax * rmax - rmin * rmin) *
                                                 amrex::Random(engine)));
}

/** \brief Convert polar coordinates (r, theta) to a Cartesian vector with
 *         z = 0: (r cos theta, r sin theta, 0). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::XDim3
PolarToCartesian (amrex::Real r, amrex::Real theta) noexcept {
    using std::cos;
    using std::sin;
    return amrex::XDim3{r * cos(theta), r * sin(theta), amrex::Real(0.0)};
}

/** \brief Wrap x into the periodic interval [0, period): fmod with a
 *         correction for negative results. */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T
WrapToPeriod (T x, T period) noexcept {
    T wrapped = std::fmod(x, period);
    if (wrapped < T(0.0)) {
        wrapped += period;
    }
    return wrapped;
}

/** \brief Signed distance between two angles on the 2 pi-periodic circle,
 *         wrapped to [-pi, pi). */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T
WrappedAngleDistance (T x, T center) noexcept {
    const T pi = amrex::Math::pi<T>();
    return WrapToPeriod(x - center + pi, T(2.0) * pi) - pi;
}

/** \brief Unnormalized Gaussian density on the 2 pi-periodic circle:
 *         exp(-d^2 / (2 sigma^2)) with the wrapped distance
 *         d = WrappedAngleDistance(x, center). */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T
WrappedGaussianDensity (T x, T center, T sigma) noexcept {
    const T normalized = WrappedAngleDistance(x, center) / sigma;
    return std::exp(T(-0.5) * normalized * normalized);
}

} // namespace Insert::Math

#endif // WARPX_INSERT_MATH_SAMPLING_H_
