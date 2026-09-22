#ifndef WARPX_INSERT_MATH_HYPERBOLIC_H_
#define WARPX_INSERT_MATH_HYPERBOLIC_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

/** \brief Numerically stable evaluations of hyperbolic functions for the
 *         spectral boundary Schur solver.
 *
 *  Both functions switch to a large-argument asymptotic form beyond the
 *  threshold 40: there exp(-2 * x) < 1e-35 is far below machine precision,
 *  so the exponential corrections can be dropped (StableCoth) or evaluated
 *  only through small, non-overflowing exponentials (SinhDecayRatio),
 *  avoiding the overflow of exp(x) / sinh(x) for large x.
 */
namespace Insert::Math {

/** \brief Evaluate `coth(x)` without overflowing for large positive `x`.
 *
 *  For x beyond the threshold, coth(x) = (1 + e^{-2x}) / (1 - e^{-2x})
 *  is replaced by its asymptotic value 1.
 *
 * \param x Positive argument.
 * \return Stable `coth(x)`.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
StableCoth (amrex::Real x) {
    constexpr amrex::Real threshold = static_cast<amrex::Real>(40.0);
    if (x > threshold) {
        return static_cast<amrex::Real>(1.0);
    }
    amrex::Real const exp_neg_2x = std::exp(-static_cast<amrex::Real>(2.0) * x);
    return (static_cast<amrex::Real>(1.0) + exp_neg_2x) /
           (static_cast<amrex::Real>(1.0) - exp_neg_2x);
}

/** \brief Compute `sinh(k(Lz-z)) / sinh(kLz)` without large-argument
 *         overflow.
 *
 *  For k*Lz beyond the threshold, the ratio is evaluated in the
 *  asymptotic form
 *      e^{-kz} (1 - e^{-2k(Lz-z)}) / (1 - e^{-2kLz}),
 *  which only involves decaying exponentials.
 *
 * \param k Transverse modal wave number.
 * \param z Distance from zmin.
 * \param lz Domain length in z.
 * \return Harmonic extension decay factor.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
SinhDecayRatio (amrex::Real k, amrex::Real z, amrex::Real lz) {
    amrex::Real const a = k * lz;
    amrex::Real const b = k * (lz - z);
    constexpr amrex::Real threshold = static_cast<amrex::Real>(40.0);
    if (a <= threshold) {
        amrex::Real const denom = std::sinh(a);
        return std::sinh(b) / denom;
    }

    amrex::Real const numerator_correction =
        static_cast<amrex::Real>(1.0) -
        std::exp(-static_cast<amrex::Real>(2.0) * b);
    amrex::Real const denominator_correction =
        static_cast<amrex::Real>(1.0) -
        std::exp(-static_cast<amrex::Real>(2.0) * a);
    return std::exp(-k * z) * numerator_correction / denominator_correction;
}

} // namespace Insert::Math

#endif // WARPX_INSERT_MATH_HYPERBOLIC_H_
