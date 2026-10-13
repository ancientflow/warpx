#ifndef WARPX_INSERT_MATH_BILINEARCELL_H_
#define WARPX_INSERT_MATH_BILINEARCELL_H_

#include <AMReX_REAL.H>

/** \brief Shared coefficient algebra for the bilinear nodal model of one
 *         r-z cell of the averaged ionization source.
 *
 *  The cell-local field is expanded as
 *
 *      s(eta, xi) = a + b * xi + c * eta + d * eta * xi,
 *
 *  where eta (radial) and xi (axial) are cell-local coordinates in [0, 1]
 *  and the nodal values are
 *      s00 = s(0, 0),  s10 = s(0, 1)  (axial-plus node),
 *      s01 = s(1, 0)  (radial-plus node),  s11 = s(1, 1).
 *
 *  The writer (IonizationSourceTable::cellRate) integrates this model to
 *  tabulate the cell rates, and the reader
 *  (IonizationSourceSampler::samplePosition) inverts the same model to
 *  sample positions within a cell. Both sides must agree on the expansion,
 *  so the coefficient algebra lives here exactly once.
 *
 *  Host-only: the tables are built and sampled at setup/injection time.
 */
namespace Insert::Math {

/** \brief Coefficients of the bilinear expansion
 *         s(eta, xi) = a + b*xi + c*eta + d*eta*xi over one cell. */
struct BilinearCoefficients
{
    amrex::Real a = amrex::Real(0.0);
    amrex::Real b = amrex::Real(0.0);
    amrex::Real c = amrex::Real(0.0);
    amrex::Real d = amrex::Real(0.0);
};

/** \brief Build the bilinear coefficients from the four nodal values. */
inline BilinearCoefficients
MakeBilinearCoefficients (
    amrex::Real s00, amrex::Real s10,
    amrex::Real s01, amrex::Real s11) noexcept
{
    return BilinearCoefficients{
        s00, s10 - s00, s01 - s00, s11 - s10 - s01 + s00};
}

/** \brief Coefficients of the r-weighted radial marginal density.
 *
 *  Averaging the bilinear field over the axial coordinate gives the radial
 *  marginal alpha + beta * eta with
 *      alpha = a + b / 2,  beta = c + d / 2.
 *  Weighting by the radius r = r0 + dr * eta (the Jacobian of the
 *  cylindrical volume element) yields the quadratic density
 *      q0 + q1 * eta + q2 * eta^2
 *      = (alpha + beta * eta) * (r0 + dr * eta).
 */
struct RadialCdfCoefficients
{
    amrex::Real q0 = amrex::Real(0.0);
    amrex::Real q1 = amrex::Real(0.0);
    amrex::Real q2 = amrex::Real(0.0);
};

/** \brief Build the r-weighted radial marginal coefficients for the cell
 *         whose radial range is [r0, r0 + dr]. */
inline RadialCdfCoefficients
MakeRadialCdfCoefficients (
    BilinearCoefficients const& coef,
    amrex::Real r0, amrex::Real dr) noexcept
{
    amrex::Real const alpha = coef.a + amrex::Real(0.5) * coef.b;
    amrex::Real const beta = coef.c + amrex::Real(0.5) * coef.d;
    return RadialCdfCoefficients{
        alpha * r0, alpha * dr + beta * r0, beta * dr};
}

/** \brief Normalization of the radial CDF: integral of the quadratic
 *         density q0 + q1 * eta + q2 * eta^2 over eta in [0, 1], i.e.
 *         q0 + q1 / 2 + q2 / 3. (The CDF itself is cubic in eta.) */
inline amrex::Real
CubicCdfIntegral (amrex::Real q0, amrex::Real q1, amrex::Real q2) noexcept
{
    return q0 + amrex::Real(0.5) * q1 +
           (amrex::Real(1.0) / amrex::Real(3.0)) * q2;
}

} // namespace Insert::Math

#endif // WARPX_INSERT_MATH_BILINEARCELL_H_
