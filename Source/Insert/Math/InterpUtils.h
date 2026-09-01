#ifndef WARPX_INSERT_MATH_INTERPUTILS_H_
#define WARPX_INSERT_MATH_INTERPUTILS_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

/** \brief Linear hat-function (first-order nodal) weight splits shared by
 *         the Insert charge deposition and interpolation code paths.
 *
 *  Two boundary strategies exist and are intentionally kept distinct:
 *    - LinearHatWeightsDiscardOutside: positions outside the node range are
 *      reported via the false return value so the caller can drop them
 *      (used by the zmin wall-charge deposition);
 *    - LinearHatWeightsClamped: the left node index is clamped into the
 *      valid cell range so boundary positions deposit into the boundary
 *      cell (used by the ECDI charge filter).
 *
 *  Deposition (scatter) and interpolation (gather) that must conserve
 *  charge have to use the SAME weight function; the ECDI charge filter
 *  relies on this by calling LinearHatWeightsClamped from both sides.
 */
namespace Insert::Math {

/** \brief Split a node-coordinate position u between the two bracketing
 *         nodes of a grid with n_nodes nodes (u in units of cells, i.e.
 *         u = (x - x0) / dx, node j at integer u = j).
 *
 *  Returns false when u is outside [0, n_nodes - 1]; the caller is
 *  expected to discard the position. On the right edge (u == n_nodes - 1)
 *  the split degenerates to the last cell with w_right = 1. The left
 *  weight is 1 - w_right.
 *
 * \param u        position in cell units
 * \param n_nodes  number of nodes (>= 2)
 * \param i_left   output: left node index, in [0, n_nodes - 2]
 * \param w_right  output: weight of node i_left + 1, in [0, 1]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
LinearHatWeightsDiscardOutside (
    amrex::Real const u, int const n_nodes,
    int& i_left, amrex::Real& w_right) noexcept
{
    if (u < amrex::Real(0.0) || u > static_cast<amrex::Real>(n_nodes - 1)) {
        return false;
    }
    i_left = static_cast<int>(amrex::Math::floor(u));
    w_right = u - static_cast<amrex::Real>(i_left);
    if (i_left >= n_nodes - 1) {
        i_left = n_nodes - 2;
        w_right = amrex::Real(1.0);
    }
    return true;
}

/** \brief Split a position x between the two bracketing nodes of a uniform
 *         grid (nodes at x = i * dx, i in [0, n_cells]), clamping the cell
 *         index into [0, n_cells - 1] instead of discarding out-of-range
 *         positions.
 *
 *  The weights are computed as (x - x_left) / denom and
 *  (x_right - x) / denom, matching the arithmetic historically used by the
 *  ECDI charge filter. Returns false for a degenerate spacing
 *  (denom < 1e-15); i_left is still set so the caller can fall back to a
 *  single-node deposit.
 *
 * \param x        position (origin at node 0)
 * \param dx       node spacing
 * \param n_cells  number of cells (= number of nodes - 1, >= 1)
 * \param i_left   output: left node index, in [0, n_cells - 1]
 * \param w_left   output: weight of node i_left
 * \param w_right  output: weight of node i_left + 1
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
LinearHatWeightsClamped (
    amrex::Real const x, amrex::Real const dx, int const n_cells,
    int& i_left, amrex::Real& w_left, amrex::Real& w_right) noexcept
{
    i_left = static_cast<int>(amrex::Math::floor(x / dx));
    i_left = std::max(0, std::min(i_left, n_cells - 1));
    int const i_right = i_left + 1;

    amrex::Real const x_left = static_cast<amrex::Real>(i_left) * dx;
    amrex::Real const x_right = static_cast<amrex::Real>(i_right) * dx;
    amrex::Real const denom = x_right - x_left;
    if (denom < amrex::Real(1.0e-15)) {
        return false;
    }

    w_right = (x - x_left) / denom;
    w_left = (x_right - x) / denom;
    return true;
}

} // namespace Insert::Math

#endif // WARPX_INSERT_MATH_INTERPUTILS_H_
