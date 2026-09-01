#ifndef WARPX_INSERT_MATH_VECTOROPS_H_
#define WARPX_INSERT_MATH_VECTOROPS_H_

#include "Utils/TextMsg.H"

#include <AMReX_Dim3.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

/** \brief Small vector operations on amrex::XDim3 shared by the Insert
 *         boundary-interaction and injection code paths.
 *
 *  All functions are templates on the arithmetic type T (amrex::Real by
 *  default). Component products and sums are evaluated in amrex::Real and
 *  narrowed to T on return, matching the arithmetic precision historically
 *  used at the call sites.
 */
namespace Insert::Math {

/** \brief Squared Euclidean norm of v. */
template <typename T = amrex::Real>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T
NormSquared (amrex::XDim3 const& v) noexcept
{
    return static_cast<T>(v.x * v.x + v.y * v.y + v.z * v.z);
}

/** \brief Euclidean norm of v. */
template <typename T = amrex::Real>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T
Norm (amrex::XDim3 const& v) noexcept
{
    using std::sqrt;
    return sqrt(NormSquared<T>(v));
}

/** \brief Normalize v into result.
 *
 *  Zero-safe variant: returns false and leaves result untouched when v is a
 *  zero vector. The caller decides how to handle the degenerate case.
 */
template <typename T = amrex::Real>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
Normalize (amrex::XDim3 const& v, amrex::XDim3& result) noexcept
{
    T const norm_sq = NormSquared<T>(v);
    if (norm_sq <= T(0.0)) {
        return false;
    }
    using std::sqrt;
    T const inv_norm = T(1.0) / sqrt(norm_sq);
    result = amrex::XDim3{
        static_cast<amrex::Real>(v.x * inv_norm),
        static_cast<amrex::Real>(v.y * inv_norm),
        static_cast<amrex::Real>(v.z * inv_norm)};
    return true;
}

/** \brief Normalize v; aborts on zero or non-finite vectors. */
template <typename T = amrex::Real>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::XDim3
NormalizeChecked (amrex::XDim3 const& v)
{
    using std::sqrt;
    T const norm = sqrt(NormSquared<T>(v));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(norm) && norm > T(0.0),
        "Cannot normalize a zero or non-finite vector.");
    return amrex::XDim3{
        static_cast<amrex::Real>(v.x / norm),
        static_cast<amrex::Real>(v.y / norm),
        static_cast<amrex::Real>(v.z / norm)};
}

/** \brief Dot product a . b. */
template <typename T = amrex::Real>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T
Dot (amrex::XDim3 const& a, amrex::XDim3 const& b) noexcept
{
    return static_cast<T>(a.x * b.x + a.y * b.y + a.z * b.z);
}

/** \brief Cross product a x b. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::XDim3
Cross (amrex::XDim3 const& a, amrex::XDim3 const& b) noexcept
{
    return amrex::XDim3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x};
}

/** \brief Build an orthonormal tangential basis (tangent1, tangent2) around
 *         a unit normal.
 *
 *  The branch avoids the near-degenerate case where the normal is almost
 *  aligned with the z axis. Note the sign convention: tangent1 is the
 *  normalized (0,0,1) x normal in the |n.z| < 0.9 branch, and the negative
 *  of the normalized (1,0,0) x normal otherwise. This matches the
 *  historical wall diffuse re-emission operator; callers that draw
 *  user-specified (non-symmetric) tangential components must not assume
 *  the opposite sign convention.
 *
 *  The behavior for a zero normal is undefined; normalize (and validate)
 *  the input first.
 */
template <typename T = amrex::Real>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
BuildOrthonormalBasis (
    amrex::XDim3 const& unit_normal,
    amrex::XDim3& tangent1, amrex::XDim3& tangent2) noexcept
{
    using std::abs;
    using std::sqrt;
    if (abs(unit_normal.z) < static_cast<T>(0.9)) {
        T const inv_tangent_norm = static_cast<T>(
            amrex::Real(1.0) / sqrt(unit_normal.x * unit_normal.x +
                                    unit_normal.y * unit_normal.y));
        tangent1 = amrex::XDim3{
            static_cast<amrex::Real>(-unit_normal.y * inv_tangent_norm),
            static_cast<amrex::Real>(unit_normal.x * inv_tangent_norm),
            amrex::Real(0.0)};
    } else {
        T const inv_tangent_norm = static_cast<T>(
            amrex::Real(1.0) / sqrt(unit_normal.y * unit_normal.y +
                                    unit_normal.z * unit_normal.z));
        tangent1 = amrex::XDim3{
            amrex::Real(0.0),
            static_cast<amrex::Real>(unit_normal.z * inv_tangent_norm),
            static_cast<amrex::Real>(-unit_normal.y * inv_tangent_norm)};
    }
    tangent2 = Cross(unit_normal, tangent1);
}

/** \brief Expand a vector given by local components in a (typically
 *         orthonormal) basis: local.x*axis1 + local.y*axis2 + local.z*axis3.
 *
 *  When (axis1, axis2, axis3) is an orthonormal basis this is the
 *  local-to-Cartesian rotation step used after sampling a velocity in a
 *  local frame (e.g. wall diffuse re-emission in the (t1, t2, n) frame, or
 *  injection velocities given in a local_normal / rotating_axis frame).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::XDim3
ExpandInBasis (
    amrex::XDim3 const& local,
    amrex::XDim3 const& axis1, amrex::XDim3 const& axis2,
    amrex::XDim3 const& axis3) noexcept
{
    return amrex::XDim3{
        local.x * axis1.x + local.y * axis2.x + local.z * axis3.x,
        local.x * axis1.y + local.y * axis2.y + local.z * axis3.y,
        local.x * axis1.z + local.y * axis2.z + local.z * axis3.z};
}

/** \brief Rotate v by the rotation that maps axis_from onto axis_to.
 *
 *  The caller provides only the vector and the two axes; the full frames
 *  are completed internally: an orthonormal tangential basis is built
 *  around each axis with BuildOrthonormalBasis, v is decomposed in the
 *  from-frame and expanded in the to-frame. The result is a proper
 *  rotation with R * axis_from = axis_to.
 *
 *  Note: this is NOT the minimal (shortest-arc) rotation between the two
 *  axes. The twist about the axis pair follows the deterministic
 *  BuildOrthonormalBasis convention, which has two consequences:
 *    - there is no singularity, not even for antiparallel axes
 *      (axis_from = -axis_to), unlike shortest-arc formulas;
 *    - components tangential to the axes acquire a convention-dependent
 *      phase. This is immaterial when the tangential components are drawn
 *      from an azimuthally symmetric distribution (e.g. wall diffuse
 *      re-emission sampled in a canonical frame, then rotated onto the
 *      wall normal); callers with user-specified tangential components
 *      should build the frames explicitly and use ExpandInBasis instead.
 *
 *  Both axes must be non-zero and finite (asserted by NormalizeChecked).
 */
template <typename T = amrex::Real>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::XDim3
RotateBetweenAxes (
    amrex::XDim3 const& v,
    amrex::XDim3 const& axis_from, amrex::XDim3 const& axis_to)
{
    amrex::XDim3 const from = NormalizeChecked<T>(axis_from);
    amrex::XDim3 const to = NormalizeChecked<T>(axis_to);

    // (tangent1, tangent2, axis) is right-handed because
    // BuildOrthonormalBasis sets tangent2 = axis x tangent1.
    amrex::XDim3 from_t1;
    amrex::XDim3 from_t2;
    amrex::XDim3 to_t1;
    amrex::XDim3 to_t2;
    BuildOrthonormalBasis<T>(from, from_t1, from_t2);
    BuildOrthonormalBasis<T>(to, to_t1, to_t2);

    amrex::XDim3 const local{
        static_cast<amrex::Real>(Dot<T>(v, from_t1)),
        static_cast<amrex::Real>(Dot<T>(v, from_t2)),
        static_cast<amrex::Real>(Dot<T>(v, from))};
    return ExpandInBasis(local, to_t1, to_t2, to);
}

/** \brief Reflect v about the plane perpendicular to unit_normal:
 *         v_out = v - 2 (v . n) n (Householder reflection).
 *
 *  The tangential component is preserved and the normal component is
 *  flipped, so |v| is unchanged. unit_normal must be a unit vector;
 *  normalize (and validate) the input first, e.g. with Normalize.
 */
template <typename T = amrex::Real>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::XDim3
Reflect (amrex::XDim3 const& v, amrex::XDim3 const& unit_normal) noexcept
{
    T const v_dot_normal = Dot<T>(v, unit_normal);
    return amrex::XDim3{
        static_cast<amrex::Real>(
            v.x - static_cast<T>(2.0) * v_dot_normal * unit_normal.x),
        static_cast<amrex::Real>(
            v.y - static_cast<T>(2.0) * v_dot_normal * unit_normal.y),
        static_cast<amrex::Real>(
            v.z - static_cast<T>(2.0) * v_dot_normal * unit_normal.z)};
}

/** \brief Rotate a vector around the z axis by theta (right-handed). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::XDim3
RotateAroundZ (amrex::XDim3 const& value, amrex::Real theta) noexcept
{
    using std::cos;
    using std::sin;
    const auto cos_theta = cos(theta);
    const auto sin_theta = sin(theta);
    return amrex::XDim3{
        value.x * cos_theta - value.y * sin_theta,
        value.x * sin_theta + value.y * cos_theta,
        value.z};
}

} // namespace Insert::Math

#endif // WARPX_INSERT_MATH_VECTOROPS_H_
