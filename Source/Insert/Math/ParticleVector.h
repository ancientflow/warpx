#pragma once

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace Insert {

/** Particle coordinates, velocities and local directions must not pass through
 * amrex::XDim3, whose components use field precision (amrex::Real).
 */
struct ParticleVector
{
    amrex::ParticleReal x;
    amrex::ParticleReal y;
    amrex::ParticleReal z;
};

namespace Math {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::ParticleReal NormSquared (ParticleVector const& v) noexcept
{
    return v.x*v.x + v.y*v.y + v.z*v.z;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool Normalize (ParticleVector const& v, ParticleVector& result) noexcept
{
    auto const norm_squared = NormSquared(v);
    if (norm_squared <= amrex::ParticleReal(0.0)) {
        return false;
    }
    auto const inverse_norm = amrex::ParticleReal(1.0) / std::sqrt(norm_squared);
    result = {v.x * inverse_norm, v.y * inverse_norm, v.z * inverse_norm};
    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::ParticleReal Dot (ParticleVector const& a, ParticleVector const& b) noexcept
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ParticleVector Cross (ParticleVector const& a, ParticleVector const& b) noexcept
{
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

/** Same tangential-axis convention as VectorOps.h, in particle precision. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void BuildOrthonormalBasis (
    ParticleVector const& normal, ParticleVector& tangent1, ParticleVector& tangent2) noexcept
{
    if (std::abs(normal.z) < amrex::ParticleReal(0.9)) {
        auto const inverse_norm = amrex::ParticleReal(1.0) /
            std::sqrt(normal.x*normal.x + normal.y*normal.y);
        tangent1 = {-normal.y * inverse_norm, normal.x * inverse_norm, amrex::ParticleReal(0.0)};
    } else {
        auto const inverse_norm = amrex::ParticleReal(1.0) /
            std::sqrt(normal.y*normal.y + normal.z*normal.z);
        tangent1 = {amrex::ParticleReal(0.0), normal.z * inverse_norm, -normal.y * inverse_norm};
    }
    tangent2 = Cross(normal, tangent1);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ParticleVector ExpandInBasis (
    ParticleVector const& local, ParticleVector const& axis1,
    ParticleVector const& axis2, ParticleVector const& axis3) noexcept
{
    return {
        local.x*axis1.x + local.y*axis2.x + local.z*axis3.x,
        local.x*axis1.y + local.y*axis2.y + local.z*axis3.y,
        local.x*axis1.z + local.y*axis2.z + local.z*axis3.z};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ParticleVector Reflect (ParticleVector const& v, ParticleVector const& unit_normal) noexcept
{
    auto const twice_projection = amrex::ParticleReal(2.0) * Dot(v, unit_normal);
    return {
        v.x - twice_projection * unit_normal.x,
        v.y - twice_projection * unit_normal.y,
        v.z - twice_projection * unit_normal.z};
}

} // namespace Math
} // namespace Insert
