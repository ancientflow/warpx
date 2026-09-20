#ifndef WARPX_INSERT_BOUNDARYMATHUTILS_H_
#define WARPX_INSERT_BOUNDARYMATHUTILS_H_

#include "Insert/Boundary/WallInteractionOperators.h"
#include "Particles/Pusher/UpdatePosition.H"

#include <AMReX_Algorithm.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_Random.H>
#include <AMReX_Dim3.H>

/** \brief Pure math operators for particle interactions with an analytic
 *         boundary.
 *
 *  All functions are templates on the geometry operator (e.g.
 *  AnalyticBoundaryGeometry::DeviceView), which must provide:
 *    - SignedValue(x): signed boundary expression, > 0 on the
 *      computational-domain side, < 0 on the solid side;
 *    - Normal(x): normal vector at x (assumed on the boundary), pointing
 *      into the computational domain.
 *
 *  Behavior selection (absorb / specular / diffuse probabilities) and
 *  particle invalidation are policy and deliberately do NOT live here; they
 *  belong to the particle-processing driver layer.
 *
 *  Velocities are relativistic proper velocities u = gamma * v, matching
 *  the WarpX particle storage convention. Position updates go through
 *  UpdatePosition and therefore require the particle mass.
 */
namespace Insert::BoundaryMath {

/** \brief Fraction-of-timestep tolerance used when bisecting for the
 *         boundary intersection point. The corresponding physical accuracy
 *         is bisection_tolerance * |v| * dt. */
inline constexpr amrex::Real bisection_tolerance = amrex::Real(1.0e-6);

/** \brief Hit test: is the particle on the solid side of the boundary? */
template <typename Boundary>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
IsOutsideDomain (Boundary const& boundary, amrex::XDim3 const& x) noexcept
{
    // A negative sign marks the solid side of the boundary.
    return boundary.SignedValue(x) < amrex::ParticleReal(0.0);
}

/** \brief Bisect along the particle trajectory for the intersection point
 *         with the boundary.
 *
 *  Starting from the end-of-step position x_end (assumed on the solid
 *  side), the particle is back-traced with its current velocity u, and the
 *  root of SignedValue = 0 is located by bisection in the fraction of dt.
 *
 *  Assumption: the start-of-step position is on the domain side, i.e. the
 *  interval brackets a root. If this does not hold (e.g. the particle was
 *  already inside the boundary at the start of the step), no root exists;
 *  the fallback returns x_end unchanged with dt_fraction_hit = 0.
 *
 * \param boundary        geometry operator
 * \param x_end           end-of-step particle position
 * \param u               particle proper velocity (assumed constant over the step)
 * \param mass            particle mass (needed for the u -> v conversion)
 * \param dt              timestep
 * \param x_hit           output: intersection point on the boundary
 * \param dt_fraction_hit output: fraction of dt between x_hit and x_end
 * \param tol             bisection tolerance in units of the dt fraction
 */
template <typename Boundary>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
BisectBoundaryIntersection (
    Boundary const& boundary, amrex::XDim3 const& x_end,
    amrex::XDim3 const& u, amrex::ParticleReal const mass,
    amrex::Real const dt, amrex::XDim3& x_hit,
    amrex::Real& dt_fraction_hit,
    amrex::Real const tol = bisection_tolerance) noexcept
{
    // Back-trace the particle along its trajectory by dt_fraction * dt and
    // evaluate the boundary expression at the trial position.
    auto const signed_value_at = [&] (amrex::Real const dt_fraction) {
        amrex::ParticleReal xt = x_end.x;
        amrex::ParticleReal yt = x_end.y;
        amrex::ParticleReal zt = x_end.z;
        UpdatePosition(xt, yt, zt, u.x, u.y, u.z, -dt_fraction * dt, mass);
        return boundary.SignedValue({xt, yt, zt});
    };

    if (signed_value_at(amrex::Real(1.0)) <= amrex::ParticleReal(0.0)) {
        // The step-start position is not on the domain side: no root is
        // bracketed. Fall back to the current position.
        x_hit = x_end;
        dt_fraction_hit = amrex::Real(0.0);
        return;
    }

    // Locate the root of SignedValue = 0 in the dt fraction by bisection,
    // then recompute the intersection point at that fraction.
    dt_fraction_hit = amrex::bisect(
        amrex::Real(0.0), amrex::Real(1.0), signed_value_at, tol);
    amrex::ParticleReal xh = x_end.x;
    amrex::ParticleReal yh = x_end.y;
    amrex::ParticleReal zh = x_end.z;
    UpdatePosition(xh, yh, zh, u.x, u.y, u.z, -dt_fraction_hit * dt, mass);
    x_hit = amrex::XDim3{xh, yh, zh};
}

/** \brief Advance a position with a constant proper velocity u over dt.
 *         The position is updated in place. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
AdvancePosition (
    amrex::XDim3& x, amrex::XDim3 const& u,
    amrex::ParticleReal const mass, amrex::Real const dt) noexcept
{
    // UpdatePosition internally converts the proper velocity u = gamma * v
    // to the physical velocity using the particle mass. XDim3 stores
    // amrex::Real while UpdatePosition works on amrex::ParticleReal, so the
    // coordinates go through local copies.
    amrex::ParticleReal xt = x.x;
    amrex::ParticleReal yt = x.y;
    amrex::ParticleReal zt = x.z;
    UpdatePosition(xt, yt, zt, u.x, u.y, u.z, dt, mass);
    x = amrex::XDim3{xt, yt, zt};
}

/** \brief Specular reflection of the proper velocity about the boundary
 *         normal: u_out = u_in - 2 (u_in . n) n. The reflection preserves
 *         |u|, so gamma is unchanged. The engine is unused but kept in the
 *         signature for interface symmetry with DiffuseVelocity. */
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::XDim3
ReflectVelocity (
    amrex::XDim3 const& normal_to_domain, amrex::XDim3 const& u_in,
    amrex::RandomEngine const& engine) noexcept
{
    // Delegate to the shared wall operator; a degenerate (zero) normal
    // leaves the velocity unchanged.
    amrex::XDim3 u_out = u_in;
    SpecularReflectionOperator{}(normal_to_domain, u_in, u_out, engine);
    return u_out;
}

/** \brief Diffuse (fully accommodating thermal wall) re-emission velocity:
 *         half-Maxwellian flux distribution along the domain-pointing
 *         normal, full Maxwellian in the tangential plane.
 *
 * \param vth thermal (proper) velocity spread. The conversion from a wall
 *            temperature, vth = sqrt(kb * T / m), requires the species mass
 *            and is the responsibility of the driver layer.
 */
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::XDim3
DiffuseVelocity (
    amrex::XDim3 const& normal_to_domain, amrex::ParticleReal const vth,
    amrex::RandomEngine const& engine) noexcept
{
    // The diffuse re-emission velocity is drawn from the wall distribution
    // and does not depend on the incident velocity; pass a zero placeholder
    // for the unused u_in argument of the shared wall operator.
    amrex::XDim3 const u_unused{
        amrex::ParticleReal(0.0), amrex::ParticleReal(0.0),
        amrex::ParticleReal(0.0)};
    amrex::XDim3 u_out = u_unused;
    DiffuseReemissionOperator{vth}(normal_to_domain, u_unused, u_out, engine);
    return u_out;
}

} // namespace Insert::BoundaryMath

#endif // WARPX_INSERT_BOUNDARYMATHUTILS_H_
