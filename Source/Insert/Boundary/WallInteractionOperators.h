#ifndef WARPX_INSERT_WALLINTERACTIONOPERATORS_H_
#define WARPX_INSERT_WALLINTERACTIONOPERATORS_H_

#include "Initialization/SampleGaussianFluxDistribution.H"

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_Random.H>
#include <AMReX_Dim3.H>

#include <cmath>

namespace Insert {

struct SpecularReflectionOperator {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator()(
        amrex::XDim3 const& normal_to_domain,
        amrex::XDim3 const& u_in, amrex::XDim3& u_out,
        amrex::RandomEngine const& engine) const noexcept {
        // A degenerate (zero) normal carries no geometric information:
        // leave the velocity unchanged.
        amrex::ParticleReal const normal_norm_sq =
            normal_to_domain.x * normal_to_domain.x
            + normal_to_domain.y * normal_to_domain.y
            + normal_to_domain.z * normal_to_domain.z;
        if (normal_norm_sq <= static_cast<amrex::ParticleReal>(0.0)) {
            u_out = u_in;
            amrex::ignore_unused(engine);
            return;
        }

        // Normalize the normal, then mirror the velocity about the boundary
        // plane: u_out = u_in - 2 (u_in . n) n.
        using std::sqrt;
        amrex::ParticleReal const inv_normal_norm =
            static_cast<amrex::ParticleReal>(1.0) / sqrt(normal_norm_sq);
        amrex::XDim3 const normal{
            normal_to_domain.x * inv_normal_norm,
            normal_to_domain.y * inv_normal_norm,
            normal_to_domain.z * inv_normal_norm};
        amrex::ParticleReal const u_dot_normal =
            u_in.x * normal.x + u_in.y * normal.y + u_in.z * normal.z;

        u_out = {
            u_in.x - static_cast<amrex::ParticleReal>(2.0) * u_dot_normal * normal.x,
            u_in.y - static_cast<amrex::ParticleReal>(2.0) * u_dot_normal * normal.y,
            u_in.z - static_cast<amrex::ParticleReal>(2.0) * u_dot_normal * normal.z};

        amrex::ignore_unused(engine);
    }
};

struct DiffuseReemissionOperator {
    amrex::ParticleReal m_vth = 0.0;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator()(
        amrex::XDim3 const& normal_to_domain,
        amrex::XDim3 const& u_in, amrex::XDim3& u_out,
        amrex::RandomEngine const& engine) const noexcept {
        // A degenerate (zero) normal carries no geometric information:
        // leave the velocity unchanged.
        amrex::ParticleReal const normal_norm_sq =
            normal_to_domain.x * normal_to_domain.x
            + normal_to_domain.y * normal_to_domain.y
            + normal_to_domain.z * normal_to_domain.z;
        if (normal_norm_sq <= static_cast<amrex::ParticleReal>(0.0)) {
            u_out = u_in;
            return;
        }

        using std::abs;
        using std::sqrt;
        amrex::ParticleReal const inv_normal_norm =
            static_cast<amrex::ParticleReal>(1.0) / sqrt(normal_norm_sq);
        amrex::XDim3 const normal{
            normal_to_domain.x * inv_normal_norm,
            normal_to_domain.y * inv_normal_norm,
            normal_to_domain.z * inv_normal_norm};

        // Build an orthonormal tangential basis (tangent_one, tangent_two)
        // around the normal. The branch avoids the near-degenerate case
        // where the normal is almost aligned with the z axis.
        amrex::XDim3 tangent_one;
        if (abs(normal.z) < static_cast<amrex::ParticleReal>(0.9)) {
            amrex::ParticleReal const inv_tangent_norm =
                static_cast<amrex::ParticleReal>(1.0)
                / sqrt(normal.x * normal.x + normal.y * normal.y);
            tangent_one = {
                -normal.y * inv_tangent_norm,
                normal.x * inv_tangent_norm,
                static_cast<amrex::ParticleReal>(0.0)};
        } else {
            amrex::ParticleReal const inv_tangent_norm =
                static_cast<amrex::ParticleReal>(1.0)
                / sqrt(normal.y * normal.y + normal.z * normal.z);
            tangent_one = {
                static_cast<amrex::ParticleReal>(0.0),
                normal.z * inv_tangent_norm,
                -normal.y * inv_tangent_norm};
        }
        amrex::XDim3 const tangent_two{
            normal.y * tangent_one.z - normal.z * tangent_one.y,
            normal.z * tangent_one.x - normal.x * tangent_one.z,
            normal.x * tangent_one.y - normal.y * tangent_one.x};

        // Sample the wall Maxwellian: full Maxwellian in both tangential
        // directions, half-Maxwellian flux distribution along the
        // domain-pointing normal.
        amrex::ParticleReal const u_tangent_one =
            amrex::RandomNormal(static_cast<amrex::ParticleReal>(0.0), m_vth, engine);
        amrex::ParticleReal const u_tangent_two =
            amrex::RandomNormal(static_cast<amrex::ParticleReal>(0.0), m_vth, engine);
        amrex::ParticleReal const u_normal =
            generateGaussianFluxDist(
                static_cast<amrex::ParticleReal>(0.0), m_vth, engine);

        // Rotate the sampled velocity from the local (t1, t2, n) frame back
        // to Cartesian components.
        u_out = {
            u_tangent_one * tangent_one.x
                + u_tangent_two * tangent_two.x + u_normal * normal.x,
            u_tangent_one * tangent_one.y
                + u_tangent_two * tangent_two.y + u_normal * normal.y,
            u_tangent_one * tangent_one.z
                + u_tangent_two * tangent_two.z + u_normal * normal.z};

        amrex::ignore_unused(u_in);
    }
};

} // namespace Insert

#endif // WARPX_INSERT_WALLINTERACTIONOPERATORS_H_
