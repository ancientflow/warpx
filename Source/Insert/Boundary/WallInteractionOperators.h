#ifndef WARPX_INSERT_WALLINTERACTIONOPERATORS_H_
#define WARPX_INSERT_WALLINTERACTIONOPERATORS_H_

#include "Initialization/SampleGaussianFluxDistribution.H"
#include "Insert/Math/VectorOps.h"

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_Random.H>
#include <AMReX_Dim3.H>

namespace Insert {

struct SpecularReflectionOperator {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator()(
        amrex::XDim3 const& normal_to_domain,
        amrex::XDim3 const& u_in, amrex::XDim3& u_out,
        amrex::RandomEngine const& engine) const noexcept {
        // A degenerate (zero) normal carries no geometric information:
        // leave the velocity unchanged.
        amrex::XDim3 normal;
        if (!Math::Normalize<amrex::ParticleReal>(normal_to_domain, normal)) {
            u_out = u_in;
            amrex::ignore_unused(engine);
            return;
        }

        // Mirror the velocity about the boundary plane:
        // u_out = u_in - 2 (u_in . n) n.
        u_out = Math::Reflect<amrex::ParticleReal>(u_in, normal);

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
        amrex::XDim3 normal;
        if (!Math::Normalize<amrex::ParticleReal>(normal_to_domain, normal)) {
            u_out = u_in;
            return;
        }

        // Build an orthonormal tangential basis (tangent_one, tangent_two)
        // around the normal.
        amrex::XDim3 tangent_one;
        amrex::XDim3 tangent_two;
        Math::BuildOrthonormalBasis<amrex::ParticleReal>(
            normal, tangent_one, tangent_two);

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
        u_out = Math::ExpandInBasis(
            amrex::XDim3{
                static_cast<amrex::Real>(u_tangent_one),
                static_cast<amrex::Real>(u_tangent_two),
                static_cast<amrex::Real>(u_normal)},
            tangent_one, tangent_two, normal);

        amrex::ignore_unused(u_in);
    }
};

} // namespace Insert

#endif // WARPX_INSERT_WALLINTERACTIONOPERATORS_H_
