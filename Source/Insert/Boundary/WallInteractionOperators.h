#ifndef WARPX_INSERT_WALLINTERACTIONOPERATORS_H_
#define WARPX_INSERT_WALLINTERACTIONOPERATORS_H_

#include "Insert/Math/ParticleVector.h"

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_Random.H>
#include <AMReX_Utility.H>

#include <cmath>

namespace Insert {

struct SpecularReflectionOperator {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator()(
        ParticleVector const& normal_to_domain,
        ParticleVector const& u_in, ParticleVector& u_out,
        amrex::RandomEngine const& engine) const noexcept {
        // A degenerate (zero) normal carries no geometric information:
        // leave the velocity unchanged.
        ParticleVector normal;
        if (!Math::Normalize(normal_to_domain, normal)) {
            u_out = u_in;
            amrex::ignore_unused(engine);
            return;
        }

        // Mirror the velocity about the boundary plane:
        // u_out = u_in - 2 (u_in . n) n.
        u_out = Math::Reflect(u_in, normal);

        amrex::ignore_unused(engine);
    }
};

struct DiffuseReemissionOperator {
    amrex::ParticleReal m_vth = 0.0;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator()(
        ParticleVector const& normal_to_domain,
        ParticleVector const& u_in, ParticleVector& u_out,
        amrex::RandomEngine const& engine) const noexcept {
        // A degenerate (zero) normal carries no geometric information:
        // leave the velocity unchanged.
        ParticleVector normal;
        if (!Math::Normalize(normal_to_domain, normal)) {
            u_out = u_in;
            return;
        }

        // Build an orthonormal tangential basis (tangent_one, tangent_two)
        // around the normal.
        ParticleVector tangent_one;
        ParticleVector tangent_two;
        Math::BuildOrthonormalBasis(
            normal, tangent_one, tangent_two);

        // Sample the wall Maxwellian: full Maxwellian in both tangential
        // directions, half-Maxwellian flux distribution along the
        // domain-pointing normal.
        amrex::ParticleReal const u_tangent_one =
            m_vth * static_cast<amrex::ParticleReal>(
                amrex::RandomNormal(amrex::Real(0.0), amrex::Real(1.0), engine));
        amrex::ParticleReal const u_tangent_two =
            m_vth * static_cast<amrex::ParticleReal>(
                amrex::RandomNormal(amrex::Real(0.0), amrex::Real(1.0), engine));
        // Exact inverse CDF of the zero-drift Gaussian flux distribution.
        // Scale in particle precision; the general flux helper uses field Real.
        // Subtract before converting so a double random draw cannot round to
        // one in single particle precision and produce log(0).
        amrex::ParticleReal const uniform = static_cast<amrex::ParticleReal>(
            amrex::Real(1.0) - amrex::Random(engine));
        amrex::ParticleReal const u_normal = m_vth *
            std::sqrt(-amrex::ParticleReal(2.0) * std::log(uniform));

        // Rotate the sampled velocity from the local (t1, t2, n) frame back
        // to Cartesian components.
        u_out = Math::ExpandInBasis(
            ParticleVector{u_tangent_one, u_tangent_two, u_normal},
            tangent_one, tangent_two, normal);

        amrex::ignore_unused(u_in);
    }
};

} // namespace Insert

#endif // WARPX_INSERT_WALLINTERACTIONOPERATORS_H_
