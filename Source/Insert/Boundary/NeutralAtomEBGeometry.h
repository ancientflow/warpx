#ifndef WARPX_INSERT_NEUTRALATOMEBGEOMETRY_H_
#define WARPX_INSERT_NEUTRALATOMEBGEOMETRY_H_

#include "Insert/Math/VectorOps.h"

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace Insert::NeutralAtomEBGeometry {

struct TruncatedConeGeometry {
    amrex::ParticleReal m_k = 0.0;
    amrex::ParticleReal m_a1 = 0.0;
    amrex::ParticleReal m_b1 = 0.0;
    amrex::ParticleReal m_axis_x = 0.0;
    amrex::ParticleReal m_axis_y = 0.0;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
GetTruncatedConeNormal (
    TruncatedConeGeometry const& cone, amrex::XDim3 const& x_hit,
    amrex::XDim3& normal_to_domain) noexcept
{
    using std::sqrt;

    amrex::ParticleReal const x_relative = x_hit.x - cone.m_axis_x;
    amrex::ParticleReal const y_relative = x_hit.y - cone.m_axis_y;
    amrex::ParticleReal const radius =
        sqrt(x_relative * x_relative + y_relative * y_relative);
    if (!(radius > cone.m_a1 && radius < cone.m_b1)) {
        normal_to_domain = {0.0, 0.0, 0.0};
        return false;
    }

    // The unnormalized normal has z component 1, so normalization cannot
    // fail; the return value is intentionally discarded.
    amrex::XDim3 const unnormalized_normal{
        static_cast<amrex::Real>(-cone.m_k * x_relative / radius),
        static_cast<amrex::Real>(-cone.m_k * y_relative / radius),
        amrex::Real(1.0)};
    static_cast<void>(Math::Normalize<amrex::ParticleReal>(
        unnormalized_normal, normal_to_domain));
    return true;
}

} // namespace Insert::NeutralAtomEBGeometry

#endif // WARPX_INSERT_NEUTRALATOMEBGEOMETRY_H_
