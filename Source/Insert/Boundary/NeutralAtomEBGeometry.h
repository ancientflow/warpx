#ifndef WARPX_INSERT_NEUTRALATOMEBGEOMETRY_H_
#define WARPX_INSERT_NEUTRALATOMEBGEOMETRY_H_

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

    amrex::ParticleReal const inv_normal_norm =
        amrex::ParticleReal(1.0) /
        sqrt(cone.m_k * cone.m_k + amrex::ParticleReal(1.0));
    normal_to_domain = {
        -cone.m_k * x_relative / radius * inv_normal_norm,
        -cone.m_k * y_relative / radius * inv_normal_norm,
        inv_normal_norm};
    return true;
}

} // namespace Insert::NeutralAtomEBGeometry

#endif // WARPX_INSERT_NEUTRALATOMEBGEOMETRY_H_
