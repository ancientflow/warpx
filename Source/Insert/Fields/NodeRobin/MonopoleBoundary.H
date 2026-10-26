#ifndef WARPX_INSERT_POISSON_MONOPOLE_BOUNDARY_H_
#define WARPX_INSERT_POISSON_MONOPOLE_BOUNDARY_H_

#include <AMReX_Box.H>
#include <AMReX_GpuContainers.H>

namespace Insert::Poisson {
// Andrews et al., PSST 31 (2022) 035022, Eq. (7), Appendix A Eq. (A.7).
// Both normal-dependent terms change sign on converting inward to outward n.
// Thus dphi/dn + (n.r)/|r|^2 * (phi-phi_infinity) = 0 also holds outward.
struct MonopoleBoundary {
    bool enabled = false;
    amrex::GpuArray<amrex::Real,3> center{0.5,0.5,0.5};
    amrex::Real phi_infinity = 0;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real coefficient(int i, int j, int k, int face,
                            amrex::Box const& domain,
                            amrex::GpuArray<amrex::Real,3> const& problo,
                            amrex::GpuArray<amrex::Real,3> const& dxinv) const noexcept
    {
        int const p[3] = {i,j,k};
        amrex::Real r[3];
        amrex::Real radius2 = 0;
        for (int d=0; d<3; ++d) {
            r[d] = problo[d] + (p[d]-domain.smallEnd(d))/dxinv[d] - center[d];
            radius2 += r[d]*r[d];
        }
        return (face%2 ? r[face/2] : -r[face/2])/radius2;
    }
};
} // namespace Insert::Poisson
#endif
