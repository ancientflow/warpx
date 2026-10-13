/*
 * AMReX Copyright (c) 2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * (2) Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * (3) Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef CHT_NODE_ROBIN_KERNELS_H_
#define CHT_NODE_ROBIN_KERNELS_H_
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>
static_assert(AMREX_SPACEDIM == 3, "CHT nodal Robin requires three dimensions");
namespace cht::detail {
using amrex::Array4;
using amrex::Real;
using amrex::GpuArray;
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_adotx_aa (int i, int j, int k, Array4<Real const> const& x,
                       Array4<Real const> const& sig, Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        return Real(0.0);
    } else {
        Real facx = Real(1.0/36.0)*dxinv[0]*dxinv[0];
        Real facy = Real(1.0/36.0)*dxinv[1]*dxinv[1];
        Real facz = Real(1.0/36.0)*dxinv[2]*dxinv[2];
        Real fxyz = facx + facy + facz;
        Real fmx2y2z = -facx + Real(2.0)*facy + Real(2.0)*facz;
        Real f2xmy2z = Real(2.0)*facx - facy + Real(2.0)*facz;
        Real f2x2ymz = Real(2.0)*facx + Real(2.0)*facy - facz;
        Real f4xm2ym2z = Real(4.0)*facx - Real(2.0)*facy - Real(2.0)*facz;
        Real fm2x4ym2z = -Real(2.0)*facx + Real(4.0)*facy - Real(2.0)*facz;
        Real fm2xm2y4z = -Real(2.0)*facx - Real(2.0)*facy + Real(4.0)*facz;
        return x(i,j,k)*Real(-4.0)*fxyz*
            (sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)
            +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ))
            + fxyz*(x(i-1,j-1,k-1)*sig(i-1,j-1,k-1)
                  + x(i+1,j-1,k-1)*sig(i  ,j-1,k-1)
                  + x(i-1,j+1,k-1)*sig(i-1,j  ,k-1)
                  + x(i+1,j+1,k-1)*sig(i  ,j  ,k-1)
                  + x(i-1,j-1,k+1)*sig(i-1,j-1,k  )
                  + x(i+1,j-1,k+1)*sig(i  ,j-1,k  )
                  + x(i-1,j+1,k+1)*sig(i-1,j  ,k  )
                  + x(i+1,j+1,k+1)*sig(i  ,j  ,k  ))
            + fmx2y2z*(x(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1))
                     + x(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1))
                     + x(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  ))
                     + x(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  )))
            + f2xmy2z*(x(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1))
                     + x(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1))
                     + x(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  ))
                     + x(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  )))
            + f2x2ymz*(x(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k))
                     + x(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k))
                     + x(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k))
                     + x(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k)))
            + f4xm2ym2z*(x(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k))
                       + x(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k)))
            + fm2x4ym2z*(x(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k))
                       + x(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k)))
            + fm2xm2y4z*(x(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1))
                       + x(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_x (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        return (w1*crse(ic,jc,kc)+w2*crse(ic+1,jc,kc))/(w1+w2);
    }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_y (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        Real w2 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        return (w1*crse(ic,jc,kc)+w2*crse(ic,jc+1,kc))/(w1+w2);
    }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_z (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        Real w2 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1*crse(ic,jc,kc)+w2*crse(ic,jc,kc+1))/(w1+w2);
    }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_face_xy (Array4<Real const> const& crse, Array4<Real const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        Real w4 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        return (w1 * aa_interp_line_y(crse,sig,i-1,j  ,k,ic  ,jc  ,kc) +
                w2 * aa_interp_line_y(crse,sig,i+1,j  ,k,ic+1,jc  ,kc) +
                w3 * aa_interp_line_x(crse,sig,i  ,j-1,k,ic  ,jc  ,kc) +
                w4 * aa_interp_line_x(crse,sig,i  ,j+1,k,ic  ,jc+1,kc)) / (w1+w2+w3+w4);
    }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_face_xz (Array4<Real const> const& crse, Array4<Real const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
        Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
        Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        Real w4 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sig,i-1,j,k  ,ic  ,jc,kc  ) +
                w2 * aa_interp_line_z(crse,sig,i+1,j,k  ,ic+1,jc,kc  ) +
                w3 * aa_interp_line_x(crse,sig,i  ,j,k-1,ic  ,jc,kc  ) +
                w4 * aa_interp_line_x(crse,sig,i  ,j,k+1,ic  ,jc,kc+1)) / (w1+w2+w3+w4);
    }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_face_yz (Array4<Real const> const& crse, Array4<Real const> const& sig,
                            int i, int j, int k, int ic, int jc, int kc) noexcept
    {
        Real w1 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
        Real w2 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
        Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
        Real w4 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
        return (w1 * aa_interp_line_z(crse,sig,i,j-1,k  ,ic,jc  ,kc  ) +
                w2 * aa_interp_line_z(crse,sig,i,j+1,k  ,ic,jc+1,kc  ) +
                w3 * aa_interp_line_y(crse,sig,i,j  ,k-1,ic,jc  ,kc  ) +
                w4 * aa_interp_line_y(crse,sig,i,j  ,k+1,ic,jc  ,kc+1)) / (w1+w2+w3+w4);
    }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_aa (int i, int j, int k, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sig,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        int kc = amrex::coarsen(k,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        bool k_is_odd = (kc*2 != k);
        if (i_is_odd && j_is_odd && k_is_odd) {
            // Fine node at center of cell
            Real w1 = sig(i-1,j-1,k-1) + sig(i-1,j,k-1) + sig(i-1,j-1,k) + sig(i-1,j,k);
            Real w2 = sig(i  ,j-1,k-1) + sig(i  ,j,k-1) + sig(i  ,j-1,k) + sig(i  ,j,k);
            Real w3 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j-1,k) + sig(i,j-1,k);
            Real w4 = sig(i-1,j  ,k-1) + sig(i,j  ,k-1) + sig(i-1,j  ,k) + sig(i,j  ,k);
            Real w5 = sig(i-1,j-1,k-1) + sig(i,j-1,k-1) + sig(i-1,j,k-1) + sig(i,j,k-1);
            Real w6 = sig(i-1,j-1,k  ) + sig(i,j-1,k  ) + sig(i-1,j,k  ) + sig(i,j,k  );
            fine(i,j,k) += (w1 * aa_interp_face_yz(crse,sig,i-1,j  ,k  ,ic  ,jc  ,kc  ) +
                            w2 * aa_interp_face_yz(crse,sig,i+1,j  ,k  ,ic+1,jc  ,kc  ) +
                            w3 * aa_interp_face_xz(crse,sig,i  ,j-1,k  ,ic  ,jc  ,kc  ) +
                            w4 * aa_interp_face_xz(crse,sig,i  ,j+1,k  ,ic  ,jc+1,kc  ) +
                            w5 * aa_interp_face_xy(crse,sig,i  ,j  ,k-1,ic  ,jc  ,kc  ) +
                            w6 * aa_interp_face_xy(crse,sig,i  ,j  ,k+1,ic  ,jc  ,kc+1))
                / (w1+w2+w3+w4+w5+w6);
        } else if (j_is_odd && k_is_odd) {
            // Node on a Y-Z face
            fine(i,j,k) += aa_interp_face_yz(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd && k_is_odd) {
            // Node on a Z-X face
            fine(i,j,k) += aa_interp_face_xz(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd && j_is_odd) {
            // Node on a X-Y face
            fine(i,j,k) += aa_interp_face_xy(crse,sig,i,j,k,ic,jc,kc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,k) += aa_interp_line_x(crse,sig,i,j,k,ic,jc,kc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,k) += aa_interp_line_y(crse,sig,i,j,k,ic,jc,kc);
        } else if (k_is_odd) {
            // Node on Z line
            fine(i,j,k) += aa_interp_line_z(crse,sig,i,j,k,ic,jc,kc);
        } else {
            // Node coincident with coarse node
            fine(i,j,k) += crse(ic,jc,kc);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction (int i, int j, int k, Array4<Real> const& crse,
                          Array4<Real const> const& fine, Array4<int const> const& msk) noexcept
{
    int ii = i*2;
    int jj = j*2;
    int kk = k*2;
    if (msk(ii,jj,kk)) {
        crse(i,j,k) = Real(0.0);
    } else {
        crse(i,j,k) = Real(1./64.)*(fine(ii-1,jj-1,kk-1)+fine(ii+1,jj-1,kk-1)
                                   +fine(ii-1,jj+1,kk-1)+fine(ii+1,jj+1,kk-1)
                                   +fine(ii-1,jj-1,kk+1)+fine(ii+1,jj-1,kk+1)
                                   +fine(ii-1,jj+1,kk+1)+fine(ii+1,jj+1,kk+1))
                    + Real(1./32.)*(fine(ii  ,jj-1,kk-1)+fine(ii  ,jj+1,kk-1)
                                   +fine(ii  ,jj-1,kk+1)+fine(ii  ,jj+1,kk+1)
                                   +fine(ii-1,jj  ,kk-1)+fine(ii+1,jj  ,kk-1)
                                   +fine(ii-1,jj  ,kk+1)+fine(ii+1,jj  ,kk+1)
                                   +fine(ii-1,jj-1,kk  )+fine(ii+1,jj-1,kk  )
                                   +fine(ii-1,jj+1,kk  )+fine(ii+1,jj+1,kk  ))
                    + Real(1./16.)*(fine(ii-1,jj,kk)+fine(ii+1,jj,kk)
                                   +fine(ii,jj-1,kk)+fine(ii,jj+1,kk)
                                   +fine(ii,jj,kk-1)+fine(ii,jj,kk+1))
                      + Real(1./8.)*fine(ii,jj,kk);
    }
}
} // namespace cht::detail
#endif
