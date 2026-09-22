#include "InsertBoundaryPhi.h"

#include "WarpX.H"

#include "Fields.H"
#include "Insert/Config/WarpXFunctionConfig.h"
#include "Insert/Config/WarpXSimulationConfig.h"
#include "Utils/WarpXConst.H"

#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#include <iostream>

namespace Insert {

void
VoltageAdjustment () {
#ifdef BENCHMARK_2D
    WarpX& warpx_instance = WarpX::GetInstance();
    auto phi_field = warpx_instance.m_fields.get(FieldType::phi_fp, 0);
    amrex::Real phisum = 0;

    static amrex::MultiFab mf(phi_field->boxArray(),
                              phi_field->DistributionMap(), phi_field->nComp(),
                              phi_field->nGrow());
    static bool DotInit = false;
    if (!DotInit) {
        mf.setVal(0);
        for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            const amrex::Box& box = mfi.tilebox();
            amrex::Array4<amrex::Real> const& arr = mf.array(mfi);
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j) {
                if (i == 491) {
                    arr(i, j, 0) = 0.48;
                } else if (i == 492) {
                    arr(i, j, 0) = 0.52;
                }
            });
        }
        DotInit = true;
    }
    phisum = amrex::MultiFab::Dot(mf, 0, *phi_field, 0, 1, 0);

    phisum /= 256;

    std::cout << "voltage adjustment: " << phisum << std::endl;

    for (amrex::MFIter mfi(*phi_field, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& box = mfi.tilebox();

        amrex::Array4<amrex::Real> const& phi = phi_field->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j) {
            phi(i, j, 0) -= phisum * i / 491.52;
        });
    }
#endif
}

void
DirichletPhiGuardSet () {
#ifndef WAVE1D
    using namespace amrex::literals;

    WarpX& warpx_instance = WarpX::GetInstance();
    auto rho = warpx_instance.m_fields.get(warpx::fields::FieldType::rho_fp, 0);
    auto phi = warpx_instance.m_fields.get(warpx::fields::FieldType::phi_fp, 0);
    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();
    const amrex::Real* dx_host = warpx_instance.Geom(0).CellSize();
    amrex::Gpu::DeviceVector<amrex::Real> dx_device(3);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, dx_host, dx_host + 3,
                     dx_device.begin());
    amrex::Real* dx = dx_device.dataPtr();

    for (amrex::MFIter mfi(*phi, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& box = mfi.validbox();
        const auto& phi_arr = phi->array(mfi);
        const auto& rho_arr = rho->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (i == domain.smallEnd(0)) {
                phi_arr(i - 1, j, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i + 1, j, k) -
                    rho_arr(i, j, k) * dx[0] * dx[0] / PhysConst::epsilon_0;
            }
            if (i == domain.bigEnd(0)) {
                phi_arr(i + 1, j, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i - 1, j, k) -
                    rho_arr(i, j, k) * dx[0] * dx[0] / PhysConst::epsilon_0;
            }
#if defined(WARPX_DIM_XZ)
            if (j == domain.smallEnd(1)) {
                phi_arr(i, j - 1, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j + 1, k) -
                    rho_arr(i, j, k) * dx[1] * dx[1] / PhysConst::epsilon_0;
            }
            if (j == domain.bigEnd(1)) {
                phi_arr(i, j + 1, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j - 1, k) -
                    rho_arr(i, j, k) * dx[1] * dx[1] / PhysConst::epsilon_0;
            }
#endif
#if defined(WARPX_DIM_3D)
            if (k == domain.smallEnd(2))
            {
                phi_arr(i, j, k - 1) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j, k + 1) -
                    rho_arr(i, j, k) * dx[2] * dx[2] / PhysConst::epsilon_0;
            }
            if (k == domain.bigEnd(2)) {
                phi_arr(i, j, k + 1) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j, k - 1) -
                    rho_arr(i, j, k) * dx[2] * dx[2] / PhysConst::epsilon_0;
            }
#endif
        });
    }
#endif
}

} // namespace Insert
