#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#include "Fields.H"
#ifdef WARPX_USE_FFT
#ifdef WARPX_DIM_RZ
#include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#else
#include "FieldSolver/SpectralSolver/SpectralSolver.H"
#endif
#endif
#include "Parallelization/GuardCellManager.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Utils/WarpXUtil.H"

#include <ablastr/utils/SignalHandling.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <Insert/WarpXFunctionConfig.h>
#include <BoundaryConditions/WarpX_PEC.H>
#include <algorithm>
#include <array>
#include <memory>
#include <ostream>
#include <vector>

using namespace amrex;
using ablastr::utils::SignalHandling;

inline void
AtomDepostiAPI (WarpXParticleContainer& pc, amrex::MultiFab& rho, const int lev) {
    for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {
        const Box& box = pti.validbox();

        auto np = pti.numParticles();
        amrex::Print() << np << " API deposit \n";
        // Extract particle data
        auto& attribs = pti.GetAttribs();
        auto& wp = attribs[PIdx::w];
        pc.DepositCharge(pti, wp, nullptr, &rho, 0, 0, np, 0, lev, lev);
    }
    // 处理边界密度
    // 尝试使用内置函数进行修改，便于进行高阶插值
    // 考虑到更多的层级，必须采用这种内置函数
    auto& warpx_instance = WarpX::GetInstance();
    /*    PEC::ApplyReflectiveBoundarytoRhofield(
            &rho, warpx_instance.field_boundary_lo,
            warpx_instance.field_boundary_hi,
       warpx_instance.particle_boundary_lo, warpx_instance.particle_boundary_hi,
       warpx_instance.Geom(lev), lev, PatchType::fine,
       warpx_instance.refRatio());*/

    amrex::Box domain = warpx_instance.Geom(lev).Domain();
    domain.surroundingNodes();
    for (MFIter mfi(rho, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const auto& rho_arr = rho[mfi].array();
        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // x direction
            if (i == domain.smallEnd(0)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
            if (i == domain.bigEnd(0)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
            // y direction
            if (j == domain.smallEnd(1)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
            if (j == domain.bigEnd(1)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
            // z direction
            if (k == domain.smallEnd(2)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
            if (k == domain.bigEnd(2)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
        });
    }
}