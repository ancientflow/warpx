#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#include "EmbeddedBoundary/Enabled.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "Fields.H"
#ifdef WARPX_USE_FFT
#ifdef WARPX_DIM_RZ
#include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#else
#include "FieldSolver/SpectralSolver/SpectralSolver.H"
#endif
#endif
#include "Fluids/MultiFluidContainer.H"
#include "Fluids/WarpXFluidContainer.H"
#include "Parallelization/GuardCellManager.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Python/callbacks.H"
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

#include "WarpXFunctionConfig.h"
#include <algorithm>
#include <array>
#include <memory>
#include <ostream>
#include <vector>


/**
 * @brief 使用WarpX的沉积函数进行密度沉积
 * @param[in] pc 粒子容器
 * @param[in] rho 密度Multifab
 
extern void AtomDepostiAPI(WarpXParticleContainer& pc, amrex::MultiFab& rho);*/

/**
 * @brief 替换粒子，本质上是先删除后添加
 * @param[in] p_delete 每个细胞删除的粒子数
 * @param[in] p_offset 偏移数组，每个细胞起始粒子序号
 * @param[in] p_indices 粒子ID数组，以p_offset中序号分组
 * @param[in] numbins 细胞数量
 * @param[in] src_pc 源粒子容器（要删除的粒子）
 * @param[in] pti 粒子迭代器
 * @param[in] rho_arr_src 源粒子密度场
 * @param[in] dst_pc 替换目标粒子容器（要添加的粒子）
 * @param[in] rho_arr_dst 替换粒子密度场
 * @param[in] inv_gap 网格边长倒数
 * @param[in] if_replace 是否替换（不替换即仅删除）
 */
template <int depos_order>
void
ReplaceParticlesEachCell (
    int* p_delete, const int* p_offset, int* p_indices, int numbins,
    WarpXParticleContainer& src_pc, WarpXParIter& pti,
    amrex::MultiFab& ground_rho, amrex::MultiFab& excitation_rho,
    WarpXParticleContainer& dst_pc, amrex::Real inv_gap,int* p_mask,
    bool if_replace = true) {

    amrex::Gpu::DeviceScalar<int> all_deleted(0);
    int* p_num = all_deleted.dataPtr();

    amrex::Real invvol = inv_gap * inv_gap * inv_gap;
    auto& ptile = src_pc.ParticlesAt(0, pti);
    auto& soa = ptile.GetStructOfArrays();

    uint64_t* const AMREX_RESTRICT idcpu = soa.GetIdCPUData().data();
    auto& soa_arr = soa.GetRealData();
    amrex::Real *px = soa_arr[PIdx::x].dataPtr(),
                *py = soa_arr[PIdx::y].dataPtr(),
                *pz = soa_arr[PIdx::z].dataPtr(),
                *pw = soa_arr[PIdx::w].dataPtr();

    amrex::Box box = pti.tilebox();
    box.grow(ground_rho.nGrowVect());
    const amrex::XDim3 xyzmin = WarpX::LowerCorner(box, 0, 0._rt);
    const amrex::Dim3 lo = lbound(box);

    const auto& rho_arr_src = ground_rho.array(pti);
    const auto& rho_arr_dst = excitation_rho.array(pti);

    amrex::ParallelForRNG(numbins, [=] AMREX_GPU_DEVICE(
                                       int ibin,
                                       const amrex::RandomEngine& engine) {
        int num = p_offset[ibin + 1] - p_offset[ibin],
            rest = amrex::min(p_delete[ibin], num);
        while (rest > 0) {
            int indices_pos = p_offset[ibin] + amrex::Random_int(num, engine);
            int pos = p_indices[indices_pos];
            auto pidw = amrex::ParticleIDWrapper{idcpu[pos]};
            if (pidw.is_valid()) {
                pidw.make_invalid();
                rest--;
                amrex::Gpu::Atomic::AddNoRet(p_num, 1);

                amrex::ParticleReal x = px[pos], y = py[pos], z = pz[pos],
                                    w = pw[pos];
                // 添加新粒子
                if (if_replace) {
                    p_mask[pos] = true;
                }
                // 处理密度变化
                w *= invvol;

                Compute_shape_factor<depos_order> const compute_shape_factor;

                amrex::Real sx[depos_order + 1] = {0._rt},
                                             sy[depos_order + 1] = {0._rt},
                                             sz[depos_order + 1] = {0._rt};

                int px = compute_shape_factor(sx, (x - xyzmin.x) * inv_gap);
                int py = compute_shape_factor(sy, (y - xyzmin.y) * inv_gap);
                int pz = compute_shape_factor(sz, (z - xyzmin.z) * inv_gap);

                for (int iz = 0; iz <= depos_order; iz++) {
                    for (int iy = 0; iy <= depos_order; iy++) {
                        for (int ix = 0; ix <= depos_order; ix++) {
                            amrex::Gpu::Atomic::AddNoRet(
                                &rho_arr_src(lo.x + px + ix, lo.y + py + iy,
                                             lo.z + pz + iz),
                                sx[ix] * sy[iy] * sz[iz] * -w);
                        }
                    }
                }

                for (int iz = 0; iz <= depos_order; iz++) {
                    for (int iy = 0; iy <= depos_order; iy++) {
                        for (int ix = 0; ix <= depos_order; ix++) {
                            amrex::Gpu::Atomic::AddNoRet(
                                &rho_arr_dst(lo.x + px + ix, lo.y + py + iy,
                                             lo.z + pz + iz),
                                sx[ix] * sy[iy] * sz[iz] * w);
                        }
                    }
                }
            }
        }
    });

    if (if_replace) {
        amrex::Print() << "excitation: replace " << all_deleted.dataValue() << " particle"
                       << src_pc.getSpeciesId() + 1 << " to particle"
                       << dst_pc.getSpeciesId() + 1 << "\n";
    } else {
        amrex::Print() << "delete " << all_deleted.dataValue() << " particle"
                       << src_pc.getSpeciesId() + 1 << "\n";
    }
}