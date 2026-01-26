#include "WarpXFunc.h"

template <int depos_order>
void
ReplaceParticlesEachCell (int* p_delete,const int* p_offset, int* p_indices,
                         int numbins, WarpXParticleContainer& src_pc,
                         WarpXParIter& pti, amrex::MultiFab& ground_rho,
                         WarpXParticleContainer& dst_pc, amrex::Real inv_gap,
                         bool if_replace) {
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
                *pvx = soa_arr[PIdx::ux].dataPtr(),
                *pvy = soa_arr[PIdx::uy].dataPtr(),
                *pvz = soa_arr[PIdx::uz].dataPtr(),
                *pw = soa_arr[PIdx::w].dataPtr();
    //amrex::Print() << "start delete"<< src_pc.S<" particles\n";

    amrex::Box box = pti.tilebox();
    box.grow(ground_rho.nGrowVect());
    const amrex::XDim3 xyzmin = WarpX::LowerCorner(box, 0, 0._rt);
    const amrex::Dim3 lo = lbound(box);

    const auto& rho_arr_src = ground_rho.array(pti);

    //统计总数
    amrex::ParallelFor(numbins, [=] AMREX_GPU_DEVICE(int ibin) {
        int num = p_offset[ibin + 1] - p_offset[ibin],
            rest = amrex::min(p_delete[ibin], num);
        amrex::Gpu::Atomic::Add(p_num, rest);
    });
    int num_of_delete = all_deleted.dataValue();
    amrex::Vector<amrex::Real> rpx(num_of_delete), rpy(num_of_delete),
        rpz(num_of_delete), rvx(num_of_delete), rvy(num_of_delete),
        rvz(num_of_delete), rpw(num_of_delete);

    amrex::Gpu::DeviceVector<amrex::Real> device_rpx(num_of_delete),
        device_rpy(num_of_delete), device_rpz(num_of_delete),
        device_rvx(num_of_delete), device_rvy(num_of_delete),
        device_rvz(num_of_delete), device_rpw(num_of_delete);

    amrex::Real *p_rpx = device_rpx.dataPtr(), *p_rpy = device_rpy.dataPtr(),
                *p_rpz = device_rpz.dataPtr(), *p_rvx = device_rvx.dataPtr(),
                *p_rvy = device_rvy.dataPtr(), *p_rvz = device_rvz.dataPtr(),
                *p_rpw = device_rpw.dataPtr();

    amrex::Gpu::DeviceScalar<int> device_now(0);
    int* p_now = device_now.dataPtr();
    amrex::ParallelForRNG(
        numbins, [=] AMREX_GPU_DEVICE(int ibin, const amrex::RandomEngine& engine) {
            int num = p_offset[ibin + 1] - p_offset[ibin],
                rest = amrex::min(p_delete[ibin], num);
            while (rest > 0) {
                int indices_pos =
                    p_offset[ibin] + amrex::Random_int(num, engine);
                int pos = p_indices[indices_pos];
                auto pidw = amrex::ParticleIDWrapper{idcpu[pos]};
                if (pidw.is_valid()) {
                    pidw.make_invalid();
                    rest--;

                    amrex::ParticleReal x = px[pos], y = py[pos], z = pz[pos],
                                 vx = pvx[pos], vy = pvy[pos], vz = pvz[pos],
                                 w = pw[pos];
                    // 添加新粒子
                    if (if_replace) {
                        int copy_now = amrex::Gpu::Atomic::Add(p_now, 1);
                        copy_now--;
                        p_rpx[copy_now] = x;
                        p_rpy[copy_now] = y;
                        p_rpz[copy_now] = z;
                        p_rvx[copy_now] = vx;
                        p_rvy[copy_now] = vy;
                        p_rvz[copy_now] = vz;
                        p_rpw[copy_now] = w;
                    }

                    // 处理密度变化
                    w *= -(invvol);

                    Compute_shape_factor<depos_order> const
                        compute_shape_factor;

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
                                    sx[ix] * sy[iy] * sz[iz] * w);
                            }
                        }
                    }

                    /*
                    int pi = static_cast<int>(x / gap),
                        pj = static_cast<int>(y / gap),
                        pk = static_cast<int>(z / gap);

                    amrex::Real pid = x / gap, pjd = y / gap, pkd = z / gap,
                                dil = pid - pi, dih = 1 - dil, djl = pjd - pj,
                                djh = 1 - djl, dkl = pkd - pk, dkh = 1 - dkl;

                    amrex::Real coff[8];
                    coff[0] = dih * djh * dkh;
                    coff[1] = dil * djh * dkh;
                    coff[2] = dih * djl * dkh;
                    coff[3] = dil * djl * dkh;
                    coff[4] = dih * djh * dkl;
                    coff[5] = dil * djh * dkl;
                    coff[6] = dih * djl * dkl;
                    coff[7] = dil * dil * dkl;

                    // 处理边界条件
                    if (pi == domain.smallEnd(0)) {
                        coff[0] *= 2.0_prt;
                        coff[2] *= 2.0_prt;
                        coff[4] *= 2.0_prt;
                        coff[6] *= 2.0_prt;
                    }
                    if ((pi + 1) == domain.bigEnd(0)) {
                        coff[1] *= 2.0_prt;
                        coff[3] *= 2.0_prt;
                        coff[5] *= 2.0_prt;
                        coff[7] *= 2.0_prt;
                    }
                    if (pj == domain.smallEnd(1)) {
                        coff[0] *= 2.0_prt;
                        coff[1] *= 2.0_prt;
                        coff[4] *= 2.0_prt;
                        coff[5] *= 2.0_prt;
                    }
                    if ((pj + 1) == domain.bigEnd(1)) {
                        coff[2] *= 2.0_prt;
                        coff[3] *= 2.0_prt;
                        coff[6] *= 2.0_prt;
                        coff[7] *= 2.0_prt;
                    }
                    if (pk == domain.smallEnd(2)) {
                        coff[0] *= 2.0_prt;
                        coff[1] *= 2.0_prt;
                        coff[2] *= 2.0_prt;
                        coff[3] *= 2.0_prt;
                    }
                    if ((pk + 1) == domain.bigEnd(2)) {
                        coff[4] *= 2.0_prt;
                        coff[5] *= 2.0_prt;
                        coff[6] *= 2.0_prt;
                        coff[7] *= 2.0_prt;
                    }

                    amrex::Gpu::Atomic::AddNoRet(&rho_arr(pi, pj, pk),
                                                 coff[0] * w);
                    amrex::Gpu::Atomic::AddNoRet(&rho_arr(pi + 1, pj, pk),
                                                 coff[1] * w);
                    amrex::Gpu::Atomic::AddNoRet(&rho_arr(pi, pj + 1, pk),
                                                 coff[2] * w);
                    amrex::Gpu::Atomic::AddNoRet(&rho_arr(pi + 1, pj + 1, pk),
                                                 coff[3] * w);
                    amrex::Gpu::Atomic::AddNoRet(&rho_arr(pi, pj, pk + 1),
                                                 coff[4] * w);
                    amrex::Gpu::Atomic::AddNoRet(&rho_arr(pi + 1, pj, pk + 1),
                                                 coff[5] * w);
                    amrex::Gpu::Atomic::AddNoRet(&rho_arr(pi, pj + 1, pk + 1),
                                                 coff[6] * w);
                    amrex::Gpu::Atomic::AddNoRet(
                        &rho_arr(pi + 1, pj + 1, pk + 1), coff[7] * w);*/
                }
            }
        });

    if (if_replace) {
        // 添加粒子，目前仅转换类型
        const amrex::Vector<amrex::Vector<int>> nattr;
        dst_pc.AddNParticles(0, all_deleted.dataValue(), rpx, rpy, rpz, rvx,
                             rvy, rvz, 1, {rpw}, 0, nattr, false);
    }

    amrex::Print() << "end delete particles, number of deleted particle: "
                   << all_deleted.dataValue() << "\n";
}