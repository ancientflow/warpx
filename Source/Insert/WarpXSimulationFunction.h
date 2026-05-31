//
// Created by ps on 2026/3/2.
//

#ifndef WARPX_WARPXSIMULATIONFUNCTION_H
#define WARPX_WARPXSIMULATIONFUNCTION_H

#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "Fields.H"
#ifdef WARPX_USE_FFT
#ifdef WARPX_DIM_RZ
#include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#else
#include "FieldSolver/SpectralSolver/SpectralSolver.H"
#endif
#endif
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Particles/ParticleCreation/FilterCopyTransform.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartUtils.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Random.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include "BackgroundCoupledDensity.h"
#include "WarpXFunctionConfig.h"
#include "WarpXSimulationConfig.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <ostream>
#include <string>
#include <vector>

template <int N = 1, typename SrcPC, typename DstPC, typename FilterFunc,
          typename TransformFunc>
amrex::Long
FilterCopyTransformParticleTiles (SrcPC& src_pc, DstPC& dst_pc,
                                  FilterFunc const& filter,
                                  TransformFunc const& transform,
                                  int lev_min = 0, int lev_max = -1) {
#ifdef AMREX_USE_GPU
    auto const* src_arena = src_pc.arena();
    auto const* dst_arena = dst_pc.arena();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_arena->isDevice() || src_arena->isManaged() ||
            src_arena == amrex::The_Pinned_Arena(),
        "FilterCopyTransformParticleTiles source particles must be in a "
        "device, managed, or AMReX pinned arena for GPU execution.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        dst_arena->isDevice() || dst_arena->isManaged() ||
            dst_arena == amrex::The_Pinned_Arena(),
        "FilterCopyTransformParticleTiles destination particles must be in a "
        "device, managed, or AMReX pinned arena for GPU execution.");
#endif

    const SmartCopyFactory copy_factory(src_pc, dst_pc);
    const auto Copy = copy_factory.getSmartCopy();

    if (lev_max < 0) {
        lev_max = src_pc.numLevels() - 1;
    }
    lev_max = std::min(lev_max, src_pc.numLevels() - 1);

    amrex::Long total_added = 0;

    for (int lev = lev_min; lev <= lev_max; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())                      \
    reduction(+ : total_added)
#endif
        for (WarpXParIter pti(src_pc, lev); pti.isValid(); ++pti) {
            auto& src_tile = pti.GetParticleTile();
            if (src_tile.numParticles() == 0) {
                continue;
            }

            auto& dst_tile = dst_pc.DefineAndReturnParticleTile(
                lev, pti.index(), pti.LocalTileIndex());
            const auto np_dst = dst_tile.numParticles();

            const auto num_added = filterCopyTransformParticles<N>(
                dst_pc, dst_tile, src_tile, np_dst, filter, Copy, transform);
            setNewParticleIDs(dst_tile, np_dst, num_added);

            total_added += static_cast<amrex::Long>(num_added);
        }
    }

    return total_added;
}

template <int N = 1, typename DstPC, typename FilterFunc,
          typename TransformFunc>
amrex::Long
FilterCopyTransformBoundaryBuffer (ParticleBoundaryBuffer& boundary_buffer,
                                   std::string const& src_species,
                                   int const boundary, DstPC& dst_pc,
                                   FilterFunc const& filter,
                                   TransformFunc const& transform,
                                   int lev_min = 0, int lev_max = -1) {
    auto* src_pc =
        boundary_buffer.getParticleBufferPointer(src_species, boundary);

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_pc != nullptr && src_pc->isDefined(),
        "The requested absorbing-boundary particle buffer is not defined. "
        "Enable the corresponding save_particles_at_* option and call this "
        "after ParticleBoundaryBuffer has gathered absorbed particles.");

    return FilterCopyTransformParticleTiles<N>(*src_pc, dst_pc, filter,
                                               transform, lev_min, lev_max);
}

struct SecondaryEmissionFilter {
    amrex::ParticleReal m_probability = 0.1;

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
    operator()(PData const& ptd, int const i,
               amrex::RandomEngine const& engine) const noexcept {
        if (!amrex::ParticleIDWrapper{ptd.m_idcpu[i]}.is_valid()) {
            return false;
        }
        return amrex::Random(engine) < m_probability;
    }
};

struct SecondaryEmissionTransform {
    int m_normal_index = 0;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_phi;
    amrex::ParticleReal m_eps = 0.0;
    amrex::ParticleReal m_vz_std = 0.0;

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    operator()(DstData& dst, SrcData const& src, int const i_src,
               int const i_dst,
               amrex::RandomEngine const& engine) const noexcept {
#if defined(WARPX_DIM_3D)
        const auto& p = src.getSuperParticle(i_src);
        amrex::ParticleReal x, y, z;
        get_particle_position(p, x, y, z);

        const amrex::ParticleReal ux_inc = src.m_rdata[PIdx::ux][i_src];
        const amrex::ParticleReal uy_inc = src.m_rdata[PIdx::uy][i_src];
        const amrex::ParticleReal uz_inc = src.m_rdata[PIdx::uz][i_src];
        const amrex::ParticleReal nz =
            src.m_runtime_rdata[m_normal_index + 2][i_src];

        const amrex::ParticleReal normal_sign = (nz >= amrex::ParticleReal(0.0))
                                                    ? amrex::ParticleReal(1.0)
                                                    : amrex::ParticleReal(-1.0);
        const amrex::ParticleReal z_boundary =
            (normal_sign > amrex::ParticleReal(0.0)) ? m_plo[2] : m_phi[2];

        amrex::ParticleReal x_emit = x;
        amrex::ParticleReal y_emit = y;
        if (uz_inc != amrex::ParticleReal(0.0)) {
            const amrex::ParticleReal dt_back = (z - z_boundary) / uz_inc;
            if (dt_back >= amrex::ParticleReal(0.0)) {
                x_emit -= dt_back * ux_inc;
                y_emit -= dt_back * uy_inc;
            }
        }

        const amrex::ParticleReal xlo = m_plo[0] + m_eps;
        const amrex::ParticleReal xhi = m_phi[0] - m_eps;
        const amrex::ParticleReal ylo = m_plo[1] + m_eps;
        const amrex::ParticleReal yhi = m_phi[1] - m_eps;
        if (x_emit < xlo) {
            x_emit = xlo;
        }
        if (x_emit > xhi) {
            x_emit = xhi;
        }
        if (y_emit < ylo) {
            y_emit = ylo;
        }
        if (y_emit > yhi) {
            y_emit = yhi;
        }

        dst.m_rdata[PIdx::x][i_dst] = x_emit;
        dst.m_rdata[PIdx::y][i_dst] = y_emit;
        dst.m_rdata[PIdx::z][i_dst] = z_boundary + normal_sign * m_eps;

        amrex::ParticleReal vz_emit =
            amrex::RandomNormal(amrex::ParticleReal(0.0), m_vz_std, engine);
        if (vz_emit < amrex::ParticleReal(0.0)) {
            vz_emit = -vz_emit;
        }

        dst.m_rdata[PIdx::ux][i_dst] = amrex::ParticleReal(0.0);
        dst.m_rdata[PIdx::uy][i_dst] = amrex::ParticleReal(0.0);
        dst.m_rdata[PIdx::uz][i_dst] = normal_sign * vz_emit;
#else
        amrex::ignore_unused(dst, src, i_src, i_dst, engine);
#endif
    }
};

#ifdef BENCHMARK_2D
// 电子数量计算
int
ElectronsCollection (WarpX& warpx_instance) {
    auto& boundary_buffer = warpx_instance.GetParticleBoundaryBuffer();
    int e_xlo =
            boundary_buffer.getNumParticlesInContainer("electrons", 0, false),
        // e_xhi =
        // boundary_buffer.getNumParticlesInContainer("electrons",1,false),
        xe_xlo =
            boundary_buffer.getNumParticlesInContainer("xe_ions", 0, false);
    // xe_xhi = boundary_buffer.getNumParticlesInContainer("xe_ions",1,false);
    return std::max(e_xlo - xe_xlo, 0);
}

// 粒子注入-模拟电离-阴极发射
void
ParticleInjection (WarpX& warpx_instance) {
    static const int nx = 512, ny = 256, nppc = 75;
    static const double lx = 0.025, ly = 0.0128, NPlasma = 5e16, s0 = 5.23e23,
                        const_dt = 5e-12;

    static const int once_injection = 1000; // 一次注入电子-离子对数量
    static double rest_macro_particles = 0; // 剩余应当注入粒子
    static const double one_step_real_particles =
        s0 * const_dt * 0.0128 * 2 / 3.1415926 *
        0.0075; // 每时步应当注入真实粒子
    static const Real global_weight =
        NPlasma / nx / ny * lx * ly / nppc; // 统一权重
    static const double one_step_macro_particles =
        one_step_real_particles / global_weight; // 每时步应当注入宏粒子

    int this_step_pair = 0,
        this_step_cathode_electron = ElectronsCollection(warpx_instance);
    rest_macro_particles += one_step_macro_particles;
    if (rest_macro_particles > once_injection) {
        rest_macro_particles -= once_injection;
        this_step_pair = once_injection;
    }

    auto& mypc = warpx_instance.GetPartContainer();
    auto& electons_pc =
        mypc.GetParticleContainer(mypc.getSpeciesID("electrons"));
    auto& xe_pc = mypc.GetParticleContainer(mypc.getSpeciesID("xe_ions"));

    const amrex::Vector<amrex::Vector<int>> nattr;
    amrex::RandomEngine uniform_engine(amrex::getRandState()),
        normal_engine(amrex::getRandState());

    int electron_size = this_step_pair + this_step_cathode_electron;
    std::cout << electron_size << " " << this_step_pair << " " << std::endl;
    amrex::Vector<amrex::Real> pxe(electron_size), pye(electron_size, 0),
        pze(electron_size), vxe(electron_size), vye(electron_size),
        vze(electron_size), we(electron_size, global_weight),
        pxxe(this_step_pair), pyxe(this_step_pair, 0), pzxe(this_step_pair),
        vxxe(this_step_pair), vyxe(this_step_pair), vzxe(this_step_pair),
        wxe(this_step_pair, global_weight);

    // 电子-离子对注入
    if (this_step_pair > 0) {
        // std::cout<<istart<<" "<<iend<<std::endl;
        for (size_t i = 0; i < this_step_pair; i++) {
            double r1 = amrex::Random(uniform_engine),
                   r2 = amrex::Random(uniform_engine);
            pxe[i] = 0.00625 + asin(2 * r1 - 1) / 3.14159 * 0.0075;
            pze[i] = 0.0128 * r2;
            pxxe[i] = pxe[i];
            pzxe[i] = pze[i];

            vxe[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vye[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vze[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vxxe[i] = amrex::RandomNormal(0, 606.34, normal_engine);
            vyxe[i] = amrex::RandomNormal(0, 606.34, normal_engine);
            vzxe[i] = amrex::RandomNormal(0, 606.34, normal_engine);
        }
    }

    // 阴极发射
    if (this_step_cathode_electron > 0) {
        for (size_t i = this_step_pair; i < electron_size; i++) {
            pxe[i] = 0.024;
            pze[i] = 0.0128 * amrex::Random(uniform_engine);
            vxe[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vye[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vze[i] = amrex::RandomNormal(0, 1326232, normal_engine);
        }
    }
    if (this_step_pair > 0) {
        std::cout << "add pairs number: " << this_step_pair << std::endl;
        xe_pc.AddNParticles(0, this_step_pair, pxxe, pyxe, pzxe, vxxe, vyxe,
                            vzxe, 1, {wxe}, 0, nattr, false);
    }

    // 阴极发射部分
    // 清除缓存
    if (electron_size > 0) {
        electons_pc.AddNParticles(0, electron_size, pxe, pye, pze, vxe, vye,
                                  vze, 1, {we}, 0, nattr, false);
        if (this_step_cathode_electron > 0) {
            auto& boundary_buffer = warpx_instance.GetParticleBoundaryBuffer();
            boundary_buffer.clearParticles();
        }
    }
}

// 电压修正
void
VoltageAdjustment (WarpX& warpx_instance) {
    auto phi_field = warpx_instance.m_fields.get(FieldType::phi_fp, 0);
    // fields.get(FieldType::phi_fp,0);
    amrex::Real phisum = 0;

    static MultiFab mf(phi_field->boxArray(), phi_field->DistributionMap(),
                       phi_field->nComp(), phi_field->nGrow());
    static bool DotInit = false;
    if (!DotInit) {
        mf.setVal(0);
        for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            const Box& box = mfi.tilebox();
            amrex::Array4<amrex::Real> const& arr = mf.array(mfi);
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j) {
                if (i == 491)
                    arr(i, j, 0) = 0.48;
                else if (i == 492)
                    arr(i, j, 0) = 0.52;
            });
        }
        DotInit = true;
    }
    phisum = amrex::MultiFab::Dot(mf, 0, *phi_field, 0, 1, 0);

    phisum /= 256;

    std::cout << "voltage adjustment: " << phisum << std::endl;

    for (amrex::MFIter mfi(*phi_field, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const Box& box = mfi.tilebox();

        amrex::Array4<amrex::Real> const& phi = phi_field->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j) {
            phi(i, j, 0) -= phisum * i / 491.52;
        });
    }
}
#endif

#ifdef HALL3D
// 3D阴极发射
void
CathodeInjection3D () {
    static const double L = 0.05; // 未缩比边长
    static double Ic, dt, l_factor, elec_weight;

    static bool ifinit = false;

    if (!ifinit) {
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("dt", dt);
        pp_mc.get("Ic", Ic);
        pp_mc.get("l_factor", l_factor);
        pp_mc.getWithParser("elec_weight", elec_weight);
        ifinit = true;
    }

    const double num_of_real_electrons =
                     Ic * dt / l_factor / l_factor / 1.60217e-19,
                 num_of_marco_particle = num_of_real_electrons / elec_weight;

    WarpX& warpx_instance = WarpX::GetInstance();

    static double electron_rest = 0;
    const int one_step_injection = static_cast<int>(num_of_marco_particle) + 1;

    electron_rest += num_of_marco_particle;

    static const double r1 = 0.016 / l_factor, r2 = 0.02 / l_factor,
                        dr = r2 - r1, sumr = r1 + r2, multr = r1 * r2,
                        z1 = 0.042 / l_factor, z2 = 0.046 / l_factor,
                        dz = z2 - z1, Pi2 = 3.1415926 * 2, sigma = 592982,
                        half_l = L / 2 / l_factor;

    double r, theta;
    static amrex::Vector<amrex::Vector<int>> nattr;
    if (electron_rest > one_step_injection) {

        electron_rest -= one_step_injection;
        // std::cout << "electron inject: " << one_step_injection << std::endl;

        amrex::Vector<amrex::ParticleReal> pz(one_step_injection),
            px(one_step_injection), py(one_step_injection),
            vx(one_step_injection), vy(one_step_injection),
            vz(one_step_injection), pw(one_step_injection, elec_weight);

        amrex::RandomEngine uniform_engine(amrex::getRandState()),
            normal_engine(amrex::getRandState());

        for (int i = 0; i < one_step_injection; i++) {
            // 柱坐标空间均匀分布
            r = amrex::Random(uniform_engine) * dr + r1;
            theta = amrex::Random(uniform_engine) * Pi2;
            pz[i] = amrex::Random(uniform_engine) * dz + z1;

            // 转化到直角坐标空间均匀分布
            r = sqrt(sumr * r - multr);
            px[i] = r * cos(theta) + half_l;
            py[i] = r * sin(theta) + half_l;

            // 速度
            vx[i] = amrex::RandomNormal(0, sigma, normal_engine);
            vy[i] = amrex::RandomNormal(0, sigma, normal_engine);
            vz[i] = amrex::RandomNormal(0, sigma, normal_engine);
        }

        auto& mypc = warpx_instance.GetPartContainer();
        auto& e_pc = mypc.GetParticleContainer(mypc.getSpeciesID("electrons"));
        e_pc.AddNParticles(0, one_step_injection, px, py, pz, vx, vy, vz, 1,
                           {pw}, 0, nattr, 0);
    }
}

// 初始注入等离子体
void
PlasmaInit () {
    const double L = 0.05; // 未缩比边长
    double Ic, dt, l_factor, elec_weight;

        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("dt", dt);
        pp_mc.get("l_factor", l_factor);
        pp_mc.getWithParser("elec_weight", elec_weight);

    const double r1 = 0.0105 / l_factor, r2 = 0.0155 / l_factor,
                        dr = r2 - r1, sumr = r1 + r2, multr = r1 * r2,
                        z1 = 0.001 / l_factor, z2 = 0.004 / l_factor,
                        dz = z2 - z1, Pi2 = 3.1415926 * 2, sigma_e = 592982,
                        sigma_xe = 1212.41, half_l = L / 2 / l_factor;

    const double volume = (r2 * r2 - r1 * r1) * 3.14159 * dz;
    const auto num_elec = static_cast<int>(volume * 1e18 / elec_weight);

    WarpX& warpx_instance = WarpX::GetInstance();

    double r, theta;
    static amrex::Vector<amrex::Vector<int>> nattr;

    amrex::Vector<amrex::ParticleReal> pz(num_elec), px(num_elec),
        py(num_elec), vx(num_elec), vy(num_elec),
        vz(num_elec), pw(num_elec, elec_weight);

    amrex::RandomEngine uniform_engine(amrex::getRandState()),
        normal_engine(amrex::getRandState());

    for (int i = 0; i < num_elec; i++) {
        // 柱坐标空间均匀分布
        r = amrex::Random(uniform_engine) * dr + r1;
        theta = amrex::Random(uniform_engine) * Pi2;
        pz[i] = amrex::Random(uniform_engine) * dz + z1;

        // 转化到直角坐标空间均匀分布
        r = sqrt(sumr * r - multr);
        px[i] = r * cos(theta) + half_l;
        py[i] = r * sin(theta) + half_l;

        // 速度
        vx[i] = amrex::RandomNormal(0, sigma_e, normal_engine);
        vy[i] = amrex::RandomNormal(0, sigma_e, normal_engine);
        vz[i] = amrex::RandomNormal(0, sigma_e, normal_engine);
    }
    auto& mypc = warpx_instance.GetPartContainer();
    auto& e_pc = mypc.GetParticleContainer(mypc.getSpeciesID("electrons"));
    e_pc.AddNParticles(0, num_elec, px, py, pz, vx, vy, vz, 1, {pw}, 0, nattr,
                       0);
    for (int i = 0; i < num_elec; i++) {
        // 速度
        vx[i] = amrex::RandomNormal(0, sigma_xe, normal_engine);
        vy[i] = amrex::RandomNormal(0, sigma_xe, normal_engine);
        vz[i] = amrex::RandomNormal(0, sigma_xe, normal_engine);
    }

    auto& xe_pc = mypc.GetParticleContainer(mypc.getSpeciesID("xe_ions"));
    xe_pc.AddNParticles(0, num_elec, px, py, pz, vx, vy, vz, 1, {pw}, 0, nattr,
                       0);
}

// 氙气注入
void
XeInjection () {
    static const double L = 0.05, // 未缩比边长
        NA = 6.03e23;             // 阿伏伽德罗常数
    static bool ifinit = false, ifhole = false;
    static int hole_num = 48;
    static double dt, l_factor, atom_weight, m_dot, Tx, Ty, Tz, vz0;
    static amrex::Vector<double> hole_x, hole_y;

    if (!ifinit) {
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("dt", dt);
        pp_mc.get("l_factor", l_factor);
        pp_mc.getWithParser("xe_weight", atom_weight);
        pp_mc.get("m_dot", m_dot);
        pp_mc.query("ifhole", ifhole);
        pp_mc.get("Tx", Tx);
        pp_mc.get("Ty", Ty);
        pp_mc.get("Tz", Tz);
        pp_mc.get("vz0", vz0);
        ifinit = true;
    }
    static const double V_per_sec = 0.06 / 60. / 1e6,
                        n_per_sec = V_per_sec * 101325 / 8.314 / 273.15 * NA,
                        n_per_step = n_per_sec * dt, mxe = 2.179e-25,
                        n_marco_per_step =
                            n_per_step / atom_weight, // 每时步应当注入宏粒子
        Pi2 = 3.1415926 * 2;

    const double n_marco_per_step_m =
        m_dot * dt / l_factor / l_factor / mxe / atom_weight;

    const unsigned one_times_inject_particle = static_cast<int>(
        (n_marco_per_step_m * 100) + 1); // 一次注入粒子数 100时步一次
    static double xe_rest = 0;

    const double rmid = (0.021 + 0.031) / 4 / l_factor,
                 rwidth = 0.001 / l_factor, r1 = rmid - rwidth,
                 r2 = rmid + rwidth, sqdr = (r2 * r2) - (r1 * r1),
                 sqr1 = r1 * r1, kb = 1.38064852e-23,
                 sigmax = std::sqrt(kb * Tx / mxe),
                 sigmay = std::sqrt(kb * Ty / mxe),
                 sigmaz = std::sqrt(kb * Tz / mxe), half_L = L / 2 / l_factor,
                 sqdr_hole = rwidth * rwidth;
    double r, theta;
    static const amrex::Vector<amrex::Vector<int>> nattr;

    static bool if_hole_init = false;
    if (!if_hole_init) {
        const double per_angle = Pi2 / hole_num;
        double angle = 0;
        for (int i = 0; i < hole_num; i++) {
            hole_x.push_back(rmid * cos(angle));
            hole_y.push_back(rmid * sin(angle));
            angle += per_angle;
        }
        if_hole_init = true;
    }

    xe_rest += n_marco_per_step_m;

    if (xe_rest > one_times_inject_particle) {
        WarpX& warpx_instance = WarpX::GetInstance();

        xe_rest -= one_times_inject_particle;

        amrex::Vector<amrex::ParticleReal> pz(one_times_inject_particle, 0),
            px(one_times_inject_particle), py(one_times_inject_particle),
            vx(one_times_inject_particle), vy(one_times_inject_particle),
            vz(one_times_inject_particle),
            pw(one_times_inject_particle, atom_weight);

        amrex::RandomEngine uniform_engine(amrex::getRandState()),
            normal_engine(amrex::getRandState());
        int hole_start = rand() % hole_num;
        for (int i = 0; i < one_times_inject_particle; i++) {

            // 柱坐标空间均匀分布
            // 转化到直角坐标空间均匀分布
            if (!ifhole) {
                // 非小孔进气
                theta = amrex::Random(uniform_engine) * Pi2;
                r = sqrt(sqdr * amrex::Random(uniform_engine) + sqr1);
                px[i] = r * cos(theta) + half_L;
                py[i] = r * sin(theta) + half_L;
            } else {
                // 小孔进气
                theta = amrex::Random(uniform_engine) * Pi2;
                r = rwidth * sqrt(amrex::Random(uniform_engine));
                // 小孔进气是r1 = 0, r2 = rwidth, 化简后得到该式
                px[i] = r * cos(theta) + half_L +
                        hole_x[(i + hole_start) % hole_num];
                py[i] = r * sin(theta) + half_L +
                        hole_y[(i + hole_start) % hole_num];
            }

            // 速度
            vx[i] = amrex::RandomNormal(0, sigmax, normal_engine);
            vy[i] = amrex::RandomNormal(0, sigmay, normal_engine);
            do {
                vz[i] = amrex::RandomNormal(0, sigmaz, normal_engine) + vz0;
            } while (vz[i] < 0);
        }

        auto& mypc = warpx_instance.GetPartContainer();
        auto& xe_pc =
            mypc.GetParticleContainer(mypc.getSpeciesID("xe_netural"));
        xe_pc.AddNParticles(0, one_times_inject_particle, px, py, pz, vx, vy,
                            vz, 1, {pw}, 0, nattr, 0);
        amrex::Print() << "Injection Xe Atom\n";
    }
}

// 阳极电压修正
void
AnodeVoltage () {
    WarpX& warpx_instance = WarpX::GetInstance();
    amrex::ParmParse pp_mc("my_constants");
    double ncell, voltage, L, l_factor;
    pp_mc.getWithParser("n_cell", ncell);
    pp_mc.get("L", L);
    pp_mc.get("l_factor", l_factor);
    pp_mc.query("voltage", voltage);
    double halfncell = ncell / 2;

    double sim_L = L / l_factor, gap = sim_L / ncell,
           index_sq_min =
               (0.021 / 2 / l_factor / gap) * (0.021 / 2 / l_factor / gap),
           index_sq_max =
               (0.031 / 2 / l_factor / gap) * (0.031 / 2 / l_factor / gap);

    auto phi_field =
        warpx_instance.m_fields.get(warpx::fields::FieldType::phi_fp, 0);
    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();
    for (amrex::MFIter mfi(*phi_field, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        amrex::Array4<amrex::Real> const& phi = phi_field->array(mfi);
        if (!domain.strictly_contains(box)) {
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k == domain.smallEnd(2)) {
                    double i_ex = i - halfncell, j_ex = j - halfncell;
                    if (i_ex * i_ex + j_ex * j_ex >= index_sq_min &&
                        i_ex * i_ex + j_ex * j_ex <= index_sq_max)
                        phi(i, j, k) = voltage;
                    else
                        phi(i, j, k) = 0;
                }
            });
        }
    }
    phi_field->FillBoundary(warpx_instance.Geom(0).periodicity());
}

// 阳极电流及全场电子计算
void
AnodeCurrentCalc () {
    static int times = 0;
    static const int gap = 10;
    times++;

    static bool ifinit = false;
    if (times == gap) {
        // 打开输出文件
        std::fstream fileout;
        if (!ifinit) {
            fileout.open("anode_current.dat", std::ios::out);
            fileout << "time\t"
                       "zmin_electron\tzmin_ion\tanode_electron\t"
                       "anode_electron_cut\n";
            ifinit = true;
        } else {
            fileout.open("anode_current.dat", std::ios::app);
        }
        times = 0;

        // 获取粒子容器信息
        WarpX& warpx_instance = WarpX::GetInstance();
        auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();

        auto& elec_zmin = mybpc.getParticleBuffer("electrons", 4);
        auto& xe_ion_zmin = mybpc.getParticleBuffer("xe_ions", 4);

        const int data_size = 4;
        amrex::Gpu::DeviceVector<amrex::ParticleReal> device_charge(data_size,
                                                                    0);
        amrex::Vector<amrex::ParticleReal> host_charge(data_size, 0);

        amrex::ParticleReal* device_ptr = device_charge.dataPtr();

        amrex::ParticleReal l_factor, half_L, L, rn1, rn2;
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("l_factor", l_factor);
        pp_mc.get("L", L);

        half_L = L / l_factor / 2;
        rn1 = 0.021 / 2 / l_factor;
        rn2 = 0.031 / 2 / l_factor;

        // 统计zmin电子权重
        for (auto pti = WarpXParIter(elec_zmin, 0); pti.isValid(); ++pti) {
            auto& arr = pti.GetStructOfArrays().GetRealData();
            auto px = arr[PIdx::x].dataPtr();
            auto py = arr[PIdx::y].dataPtr();
            auto pw = arr[PIdx::w].dataPtr();
            auto pz = arr[PIdx::z].dataPtr();
            auto pvx = arr[PIdx::ux].dataPtr();
            auto pvy = arr[PIdx::uy].dataPtr();
            auto pvz = arr[PIdx::uz].dataPtr();
            int np = pti.numParticles();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) {
                amrex::ParticleReal x = px[ip], y = py[ip], w = pw[ip],
                                    z = pz[ip], vx = pvx[ip], vy = pvy[ip],
                                    vz = pvz[ip];
                amrex::ParticleReal dt = z / vz;
                x -= half_L;
                y -= half_L;
                amrex::ParticleReal r = sqrt(x * x + y * y);
                if (r >= rn1 && r <= rn2) {
                    amrex::Gpu::Atomic::Add(&device_ptr[data_size - 2], w);
                }
                x -= dt * vx;
                y -= dt * vy;
                r = sqrt(x * x + y * y);
                if (r >= rn1 && r <= rn2) {
                    amrex::Gpu::Atomic::Add(&device_ptr[data_size - 1], w);
                }
                amrex::Gpu::Atomic::Add(&device_ptr[0], w);
            });
        }

        // 统计zmin 氙离子权重
        for (auto pti = WarpXParIter(xe_ion_zmin, 0); pti.isValid(); ++pti) {
            auto& arr = pti.GetStructOfArrays().GetRealData();
            auto px = arr[PIdx::x].dataPtr();
            auto py = arr[PIdx::y].dataPtr();
            auto pw = arr[PIdx::w].dataPtr();
            auto pz = arr[PIdx::z].dataPtr();
            auto pvx = arr[PIdx::ux].dataPtr();
            auto pvy = arr[PIdx::uy].dataPtr();
            auto pvz = arr[PIdx::uz].dataPtr();
            int np = pti.numParticles();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) {
                amrex::ParticleReal x = px[ip], y = py[ip], w = pw[ip],
                                    z = pz[ip], vx = pvx[ip], vy = pvy[ip],
                                    vz = pvz[ip];
                amrex::ParticleReal dt = z / vz;
                x -= half_L;
                y -= half_L;
                x -= dt * vx;
                y -= dt * vy;
                amrex::ParticleReal r = sqrt(x * x + y * y);
                if (r >= rn1 && r <= rn2) {
                    amrex::Gpu::Atomic::Add(&device_ptr[1], w);
                }
            });
        }

        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_charge.begin(),
                         device_charge.end(), host_charge.begin());
        // 写入文件 清理缓存
        fileout << warpx_instance.gett_new(0);
        for (int i = 0; i < host_charge.size(); i++) {
            fileout << "\t" << host_charge[i];
        }
        fileout << "\n";
        fileout.close();
        mybpc.clearParticles();
    }
}

inline void
SecondaryEmission () {
    static int times = 0;
    static const int gap = 10;
    times++;

    if (times == gap) {
        times = 0;

        WarpX& warpx_instance = WarpX::GetInstance();

        // 获取边界电子容器
        auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();
        auto& elec_zmin = mybpc.getParticleBuffer("electrons", 4);
        const int normal_index = elec_zmin.GetRealCompIndex("nx") -
                                 WarpXParticleContainer::NArrayReal;

        // 获取电子容器
        auto& mypc = warpx_instance.GetPartContainer();
        auto& elec_pc = mypc.GetParticleContainerFromName("electrons");

        const auto plo = warpx_instance.Geom(0).ProbLoArray();
        const auto phi = warpx_instance.Geom(0).ProbHiArray();
        const auto dx = warpx_instance.Geom(0).CellSizeArray();
        amrex::Real min_dx = dx[0];
        for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
            min_dx = std::min(min_dx, dx[idim]);
        }

        constexpr amrex::ParticleReal emission_energy_eV =
            amrex::ParticleReal(2.0);
        const amrex::ParticleReal vz_std =
            static_cast<amrex::ParticleReal>(std::sqrt(
                2.0 * emission_energy_eV * PhysConst::q_e / PhysConst::m_e));

        const SecondaryEmissionFilter filter{amrex::ParticleReal(0.5)};
        const SecondaryEmissionTransform transform{
            normal_index, plo, phi, amrex::ParticleReal(0.1 * min_dx), vz_std};

        int num_added = FilterCopyTransformBoundaryBuffer(
            mybpc, "electrons", 4, elec_pc, filter, transform);
        amrex::Print() << "Emission electron: " << num_added << "\n";
    }
}

#endif

#ifdef HALL3D_INIT
// 初始化氙气注入
void
XeFastInjection () {
    double atom_weight, l_factor, m_dot, Tx, Ty, Tz, vz0;
    bool ifhole = false;
    static amrex::Vector<double> hole_x, hole_y;
    static int hole_num = 48;

    amrex::ParmParse pp_mc("my_constants");
    pp_mc.get("l_factor", l_factor);
    pp_mc.getWithParser("xe_weight", atom_weight);
    pp_mc.get("m_dot", m_dot);
    pp_mc.query("ifhole", ifhole);
    pp_mc.get("Tx", Tx);
    pp_mc.get("Ty", Ty);
    pp_mc.get("Tz", Tz);
    pp_mc.get("vz0", vz0);

    static const double L = 0.05, // 未缩比边长
        NA = 6.03e23;             // 阿伏伽德罗常数
    const double V_per_sec = 0.06 / 60. / 1e6, dt = 5.6e-10,
                 n_per_sec = V_per_sec * 101325 / 8.314 / 273.15 * NA,
                 n_per_step = n_per_sec * dt,
                 n_marco_per_step =
                     n_per_step / atom_weight, // 每时步应当注入宏粒子;
        mxe = 2.179e-25, Pi2 = 3.1415926 * 2;
    // 质量流量
    const double n_marco_per_step_m =
        m_dot * dt / mxe / l_factor / l_factor / atom_weight;

    const double one_times_inject_particle =
        static_cast<int>(n_marco_per_step_m + 1); // 一次注入粒子数
    static double xe_rest = 0;
    static const double rmid = (0.021 + 0.031) / 4 / l_factor,
                        rwidth = 0.001 / l_factor, r1 = rmid - rwidth,
                        r2 = rmid + rwidth, sqdr = r2 * r2 - r1 * r1,
                        sqr1 = r1 * r1, kb = 1.38064852e-23,
                        sigmax = std::sqrt(kb * Tx / mxe),
                        sigmay = std::sqrt(kb * Ty / mxe),
                        sigmaz = std::sqrt(kb * Tz / mxe),
                        half_L = L / 2 / l_factor;
    double r, theta;
    static const amrex::Vector<amrex::Vector<int>> nattr;

    static bool if_hole_init = false;
    if (!if_hole_init) {
        double per_angle = Pi2 / hole_num, angle = 0;
        for (int i = 0; i < hole_num; i++) {
            hole_x.push_back(rmid * cos(angle));
            hole_y.push_back(rmid * sin(angle));
            angle += per_angle;
        }
        if_hole_init = true;
    }

    xe_rest += n_marco_per_step_m;
    if (xe_rest > one_times_inject_particle) {
        WarpX& warpx_instance = WarpX::GetInstance();

        xe_rest -= one_times_inject_particle;

        amrex::Vector<amrex::ParticleReal> pz(one_times_inject_particle, 0),
            px(one_times_inject_particle), py(one_times_inject_particle),
            vx(one_times_inject_particle), vy(one_times_inject_particle),
            vz(one_times_inject_particle, vz0),
            pw(one_times_inject_particle, atom_weight);

        amrex::RandomEngine uniform_engine(amrex::getRandState()),
            normal_engine(amrex::getRandState());

        int hole_start = rand() % hole_num;
        for (int i = 0; i < one_times_inject_particle; i++) {
            // 柱坐标空间均匀分布
            // 转化到直角坐标空间均匀分布

            theta = amrex::Random(uniform_engine) * Pi2;
            r = sqrt(sqdr * amrex::Random(uniform_engine) + sqr1);

            if (!ifhole) {
                // 非小孔进气
                theta = amrex::Random(uniform_engine) * Pi2;
                r = sqrt(sqdr * amrex::Random(uniform_engine) + sqr1);
                px[i] = r * cos(theta) + half_L;
                py[i] = r * sin(theta) + half_L;
            } else {
                // 小孔进气
                theta = amrex::Random(uniform_engine) * Pi2;
                r = rwidth * sqrt(amrex::Random(
                                 uniform_engine)); // 小孔进气是r1 = 0, r2 =
                                                   // rwidth, 化简后得到该式
                px[i] = r * cos(theta) + half_L +
                        hole_x[(i + hole_start) % hole_num];
                py[i] = r * sin(theta) + half_L +
                        hole_y[(i + hole_start) % hole_num];
            }

            // 速度
            vx[i] = amrex::RandomNormal(0, sigmax, normal_engine);
            vy[i] = amrex::RandomNormal(0, sigmay, normal_engine);
            do {
                vz[i] = amrex::RandomNormal(0, sigmaz, normal_engine) + vz0;
            } while (vz[i] < 0);
        }

        auto& mypc = warpx_instance.GetPartContainer();
        auto& xe_pc =
            mypc.GetParticleContainer(mypc.getSpeciesID("xe_netural"));
        xe_pc.AddNParticles(0, one_times_inject_particle, px, py, pz, vx, vy,
                            vz, 1, {pw}, 0, nattr, 0);
        amrex::Print() << "Injection Xe Atom\n";
    }
}
#endif

#endif // WARPX_WARPXSIMULATIONFUNCTION_H
