#include "ZmaxRadialExitStats.h"

#include "Insert/Utils/InsertUtils.h"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#if defined(WARPX_DIM_3D)
using namespace amrex::literals;

namespace {

constexpr int zhi_boundary = 5;

/** \brief Static configuration, parsed once from my_constants.zmax_radial_*. */
struct ZmaxRadialExitConfig {
    bool enabled = false;
    std::string species = "xe_neutral";
    int start_step = 0;
    amrex::Real r_max = 0.0_rt;
    amrex::Real v_max = 1000.0_rt;
    int v_bins = 100;
    std::string dir = "zmax_radial";
    amrex::Real dx = 0.0_rt;
    amrex::Real zmax = 0.0_rt;
    int n_r = 0;
};

ZmaxRadialExitConfig
ReadZmaxRadialExitConfig ()
{
    ZmaxRadialExitConfig cfg;

    amrex::ParmParse pp_mc("my_constants");
    bool enabled = false;
    pp_mc.query("zmax_radial_exit_diag", enabled);
    cfg.enabled = enabled;
    if (!cfg.enabled) { return cfg; }

    pp_mc.query("zmax_radial_species", cfg.species);
    pp_mc.query("zmax_radial_start_step", cfg.start_step);
    pp_mc.query("zmax_radial_v_max", cfg.v_max);
    pp_mc.query("zmax_radial_v_bins", cfg.v_bins);
    pp_mc.query("zmax_radial_dir", cfg.dir);

    amrex::Geometry const& geom = WarpX::GetInstance().Geom(0);
    amrex::Real const hi_x = geom.ProbHi(0);
    amrex::Real const hi_y = geom.ProbHi(1);
    amrex::Real const lo_x = geom.ProbLo(0);
    amrex::Real const lo_y = geom.ProbLo(1);
    // Default: radius of the domain corner, so the whole exit plane is covered.
    amrex::Real const r_corner = std::sqrt(
        std::max(hi_x * hi_x, lo_x * lo_x) + std::max(hi_y * hi_y, lo_y * lo_y));
    cfg.r_max = r_corner;
    pp_mc.query("zmax_radial_r_max", cfg.r_max);

    cfg.dx = geom.CellSize(0);
    cfg.zmax = geom.ProbHi(2);

    // The input file must guarantee that r_max is an integer multiple of dx.
    amrex::Real const n_r_real = cfg.r_max / cfg.dx;
    cfg.n_r = static_cast<int>(std::llround(n_r_real));
    cfg.n_r = std::max(cfg.n_r, 1);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::abs(cfg.n_r * cfg.dx - cfg.r_max) < 1e-6_rt * cfg.dx,
        "my_constants.zmax_radial_r_max must be an integer multiple of the "
        "grid spacing dx (radial bins follow the grid).");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        cfg.v_max > 0.0_rt && cfg.v_bins > 0,
        "my_constants.zmax_radial_v_max and zmax_radial_v_bins must be "
        "positive.");

    amrex::Print() << "zmax_radial_exit_diag: species=" << cfg.species
                   << " start_step=" << cfg.start_step
                   << " r_max=" << cfg.r_max << " n_r=" << cfg.n_r
                   << " dx=" << cfg.dx << " v_max=" << cfg.v_max
                   << " v_bins=" << cfg.v_bins << "\n";
    return cfg;
}

/** \brief Host-side accumulators, persistent across steps (only the values on
 *         the IO process are meaningful, the output is written there). */
struct ZmaxRadialExitAccum {
    std::vector<amrex::Real> weight;      // per radial bin, total weight
    std::vector<amrex::Real> num_macro;   // per radial bin, macro count
    std::vector<amrex::Real> vhist;       // [n_r][3][v_bins] flat, weights
    amrex::Real overflow_weight = 0.0_rt; // weight with r >= r_max
    amrex::Real accum_time = 0.0_rt;      // physical time covered [s]
    int last_clear_step = -1;             // step at which buffer last cleared
    bool started = false;
};

} // namespace
#endif

namespace Insert {

void
ZmaxRadialExitStatsCalc ()
{
#if defined(WARPX_DIM_3D)
    static ZmaxRadialExitConfig const cfg = ReadZmaxRadialExitConfig();
    if (!cfg.enabled) { return; }

    WarpX& warpx_instance = WarpX::GetInstance();
    int const step = warpx_instance.getistep(0);
    bool const is_final = (step == warpx_instance.maxStep() - 1);

    static ZmaxRadialExitAccum accum;
    if (accum.weight.empty()) {
        accum.weight.assign(cfg.n_r, 0.0_rt);
        accum.num_macro.assign(cfg.n_r, 0.0_rt);
        accum.vhist.assign(
            static_cast<size_t>(cfg.n_r) * 3 * cfg.v_bins, 0.0_rt);
    }

    if (DoBoundaryParticleDiag(step)) {
        if (step >= cfg.start_step) {
            accum.started = true;

            // Accumulate the zhi buffer content into device arrays, then
            // reduce to the IO process and merge into the host accumulators.
            int const n_hist =
                cfg.n_r * 3 * cfg.v_bins;
            amrex::Gpu::DeviceVector<amrex::Real> d_weight(cfg.n_r, 0.0_rt);
            amrex::Gpu::DeviceVector<amrex::Real> d_num(cfg.n_r, 0.0_rt);
            amrex::Gpu::DeviceVector<amrex::Real> d_vhist(n_hist, 0.0_rt);
            amrex::Gpu::DeviceVector<amrex::Real> d_overflow(1, 0.0_rt);

            auto* const AMREX_RESTRICT weight_ptr = d_weight.dataPtr();
            auto* const AMREX_RESTRICT num_ptr = d_num.dataPtr();
            auto* const AMREX_RESTRICT vhist_ptr = d_vhist.dataPtr();
            auto* const AMREX_RESTRICT overflow_ptr = d_overflow.dataPtr();

            auto& buffer = warpx_instance.GetParticleBoundaryBuffer();
            auto* scraped =
                buffer.getParticleBufferPointer(cfg.species, zhi_boundary);

            amrex::ParticleReal const zmax =
                static_cast<amrex::ParticleReal>(cfg.zmax);
            amrex::ParticleReal const dx =
                static_cast<amrex::ParticleReal>(cfg.dx);
            amrex::ParticleReal const v_min =
                static_cast<amrex::ParticleReal>(-cfg.v_max);
            amrex::ParticleReal const inv_dv =
                static_cast<amrex::ParticleReal>(
                    cfg.v_bins / (2.0_rt * cfg.v_max));
            int const n_r = cfg.n_r;
            int const v_bins = cfg.v_bins;

            if (scraped != nullptr && scraped->isDefined()) {
                for (auto pti = WarpXParIter(*scraped, 0); pti.isValid();
                     ++pti) {
                    auto& arr = pti.GetStructOfArrays().GetRealData();
                    auto const* const AMREX_RESTRICT px =
                        arr[PIdx::x].dataPtr();
                    auto const* const AMREX_RESTRICT py =
                        arr[PIdx::y].dataPtr();
                    auto const* const AMREX_RESTRICT pz =
                        arr[PIdx::z].dataPtr();
                    auto const* const AMREX_RESTRICT pw =
                        arr[PIdx::w].dataPtr();
                    auto const* const AMREX_RESTRICT pvx =
                        arr[PIdx::ux].dataPtr();
                    auto const* const AMREX_RESTRICT pvy =
                        arr[PIdx::uy].dataPtr();
                    auto const* const AMREX_RESTRICT pvz =
                        arr[PIdx::uz].dataPtr();
                    int const np = pti.numParticles();

                    amrex::ParallelFor(
                        np, [=] AMREX_GPU_DEVICE(long ip) noexcept {
                            amrex::ParticleReal const vx = pvx[ip];
                            amrex::ParticleReal const vy = pvy[ip];
                            amrex::ParticleReal const vz = pvz[ip];
                            amrex::ParticleReal x_hit, y_hit;
                            if (!BacktraceParticleToZPlane(
                                    px[ip], py[ip], pz[ip], vx, vy, vz, zmax,
                                    x_hit, y_hit)) {
                                return;
                            }
                            amrex::ParticleReal const r =
                                std::sqrt(x_hit * x_hit + y_hit * y_hit);
                            int const ir = static_cast<int>(r / dx);
                            amrex::ParticleReal const w = pw[ip];
                            if (ir < 0 || ir >= n_r) {
                                amrex::Gpu::Atomic::Add(
                                    overflow_ptr, static_cast<amrex::Real>(w));
                                return;
                            }
                            amrex::Gpu::Atomic::Add(
                                &weight_ptr[ir], static_cast<amrex::Real>(w));
                            amrex::Gpu::Atomic::Add(&num_ptr[ir], 1.0_rt);
                            amrex::ParticleReal const vcomp[3] = {vx, vy, vz};
                            for (int comp = 0; comp < 3; ++comp) {
                                int const iv = static_cast<int>(
                                    (vcomp[comp] - v_min) * inv_dv);
                                if (iv >= 0 && iv < v_bins) {
                                    amrex::Gpu::Atomic::Add(
                                        &vhist_ptr[(ir * 3 + comp) * v_bins +
                                                   iv],
                                        static_cast<amrex::Real>(w));
                                }
                            }
                        });
                }
            }

            // Device -> host, reduce to the IO process, merge.
            amrex::Vector<amrex::Real> h_weight(cfg.n_r);
            amrex::Vector<amrex::Real> h_num(cfg.n_r);
            amrex::Vector<amrex::Real> h_vhist(n_hist);
            amrex::Vector<amrex::Real> h_overflow(1);
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_weight.begin(),
                             d_weight.end(), h_weight.begin());
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_num.begin(),
                             d_num.end(), h_num.begin());
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_vhist.begin(),
                             d_vhist.end(), h_vhist.begin());
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_overflow.begin(),
                             d_overflow.end(), h_overflow.begin());
            amrex::ParallelDescriptor::ReduceRealSum(
                h_weight.data(), static_cast<int>(h_weight.size()),
                amrex::ParallelDescriptor::IOProcessorNumber());
            amrex::ParallelDescriptor::ReduceRealSum(
                h_num.data(), static_cast<int>(h_num.size()),
                amrex::ParallelDescriptor::IOProcessorNumber());
            amrex::ParallelDescriptor::ReduceRealSum(
                h_vhist.data(), static_cast<int>(h_vhist.size()),
                amrex::ParallelDescriptor::IOProcessorNumber());
            amrex::ParallelDescriptor::ReduceRealSum(
                h_overflow.data(), static_cast<int>(h_overflow.size()),
                amrex::ParallelDescriptor::IOProcessorNumber());

            if (amrex::ParallelDescriptor::IOProcessor()) {
                for (int ir = 0; ir < cfg.n_r; ++ir) {
                    accum.weight[ir] += h_weight[ir];
                    accum.num_macro[ir] += h_num[ir];
                }
                for (int i = 0; i < n_hist; ++i) {
                    accum.vhist[i] += h_vhist[i];
                }
                accum.overflow_weight += h_overflow[0];
            }
        }

        // The buffer covers all steps since the previous clear; keep the
        // covered physical time in sync and prevent unbounded growth.
        if (amrex::ParallelDescriptor::IOProcessor()) {
            if (step >= cfg.start_step) {
                accum.accum_time += static_cast<amrex::Real>(
                    step - accum.last_clear_step) * warpx_instance.getdt(0);
            }
        }
        accum.last_clear_step = step;
        warpx_instance.GetParticleBoundaryBuffer().clearParticles(
            zhi_boundary);
    }

    if (is_final && accum.started &&
        amrex::ParallelDescriptor::IOProcessor()) {
        CreateDirectoryTree(cfg.dir);

        // Density distribution: exiting number flux per radial annulus.
        {
            std::fstream density_file(
                PathJoin(cfg.dir, "density.dat"), std::ios::out);
            density_file << "# zmax exit radial density distribution\n"
                         << "# species " << cfg.species << "\n"
                         << "# start_step " << cfg.start_step << "\n"
                         << "# end_step " << step << "\n"
                         << "# accum_time_s " << accum.accum_time << "\n"
                         << "# dx_m " << cfg.dx << "\n"
                         << "# r_max_m " << cfg.r_max << "\n"
                         << "# overflow_weight (r >= r_max) "
                         << accum.overflow_weight << "\n"
                         << "r_lo_m\tr_hi_m\tr_center_m\tweight\tnum_macro\t"
                            "flux_per_m2_s\n";
            for (int ir = 0; ir < cfg.n_r; ++ir) {
                amrex::Real const r_lo = cfg.dx * ir;
                amrex::Real const r_hi = r_lo + cfg.dx;
                amrex::Real const area =
                    MathConst::pi * (r_hi * r_hi - r_lo * r_lo);
                amrex::Real const flux =
                    accum.accum_time > 0.0_rt
                        ? accum.weight[ir] / (area * accum.accum_time)
                        : 0.0_rt;
                density_file << r_lo << "\t" << r_hi << "\t"
                             << 0.5_rt * (r_lo + r_hi) << "\t"
                             << accum.weight[ir] << "\t" << accum.num_macro[ir]
                             << "\t" << flux << "\n";
            }
        }

        // Velocity distributions: per radial bin, vx/vy/vz histograms.
        {
            std::fstream vdist_file(
                PathJoin(cfg.dir, "vdist.dat"), std::ios::out);
            vdist_file << "# zmax exit velocity distributions per radial bin\n"
                       << "# species " << cfg.species << "\n"
                       << "# start_step " << cfg.start_step << "\n"
                       << "# end_step " << step << "\n"
                       << "# accum_time_s " << accum.accum_time << "\n"
                       << "# v range [m/s]: [" << -cfg.v_max << ", "
                       << cfg.v_max << "], bins " << cfg.v_bins << "\n"
                       << "# ux/uy/uz are treated as vx/vy/vz (gamma ~ 1)\n"
                       << "r_center_m\tv_lo\tv_hi\tv_center\tw_vx\tw_vy\tw_vz\n";
            amrex::Real const dv =
                2.0_rt * cfg.v_max / static_cast<amrex::Real>(cfg.v_bins);
            for (int ir = 0; ir < cfg.n_r; ++ir) {
                amrex::Real const r_center = cfg.dx * (ir + 0.5_rt);
                for (int iv = 0; iv < cfg.v_bins; ++iv) {
                    amrex::Real const v_lo = -cfg.v_max + dv * iv;
                    int const base = (ir * 3) * cfg.v_bins + iv;
                    vdist_file << r_center << "\t" << v_lo << "\t"
                               << v_lo + dv << "\t" << v_lo + 0.5_rt * dv
                               << "\t" << accum.vhist[base] << "\t"
                               << accum.vhist[base + cfg.v_bins] << "\t"
                               << accum.vhist[base + 2 * cfg.v_bins] << "\n";
                }
                vdist_file << "\n";
            }
        }

        amrex::Print() << "zmax_radial_exit_diag: wrote " << cfg.dir
                       << "/density.dat and vdist.dat (step " << step
                       << ")\n";
    }
#endif
}

} // namespace Insert
