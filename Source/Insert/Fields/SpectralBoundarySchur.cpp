#include "SpectralBoundarySchur.h"

#include "Insert/Config/WarpXSimulationConfig.h"

#ifdef HALL3D

#include "EmbeddedBoundary/Enabled.H"
#include "Fields.H"
#include "Insert/Boundary/ZMinWallCharge.h"
#include "Insert/Diagnostics/InsertRuntimeDiagnostics.h"
#include "Insert/Utils/InsertUtils.h"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_Math.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

namespace {

struct HallAnodeRingConfig {
    amrex::Real center_x;
    amrex::Real center_y;
    amrex::Real r_min_sq;
    amrex::Real r_max_sq;
};

/**
 * Read the HALL3D anode-ring geometry from the existing Insert constants.
 *
 * @param geom Level geometry used to place the ring center in physical space.
 * @return Ring center and squared radius bounds.
 */
HallAnodeRingConfig
ReadHallAnodeRingConfig (amrex::Geometry const& geom) {
    amrex::ParmParse pp_mc("my_constants");
    auto l_factor = static_cast<amrex::Real>(1.0);
    pp_mc.get("l_factor", l_factor);

    amrex::Real const r_min = static_cast<amrex::Real>(0.021) /
                              (static_cast<amrex::Real>(2.0) * l_factor);
    amrex::Real const r_max = static_cast<amrex::Real>(0.031) /
                              (static_cast<amrex::Real>(2.0) * l_factor);

    return HallAnodeRingConfig{
        static_cast<amrex::Real>(0.5) * (geom.ProbLo(0) + geom.ProbHi(0)),
        static_cast<amrex::Real>(0.5) * (geom.ProbLo(1) + geom.ProbHi(1)),
        r_min * r_min, r_max * r_max};
}

/**
 * Test whether a nodal point belongs to the zmin HALL3D anode ring.
 *
 * @param i x node index.
 * @param j y node index.
 * @param k z node index.
 * @param zlo zmin node index.
 * @param xlo xmin node index.
 * @param ylo ymin node index.
 * @param problo_x Physical x lower bound.
 * @param problo_y Physical y lower bound.
 * @param dx x cell size.
 * @param dy y cell size.
 * @param config Anode ring geometry.
 * @return true if the node is on the anode ring.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
IsHallAnodeRingNode (int i, int j, int k, int zlo, int xlo, int ylo,
                     amrex::Real problo_x, amrex::Real problo_y, amrex::Real dx,
                     amrex::Real dy, HallAnodeRingConfig const& config) {
    amrex::Real const x = problo_x + (i - xlo) * dx;
    amrex::Real const y = problo_y + (j - ylo) * dy;
    amrex::Real const r_sq = (x - config.center_x) * (x - config.center_x) +
                             (y - config.center_y) * (y - config.center_y);

    return k == zlo && r_sq >= config.r_min_sq && r_sq <= config.r_max_sq;
}

struct SchurConfig {
    bool enabled = false;
    int max_iter = 30;
    amrex::Real rel_tol = static_cast<amrex::Real>(1.0e-6);
    amrex::Real abs_tol = static_cast<amrex::Real>(0.0);
    std::string face = "zmin";
};

struct SchurState {
    bool config_read = false;
    SchurConfig config;
    std::unique_ptr<amrex::MultiFab> phi_corr;
    std::unique_ptr<amrex::MultiFab> neumann_mask;
    std::unique_ptr<amrex::MultiFab> u_face;
    std::unique_ptr<amrex::FFT::R2X<amrex::Real>> r2x;
    amrex::Box r2x_domain;
    int nx_internal = 0;
    int ny_internal = 0;
    bool fft_cleanup_registered = false;
};

/**
 * Access the module-local persistent state.
 *
 * @return Mutable Schur state singleton.
 */
SchurState&
State () {
    static SchurState state;
    return state;
}

/**
 * Read Schur runtime parameters once from `insert.schur_boundary`.
 *
 * @return Cached Schur configuration.
 */
SchurConfig const&
ReadConfig () {
    auto& state = State();
    if (!state.config_read) {
        amrex::ParmParse pp("insert.schur_boundary");
        pp.query("enabled", state.config.enabled);
        pp.query("max_iter", state.config.max_iter);
        pp.query("rel_tol", state.config.rel_tol);
        pp.query("abs_tol", state.config.abs_tol);
        pp.query("face", state.config.face);
        state.config_read = true;
    }
    return state.config;
}

/**
 * Evaluate `coth(x)` without overflowing for large positive `x`.
 *
 * @param x Positive argument.
 * @return Stable `coth(x)`.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
StableCoth (amrex::Real x) {
    constexpr amrex::Real threshold = static_cast<amrex::Real>(40.0);
    if (x > threshold) {
        return static_cast<amrex::Real>(1.0);
    }
    amrex::Real const exp_neg_2x = std::exp(-static_cast<amrex::Real>(2.0) * x);
    return (static_cast<amrex::Real>(1.0) + exp_neg_2x) /
           (static_cast<amrex::Real>(1.0) - exp_neg_2x);
}

#if defined(WARPX_DIM_3D)

/**
 * Select retained modal count for a normalized z depth.
 *
 * @param full_modes Full number of modes in one transverse direction.
 * @param z_fraction Depth from zmin normalized by `Lz`.
 * @return Number of modes retained in that direction.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
LayerModes (int full_modes, amrex::Real z_fraction) {
    int divisor = 1;
    if (z_fraction >= static_cast<amrex::Real>(0.5)) {
        divisor = 16;
    } else if (z_fraction >= static_cast<amrex::Real>(0.25)) {
        divisor = 8;
    } else if (z_fraction >= static_cast<amrex::Real>(0.125)) {
        divisor = 4;
    } else if (z_fraction >= static_cast<amrex::Real>(0.0625)) {
        divisor = 2;
    }
    return std::max(1, full_modes / divisor);
}

/**
 * Build a zero-based slab box for one zmin-face field.
 *
 * @param nx Number of x face entries.
 * @param ny Number of y face entries.
 * @return Box `[0,nx-1] x [0,ny-1] x [0,0]`.
 */
amrex::Box
MakeFaceSlabBox (int nx, int ny) {
    return amrex::Box(amrex::IntVect(0, 0, 0),
                      amrex::IntVect(nx - 1, ny - 1, 0));
}

/**
 * Build a one-component, zero-grow MultiFab on a zero-based face slab.
 *
 * @param nx Number of x face entries.
 * @param ny Number of y face entries.
 * @return Newly allocated face MultiFab.
 */
std::unique_ptr<amrex::MultiFab>
MakeFaceSlabMultiFab (int nx, int ny) {
    amrex::BoxArray const ba(MakeFaceSlabBox(nx, ny));
    amrex::DistributionMapping const dm(ba);
    return std::make_unique<amrex::MultiFab>(ba, dm, 1, 0);
}

/**
 * Register AMReX-finalize cleanup for cached FFT plans.
 */
void
RegisterFFTCleanup () {
    auto& state = State();
    if (!state.fft_cleanup_registered) {
        amrex::ExecOnFinalize([]() {
            State().r2x.reset();
        });
        state.fft_cleanup_registered = true;
    }
}

/**
 * Prepare the cached AMReX odd/odd real-to-real FFT plan.
 *
 * @param nx Number of x interior nodes/modes.
 * @param ny Number of y interior nodes/modes.
 * @return Cached R2X plan for the requested face.
 */
amrex::FFT::R2X<amrex::Real>&
PrepareSchurFFT (int nx, int ny) {
    auto& state = State();
    amrex::Box const domain = MakeFaceSlabBox(nx, ny);
    if (state.r2x == nullptr || !(state.r2x_domain == domain)) {
        using amrex::FFT::Boundary;
        amrex::Array<std::pair<Boundary, Boundary>, AMREX_SPACEDIM> bc;
        bc[0] = {Boundary::odd, Boundary::odd};
        bc[1] = {Boundary::odd, Boundary::odd};
        bc[2] = {Boundary::periodic, Boundary::periodic};

        amrex::FFT::Info info{};
        info.setNumProcs(1);
        state.r2x = std::make_unique<amrex::FFT::R2X<amrex::Real>>(domain, bc,
                                                                   info);
        state.r2x_domain = domain;
        RegisterFFTCleanup();
    }

    return *state.r2x;
}

/**
 * Compute a one-component MultiFab dot product.
 */
amrex::Real
DotMF (amrex::MultiFab const& a, amrex::MultiFab const& b) {
    return amrex::MultiFab::Dot(a, 0, b, 0, 1, 0);
}

/**
 * Apply a real-valued mask to one face MultiFab.
 */
void
ApplyMaskMF (amrex::MultiFab const& in, amrex::MultiFab& out,
             amrex::MultiFab const& mask) {
    for (amrex::MFIter mfi(out, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& box = mfi.tilebox();
        amrex::Array4<amrex::Real const> const& in_arr = in.const_array(mfi);
        amrex::Array4<amrex::Real const> const& mask_arr = mask.const_array(mfi);
        amrex::Array4<amrex::Real> const& out_arr = out.array(mfi);
        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                out_arr(i, j, k) = mask_arr(i, j, k) * in_arr(i, j, k);
            });
    }
}

/**
 * Update the CG direction `p = r + beta*p`.
 */
void
UpdateCGDirectionMF (amrex::MultiFab& p, amrex::MultiFab const& r,
                     amrex::Real beta) {
    for (amrex::MFIter mfi(p, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& box = mfi.tilebox();
        amrex::Array4<amrex::Real> const& p_arr = p.array(mfi);
        amrex::Array4<amrex::Real const> const& r_arr = r.const_array(mfi);
        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                p_arr(i, j, k) = r_arr(i, j, k) + beta * p_arr(i, j, k);
            });
    }
}

/**
 * Build the Neumann mask and full-face Schur right-hand side.
 *
 * WarpX already applies the anode Dirichlet potential and homogeneous Neumann
 * condition on the remaining ceramic surface. With the zmin outward normal
 * n=-z, the accumulated wall charge density sets dphi/dn = sigma_s/epsilon0.
 * The Schur RHS is therefore just this target nonhomogeneous normal derivative.
 *
 * @param neumann_mask Output interior-face mask as 0/1 reals.
 * @param rhs Output full interior RHS with zero Dirichlet/anode entries.
 * @param sigma_s_mf Full-face surface charge density slab.
 * @param geom Level geometry.
 */
void
BuildMaskAndRhsMF (
    amrex::MultiFab& neumann_mask, amrex::MultiFab& rhs,
    amrex::MultiFab const& sigma_s_mf, amrex::Geometry const& geom) {
    amrex::Box domain = geom.Domain();
    domain.surroundingNodes();
    int const xlo = domain.smallEnd(0);
    int const ylo = domain.smallEnd(1);
    int const zlo = domain.smallEnd(2);
    amrex::Real const problo_x = geom.ProbLo(0);
    amrex::Real const problo_y = geom.ProbLo(1);
    amrex::Real const dx = geom.CellSize(0);
    amrex::Real const dy = geom.CellSize(1);

    auto const anode_config = ReadHallAnodeRingConfig(geom);

    amrex::Array4<amrex::Real const> sigma_arr;
    for (amrex::MFIter mfi(sigma_s_mf); mfi.isValid(); ++mfi) {
        sigma_arr = sigma_s_mf.const_array(mfi);
    }

    // Only interior transverse nodes enter the sine basis. Dirichlet/anode
    // nodes are represented by NeumannMask = false and explicit zeros in
    // RHS/u_face.
    amrex::Real const inv_epsilon0 =
        static_cast<amrex::Real>(1.0) / PhysConst::epsilon_0;

    for (amrex::MFIter mfi(rhs, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& box = mfi.tilebox();
        amrex::Array4<amrex::Real> const& mask_arr = neumann_mask.array(mfi);
        amrex::Array4<amrex::Real> const& rhs_arr = rhs.array(mfi);
        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::ignore_unused(k);
                int const gj = ylo + 1 + j;
                int const gi = xlo + 1 + i;
                bool const is_neumann =
                    !IsHallAnodeRingNode(gi, gj, zlo, zlo, xlo, ylo, problo_x,
                                          problo_y, dx, dy, anode_config);

                if (is_neumann) {
                    amrex::Real const sigma = sigma_arr(gi - xlo, gj - ylo, 0);
                    // With zmin outward normal n=-z, dphi/dn = sigma_s/epsilon0.
                    mask_arr(i, j, 0) = static_cast<amrex::Real>(1.0);
                    rhs_arr(i, j, 0) = sigma * inv_epsilon0;
                } else {
                    mask_arr(i, j, 0) = static_cast<amrex::Real>(0.0);
                    rhs_arr(i, j, 0) = static_cast<amrex::Real>(0.0);
                }
            });
    }
}

/**
 * Compute `sinh(k(Lz-z)) / sinh(kLz)` without large-argument overflow.
 *
 * @param k Transverse modal wave number.
 * @param z Distance from zmin.
 * @param lz Domain length in z.
 * @return Harmonic extension decay factor.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
SinhDecayRatio (amrex::Real k, amrex::Real z, amrex::Real lz) {
    amrex::Real const a = k * lz;
    amrex::Real const b = k * (lz - z);
    constexpr amrex::Real threshold = static_cast<amrex::Real>(40.0);
    if (a <= threshold) {
        amrex::Real const denom = std::sinh(a);
        return std::sinh(b) / denom;
    }

    amrex::Real const numerator_correction =
        static_cast<amrex::Real>(1.0) -
        std::exp(-static_cast<amrex::Real>(2.0) * b);
    amrex::Real const denominator_correction =
        static_cast<amrex::Real>(1.0) -
        std::exp(-static_cast<amrex::Real>(2.0) * a);
    return std::exp(-k * z) * numerator_correction / denominator_correction;
}

/**
 * Apply the masked Schur operator to a full-face MultiFab.
 *
 * @param v Full interior face vector with zero Dirichlet/anode entries.
 * @param av Output full interior face vector with zero Dirichlet/anode entries.
 * @param neumann_mask Full interior-face mask.
 * @param r2x AMReX real-to-real FFT plan.
 * @param geom Level geometry that defines physical lengths.
 */
void
ApplySchurOperatorMF (amrex::MultiFab const& v, amrex::MultiFab& av,
                      amrex::MultiFab const& neumann_mask,
                      amrex::FFT::R2X<amrex::Real>& r2x,
                      amrex::Geometry const& geom) {
    amrex::MultiFab face(v.boxArray(), v.DistributionMap(), 1, 0);
    ApplyMaskMF(v, face, neumann_mask);

    amrex::Real const lx = geom.ProbHi(0) - geom.ProbLo(0);
    amrex::Real const ly = geom.ProbHi(1) - geom.ProbLo(1);
    amrex::Real const lz = geom.ProbHi(2) - geom.ProbLo(2);
    amrex::Real const scale = r2x.scalingFactor();

    r2x.forwardThenBackward(
        face, av,
        [=] AMREX_GPU_DEVICE (int m, int n, int, auto& spectral_data)
        {
            amrex::Real const kx =
                static_cast<amrex::Real>(m + 1) * MathConst::pi / lx;
            amrex::Real const ky =
                static_cast<amrex::Real>(n + 1) * MathConst::pi / ly;
            amrex::Real const k = std::sqrt(kx * kx + ky * ky);
            spectral_data *= k * StableCoth(k * lz) * scale;
        });

    ApplyMaskMF(av, av, neumann_mask);
}

/**
 * Solve the masked full-face Schur system with unpreconditioned conjugate
 * gradients.
 *
 * @param x Output full interior boundary trace.
 * @param b Full interior right-hand side with zero Dirichlet/anode entries.
 * @param neumann_mask Full interior-face mask.
 * @param r2x AMReX real-to-real FFT plan.
 * @param geom Level geometry.
 * @param config Iteration limits and tolerances.
 */
void
SolveCGMF (amrex::MultiFab& x, amrex::MultiFab const& b,
           amrex::MultiFab const& neumann_mask,
           amrex::FFT::R2X<amrex::Real>& r2x, amrex::Geometry const& geom,
           SchurConfig const& config) {
    x.setVal(static_cast<amrex::Real>(0.0));

    amrex::MultiFab r(b.boxArray(), b.DistributionMap(), 1, 0);
    amrex::MultiFab p(b.boxArray(), b.DistributionMap(), 1, 0);
    amrex::MultiFab ap(b.boxArray(), b.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(r, b, 0, 0, 1, 0);
    amrex::MultiFab::Copy(p, b, 0, 0, 1, 0);

    amrex::Real rs_old = DotMF(r, r);
    amrex::Real const b_norm = std::sqrt(DotMF(b, b));
    amrex::Real const tolerance =
        std::max(config.abs_tol, config.rel_tol * b_norm);
    amrex::Real const tolerance_sq = tolerance * tolerance;

    if (rs_old <= tolerance_sq) {
        return;
    }

    for (int iter = 0; iter < config.max_iter; ++iter) {
        ApplySchurOperatorMF(p, ap, neumann_mask, r2x, geom);
        amrex::Real const denom = DotMF(p, ap);
        if (denom == static_cast<amrex::Real>(0.0)) {
            amrex::Abort("insert.schur_boundary CG encountered zero denominator");
        }

        amrex::Real const alpha = rs_old / denom;
        amrex::MultiFab::Saxpy(x, alpha, p, 0, 0, 1, 0);
        amrex::MultiFab::Saxpy(r, -alpha, ap, 0, 0, 1, 0);

        amrex::Real const rs_new = DotMF(r, r);
        if (rs_new <= tolerance_sq) {
            break;
        }

        amrex::Real const beta = rs_new / rs_old;
        UpdateCGDirectionMF(p, r, beta);
        rs_old = rs_new;
    }
}

/**
 * Reconstruct `phi_corr` on the volume using layered AMReX R2X transforms.
 *
 * @param geom Level geometry.
 * @param u_face Full boundary trace from the Schur solve.
 * @param r2x AMReX real-to-real FFT plan.
 * @param nx Number of x interior nodes/modes.
 * @param ny Number of y interior nodes/modes.
 */
void
ReconstructPhiCorrection (amrex::Geometry const& geom,
                          amrex::MultiFab const& u_face,
                          amrex::FFT::R2X<amrex::Real>& r2x, int nx, int ny) {
    // The correction MultiFab follows rho_fp level 0 layout so later field
    // updates can use WarpX-owned decomposition and guard-cell shape.
    auto& state = State();
    WarpX& warpx = WarpX::GetInstance();
    auto const& rho_fp = *warpx.m_fields.get(warpx::fields::FieldType::rho_fp, 0);

    state.phi_corr = std::make_unique<amrex::MultiFab>(
        rho_fp.boxArray(), rho_fp.DistributionMap(), 1, rho_fp.nGrowVect());
    state.phi_corr->setVal(static_cast<amrex::Real>(0.0));

    amrex::Box domain = geom.Domain();
    domain.surroundingNodes();
    int const xlo = domain.smallEnd(0);
    int const ylo = domain.smallEnd(1);
    int const zlo = domain.smallEnd(2);
    int const zhi = domain.bigEnd(2);
    int const nz_nodes = zhi - zlo + 1;
    amrex::Real const lz = geom.ProbHi(2) - geom.ProbLo(2);
    amrex::Real const dz = geom.CellSize(2);
    amrex::Real const lx = geom.ProbHi(0) - geom.ProbLo(0);
    amrex::Real const ly = geom.ProbHi(1) - geom.ProbLo(1);

    amrex::MultiFab layer(u_face.boxArray(), u_face.DistributionMap(), 1, 0);
    amrex::Real const scale = r2x.scalingFactor();

    for (int kk = 0; kk < nz_nodes; ++kk) {
        amrex::Real const z = static_cast<amrex::Real>(kk) * dz;
        amrex::Real const z_fraction = z / lz;
        int const mx = LayerModes(nx, z_fraction);
        int const my = LayerModes(ny, z_fraction);

        r2x.forwardThenBackward(
            u_face, layer,
            [=] AMREX_GPU_DEVICE (int m, int n, int, auto& spectral_data)
            {
                if (m >= mx || n >= my) {
                    spectral_data *= static_cast<amrex::Real>(0.0);
                    return;
                }
                amrex::Real const kx =
                    static_cast<amrex::Real>(m + 1) * MathConst::pi / lx;
                amrex::Real const ky =
                    static_cast<amrex::Real>(n + 1) * MathConst::pi / ly;
                amrex::Real const kmn = std::sqrt(kx * kx + ky * ky);
                amrex::Real const r = SinhDecayRatio(kmn, z, lz);
                spectral_data *= r * scale;
            });

        amrex::Array4<amrex::Real const> layer_arr;
        for (amrex::MFIter mfi(layer); mfi.isValid(); ++mfi) {
            layer_arr = layer.const_array(mfi);
        }

        // Scatter the reconstructed z layer into whatever boxes own that slab.
        amrex::Box const slab(
            amrex::IntVect(xlo + 1, ylo + 1, zlo + kk),
            amrex::IntVect(domain.bigEnd(0) - 1, domain.bigEnd(1) - 1, zlo + kk));
        for (amrex::MFIter mfi(*state.phi_corr, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            amrex::Box const write_box = mfi.validbox() & slab;
            if (!write_box.ok()) {
                continue;
            }
            amrex::Array4<amrex::Real> const& arr = state.phi_corr->array(mfi);
            amrex::ParallelFor(
                write_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    int const ii = i - xlo - 1;
                    int const jj = j - ylo - 1;
                    arr(i, j, k) = layer_arr(ii, jj, 0);
                });
        }
    }

    amrex::Gpu::streamSynchronize();
}

/**
 * Run the full MultiFab Schur solve and volume reconstruction.
 *
 * @param geom Level geometry.
 * @param sigma_s_mf Surface charge density on the full zmin face.
 */
void
SolveAndReconstructMF (amrex::Geometry const& geom,
                       amrex::MultiFab const& sigma_s_mf) {
    SchurConfig const& config = ReadConfig();

    amrex::Box const sigma_box = sigma_s_mf.boxArray().minimalBox();
    int const nx_face = sigma_box.length(0);
    int const ny_face = sigma_box.length(1);
    int const nx_internal = nx_face - 2;
    int const ny_internal = ny_face - 2;

    auto& state = State();
    amrex::Box const internal_box = MakeFaceSlabBox(nx_internal, ny_internal);
    if (state.neumann_mask == nullptr ||
        !(state.neumann_mask->boxArray().minimalBox() == internal_box))
    {
        state.neumann_mask = MakeFaceSlabMultiFab(nx_internal, ny_internal);
        state.u_face = MakeFaceSlabMultiFab(nx_internal, ny_internal);
    }

    state.nx_internal = nx_internal;
    state.ny_internal = ny_internal;
    amrex::FFT::R2X<amrex::Real>& r2x = PrepareSchurFFT(nx_internal,
                                                        ny_internal);

    amrex::MultiFab rhs(state.neumann_mask->boxArray(),
                        state.neumann_mask->DistributionMap(), 1, 0);
    BuildMaskAndRhsMF(*state.neumann_mask, rhs, sigma_s_mf, geom);
    SolveCGMF(*state.u_face, rhs, *state.neumann_mask, r2x, geom, config);
    ReconstructPhiCorrection(geom, *state.u_face, r2x, nx_internal,
                             ny_internal);
    amrex::Gpu::streamSynchronize();
}

#endif

void
AddPhiCorrectionToPotential (
    ablastr::fields::MultiLevelScalarField const& phi)
{
    auto const* const phi_corr = State().phi_corr.get();
    if (phi_corr == nullptr) {
        return;
    }

    auto& phi_mf = *phi[0];
    if (!phi_mf.boxArray().CellEqual(phi_corr->boxArray()) ||
        !(phi_mf.DistributionMap() == phi_corr->DistributionMap()))
    {
        amrex::Abort("insert.schur_boundary phi_corr layout does not match phi");
    }

    amrex::MultiFab::Saxpy(phi_mf, static_cast<amrex::Real>(1.0), *phi_corr,
                           0, 0, 1, 0);
}

} // namespace

namespace Insert {

/**
 * Check whether the Schur boundary correction is enabled.
 *
 * @return true when `insert.schur_boundary.enabled` is nonzero.
 */
bool
SpectralBoundarySchur::Enabled () {
    return ReadConfig().enabled;
}

/**
 * Apply the zmin Schur correction from accumulated wall-charge density.
 *
 * @param phi Multi-level nodal potential field.
 */
void
SpectralBoundarySchur::ApplyZMinCorrection (
    ablastr::fields::MultiLevelScalarField const& phi) {
    SchurConfig const& config = ReadConfig();
    if (!config.enabled) {
        return;
    }

#if defined(WARPX_DIM_3D)
    if (config.face != "zmin") {
        amrex::Abort("insert.schur_boundary only supports face = zmin");
    }
    if (phi.size() != 1) {
        amrex::Abort(
            "insert.schur_boundary only supports single-level level 0");
    }
    if (EB::enabled()) {
        amrex::Abort("insert.schur_boundary does not support EB");
    }
    if (amrex::ParallelDescriptor::NProcs() != 1) {
        amrex::Abort(
            "insert.schur_boundary requires one MPI rank");
    }

    WarpX& warpx = WarpX::GetInstance();
    ZMinWallChargeGrid const grid = MakeZMinWallChargeGrid(warpx.Geom(0));
    InitializeAccumulatedZMinWallChargeDensity(grid);

    auto& state = State();
    int const step = warpx.getistep(0);
    // Zmin wall charge is accumulated in AfterDiagnostics() after istep has
    // advanced.  The next field solve sees that same step value, so the Schur
    // update follows the Hall diagnostic cadence directly.
    if (state.phi_corr == nullptr ||
        step % BoundaryParticleDiagInterval() == 0)
    {
        SolveAndReconstruct(warpx.Geom(0), *g_accumulated_wall_charge_density);
    }

    AddPhiCorrectionToPotential(phi);
#else
    amrex::Abort("insert.schur_boundary requires WarpX_DIMS=3");
#endif
}

/**
 * Solve the zmin Schur boundary problem from an externally provided surface
 * charge.
 *
 * @param geom Level-0 geometry.
 * @param sigma_s_mf Surface charge density on the full zmin face.
 */
void
SpectralBoundarySchur::SolveAndReconstruct (
    amrex::Geometry const& geom,
    amrex::MultiFab const& sigma_s_mf) {
#if defined(WARPX_DIM_3D)
    if (amrex::ParallelDescriptor::NProcs() != 1) {
        amrex::Abort(
            "insert.schur_boundary requires one MPI rank");
    }
    if (sigma_s_mf.nComp() != 1) {
        amrex::Abort("insert.schur_boundary sigma_s_mf must have one component");
    }
    if (sigma_s_mf.boxArray().size() != 1) {
        amrex::Abort("insert.schur_boundary currently requires one wall-charge box");
    }

    amrex::Box domain = geom.Domain();
    domain.surroundingNodes();
    amrex::Box const sigma_box = sigma_s_mf.boxArray().minimalBox();
    amrex::Box const expected_box =
        amrex::Box(amrex::IntVect(0, 0, 0),
                   amrex::IntVect(domain.length(0) - 1,
                                  domain.length(1) - 1, 0));
    if (!(sigma_box == expected_box)) {
        amrex::Abort("insert.schur_boundary sigma_s_mf shape does not match "
                     "zmin face");
    }
    if (domain.length(0) < 3 || domain.length(1) < 3) {
        amrex::Abort("insert.schur_boundary requires at least one interior "
                     "transverse node");
    }

    SolveAndReconstructMF(geom, sigma_s_mf);
#else
    (void)geom;
    (void)sigma_s_mf;
    amrex::Abort("insert.schur_boundary requires WarpX_DIMS=3");
#endif
}

/**
 * Access the latest reconstructed correction potential.
 *
 * @return Pointer to `phi_corr`, or null if volume reconstruction was not run.
 */
amrex::MultiFab const*
SpectralBoundarySchur::PhiCorrection () {
    return State().phi_corr.get();
}

/**
 * Narrow Insert-level entry point used by the electrostatic solver.
 *
 * @param phi Multi-level nodal potential field.
 */
void
ApplyElectrostaticBoundaryCorrection (
    ablastr::fields::MultiLevelScalarField const& phi) {
    SpectralBoundarySchur::ApplyZMinCorrection(phi);
}

} // namespace Insert

#endif // HALL3D
