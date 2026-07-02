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

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

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
    amrex::Gpu::DeviceVector<int> d_neumann_mask;
    amrex::Gpu::DeviceVector<amrex::Real> d_u_face;
    amrex::Gpu::DeviceVector<amrex::Real> d_u_hat;
    int nx_internal = 0;
    int ny_internal = 0;
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
 * Convert an `(i,j)` index into the module's x-fast 1D layout.
 *
 * @param i x-like index.
 * @param j y-like index.
 * @param nx Extent in the x-like direction.
 * @return Flattened x-fast index.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
Index2D (int i, int j, int nx) {
    // x is the fast index for all dense face and mode arrays in this module.
    return i + nx * j;
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
 * Fill a device vector with a scalar.
 *
 * @param values Vector to resize and fill.
 * @param n Number of entries.
 * @param value Fill value.
 */
void
SetDeviceVector (amrex::Gpu::DeviceVector<amrex::Real>& values, std::size_t n,
    amrex::Real value) {
    values.resize(n);

    amrex::Real* values_ptr = values.dataPtr();
    amrex::ParallelFor(static_cast<int>(n), [=] AMREX_GPU_DEVICE(int i) {
        values_ptr[i] = value;
    });
}

/**
 * Build an orthonormal DST-I sine basis for interior nodes on device.
 *
 * @param basis Output x-fast dense basis matrix `basis[node, mode]`.
 * @param n Number of interior nodes and modes.
 */
void
BuildSineBasisDevice (amrex::Gpu::DeviceVector<amrex::Real>& basis, int n) {
    basis.resize(static_cast<std::size_t>(n) * n);
    amrex::Real* basis_ptr = basis.dataPtr();
    amrex::Real const norm = std::sqrt(static_cast<amrex::Real>(2.0) /
                                       static_cast<amrex::Real>(n + 1));

    amrex::Box const box(amrex::IntVect(0, 0, 0),
                         amrex::IntVect(n - 1, n - 1, 0));
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int m, int) {
        amrex::Real const angle = MathConst::pi *
                                  static_cast<amrex::Real>((i + 1) * (m + 1)) /
                                  static_cast<amrex::Real>(n + 1);
        basis_ptr[Index2D(i, m, n)] = norm * std::sin(angle);
    });
}

/**
 * Apply the 2D orthonormal DST-I transform on device.
 *
 * @param face Input physical-space face values, x-fast.
 * @param face_hat Output spectral coefficients, x-mode-fast.
 * @param sx x-direction sine basis.
 * @param sy y-direction sine basis.
 * @param tmp Workspace with `nx*ny` entries.
 * @param nx Number of x interior nodes/modes.
 * @param ny Number of y interior nodes/modes.
 */
void
ApplyDST2Device (amrex::Gpu::DeviceVector<amrex::Real> const& face,
                 amrex::Gpu::DeviceVector<amrex::Real>& face_hat,
                 amrex::Gpu::DeviceVector<amrex::Real> const& sx,
                 amrex::Gpu::DeviceVector<amrex::Real> const& sy,
                 amrex::Gpu::DeviceVector<amrex::Real>& tmp, int nx, int ny) {
    std::size_t const n_total = static_cast<std::size_t>(nx) * ny;
    tmp.resize(n_total);
    face_hat.resize(n_total);

    amrex::Real const* face_ptr = face.dataPtr();
    amrex::Real const* sx_ptr = sx.dataPtr();
    amrex::Real const* sy_ptr = sy.dataPtr();
    amrex::Real* tmp_ptr = tmp.dataPtr();
    amrex::Real* face_hat_ptr = face_hat.dataPtr();

    amrex::Box const box(amrex::IntVect(0, 0, 0),
                         amrex::IntVect(nx - 1, ny - 1, 0));

    // face_hat = Sx^T * face * Sy, implemented as two dense matrix
    // contractions.
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int m, int j, int) {
        amrex::Real sum = static_cast<amrex::Real>(0.0);
        for (int i = 0; i < nx; ++i) {
            sum += sx_ptr[Index2D(i, m, nx)] * face_ptr[Index2D(i, j, nx)];
        }
        tmp_ptr[Index2D(m, j, nx)] = sum;
    });

    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int m, int n, int) {
        amrex::Real sum = static_cast<amrex::Real>(0.0);
        for (int j = 0; j < ny; ++j) {
            sum += tmp_ptr[Index2D(m, j, nx)] * sy_ptr[Index2D(j, n, ny)];
        }
        face_hat_ptr[Index2D(m, n, nx)] = sum;
    });
}

/**
 * Apply the inverse 2D orthonormal DST-I transform on device.
 *
 * @param face_hat Input spectral coefficients, x-mode-fast.
 * @param face Output physical-space face values, x-fast.
 * @param sx x-direction sine basis.
 * @param sy y-direction sine basis.
 * @param tmp Workspace with `nx*ny` entries.
 * @param nx Number of x interior nodes/modes.
 * @param ny Number of y interior nodes/modes.
 */
void
ApplyIDST2Device (amrex::Gpu::DeviceVector<amrex::Real> const& face_hat,
                  amrex::Gpu::DeviceVector<amrex::Real>& face,
                  amrex::Gpu::DeviceVector<amrex::Real> const& sx,
                  amrex::Gpu::DeviceVector<amrex::Real> const& sy,
                  amrex::Gpu::DeviceVector<amrex::Real>& tmp, int nx, int ny) {
    std::size_t const n_total = static_cast<std::size_t>(nx) * ny;
    tmp.resize(n_total);
    face.resize(n_total);

    amrex::Real const* face_hat_ptr = face_hat.dataPtr();
    amrex::Real const* sx_ptr = sx.dataPtr();
    amrex::Real const* sy_ptr = sy.dataPtr();
    amrex::Real* tmp_ptr = tmp.dataPtr();
    amrex::Real* face_ptr = face.dataPtr();

    amrex::Box const box(amrex::IntVect(0, 0, 0),
                         amrex::IntVect(nx - 1, ny - 1, 0));

    // face = Sx * face_hat * Sy^T. With the normalized basis this is the
    // inverse.
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int n, int) {
        amrex::Real sum = static_cast<amrex::Real>(0.0);
        for (int m = 0; m < nx; ++m) {
            sum += sx_ptr[Index2D(i, m, nx)] * face_hat_ptr[Index2D(m, n, nx)];
        }
        tmp_ptr[Index2D(i, n, nx)] = sum;
    });

    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int) {
        amrex::Real sum = static_cast<amrex::Real>(0.0);
        for (int n = 0; n < ny; ++n) {
            sum += tmp_ptr[Index2D(i, n, nx)] * sy_ptr[Index2D(j, n, ny)];
        }
        face_ptr[Index2D(i, j, nx)] = sum;
    });
}

/**
 * Compute a device dot product.
 *
 * @param a First vector.
 * @param b Second vector.
 * @return Sum of elementwise products.
 */
amrex::Real
DotDevice (amrex::Gpu::DeviceVector<amrex::Real> const& a,
           amrex::Gpu::DeviceVector<amrex::Real> const& b) {
    amrex::Real const* a_ptr = a.dataPtr();
    amrex::Real const* b_ptr = b.dataPtr();
    return amrex::Reduce::Sum<amrex::Real>(
        static_cast<int>(a.size()),
        [=] AMREX_GPU_DEVICE(int i) -> amrex::Real {
            return a_ptr[i] * b_ptr[i];
        });
}

/**
 * Add a scaled vector into another device vector in place.
 *
 * @param y Vector updated as `y += alpha * x`.
 * @param x Input vector.
 * @param alpha Scale factor.
 */
void
AddScaledDevice (amrex::Gpu::DeviceVector<amrex::Real>& y,
                 amrex::Gpu::DeviceVector<amrex::Real> const& x,
                 amrex::Real alpha) {
    amrex::Real* y_ptr = y.dataPtr();
    amrex::Real const* x_ptr = x.dataPtr();
    amrex::ParallelFor(
        static_cast<int>(y.size()),
        [=] AMREX_GPU_DEVICE(int i) { y_ptr[i] += alpha * x_ptr[i]; });
}

/**
 * Update the CG search direction on device.
 *
 * @param p Search direction updated as `p = r + beta*p`.
 * @param r Residual vector.
 * @param beta CG beta coefficient.
 */
void
UpdateCGDirectionDevice (amrex::Gpu::DeviceVector<amrex::Real>& p,
                         amrex::Gpu::DeviceVector<amrex::Real> const& r,
                         amrex::Real beta) {
    amrex::Real* p_ptr = p.dataPtr();
    amrex::Real const* r_ptr = r.dataPtr();
    amrex::ParallelFor(
        static_cast<int>(p.size()),
        [=] AMREX_GPU_DEVICE(int i) { p_ptr[i] = r_ptr[i] + beta * p_ptr[i]; });
}

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
 * Build the spectral Dirichlet-to-Neumann eigenvalue table on device.
 *
 * @param lambda Output `lambda[m,n] = k_mn coth(k_mn Lz)`.
 * @param nx Number of x modes.
 * @param ny Number of y modes.
 * @param geom Level geometry that defines physical lengths.
 */
void
BuildLambdaDevice (amrex::Gpu::DeviceVector<amrex::Real>& lambda, int nx,
                   int ny, amrex::Geometry const& geom) {
    lambda.resize(static_cast<std::size_t>(nx) * ny);
    amrex::Real const lx = geom.ProbHi(0) - geom.ProbLo(0);
    amrex::Real const ly = geom.ProbHi(1) - geom.ProbLo(1);
    amrex::Real const lz = geom.ProbHi(2) - geom.ProbLo(2);
    amrex::Real* lambda_ptr = lambda.dataPtr();

    amrex::Box const box(amrex::IntVect(0, 0, 0),
                         amrex::IntVect(nx - 1, ny - 1, 0));
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int m, int n, int) {
        amrex::Real const kx =
            static_cast<amrex::Real>(m + 1) * MathConst::pi / lx;
        amrex::Real const ky =
            static_cast<amrex::Real>(n + 1) * MathConst::pi / ly;
        amrex::Real const k = std::sqrt(kx * kx + ky * ky);
        lambda_ptr[Index2D(m, n, nx)] = k * StableCoth(k * lz);
    });
}

/**
 * Apply the masked Schur operator to a full-face device vector.
 *
 * @param v Full interior face vector with zero Dirichlet/anode entries.
 * @param av Output full interior face vector with zero Dirichlet/anode entries.
 * @param neumann_mask Full interior-face mask.
 * @param lambda Spectral DtN eigenvalues.
 * @param sx x-direction sine basis.
 * @param sy y-direction sine basis.
 * @param face Workspace for masked physical-space data.
 * @param face_hat Workspace for spectral data.
 * @param tmp Workspace for dense matrix contractions.
 * @param nx Number of x interior nodes/modes.
 * @param ny Number of y interior nodes/modes.
 */
void
ApplySchurOperatorDevice (amrex::Gpu::DeviceVector<amrex::Real> const& v,
                          amrex::Gpu::DeviceVector<amrex::Real>& av,
                          amrex::Gpu::DeviceVector<int> const& neumann_mask,
                          amrex::Gpu::DeviceVector<amrex::Real> const& lambda,
                          amrex::Gpu::DeviceVector<amrex::Real> const& sx,
                          amrex::Gpu::DeviceVector<amrex::Real> const& sy,
                          amrex::Gpu::DeviceVector<amrex::Real>& face,
                          amrex::Gpu::DeviceVector<amrex::Real>& face_hat,
                          amrex::Gpu::DeviceVector<amrex::Real>& tmp, int nx,
                          int ny) {
    int const n_total = nx * ny;
    std::size_t const n_total_size = static_cast<std::size_t>(n_total);
    face.resize(n_total_size);
    av.resize(n_total_size);
    amrex::Real const* v_ptr = v.dataPtr();
    int const* mask_ptr = neumann_mask.dataPtr();
    amrex::Real* face_ptr = face.dataPtr();

    // Apply P_N Lambda P_N^T on the full face. Dirichlet/anode entries are kept
    // as explicit zeros, avoiding GPU-side variable-length compaction.
    amrex::ParallelFor(n_total, [=] AMREX_GPU_DEVICE(int i) {
        face_ptr[i] =
            mask_ptr[i] != 0 ? v_ptr[i] : static_cast<amrex::Real>(0.0);
    });

    ApplyDST2Device(face, face_hat, sx, sy, tmp, nx, ny);

    amrex::Real* face_hat_ptr = face_hat.dataPtr();
    amrex::Real const* lambda_ptr = lambda.dataPtr();
    amrex::ParallelFor(n_total, [=] AMREX_GPU_DEVICE(int i) {
        face_hat_ptr[i] *= lambda_ptr[i];
    });

    ApplyIDST2Device(face_hat, face, sx, sy, tmp, nx, ny);

    face_ptr = face.dataPtr();
    amrex::Real* av_ptr = av.dataPtr();
    amrex::ParallelFor(n_total, [=] AMREX_GPU_DEVICE(int i) {
        av_ptr[i] =
            mask_ptr[i] != 0 ? face_ptr[i] : static_cast<amrex::Real>(0.0);
    });
}

/**
 * Solve the masked full-face Schur system with unpreconditioned conjugate
 * gradients.
 *
 * @param x Output full interior boundary trace.
 * @param b Full interior right-hand side with zero Dirichlet/anode entries.
 * @param neumann_mask Full interior-face mask.
 * @param lambda Spectral DtN eigenvalues.
 * @param sx x-direction sine basis.
 * @param sy y-direction sine basis.
 * @param nx Number of x interior nodes/modes.
 * @param ny Number of y interior nodes/modes.
 * @param config Iteration limits and tolerances.
 */
void
SolveCGDevice (amrex::Gpu::DeviceVector<amrex::Real>& x,
               amrex::Gpu::DeviceVector<amrex::Real> const& b,
               amrex::Gpu::DeviceVector<int> const& neumann_mask,
               amrex::Gpu::DeviceVector<amrex::Real> const& lambda,
               amrex::Gpu::DeviceVector<amrex::Real> const& sx,
               amrex::Gpu::DeviceVector<amrex::Real> const& sy, int nx, int ny,
               SchurConfig const& config) {
    SetDeviceVector(x, b.size(), static_cast<amrex::Real>(0.0));

    amrex::Gpu::DeviceVector<amrex::Real> r(b.size());
    amrex::Gpu::DeviceVector<amrex::Real> p(b.size());
    amrex::Gpu::DeviceVector<amrex::Real> ap(b.size());
    amrex::Gpu::DeviceVector<amrex::Real> face;
    amrex::Gpu::DeviceVector<amrex::Real> face_hat;
    amrex::Gpu::DeviceVector<amrex::Real> tmp;
    amrex::Gpu::copy(amrex::Gpu::deviceToDevice, b.begin(), b.end(), r.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToDevice, b.begin(), b.end(), p.begin());

    amrex::Real rs_old = DotDevice(r, r);
    amrex::Real const b_norm = std::sqrt(DotDevice(b, b));
    amrex::Real const tolerance =
        std::max(config.abs_tol, config.rel_tol * b_norm);
    amrex::Real const tolerance_sq = tolerance * tolerance;

    if (rs_old <= tolerance_sq) {
        amrex::Gpu::streamSynchronize();
        return;
    }

    for (int iter = 0; iter < config.max_iter; ++iter) {
        ApplySchurOperatorDevice(p, ap, neumann_mask, lambda, sx, sy, face,
                                 face_hat, tmp, nx, ny);
        amrex::Real const denom = DotDevice(p, ap);
        amrex::Real const alpha = rs_old / denom;
        AddScaledDevice(x, p, alpha);
        AddScaledDevice(r, ap, -alpha);

        amrex::Real const rs_new = DotDevice(r, r);
        if (rs_new <= tolerance_sq) {
            break;
        }

        amrex::Real const beta = rs_new / rs_old;
        UpdateCGDirectionDevice(p, r, beta);
        rs_old = rs_new;
    }
    amrex::Gpu::streamSynchronize();
}

/**
 * Build the Neumann mask and full-face Schur right-hand side on device.
 *
 * WarpX already applies the anode Dirichlet potential and homogeneous Neumann
 * condition on the remaining ceramic surface. With the zmin outward normal
 * n=-z, the accumulated wall charge density sets dphi/dn = sigma_s/epsilon0.
 * The Schur RHS is therefore just this target nonhomogeneous normal derivative.
 *
 * @param neumann_mask Output interior-face mask.
 * @param rhs Output full interior RHS with zero Dirichlet/anode entries.
 * @param sigma_s_device Device full-face surface charge density.
 * @param geom Level geometry.
 * @param nx_internal Number of x interior nodes.
 * @param ny_internal Number of y interior nodes.
 * @param nx_face Full zmin face node count in x.
 */
void
BuildMaskAndRhsDevice (
    amrex::Gpu::DeviceVector<int>& neumann_mask,
    amrex::Gpu::DeviceVector<amrex::Real>& rhs,
    amrex::Gpu::DeviceVector<amrex::Real> const& sigma_s_device,
    amrex::Geometry const& geom, int nx_internal, int ny_internal,
    int nx_face) {
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

    std::size_t const n_total =
        static_cast<std::size_t>(nx_internal) * ny_internal;
    neumann_mask.resize(n_total);
    rhs.resize(n_total);

    // Only interior transverse nodes enter the sine basis. Dirichlet/anode
    // nodes are represented by NeumannMask = false and explicit zeros in
    // RHS/u_face.
    int* mask_ptr = neumann_mask.dataPtr();
    amrex::Real* rhs_ptr = rhs.dataPtr();
    amrex::Real const* sigma_ptr = sigma_s_device.dataPtr();
    amrex::Real const inv_epsilon0 =
        static_cast<amrex::Real>(1.0) / PhysConst::epsilon_0;

    amrex::Box const box(amrex::IntVect(0, 0, 0),
                         amrex::IntVect(nx_internal - 1, ny_internal - 1, 0));
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int) {
        int const gj = ylo + 1 + j;
        int const gi = xlo + 1 + i;
        bool is_neumann = true;
        if (IsHallAnodeRingNode(gi, gj, zlo, zlo, xlo, ylo, problo_x, problo_y,
                                dx, dy, anode_config)) {
            is_neumann = false;
        }

        int const idx = Index2D(i, j, nx_internal);
        if (is_neumann) {
            amrex::Real const sigma =
                sigma_ptr[Index2D(gi - xlo, gj - ylo, nx_face)];
            // With zmin outward normal n=-z, dphi/dn = sigma_s/epsilon0.
            mask_ptr[idx] = 1;
            rhs_ptr[idx] = sigma * inv_epsilon0;
        } else {
            mask_ptr[idx] = 0;
            rhs_ptr[idx] = static_cast<amrex::Real>(0.0);
        }
    });
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
 * Reconstruct `phi_corr` on the volume using layered GPU matrix contractions.
 *
 * @param geom Level geometry.
 * @param u_hat Full boundary spectrum from the solved boundary trace.
 * @param sx x-direction sine basis.
 * @param sy y-direction sine basis.
 * @param nx Number of x interior nodes/modes.
 * @param ny Number of y interior nodes/modes.
 */
void
ReconstructPhiCorrection (amrex::Geometry const& geom,
                          amrex::Gpu::DeviceVector<amrex::Real> const& u_hat,
                          amrex::Gpu::DeviceVector<amrex::Real> const& sx,
                          amrex::Gpu::DeviceVector<amrex::Real> const& sy,
                          int nx, int ny) {
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

    amrex::Real const* u_hat_ptr = u_hat.dataPtr();
    amrex::Real const* sx_ptr = sx.dataPtr();
    amrex::Real const* sy_ptr = sy.dataPtr();

    amrex::Gpu::DeviceVector<amrex::Real> d_a(static_cast<std::size_t>(nx) * ny,
                                              static_cast<amrex::Real>(0.0));
    amrex::Gpu::DeviceVector<amrex::Real> d_tmp(
        static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));
    amrex::Gpu::DeviceVector<amrex::Real> d_layer(
        static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));
    amrex::Real* a_ptr = d_a.dataPtr();
    amrex::Real* tmp_ptr = d_tmp.dataPtr();
    amrex::Real* layer_ptr = d_layer.dataPtr();

    for (int kk = 0; kk < nz_nodes; ++kk) {
        amrex::Real const z = static_cast<amrex::Real>(kk) * dz;
        amrex::Real const z_fraction = z / lz;
        int const mx = LayerModes(nx, z_fraction);
        int const my = LayerModes(ny, z_fraction);

        // Step 1: build the retained spectral coefficients A(m,n,z).
        amrex::Box const spectral_box(amrex::IntVect(0, 0, 0),
                                      amrex::IntVect(mx - 1, my - 1, 0));
        amrex::ParallelFor(
            spectral_box, [=] AMREX_GPU_DEVICE(int m, int n, int) {
                amrex::Real const kx =
                    static_cast<amrex::Real>(m + 1) * MathConst::pi / lx;
                amrex::Real const ky =
                    static_cast<amrex::Real>(n + 1) * MathConst::pi / ly;
                amrex::Real const kmn = std::sqrt(kx * kx + ky * ky);
                amrex::Real const r = SinhDecayRatio(kmn, z, lz);
                a_ptr[Index2D(m, n, nx)] = r * u_hat_ptr[Index2D(m, n, nx)];
            });

        // Step 2: tmp(i,n,z) = sum_m Sx(i,m) A(m,n,z).
        amrex::Box const tmp_box(amrex::IntVect(0, 0, 0),
                                 amrex::IntVect(nx - 1, my - 1, 0));
        amrex::ParallelFor(tmp_box, [=] AMREX_GPU_DEVICE(int i, int n, int) {
            amrex::Real sum = static_cast<amrex::Real>(0.0);
            for (int m = 0; m < mx; ++m) {
                sum += sx_ptr[Index2D(i, m, nx)] * a_ptr[Index2D(m, n, nx)];
            }
            tmp_ptr[Index2D(i, n, nx)] = sum;
        });

        // Step 3: layer(i,j,z) = sum_n tmp(i,n,z) Sy(j,n).
        amrex::Box const layer_box(amrex::IntVect(0, 0, 0),
                                   amrex::IntVect(nx - 1, ny - 1, 0));
        amrex::ParallelFor(layer_box, [=] AMREX_GPU_DEVICE(int i, int j, int) {
            amrex::Real sum = static_cast<amrex::Real>(0.0);
            for (int n = 0; n < my; ++n) {
                sum += tmp_ptr[Index2D(i, n, nx)] * sy_ptr[Index2D(j, n, ny)];
            }
            layer_ptr[Index2D(i, j, nx)] = sum;
        });

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
                    arr(i, j, k) = layer_ptr[Index2D(ii, jj, nx)];
                });
        }
    }

    amrex::Gpu::streamSynchronize();
}

/**
 * Run the full device-side Schur solve and optional device-side reconstruction.
 *
 * @param geom Level geometry.
 * @param sigma_s_device Device surface charge density in full-face x-fast
 * layout.
 * @param nx_face Full zmin face node count in x.
 * @param ny_face Full zmin face node count in y.
 */
void
SolveAndReconstructDevice (
    amrex::Geometry const& geom,
    amrex::Gpu::DeviceVector<amrex::Real> const& sigma_s_device, int nx_face,
    int ny_face) {
    SchurConfig const& config = ReadConfig();

    int const nx_internal = nx_face - 2;
    int const ny_internal = ny_face - 2;

    amrex::Gpu::DeviceVector<amrex::Real> sx;
    amrex::Gpu::DeviceVector<amrex::Real> sy;
    amrex::Gpu::DeviceVector<amrex::Real> lambda;
    amrex::Gpu::DeviceVector<amrex::Real> rhs;
    amrex::Gpu::DeviceVector<amrex::Real> tmp;
    BuildSineBasisDevice(sx, nx_internal);
    BuildSineBasisDevice(sy, ny_internal);
    BuildLambdaDevice(lambda, nx_internal, ny_internal, geom);

    auto& state = State();
    state.nx_internal = nx_internal;
    state.ny_internal = ny_internal;
    BuildMaskAndRhsDevice(state.d_neumann_mask, rhs, sigma_s_device, geom,
                          nx_internal, ny_internal, nx_face);

    SolveCGDevice(state.d_u_face, rhs, state.d_neumann_mask, lambda, sx, sy,
                  nx_internal, ny_internal, config);

    ApplyDST2Device(state.d_u_face, state.d_u_hat, sx, sy, tmp, nx_internal,
                    ny_internal);
    ReconstructPhiCorrection(geom, state.d_u_hat, sx, sy, nx_internal,
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
    amrex::Box domain = warpx.Geom(0).Domain();
    domain.surroundingNodes();
    int const nx_face = domain.length(0);
    int const ny_face = domain.length(1);

    ZMinWallChargeGrid const grid = MakeZMinWallChargeGrid(warpx.Geom(0));
    InitializeAccumulatedZMinWallChargeDensity(ZMinWallChargeSize(grid));

    auto& state = State();
    int const step = warpx.getistep(0);
    // Zmin wall charge is accumulated in AfterDiagnostics() after istep has
    // advanced.  The next field solve sees that same step value, so the Schur
    // update follows the Hall diagnostic cadence directly.
    if (state.phi_corr == nullptr ||
        step % BoundaryParticleDiagInterval() == 0)
    {
        SolveAndReconstruct(warpx.Geom(0), g_accumulated_wall_charge_density,
                            nx_face, ny_face);
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
 * @param sigma_s_device Surface charge density in full-face x-fast layout.
 * @param nx_face Full zmin face node count in x.
 * @param ny_face Full zmin face node count in y.
 */
void
SpectralBoundarySchur::SolveAndReconstruct (
    amrex::Geometry const& geom,
    amrex::Gpu::DeviceVector<amrex::Real> const& sigma_s_device, int nx_face,
    int ny_face) {
#if defined(WARPX_DIM_3D)
    if (amrex::ParallelDescriptor::NProcs() != 1) {
        amrex::Abort(
            "insert.schur_boundary requires one MPI rank");
    }
    if (sigma_s_device.size() != static_cast<std::size_t>(nx_face) * ny_face) {
        amrex::Abort("insert.schur_boundary sigma_s_device size does not match "
                     "face size");
    }

    SolveAndReconstructDevice(geom, sigma_s_device, nx_face, ny_face);
#else
    (void)geom;
    (void)sigma_s_device;
    (void)nx_face;
    (void)ny_face;
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
