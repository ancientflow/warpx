#include "SpectralBoundarySchur.h"

#include "EmbeddedBoundary/Enabled.H"
#include "Fields.H"
#include "Insert/Boundary/ZMinWallCharge.h"
#include "Insert/Config/WarpXSimulationConfig.h"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {

#ifdef HALL3D
struct HallAnodeRingConfig
{
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
ReadHallAnodeRingConfig (amrex::Geometry const& geom)
{
    amrex::ParmParse pp_mc("my_constants");
    auto l_factor = static_cast<amrex::Real>(1.0);
    pp_mc.get("l_factor", l_factor);

    amrex::Real const r_min =
        static_cast<amrex::Real>(0.021) / (static_cast<amrex::Real>(2.0) * l_factor);
    amrex::Real const r_max =
        static_cast<amrex::Real>(0.031) / (static_cast<amrex::Real>(2.0) * l_factor);

    return HallAnodeRingConfig{
        static_cast<amrex::Real>(0.5) * (geom.ProbLo(0) + geom.ProbHi(0)),
        static_cast<amrex::Real>(0.5) * (geom.ProbLo(1) + geom.ProbHi(1)),
        r_min * r_min,
        r_max * r_max
    };
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
bool
IsHallAnodeRingNode (
    int i, int j, int k,
    int zlo, int xlo, int ylo,
    amrex::Real problo_x, amrex::Real problo_y,
    amrex::Real dx, amrex::Real dy,
    HallAnodeRingConfig const& config)
{
    amrex::Real const x = problo_x + (i - xlo) * dx;
    amrex::Real const y = problo_y + (j - ylo) * dy;
    amrex::Real const r_sq =
        (x - config.center_x) * (x - config.center_x) +
        (y - config.center_y) * (y - config.center_y);

    return k == zlo && r_sq >= config.r_min_sq && r_sq <= config.r_max_sq;
}
#endif

bool
UseZMinWallChargeSource (std::string const& source)
{
    return source == "zmin_wall_charge" ||
           source == "zmin_wall_charge_density" ||
           source == "zmin_wall_charge_deposit";
}

bool
UseConstantSigmaSource (std::string const& source)
{
    return source == "sigma_s_device_vector" ||
           source == "constant" ||
           source == "constant_sigma_s";
}

struct SchurConfig
{
    bool enabled = false;
    bool rebuild_volume_field = false;
    int max_iter = 100;
    amrex::Real rel_tol = static_cast<amrex::Real>(1.0e-10);
    amrex::Real abs_tol = static_cast<amrex::Real>(0.0);
    amrex::Real sigma_s = static_cast<amrex::Real>(0.0);
    std::string face = "zmin";
    std::string backend = "cpu_serial";
    std::string source = "sigma_s_device_vector";
};

struct SchurState
{
    bool config_read = false;
    SchurConfig config;
    std::unique_ptr<amrex::MultiFab> phi_corr;
    // Boundary solve products are kept on host for the current cpu_serial backend.
    std::vector<amrex::Real> u_N;
    std::vector<amrex::Real> u_face;
    std::vector<amrex::Real> u_hat;
    int nx_internal = 0;
    int ny_internal = 0;
};

/**
 * Access the module-local persistent state.
 *
 * @return Mutable Schur state singleton.
 */
SchurState&
State ()
{
    static SchurState state;
    return state;
}

/**
 * Read Schur runtime parameters once from `insert.schur_boundary`.
 *
 * @return Cached Schur configuration.
 */
SchurConfig const&
ReadConfig ()
{
    auto& state = State();
    if (!state.config_read) {
        amrex::ParmParse pp("insert.schur_boundary");
        pp.query("enabled", state.config.enabled);
        pp.query("rebuild_volume_field", state.config.rebuild_volume_field);
        pp.query("max_iter", state.config.max_iter);
        pp.query("rel_tol", state.config.rel_tol);
        pp.query("abs_tol", state.config.abs_tol);
        pp.query("sigma_s", state.config.sigma_s);
        pp.query("face", state.config.face);
        pp.query("backend", state.config.backend);
        pp.query("source", state.config.source);
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int
Index2D (int i, int j, int nx)
{
    // x is the fast index for all dense face and mode arrays in this module.
    return i + nx * j;
}

/**
 * Evaluate `coth(x)` without overflowing for large positive `x`.
 *
 * @param x Positive argument.
 * @return Stable `coth(x)`.
 */
amrex::Real
StableCoth (amrex::Real x)
{
    constexpr amrex::Real threshold = static_cast<amrex::Real>(40.0);
    if (x > threshold) {
        return static_cast<amrex::Real>(1.0);
    }
    amrex::Real const exp_neg_2x = std::exp(-static_cast<amrex::Real>(2.0) * x);
    return (static_cast<amrex::Real>(1.0) + exp_neg_2x) /
           (static_cast<amrex::Real>(1.0) - exp_neg_2x);
}

/**
 * Build an orthonormal DST-I sine basis for interior nodes.
 *
 * @param basis Output x-fast dense basis matrix `basis[node, mode]`.
 * @param n Number of interior nodes and modes.
 */
void
BuildSineBasis (
    std::vector<amrex::Real>& basis,
    int n)
{
    // Orthonormal DST-I basis for interior nodes; forward and inverse are transposes.
    basis.assign(static_cast<std::size_t>(n) * n, static_cast<amrex::Real>(0.0));
    amrex::Real const norm = std::sqrt(
        static_cast<amrex::Real>(2.0) / static_cast<amrex::Real>(n + 1));

    for (int i = 0; i < n; ++i) {
        for (int m = 0; m < n; ++m) {
            amrex::Real const angle =
                MathConst::pi * static_cast<amrex::Real>((i + 1) * (m + 1)) /
                static_cast<amrex::Real>(n + 1);
            basis[Index2D(i, m, n)] = norm * std::sin(angle);
        }
    }
}

/**
 * Apply the 2D orthonormal DST-I transform on a face array.
 *
 * @param face Input physical-space face values, x-fast.
 * @param face_hat Output spectral coefficients, x-mode-fast.
 * @param sx x-direction sine basis.
 * @param sy y-direction sine basis.
 * @param nx Number of x interior nodes/modes.
 * @param ny Number of y interior nodes/modes.
 */
void
ApplyDST2 (
    std::vector<amrex::Real> const& face,
    std::vector<amrex::Real>& face_hat,
    std::vector<amrex::Real> const& sx,
    std::vector<amrex::Real> const& sy,
    int nx,
    int ny)
{
    // face_hat = Sx^T * face * Sy, implemented as two dense matrix contractions.
    std::vector<amrex::Real> tmp(
        static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));
    face_hat.assign(static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));

    for (int m = 0; m < nx; ++m) {
        for (int j = 0; j < ny; ++j) {
            amrex::Real sum = static_cast<amrex::Real>(0.0);
            for (int i = 0; i < nx; ++i) {
                sum += sx[Index2D(i, m, nx)] * face[Index2D(i, j, nx)];
            }
            tmp[Index2D(m, j, nx)] = sum;
        }
    }

    for (int m = 0; m < nx; ++m) {
        for (int n = 0; n < ny; ++n) {
            amrex::Real sum = static_cast<amrex::Real>(0.0);
            for (int j = 0; j < ny; ++j) {
                sum += tmp[Index2D(m, j, nx)] * sy[Index2D(j, n, ny)];
            }
            face_hat[Index2D(m, n, nx)] = sum;
        }
    }
}

/**
 * Apply the inverse 2D orthonormal DST-I transform on a face array.
 *
 * @param face_hat Input spectral coefficients, x-mode-fast.
 * @param face Output physical-space face values, x-fast.
 * @param sx x-direction sine basis.
 * @param sy y-direction sine basis.
 * @param nx Number of x interior nodes/modes.
 * @param ny Number of y interior nodes/modes.
 */
void
ApplyIDST2 (
    std::vector<amrex::Real> const& face_hat,
    std::vector<amrex::Real>& face,
    std::vector<amrex::Real> const& sx,
    std::vector<amrex::Real> const& sy,
    int nx,
    int ny)
{
    // face = Sx * face_hat * Sy^T. With the normalized basis this is the inverse.
    std::vector<amrex::Real> tmp(
        static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));
    face.assign(static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));

    for (int i = 0; i < nx; ++i) {
        for (int n = 0; n < ny; ++n) {
            amrex::Real sum = static_cast<amrex::Real>(0.0);
            for (int m = 0; m < nx; ++m) {
                sum += sx[Index2D(i, m, nx)] * face_hat[Index2D(m, n, nx)];
            }
            tmp[Index2D(i, n, nx)] = sum;
        }
    }

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            amrex::Real sum = static_cast<amrex::Real>(0.0);
            for (int n = 0; n < ny; ++n) {
                sum += tmp[Index2D(i, n, nx)] * sy[Index2D(j, n, ny)];
            }
            face[Index2D(i, j, nx)] = sum;
        }
    }
}

/**
 * Compute a host dot product.
 *
 * @param a First vector.
 * @param b Second vector.
 * @return Sum of elementwise products.
 */
amrex::Real
Dot (
    std::vector<amrex::Real> const& a,
    std::vector<amrex::Real> const& b)
{
    AMREX_ALWAYS_ASSERT(a.size() == b.size());
    amrex::Real result = static_cast<amrex::Real>(0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

/**
 * Add a scaled vector into another vector in place.
 *
 * @param y Vector updated as `y += alpha * x`.
 * @param x Input vector.
 * @param alpha Scale factor.
 */
void
AddScaled (
    std::vector<amrex::Real>& y,
    std::vector<amrex::Real> const& x,
    amrex::Real alpha)
{
    AMREX_ALWAYS_ASSERT(y.size() == x.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        y[i] += alpha * x[i];
    }
}

/**
 * Select retained modal count for a normalized z depth.
 *
 * @param full_modes Full number of modes in one transverse direction.
 * @param z_fraction Depth from zmin normalized by `Lz`.
 * @return Number of modes retained in that direction.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int
LayerModes (int full_modes, amrex::Real z_fraction)
{
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

#if defined(WARPX_DIM_3D)

/**
 * Build the spectral Dirichlet-to-Neumann eigenvalue table.
 *
 * @param lambda Output `lambda[m,n] = k_mn coth(k_mn Lz)`.
 * @param nx Number of x modes.
 * @param ny Number of y modes.
 * @param geom Level geometry that defines physical lengths.
 */
void
BuildLambda (
    std::vector<amrex::Real>& lambda,
    int nx,
    int ny,
    amrex::Geometry const& geom)
{
    lambda.assign(static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));
    amrex::Real const lx = geom.ProbHi(0) - geom.ProbLo(0);
    amrex::Real const ly = geom.ProbHi(1) - geom.ProbLo(1);
    amrex::Real const lz = geom.ProbHi(2) - geom.ProbLo(2);

    for (int m = 0; m < nx; ++m) {
        amrex::Real const kx =
            static_cast<amrex::Real>(m + 1) * MathConst::pi / lx;
        for (int n = 0; n < ny; ++n) {
            amrex::Real const ky =
                static_cast<amrex::Real>(n + 1) * MathConst::pi / ly;
            amrex::Real const k = std::sqrt(kx * kx + ky * ky);
            lambda[Index2D(m, n, nx)] = k * StableCoth(k * lz);
        }
    }
}

/**
 * Apply the compact masked Schur operator to a trial vector.
 *
 * @param v Compact vector over `NeumannMask == true` nodes.
 * @param av Output compact operator result.
 * @param neumann_mask Full interior-face mask.
 * @param lambda Spectral DtN eigenvalues.
 * @param sx x-direction sine basis.
 * @param sy y-direction sine basis.
 * @param nx Number of x interior nodes/modes.
 * @param ny Number of y interior nodes/modes.
 */
void
ApplySchurOperator (
    std::vector<amrex::Real> const& v,
    std::vector<amrex::Real>& av,
    std::vector<char> const& neumann_mask,
    std::vector<amrex::Real> const& lambda,
    std::vector<amrex::Real> const& sx,
    std::vector<amrex::Real> const& sy,
    int nx,
    int ny)
{
    // Apply P_N Lambda P_N^T: expand compact Neumann unknowns to the full face,
    // apply the full DtN operator spectrally, then restrict back to Gamma_N.
    std::vector<amrex::Real> face(
        static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));

    int p = 0;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (neumann_mask[Index2D(i, j, nx)] != 0) {
                face[Index2D(i, j, nx)] = v[p++];
            }
        }
    }

    std::vector<amrex::Real> face_hat;
    ApplyDST2(face, face_hat, sx, sy, nx, ny);
    for (std::size_t i = 0; i < face_hat.size(); ++i) {
        face_hat[i] *= lambda[i];
    }

    std::vector<amrex::Real> q_face;
    ApplyIDST2(face_hat, q_face, sx, sy, nx, ny);

    av.assign(v.size(), static_cast<amrex::Real>(0.0));
    p = 0;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (neumann_mask[Index2D(i, j, nx)] != 0) {
                av[p++] = q_face[Index2D(i, j, nx)];
            }
        }
    }
}

/**
 * Solve the compact Schur system with unpreconditioned conjugate gradients.
 *
 * @param x Output compact boundary unknown vector `u_N`.
 * @param b Compact right-hand side.
 * @param neumann_mask Full interior-face mask.
 * @param lambda Spectral DtN eigenvalues.
 * @param sx x-direction sine basis.
 * @param sy y-direction sine basis.
 * @param nx Number of x interior nodes/modes.
 * @param ny Number of y interior nodes/modes.
 * @param config Iteration limits and tolerances.
 */
void
SolveCG (
    std::vector<amrex::Real>& x,
    std::vector<amrex::Real> const& b,
    std::vector<char> const& neumann_mask,
    std::vector<amrex::Real> const& lambda,
    std::vector<amrex::Real> const& sx,
    std::vector<amrex::Real> const& sy,
    int nx,
    int ny,
    SchurConfig const& config)
{
    x.assign(b.size(), static_cast<amrex::Real>(0.0));
    if (b.empty()) {
        return;
    }

    std::vector<amrex::Real> r = b;
    std::vector<amrex::Real> p = r;
    std::vector<amrex::Real> ap(b.size(), static_cast<amrex::Real>(0.0));

    amrex::Real rs_old = Dot(r, r);
    amrex::Real const b_norm = std::sqrt(
        std::max(Dot(b, b), std::numeric_limits<amrex::Real>::min()));
    amrex::Real const tolerance = std::max(config.abs_tol, config.rel_tol * b_norm);
    amrex::Real const tolerance_sq = tolerance * tolerance;

    if (rs_old <= tolerance_sq) {
        return;
    }

    for (int iter = 0; iter < config.max_iter; ++iter) {
        ApplySchurOperator(p, ap, neumann_mask, lambda, sx, sy, nx, ny);
        amrex::Real const denom = Dot(p, ap);
        if (std::abs(denom) <= std::numeric_limits<amrex::Real>::min()) {
            break;
        }

        amrex::Real const alpha = rs_old / denom;
        AddScaled(x, p, alpha);
        AddScaled(r, ap, -alpha);

        amrex::Real const rs_new = Dot(r, r);
        if (rs_new <= tolerance_sq) {
            break;
        }

        amrex::Real const beta = rs_new / rs_old;
        for (std::size_t i = 0; i < p.size(); ++i) {
            p[i] = r[i] + beta * p[i];
        }
        rs_old = rs_new;
    }
}

/**
 * Build the Neumann mask and compact Schur right-hand side.
 *
 * WarpX already applies the anode Dirichlet potential and homogeneous Neumann
 * condition on the remaining ceramic surface. The Schur RHS is therefore just
 * the target nonhomogeneous normal derivative from `sigma_s`.
 *
 * @param neumann_mask Output interior-face bool mask stored as `char`.
 * @param rhs Output compact RHS over true mask entries.
 * @param sigma_s Full-face surface charge density. Empty means use constant config value.
 * @param geom Level geometry.
 * @param nx_internal Number of x interior nodes.
 * @param ny_internal Number of y interior nodes.
 * @param nx_face Full zmin face node count in x.
 * @param config Schur configuration.
 */
void
BuildMaskAndRhs (
    std::vector<char>& neumann_mask,
    std::vector<amrex::Real>& rhs,
    std::vector<amrex::Real> const& sigma_s,
    amrex::Geometry const& geom,
    int nx_internal,
    int ny_internal,
    int nx_face,
    SchurConfig const& config)
{
    amrex::Box domain = geom.Domain();
    domain.surroundingNodes();
    int const xlo = domain.smallEnd(0);
    int const ylo = domain.smallEnd(1);
    int const zlo = domain.smallEnd(2);
    amrex::Real const problo_x = geom.ProbLo(0);
    amrex::Real const problo_y = geom.ProbLo(1);
    amrex::Real const dx = geom.CellSize(0);
    amrex::Real const dy = geom.CellSize(1);

#ifdef HALL3D
    auto const anode_config = ReadHallAnodeRingConfig(geom);
#endif

    neumann_mask.assign(
        static_cast<std::size_t>(nx_internal) * ny_internal, static_cast<char>(0));
    rhs.clear();
    rhs.reserve(static_cast<std::size_t>(nx_internal) * ny_internal);

    // Only interior transverse nodes enter the sine basis. Dirichlet/anode nodes
    // are represented by NeumannMask = false and omitted from compact RHS/u_N.
    for (int j = 0; j < ny_internal; ++j) {
        int const gj = ylo + 1 + j;
        for (int i = 0; i < nx_internal; ++i) {
            int const gi = xlo + 1 + i;

            bool is_neumann = true;
#ifdef HALL3D
            if (IsHallAnodeRingNode(
                    gi, gj, zlo, zlo, xlo, ylo, problo_x, problo_y, dx, dy, anode_config))
            {
                is_neumann = false;
            }
#endif

            if (is_neumann) {
                neumann_mask[Index2D(i, j, nx_internal)] = static_cast<char>(1);
                amrex::Real const sigma =
                    sigma_s.empty() ? config.sigma_s :
                                      sigma_s[Index2D(gi - xlo, gj - ylo, nx_face)];
                // sigma_s = -epsilon0 * dphi/dn, so dphi/dn = -sigma_s/epsilon0.
                rhs.push_back(-sigma / PhysConst::epsilon_0);
            }
        }
    }
}

/**
 * Expand compact boundary unknowns into a full interior face array.
 *
 * @param u_face Output full interior face values.
 * @param u_N Compact boundary unknowns.
 * @param neumann_mask Full interior-face mask.
 * @param nx Number of x interior nodes.
 * @param ny Number of y interior nodes.
 */
void
ExpandBoundaryUnknowns (
    std::vector<amrex::Real>& u_face,
    std::vector<amrex::Real> const& u_N,
    std::vector<char> const& neumann_mask,
    int nx,
    int ny)
{
    u_face.assign(static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));
    int p = 0;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (neumann_mask[Index2D(i, j, nx)] != 0) {
                u_face[Index2D(i, j, nx)] = u_N[p++];
            }
        }
    }
}

/**
 * Copy surface charge density from device memory to a host vector.
 *
 * @param sigma_s_device Device vector with x-fast full-face layout.
 * @return Host copy of the input surface charge density.
 */
std::vector<amrex::Real>
BuildHostSigma (
    amrex::Gpu::DeviceVector<amrex::Real> const& sigma_s_device)
{
    amrex::Gpu::HostVector<amrex::Real> host_sigma(sigma_s_device.size());
    if (!host_sigma.empty()) {
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost,
            sigma_s_device.begin(),
            sigma_s_device.end(),
            host_sigma.begin());
        amrex::Gpu::streamSynchronize();
    }
    return std::vector<amrex::Real>(host_sigma.begin(), host_sigma.end());
}

/**
 * Compute `sinh(k(Lz-z)) / sinh(kLz)` without large-argument overflow.
 *
 * @param k Transverse modal wave number.
 * @param z Distance from zmin.
 * @param lz Domain length in z.
 * @return Harmonic extension decay factor.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
SinhDecayRatio (
    amrex::Real k,
    amrex::Real z,
    amrex::Real lz)
{
    amrex::Real const a = k * lz;
    amrex::Real const b = k * (lz - z);
    constexpr amrex::Real threshold = static_cast<amrex::Real>(40.0);
    if (a <= threshold) {
        amrex::Real const denom = std::sinh(a);
        return denom == static_cast<amrex::Real>(0.0) ?
            static_cast<amrex::Real>(0.0) : std::sinh(b) / denom;
    }

    amrex::Real const numerator_correction =
        static_cast<amrex::Real>(1.0) - std::exp(
            -static_cast<amrex::Real>(2.0) * b);
    amrex::Real const denominator_correction =
        static_cast<amrex::Real>(1.0) - std::exp(
            -static_cast<amrex::Real>(2.0) * a);
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
ReconstructPhiCorrection (
    amrex::Geometry const& geom,
    std::vector<amrex::Real> const& u_hat,
    std::vector<amrex::Real> const& sx,
    std::vector<amrex::Real> const& sy,
    int nx,
    int ny)
{
    // The correction MultiFab follows rho_fp level 0 layout so later field updates
    // can use WarpX-owned decomposition and guard-cell shape.
    auto& state = State();
    WarpX& warpx = WarpX::GetInstance();
    auto const* rho_fp = warpx.m_fields.get(warpx::fields::FieldType::rho_fp, 0);
    AMREX_ALWAYS_ASSERT(rho_fp != nullptr);

    state.phi_corr = std::make_unique<amrex::MultiFab>(
        rho_fp->boxArray(), rho_fp->DistributionMap(), 1, rho_fp->nGrowVect());
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

    amrex::Gpu::DeviceVector<amrex::Real> d_u_hat(u_hat.size());
    amrex::Gpu::DeviceVector<amrex::Real> d_sx(sx.size());
    amrex::Gpu::DeviceVector<amrex::Real> d_sy(sy.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, u_hat.begin(), u_hat.end(), d_u_hat.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, sx.begin(), sx.end(), d_sx.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, sy.begin(), sy.end(), d_sy.begin());

    amrex::Real const* u_hat_ptr = d_u_hat.dataPtr();
    amrex::Real const* sx_ptr = d_sx.dataPtr();
    amrex::Real const* sy_ptr = d_sy.dataPtr();

    amrex::Gpu::DeviceVector<amrex::Real> d_a(
        static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));
    amrex::Gpu::DeviceVector<amrex::Real> d_tmp(
        static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));
    amrex::Gpu::DeviceVector<amrex::Real> d_layer(
        static_cast<std::size_t>(nx) * ny, static_cast<amrex::Real>(0.0));
    amrex::Real* a_ptr = d_a.dataPtr();
    amrex::Real* tmp_ptr = d_tmp.dataPtr();
    amrex::Real* layer_ptr = d_layer.dataPtr();

    for (int kk = 0; kk < nz_nodes; ++kk) {
        amrex::Real const z = static_cast<amrex::Real>(kk) * dz;
        amrex::Real const z_fraction =
            lz > static_cast<amrex::Real>(0.0) ?
            z / lz : static_cast<amrex::Real>(0.0);
        int const mx = LayerModes(nx, z_fraction);
        int const my = LayerModes(ny, z_fraction);

        // Step 1: build the retained spectral coefficients A(m,n,z).
        amrex::Box const spectral_box(
            amrex::IntVect(0, 0, 0),
            amrex::IntVect(mx - 1, my - 1, 0));
        amrex::ParallelFor(
            spectral_box,
            [=] AMREX_GPU_DEVICE (int m, int n, int)
            {
                amrex::Real const kx =
                    static_cast<amrex::Real>(m + 1) * MathConst::pi / lx;
                amrex::Real const ky =
                    static_cast<amrex::Real>(n + 1) * MathConst::pi / ly;
                amrex::Real const kmn = std::sqrt(kx * kx + ky * ky);
                amrex::Real const r = SinhDecayRatio(kmn, z, lz);
                a_ptr[Index2D(m, n, nx)] = r * u_hat_ptr[Index2D(m, n, nx)];
            });

        // Step 2: tmp(i,n,z) = sum_m Sx(i,m) A(m,n,z).
        amrex::Box const tmp_box(
            amrex::IntVect(0, 0, 0),
            amrex::IntVect(nx - 1, my - 1, 0));
        amrex::ParallelFor(
            tmp_box,
            [=] AMREX_GPU_DEVICE (int i, int n, int)
            {
                amrex::Real sum = static_cast<amrex::Real>(0.0);
                for (int m = 0; m < mx; ++m) {
                    sum += sx_ptr[Index2D(i, m, nx)] * a_ptr[Index2D(m, n, nx)];
                }
                tmp_ptr[Index2D(i, n, nx)] = sum;
            });

        // Step 3: layer(i,j,z) = sum_n tmp(i,n,z) Sy(j,n).
        amrex::Box const layer_box(
            amrex::IntVect(0, 0, 0),
            amrex::IntVect(nx - 1, ny - 1, 0));
        amrex::ParallelFor(
            layer_box,
            [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                amrex::Real sum = static_cast<amrex::Real>(0.0);
                for (int n = 0; n < my; ++n) {
                    sum += tmp_ptr[Index2D(i, n, nx)] * sy_ptr[Index2D(j, n, ny)];
                }
                layer_ptr[Index2D(i, j, nx)] = sum;
            });

        // Scatter the reconstructed z layer into whatever boxes own that slab.
        amrex::Box const slab(
            amrex::IntVect(domain.smallEnd(0), domain.smallEnd(1), zlo + kk),
            amrex::IntVect(domain.bigEnd(0), domain.bigEnd(1), zlo + kk));
        for (amrex::MFIter mfi(*state.phi_corr, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box const write_box = mfi.validbox() & slab;
            if (!write_box.ok()) {
                continue;
            }
            amrex::Array4<amrex::Real> const& arr = state.phi_corr->array(mfi);
            amrex::ParallelFor(
                write_box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    int const ii = i - xlo - 1;
                    int const jj = j - ylo - 1;
                    if (ii >= 0 && ii < nx && jj >= 0 && jj < ny) {
                        arr(i, j, k) = layer_ptr[Index2D(ii, jj, nx)];
                    } else {
                        arr(i, j, k) = static_cast<amrex::Real>(0.0);
                    }
                });
        }
    }

    amrex::Gpu::streamSynchronize();
}

/**
 * Run the full host-side Schur solve and optional device-side reconstruction.
 *
 * @param geom Level geometry.
 * @param sigma_s Host surface charge density in full-face x-fast layout.
 * @param nx_face Full zmin face node count in x.
 * @param ny_face Full zmin face node count in y.
 */
void
SolveAndReconstructHost (
    amrex::Geometry const& geom,
    std::vector<amrex::Real> const& sigma_s,
    int nx_face,
    int ny_face)
{
    SchurConfig const& config = ReadConfig();

    amrex::Box domain = geom.Domain();
    domain.surroundingNodes();
    int const expected_nx_face = domain.length(0);
    int const expected_ny_face = domain.length(1);
    AMREX_ALWAYS_ASSERT(nx_face == expected_nx_face);
    AMREX_ALWAYS_ASSERT(ny_face == expected_ny_face);

    int const nx_internal = nx_face - 2;
    int const ny_internal = ny_face - 2;
    AMREX_ALWAYS_ASSERT(nx_internal > 0 && ny_internal > 0);

    std::vector<amrex::Real> sx;
    std::vector<amrex::Real> sy;
    std::vector<amrex::Real> lambda;
    BuildSineBasis(sx, nx_internal);
    BuildSineBasis(sy, ny_internal);
    BuildLambda(lambda, nx_internal, ny_internal, geom);

    std::vector<char> neumann_mask;
    std::vector<amrex::Real> rhs;
    BuildMaskAndRhs(
        neumann_mask, rhs, sigma_s, geom,
        nx_internal, ny_internal, nx_face, config);

    auto& state = State();
    state.nx_internal = nx_internal;
    state.ny_internal = ny_internal;
    SolveCG(
        state.u_N, rhs, neumann_mask, lambda, sx, sy,
        nx_internal, ny_internal, config);

    ExpandBoundaryUnknowns(
        state.u_face, state.u_N, neumann_mask, nx_internal, ny_internal);
    ApplyDST2(state.u_face, state.u_hat, sx, sy, nx_internal, ny_internal);

    if (config.rebuild_volume_field) {
        ReconstructPhiCorrection(geom, state.u_hat, sx, sy, nx_internal, ny_internal);
    }
}

#endif

} // namespace

namespace Insert {

/**
 * Check whether the Schur boundary correction is enabled.
 *
 * @return true when `insert.schur_boundary.enabled` is nonzero.
 */
bool
SpectralBoundarySchur::Enabled ()
{
    return ReadConfig().enabled;
}

/**
 * Apply the configured zmin Schur correction path to WarpX level-0 fields.
 *
 * This convenience entry creates the configured `sigma_s` device vector from
 * runtime input. External wall-charge deposition can call `SolveAndReconstruct`
 * directly.
 *
 * @param phi Multi-level nodal potential field.
 */
void
SpectralBoundarySchur::ApplyZMinCorrection (
    ablastr::fields::MultiLevelScalarField const& phi)
{
    SchurConfig const& config = ReadConfig();
    if (!config.enabled) {
        return;
    }

#if defined(WARPX_DIM_3D)
    if (config.face != "zmin") {
        amrex::Abort("insert.schur_boundary only supports face = zmin");
    }
    if (config.backend != "cpu_serial") {
        amrex::Abort("insert.schur_boundary only supports backend = cpu_serial");
    }
    if (phi.size() != 1) {
        amrex::Abort("insert.schur_boundary only supports single-level level 0");
    }
    if (EB::enabled()) {
        amrex::Abort("insert.schur_boundary does not support EB");
    }
    if (amrex::ParallelDescriptor::NProcs() != 1) {
        amrex::Abort("insert.schur_boundary cpu_serial backend requires one MPI rank");
    }

    WarpX& warpx = WarpX::GetInstance();
    amrex::Box domain = warpx.Geom(0).Domain();
    domain.surroundingNodes();
    int const nx_face = domain.length(0);
    int const ny_face = domain.length(1);

    if (UseZMinWallChargeSource(config.source)) {
        ZMinWallChargeGrid const grid = MakeZMinWallChargeGrid(warpx.Geom(0));
        if (grid.nx != nx_face || grid.ny != ny_face) {
            amrex::Abort("insert.schur_boundary zmin wall-charge grid size mismatch");
        }
        auto sigma_s_device = DepositZMinWallChargeDensity(warpx, grid);
        SolveAndReconstruct(warpx.Geom(0), sigma_s_device, nx_face, ny_face);
    } else if (UseConstantSigmaSource(config.source)) {
        amrex::Gpu::DeviceVector<amrex::Real> sigma_s_device(
            static_cast<std::size_t>(nx_face) * ny_face, config.sigma_s);
        SolveAndReconstruct(warpx.Geom(0), sigma_s_device, nx_face, ny_face);
    } else {
        amrex::Abort("unknown insert.schur_boundary.source");
    }
#else
    amrex::Abort("insert.schur_boundary requires WarpX_DIMS=3");
#endif
}

/**
 * Solve the zmin Schur boundary problem from an externally provided surface charge.
 *
 * @param geom Level-0 geometry.
 * @param sigma_s_device Surface charge density in full-face x-fast layout.
 * @param nx_face Full zmin face node count in x.
 * @param ny_face Full zmin face node count in y.
 */
void
SpectralBoundarySchur::SolveAndReconstruct (
    amrex::Geometry const& geom,
    amrex::Gpu::DeviceVector<amrex::Real> const& sigma_s_device,
    int nx_face,
    int ny_face)
{
#if defined(WARPX_DIM_3D)
    if (amrex::ParallelDescriptor::NProcs() != 1) {
        amrex::Abort("insert.schur_boundary cpu_serial backend requires one MPI rank");
    }
    if (sigma_s_device.size() != static_cast<std::size_t>(nx_face) * ny_face) {
        amrex::Abort("insert.schur_boundary sigma_s_device size does not match face size");
    }

    std::vector<amrex::Real> const sigma_s = BuildHostSigma(sigma_s_device);
    SolveAndReconstructHost(geom, sigma_s, nx_face, ny_face);
#else
    (void)geom;
    (void)sigma_s_device;
    (void)nx_face;
    (void)ny_face;
    amrex::Abort("insert.schur_boundary requires WarpX_DIMS=3");
#endif
}

/**
 * Access the latest compact boundary unknown vector.
 *
 * @return Compact `u_N` over `NeumannMask == true` nodes.
 */
std::vector<amrex::Real> const&
SpectralBoundarySchur::BoundaryUnknowns ()
{
    return State().u_N;
}

/**
 * Access the latest full interior zmin boundary trace.
 *
 * @return Full interior-face `u_face`, x-fast.
 */
std::vector<amrex::Real> const&
SpectralBoundarySchur::BoundaryFace ()
{
    return State().u_face;
}

/**
 * Access the latest full boundary sine spectrum.
 *
 * @return Full `u_hat[m,n]`, x-mode-fast.
 */
std::vector<amrex::Real> const&
SpectralBoundarySchur::BoundarySpectrum ()
{
    return State().u_hat;
}

/**
 * Access the latest reconstructed correction potential.
 *
 * @return Pointer to `phi_corr`, or null if volume reconstruction was not run.
 */
amrex::MultiFab const*
SpectralBoundarySchur::PhiCorrection ()
{
    return State().phi_corr.get();
}

/**
 * Narrow Insert-level entry point used by the electrostatic solver.
 *
 * @param phi Multi-level nodal potential field.
 */
void
ApplyElectrostaticBoundaryCorrection (
    ablastr::fields::MultiLevelScalarField const& phi)
{
    SpectralBoundarySchur::ApplyZMinCorrection(phi);
}

} // namespace Insert
