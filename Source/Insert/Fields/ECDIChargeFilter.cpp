#include "ECDIChargeFilter.h"

#include "WarpX.H"

#include <ablastr/profiler/ProfilerWrapper.H>

#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace
{

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
NodeControlVolume (
    int const i, int const j, int const k,
    int const ilo, int const ihi,
    int const jlo, int const jhi,
    int const klo, int const khi,
    amrex::Real const dV) noexcept
{
    // Nodal source terms represent half/quarter/eighth control volumes on physical boundaries.
    amrex::Real omega = dV;
    if (i == ilo || i == ihi) {
        omega *= amrex::Real(0.5);
    }
    if (j == jlo || j == jhi) {
        omega *= amrex::Real(0.5);
    }
    if (k == klo || k == khi) {
        omega *= amrex::Real(0.5);
    }
    return omega;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
AddWeightedChargeToBins (
    amrex::Real* const AMREX_RESTRICT q,
    amrex::Real* const AMREX_RESTRICT d,
    int const nr,
    amrex::Real const dr,
    amrex::Real const rmax,
    int const k_bin,
    amrex::Real const r,
    amrex::Real const charge,
    amrex::Real const volume) noexcept
{
    long const base = static_cast<long>(k_bin) * static_cast<long>(nr + 1);

    // Safety path for degenerate radial grids and for explicitly under-resolved user input.
    if (nr == 0 || r >= rmax) {
        int const ir = nr;
        amrex::HostDevice::Atomic::Add(&q[base + ir], charge * volume);
        amrex::HostDevice::Atomic::Add(&d[base + ir], volume);
        return;
    }

    // Linear hat-function weights: each node contributes to at most two radial bins.
    int i_left = static_cast<int>(amrex::Math::floor(r / dr));
    i_left = amrex::max(0, amrex::min(i_left, nr - 1));
    int const i_right = i_left + 1;

    amrex::Real const r_left = static_cast<amrex::Real>(i_left) * dr;
    amrex::Real const r_right = static_cast<amrex::Real>(i_right) * dr;
    amrex::Real const denom = r_right - r_left;

    if (denom < amrex::Real(1.0e-15)) {
        amrex::HostDevice::Atomic::Add(&q[base + i_left], charge * volume);
        amrex::HostDevice::Atomic::Add(&d[base + i_left], volume);
        return;
    }

    amrex::Real const w_right = (r - r_left) / denom;
    amrex::Real const w_left = (r_right - r) / denom;

    amrex::HostDevice::Atomic::Add(&q[base + i_left], charge * volume * w_left);
    amrex::HostDevice::Atomic::Add(&d[base + i_left], volume * w_left);
    amrex::HostDevice::Atomic::Add(&q[base + i_right], charge * volume * w_right);
    amrex::HostDevice::Atomic::Add(&d[base + i_right], volume * w_right);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
InterpolateFromBins (
    amrex::Real const* const AMREX_RESTRICT bar_rho,
    int const nr,
    amrex::Real const dr,
    amrex::Real const rmax,
    int const k_bin,
    amrex::Real const r) noexcept
{
    long const base = static_cast<long>(k_bin) * static_cast<long>(nr + 1);

    // Use exactly the same radial weights as projection to preserve charge.
    if (nr == 0 || r >= rmax) {
        return bar_rho[base + nr];
    }

    int i_left = static_cast<int>(amrex::Math::floor(r / dr));
    i_left = amrex::max(0, amrex::min(i_left, nr - 1));
    int const i_right = i_left + 1;

    amrex::Real const r_left = static_cast<amrex::Real>(i_left) * dr;
    amrex::Real const r_right = static_cast<amrex::Real>(i_right) * dr;
    amrex::Real const denom = r_right - r_left;

    if (denom < amrex::Real(1.0e-15)) {
        return bar_rho[base + i_left];
    }

    amrex::Real const w_right = (r - r_left) / denom;
    amrex::Real const w_left = (r_right - r) / denom;
    return bar_rho[base + i_left] * w_left + bar_rho[base + i_right] * w_right;
}

amrex::Real
MaxCornerRadius (
    amrex::Geometry const& geom,
    std::array<amrex::Real, 2> const& center)
{
    amrex::Box nodal_domain = geom.Domain();
    nodal_domain.surroundingNodes();

    int const ilo = nodal_domain.smallEnd(0);
    int const ihi = nodal_domain.bigEnd(0);
    int const jlo = nodal_domain.smallEnd(1);
    int const jhi = nodal_domain.bigEnd(1);
    int const cell_ilo = geom.Domain().smallEnd(0);
    int const cell_jlo = geom.Domain().smallEnd(1);
    amrex::Real const dx = geom.CellSize(0);
    amrex::Real const dy = geom.CellSize(1);
    amrex::Real const xlo = geom.ProbLo(0);
    amrex::Real const ylo = geom.ProbLo(1);

    // The farthest transverse nodal point is always one of the four nodal corners.
    amrex::Real rmax = amrex::Real(0.0);
    for (int const i : {ilo, ihi}) {
        amrex::Real const x = xlo + static_cast<amrex::Real>(i - cell_ilo) * dx - center[0];
        for (int const j : {jlo, jhi}) {
            amrex::Real const y = ylo + static_cast<amrex::Real>(j - cell_jlo) * dy - center[1];
            rmax = std::max(rmax, std::sqrt(x*x + y*y));
        }
    }
    return rmax;
}

void
PrintDiagnostics (Insert::ECDIChargeFilterDiagnostics const& diag)
{
    amrex::Print() << "insert_ecdi_filter_charge_before = "
                   << diag.charge_before << "\n";
    amrex::Print() << "insert_ecdi_filter_charge_after = "
                   << diag.charge_after << "\n";
    amrex::Print() << "insert_ecdi_filter_charge_abs_sum = "
                   << diag.charge_abs_sum << "\n";
    amrex::Print() << "insert_ecdi_filter_charge_relative_error = "
                   << diag.charge_relative_error << "\n";
    amrex::Print() << "insert_ecdi_filter_rho_l2_ratio = "
                   << diag.rho_l2_ratio << "\n";
    amrex::Print() << "insert_ecdi_filter_max_abs_difference = "
                   << diag.max_abs_difference << "\n";
    amrex::Print() << "insert_ecdi_filter_max_relative_difference = "
                   << diag.max_relative_difference << "\n";
    amrex::Print() << "insert_ecdi_filter_nr = " << diag.nr << "\n";
    amrex::Print() << "insert_ecdi_filter_dr = " << diag.dr << "\n";
    amrex::Print() << "insert_ecdi_filter_r_domain_max = "
                   << diag.r_domain_max << "\n";
    amrex::Print() << "insert_ecdi_filter_empty_bins = "
                   << diag.empty_bins << "\n";
}

} // namespace

namespace Insert
{

ECDIChargeFilterOptions
ReadECDIChargeFilterOptions ()
{
    ECDIChargeFilterOptions options;

    amrex::ParmParse const pp("insert.ecdi_control");
    int enabled = options.enabled ? 1 : 0;
    pp.query("enabled", enabled);
    options.enabled = (enabled != 0);
    pp.query("axis", options.axis);
    pp.query("nr", options.nr);
    pp.query("dr", options.dr);
    pp.query("filter_interval", options.filter_interval);

    int diagnostics = options.diagnostics ? 1 : 0;
    pp.query("diagnostics", diagnostics);
    options.diagnostics = (diagnostics != 0);

    if (pp.countval("center") > 0) {
        amrex::Vector<amrex::Real> center;
        pp.getarr("center", center);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            center.size() == 2,
            "insert.ecdi_control.center must contain exactly two values.");
        options.center = {center[0], center[1]};
    }

    return options;
}

void
ApplyECDIChargeFilter (
    amrex::MultiFab& rho,
    amrex::Geometry const& geom,
    ECDIChargeFilterOptions options,
    ECDIChargeFilterDiagnostics* diagnostics)
{
    ABLASTR_PROFILE("Insert::ApplyECDIChargeFilter");

#if !defined(WARPX_DIM_3D)
    amrex::ignore_unused(rho, geom, options, diagnostics);
    WARPX_ABORT_WITH_MESSAGE("Insert ECDI charge filter currently supports only WARPX_DIM_3D.");
#else
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        options.axis == "z",
        "insert.ecdi_control.axis currently supports only z.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        rho.nComp() == 1,
        "Insert ECDI charge filter currently supports only one rho component.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        rho.is_nodal(0) && rho.is_nodal(1) && rho.is_nodal(2),
        "Insert ECDI charge filter requires a fully nodal rho MultiFab.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        rho.boxArray().size() == 1,
        "Insert ECDI charge filter currently assumes a single rho box.");

    // Work with the nodal domain explicitly; geom.Domain() itself is cell-centered indexing.
    amrex::Box nodal_domain = geom.Domain();
    nodal_domain.surroundingNodes();

    int const ilo = nodal_domain.smallEnd(0);
    int const ihi = nodal_domain.bigEnd(0);
    int const jlo = nodal_domain.smallEnd(1);
    int const jhi = nodal_domain.bigEnd(1);
    int const klo = nodal_domain.smallEnd(2);
    int const khi = nodal_domain.bigEnd(2);
    int const cell_ilo = geom.Domain().smallEnd(0);
    int const cell_jlo = geom.Domain().smallEnd(1);
    int const cell_klo = geom.Domain().smallEnd(2);

    amrex::Real const dx = geom.CellSize(0);
    amrex::Real const dy = geom.CellSize(1);
    amrex::Real const dz = geom.CellSize(2);
    amrex::Real const xlo = geom.ProbLo(0);
    amrex::Real const ylo = geom.ProbLo(1);
    amrex::Real const                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             dV = dx * dy * dz;

    // Auto radial resolution: choose square-cell-scale bins unless the user overrides dr.
    if (options.dr < amrex::Real(0.0)) {
        options.dr = std::min(dx, dy);
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        options.dr > amrex::Real(0.0),
        "insert.ecdi_control.dr must be positive or negative for automatic selection.");

    amrex::Real const r_domain_max = MaxCornerRadius(geom, options.center);
    if (options.nr < 0) {
        options.nr = static_cast<int>(std::ceil(r_domain_max / options.dr));
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        options.nr >= 0,
        "insert.ecdi_control.nr must be non-negative or -1 for automatic selection.");

    amrex::Real const rmax = static_cast<amrex::Real>(options.nr) * options.dr;
    amrex::Real const coverage_tol = amrex::Real(64.0)
                                     * std::numeric_limits<amrex::Real>::epsilon()
                                     * std::max(amrex::Real(1.0), r_domain_max);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        rmax + coverage_tol >= r_domain_max,
        "insert.ecdi_control.nr * insert.ecdi_control.dr must cover the farthest nodal "
        "point in the transverse domain.");

    int const n_nodes_z = nodal_domain.length(2);
    long const n_bins = static_cast<long>(options.nr + 1) * static_cast<long>(n_nodes_z);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        n_bins > 0,
        "Insert ECDI charge filter requires at least one cylindrical bin.");

    amrex::Gpu::DeviceVector<amrex::Real> d_Q(n_bins, amrex::Real(0.0));
    amrex::Gpu::DeviceVector<amrex::Real> d_D(n_bins, amrex::Real(0.0));
    amrex::Gpu::DeviceVector<amrex::Real> d_bar_rho(n_bins, amrex::Real(0.0));
    amrex::Gpu::DeviceVector<amrex::Real> d_diag(7, amrex::Real(0.0));
    amrex::Gpu::DeviceVector<int> d_empty_bins(1, 0);

    amrex::Real* const AMREX_RESTRICT q_ptr = d_Q.dataPtr();
    amrex::Real* const AMREX_RESTRICT d_ptr = d_D.dataPtr();
    amrex::Real* const AMREX_RESTRICT bar_ptr = d_bar_rho.dataPtr();
    int* const empty_bins_ptr = (diagnostics != nullptr) ? d_empty_bins.dataPtr() : nullptr;

    int const nr = options.nr;
    amrex::Real const dr = options.dr;
    amrex::Real const center_x = options.center[0];
    amrex::Real const center_y = options.center[1];

    // Projection: accumulate charge Q and discrete volume D on (r,z) bins.
    for (amrex::MFIter mfi(rho, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& tbx = mfi.tilebox();
        amrex::Array4<amrex::Real const> const& rho_arr = rho.const_array(mfi);

        if (diagnostics != nullptr) {
            amrex::Real* const AMREX_RESTRICT diag_ptr = d_diag.dataPtr();
            amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::Real const x = xlo + static_cast<amrex::Real>(i - cell_ilo) * dx - center_x;
                amrex::Real const y = ylo + static_cast<amrex::Real>(j - cell_jlo) * dy - center_y;
                amrex::Real const r = std::sqrt(x*x + y*y);
                amrex::Real const volume = NodeControlVolume(
                    i, j, k, ilo, ihi, jlo, jhi, klo, khi, dV);
                amrex::Real const rho_val = rho_arr(i, j, k, 0);

                AddWeightedChargeToBins(
                    q_ptr, d_ptr, nr, dr, rmax, k - cell_klo, r, rho_val, volume);

                amrex::HostDevice::Atomic::Add(&diag_ptr[0], rho_val * volume);
                amrex::HostDevice::Atomic::Add(&diag_ptr[1], amrex::Math::abs(rho_val) * volume);
                amrex::HostDevice::Atomic::Add(&diag_ptr[2], rho_val * rho_val * volume);
            });
        } else {
            amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::Real const x = xlo + static_cast<amrex::Real>(i - cell_ilo) * dx - center_x;
                amrex::Real const y = ylo + static_cast<amrex::Real>(j - cell_jlo) * dy - center_y;
                amrex::Real const r = std::sqrt(x*x + y*y);
                amrex::Real const volume = NodeControlVolume(
                    i, j, k, ilo, ihi, jlo, jhi, klo, khi, dV);
                amrex::Real const rho_val = rho_arr(i, j, k, 0);

                AddWeightedChargeToBins(
                    q_ptr, d_ptr, nr, dr, rmax, k - cell_klo, r, rho_val, volume);
            });
        }
    }

    // Axisymmetric mean density: bar_rho = Q / D, using discrete volume D (on GPU).
    amrex::ParallelFor(n_bins, [=] AMREX_GPU_DEVICE (long ibin)
    {
        if (d_ptr[ibin] > amrex::Real(0.0)) {
            bar_ptr[ibin] = q_ptr[ibin] / d_ptr[ibin];
        } else {
            bar_ptr[ibin] = amrex::Real(0.0);
            if (empty_bins_ptr != nullptr) {
                amrex::Gpu::Atomic::Add(empty_bins_ptr, 1);
            }
        }
    });

    // Copy empty-bin count only when diagnostics are requested.
    int empty_bins = 0;
    if (diagnostics != nullptr) {
        amrex::Vector<int> h_empty_bins(1, 0);
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, d_empty_bins.begin(), d_empty_bins.end(),
            h_empty_bins.begin());
        empty_bins = h_empty_bins[0];
    }

    // Backfill: overwrite valid rho nodes with the m=0 value interpolated from the bins.
    for (amrex::MFIter mfi(rho, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& tbx = mfi.tilebox();
        amrex::Array4<amrex::Real> const& rho_arr = rho.array(mfi);

        if (diagnostics != nullptr) {
            amrex::Real* const AMREX_RESTRICT diag_ptr = d_diag.dataPtr();
            amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::Real const x = xlo + static_cast<amrex::Real>(i - cell_ilo) * dx - center_x;
                amrex::Real const y = ylo + static_cast<amrex::Real>(j - cell_jlo) * dy - center_y;
                amrex::Real const r = std::sqrt(x*x + y*y);
                amrex::Real const volume = NodeControlVolume(
                    i, j, k, ilo, ihi, jlo, jhi, klo, khi, dV);

                amrex::Real const rho_old = rho_arr(i, j, k, 0);
                amrex::Real const rho_new = InterpolateFromBins(
                    bar_ptr, nr, dr, rmax, k - cell_klo, r);
                amrex::Real const diff = rho_old - rho_new;
                amrex::Real const abs_diff = amrex::Math::abs(diff);
                amrex::Real const abs_old = amrex::Math::abs(rho_old);
                amrex::Real const rel_denom =
                    (abs_old > amrex::Real(1.0e-30)) ? abs_old : amrex::Real(1.0e-30);

                rho_arr(i, j, k, 0) = rho_new;

                amrex::HostDevice::Atomic::Add(&diag_ptr[3], rho_new * volume);
                amrex::HostDevice::Atomic::Add(&diag_ptr[4], diff * diff * volume);
                amrex::Gpu::Atomic::Max(&diag_ptr[5], abs_diff);
                amrex::Gpu::Atomic::Max(&diag_ptr[6], abs_diff / rel_denom);
            });
        } else {
            amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::Real const x = xlo + static_cast<amrex::Real>(i - cell_ilo) * dx - center_x;
                amrex::Real const y = ylo + static_cast<amrex::Real>(j - cell_jlo) * dy - center_y;
                amrex::Real const r = std::sqrt(x*x + y*y);

                rho_arr(i, j, k, 0) = InterpolateFromBins(
                    bar_ptr, nr, dr, rmax, k - cell_klo, r);
            });
        }
    }

    // Diagnostics are local by design: this first implementation has no MPI reduction.
    if (diagnostics != nullptr) {
        amrex::Vector<amrex::Real> h_diag(d_diag.size(), amrex::Real(0.0));
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_diag.begin(), d_diag.end(), h_diag.begin());

        diagnostics->charge_before = h_diag[0];
        diagnostics->charge_abs_sum = h_diag[1];
        diagnostics->charge_after = h_diag[3];
        diagnostics->charge_relative_error =
            (h_diag[1] > amrex::Real(0.0))
                ? amrex::Math::abs(h_diag[3] - h_diag[0]) / h_diag[1]
                : amrex::Real(0.0);
        diagnostics->rho_l2_ratio =
            (h_diag[2] > amrex::Real(0.0))
                ? std::sqrt(h_diag[4] / h_diag[2])
                : amrex::Real(0.0);
        diagnostics->max_abs_difference = h_diag[5];
        diagnostics->max_relative_difference = h_diag[6];
        diagnostics->r_domain_max = r_domain_max;
        diagnostics->dr = dr;
        diagnostics->nr = nr;
        diagnostics->empty_bins = empty_bins;
    }
#endif
}

void
FilterRhoForECDIControl (
    ablastr::fields::MultiLevelScalarField const& rho_fp,
    int const max_level)
{
    ECDIChargeFilterOptions const options = ReadECDIChargeFilterOptions();
    if (!options.enabled) {
        return;
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        max_level == 0,
        "Insert ECDI charge filter currently supports only max_level = 0.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::ParallelDescriptor::NProcs() == 1,
        "Insert ECDI charge filter currently supports only one MPI rank.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        options.filter_interval > 0,
        "insert.ecdi_control.filter_interval must be positive.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        rho_fp.size() == 1 && rho_fp[0] != nullptr,
        "Insert ECDI charge filter requires rho_fp[0].");

    WarpX& warpx = WarpX::GetInstance();
    if (warpx.getistep(0) % options.filter_interval != 0) {
        return;
    }

    ECDIChargeFilterDiagnostics diagnostics;
    ECDIChargeFilterDiagnostics* diagnostics_ptr = options.diagnostics ? &diagnostics : nullptr;
    ApplyECDIChargeFilter(*rho_fp[0], warpx.Geom(0), options, diagnostics_ptr);

    // Refresh guard cells and physical boundaries after modifying valid-domain rho.
    rho_fp[0]->FillBoundary(warpx.Geom(0).periodicity());
    warpx.ApplyRhofieldBoundary(0, rho_fp[0], PatchType::fine);

    if (options.diagnostics) {
        PrintDiagnostics(diagnostics);
    }
}

} // namespace Insert
