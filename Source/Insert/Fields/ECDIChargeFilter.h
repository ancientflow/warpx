#pragma once

#include "ablastr/fields/MultiFabRegister.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <array>
#include <string>

namespace Insert
{

struct ECDIChargeFilterOptions
{
    bool enabled = false;
    std::string axis = "z";
    // Transverse coordinates of the cylindrical averaging axis.
    std::array<amrex::Real, 2> center = {amrex::Real(0.0), amrex::Real(0.0)};
    // Negative values request automatic selection from the nodal domain.
    int nr = -1;
    amrex::Real dr = amrex::Real(-1.0);
    int filter_interval = 1;
    bool diagnostics = true;
};

struct ECDIChargeFilterDiagnostics
{
    // Charge and norm diagnostics use the same nodal control volumes as the filter.
    amrex::Real charge_before = amrex::Real(0.0);
    amrex::Real charge_after = amrex::Real(0.0);
    amrex::Real charge_abs_sum = amrex::Real(0.0);
    amrex::Real charge_relative_error = amrex::Real(0.0);
    amrex::Real rho_l2_ratio = amrex::Real(0.0);
    amrex::Real max_abs_difference = amrex::Real(0.0);
    amrex::Real max_relative_difference = amrex::Real(0.0);
    amrex::Real r_domain_max = amrex::Real(0.0);
    amrex::Real dr = amrex::Real(0.0);
    int nr = 0;
    int empty_bins = 0;
};

ECDIChargeFilterOptions ReadECDIChargeFilterOptions ();

/**
 * Apply the ECDI control filter to one fully nodal rho MultiFab.
 *
 * Current implementation scope follows ECDI_FILTER_ALGORITHM.md:
 * 3D, z-axis averaging, single component, single box, no AMR/MPI reduction.
 */
void ApplyECDIChargeFilter (
    amrex::MultiFab& rho,
    amrex::Geometry const& geom,
    ECDIChargeFilterOptions options,
    ECDIChargeFilterDiagnostics* diagnostics = nullptr);

/**
 * Runtime-parameter wrapper for future insertion into the electrostatic solve path.
 *
 * This function is intentionally not called from WarpX field-solver code yet.
 */
void FilterRhoForECDIControl (
    ablastr::fields::MultiLevelScalarField const& rho_fp,
    int max_level);

} // namespace Insert
