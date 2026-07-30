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
    // Bin geometry used by the filter.
    amrex::Real r_domain_max = amrex::Real(0.0);
    amrex::Real dr = amrex::Real(0.0);
    int nr = 0;
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
 * Runtime-parameter wrapper for the electrostatic solve path.
 *
 * If enabled, applies the ECDI charge filter to rho_fp[0] and refreshes guard
 * cells and physical boundaries. Does nothing when disabled.
 */
void FilterRhoForECDIControl (
    ablastr::fields::MultiLevelScalarField const& rho_fp,
    int max_level);

} // namespace Insert
