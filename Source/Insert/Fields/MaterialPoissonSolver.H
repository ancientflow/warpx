#ifndef WARPX_INSERT_MATERIAL_POISSON_SOLVER_H_
#define WARPX_INSERT_MATERIAL_POISSON_SOLVER_H_

#include "MaterialPoissonCache.H"

#include <ablastr/fields/MLMGOptions.H>
#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX_Array.H>
#include <AMReX_LO_BCTYPES.H>

namespace Insert {

/** Adapter for the existing material solver, borrowed for one computePhi call.
 *
 * The caller owns all level data and the persistent cache. This adapter owns no
 * time history and performs no interlevel interpolation. Stage 0 retains the
 * native material operator; replacing it is a separate implementation step.
 */
class MaterialPoissonSolver
{
public:
    MaterialPoissonSolver (
        MaterialPoissonCache& cache,
        ablastr::fields::ConstMultiLevelScalarField const& epsilon,
        amrex::Vector<amrex::iMultiFab const*> const& masks,
        ablastr::fields::ConstMultiLevelScalarField const& voltage,
        amrex::Vector<amrex::Geometry> const& geom,
        amrex::Vector<amrex::BoxArray> const& grids,
        amrex::Vector<amrex::DistributionMapping> const& dmap,
        ablastr::utils::enums::GridType grid_type,
        amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> const& lobc,
        amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> const& hibc);

    /** RHS is already scaled by -1/epsilon0; return the cached MLMG after solve. */
    amrex::MLMG& operator() (
        int lev, amrex::LPInfo const& base_info, amrex::MultiFab& phi,
        amrex::MultiFab& rhs, ablastr::fields::MLMGOptions const& options,
        amrex::Real absolute_tolerance) const;

private:
    MaterialPoissonCache& m_cache;
    ablastr::fields::ConstMultiLevelScalarField const& m_epsilon;
    amrex::Vector<amrex::iMultiFab const*> const& m_masks;
    ablastr::fields::ConstMultiLevelScalarField const& m_voltage;
    amrex::Vector<amrex::Geometry> const& m_geom;
    amrex::Vector<amrex::BoxArray> const& m_grids;
    amrex::Vector<amrex::DistributionMapping> const& m_dmap;
    ablastr::utils::enums::GridType m_grid_type;
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> m_lobc;
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> m_hibc;
};

} // namespace Insert

#endif
