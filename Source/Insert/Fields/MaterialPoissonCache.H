#ifndef WARPX_INSERT_MATERIALPOISSONCACHE_H_
#define WARPX_INSERT_MATERIALPOISSONCACHE_H_

#include <ablastr/utils/Enums.H>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_Vector.H>

#include <memory>

namespace Insert {

/** One immutable material operator and its solver on a static MR level.
 * Declaration order ensures MLMG is destroyed before the operator it references.
 * Fixed potentials and RHS are supplied afresh on every solve.
 */
struct MaterialPoissonLevelCache
{
    amrex::Geometry m_geometry;
    amrex::BoxArray m_grids;
    amrex::DistributionMapping m_distribution;
    amrex::MultiFab const* m_epsilon = nullptr;
    amrex::iMultiFab const* m_mask = nullptr;
    ablastr::utils::enums::GridType m_grid_type = ablastr::utils::enums::GridType::Staggered;
    std::unique_ptr<amrex::MLNodeLaplacian> m_linop;
    std::unique_ptr<amrex::MLMG> m_mlmg;
};

using MaterialPoissonCache = amrex::Vector<std::unique_ptr<MaterialPoissonLevelCache>>;

} // namespace Insert

#endif // WARPX_INSERT_MATERIALPOISSONCACHE_H_
