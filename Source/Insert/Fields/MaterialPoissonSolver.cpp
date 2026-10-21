#include "MaterialPoissonSolver.H"

#include "Utils/TextMsg.H"

#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <algorithm>
#include <memory>

namespace Insert {

MaterialPoissonSolver::MaterialPoissonSolver (
    MaterialPoissonCache& cache,
    ablastr::fields::ConstMultiLevelScalarField const& epsilon,
    amrex::Vector<amrex::iMultiFab const*> const& masks,
    ablastr::fields::ConstMultiLevelScalarField const& voltage,
    amrex::Vector<amrex::Geometry> const& geom,
    amrex::Vector<amrex::BoxArray> const& grids,
    amrex::Vector<amrex::DistributionMapping> const& dmap,
    ablastr::utils::enums::GridType const grid_type,
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> const& lobc,
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> const& hibc)
    : m_cache(cache), m_epsilon(epsilon), m_masks(masks), m_voltage(voltage),
      m_geom(geom), m_grids(grids), m_dmap(dmap), m_grid_type(grid_type),
      m_lobc(lobc), m_hibc(hibc)
{
    auto const nlev = epsilon.size();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        nlev > 0 && masks.size() == nlev && voltage.size() == nlev &&
        geom.size() >= nlev && grids.size() >= nlev && dmap.size() >= nlev,
        "Material Poisson data must cover every solver level.");
    if (m_cache.empty()) { m_cache.resize(nlev); }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_cache.size() == nlev,
        "Material Poisson cache level count changed after initialization.");
}

amrex::MLMG& MaterialPoissonSolver::operator() (
    int const lev, amrex::LPInfo const& base_info, amrex::MultiFab& phi,
    amrex::MultiFab& rhs, ablastr::fields::MLMGOptions const& options,
    amrex::Real const absolute_tolerance) const
{
    auto const* epsilon = m_epsilon.at(lev);
    auto const* mask = m_masks.at(lev);
    auto const* voltage = m_voltage.at(lev);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(epsilon && mask && voltage,
        "Material coefficients, masks and fixed potentials must not be null.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        epsilon->ixType().cellCentered() && epsilon->nComp() == 1 &&
        epsilon->boxArray().CellEqual(m_grids[lev]) &&
        epsilon->DistributionMap() == m_dmap[lev],
        "Relative permittivity must be a scalar cell-centered field on the solver layout.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        phi.nComp() == 1 && rhs.nComp() == 1 && voltage->nComp() == 1 && mask->nComp() == 1 &&
        voltage->boxArray() == phi.boxArray() && mask->boxArray() == phi.boxArray() &&
        rhs.boxArray() == phi.boxArray() && phi.boxArray().CellEqual(m_grids[lev]) &&
        voltage->DistributionMap() == phi.DistributionMap() &&
        mask->DistributionMap() == phi.DistributionMap() &&
        rhs.DistributionMap() == phi.DistributionMap() && phi.DistributionMap() == m_dmap[lev],
        "Material sources, potentials and masks must share the solver layout.");

    auto& cached = m_cache.at(lev);
    if (!cached) {
        auto info = base_info;
        // Preserve the existing hierarchy and smoothing for the interface-only refactor.
        info.setMaxCoarseningLevel(options.max_coarsening_level >= 0
            ? std::min(4, options.max_coarsening_level) : 4);
        cached = std::make_unique<MaterialPoissonLevelCache>();
        cached->m_geometry = m_geom[lev];
        cached->m_grids = m_grids[lev];
        cached->m_distribution = m_dmap[lev];
        cached->m_epsilon = epsilon;
        cached->m_mask = mask;
        cached->m_grid_type = m_grid_type;
        cached->m_linop = std::make_unique<amrex::MLNodeLaplacian>(
            amrex::Vector<amrex::Geometry>{m_geom[lev]},
            amrex::Vector<amrex::BoxArray>{m_grids[lev]},
            amrex::Vector<amrex::DistributionMapping>{m_dmap[lev]}, info);
        cached->m_linop->setSigma(0, *epsilon);
        cached->m_linop->setDomainBC(m_lobc, m_hibc);
        cached->m_linop->setOversetMask(0, *mask);
        cached->m_mlmg = std::make_unique<amrex::MLMG>(*cached->m_linop);
    } else {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            cached->m_linop && cached->m_mlmg && cached->m_grids == m_grids[lev] &&
            cached->m_distribution == m_dmap[lev] && cached->m_epsilon == epsilon &&
            cached->m_mask == mask && cached->m_grid_type == m_grid_type &&
            cached->m_geometry.Domain() == m_geom[lev].Domain(),
            "Material Poisson cache operator, material or layout changed.");
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                cached->m_geometry.ProbLo(d) == m_geom[lev].ProbLo(d) &&
                cached->m_geometry.ProbHi(d) == m_geom[lev].ProbHi(d) &&
                cached->m_geometry.isPeriodic(d) == m_geom[lev].isPeriodic(d),
                "Material Poisson cache geometry changed.");
        }
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            cached->m_linop->m_lobc_orig[0] == m_lobc &&
            cached->m_linop->m_hibc_orig[0] == m_hibc,
            "Material Poisson cache boundary types changed.");
    }

    auto& mlmg = *cached->m_mlmg;
    ablastr::fields::applyMLMGOptions(mlmg, options);
    mlmg.setPreSmooth(1);
    mlmg.setPostSmooth(2);
    mlmg.setConvergenceNormType(amrex::MLMGNormType::greater);
    mlmg.setNoGpuSync(true);
    if (m_grid_type == ablastr::utils::enums::GridType::Collocated) {
        mlmg.setFinalFillBC(true);
    }
    rhs.OverrideSync(m_geom[lev].periodicity());

    // Restore prescribed potentials after the caller's coarse-to-fine interpolation.
    for (amrex::MFIter mfi(phi, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const p = phi.array(mfi);
        auto const v = voltage->const_array(mfi);
        auto const m = mask->const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (m(i,j,k) == 0) { p(i,j,k) = v(i,j,k); }
            });
    }
    mlmg.solve({&phi}, {&rhs}, options.relative_tolerance, absolute_tolerance);
    return mlmg;
}

} // namespace Insert
