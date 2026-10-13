/* Copyright 2019-2022 Axel Huebl, Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_POISSON_SOLVER_H
#define ABLASTR_POISSON_SOLVER_H

#include <ablastr/constant.H>
#include <ablastr/utils/Communication.H>
#include <ablastr/utils/Enums.H>
#include <ablastr/utils/TextMsg.H>
#include <ablastr/warn_manager/WarnManager.H>
#include <ablastr/math/fft/AnyFFT.H>
#include <ablastr/fields/Interpolate.H>
#include <ablastr/fields/MLMGOptions.H>
#include <ablastr/fields/MultiFabRegister.H>
#include <ablastr/profiler/ProfilerWrapper.H>

#include "Insert/Config/HallElectrostaticConfig.H"
#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
#include "Insert/Fields/MaterialPoissonCache.H"
#endif

#if defined(ABLASTR_USE_FFT) && defined(WARPX_DIM_3D)
#include <ablastr/fields/IntegratedGreenFunctionSolver.H>
#endif

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MFIter.H>
#include <AMReX_MFInterp_C.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLLinOp.H>
#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
#   include <AMReX_MLNodeLaplacian.H>
#endif
#include <AMReX_MLNodeTensorLaplacian.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>
#include <AMReX_MLEBNodeFDLaplacian.H>
#ifdef AMREX_USE_EB
#   include <AMReX_EBFabFactory.H>
#endif

#include <algorithm>
#include <array>
#include <memory>
#include <optional>
#include <stdexcept>


namespace ablastr::fields {

/** Compute the L-infinity norm of the charge density `rho` across all MR levels
 * to determine if `rho` is zero everywhere
 *
 * \param[in] rho The charge density a given species
 * \param[in] finest_level The most refined mesh refinement level
 */
inline amrex::Real getMaxNormRho (
    ablastr::fields::ConstMultiLevelScalarField const& rho,
    int finest_level)
{
    amrex::Real max_norm_rho = 0.0;
    for (int lev=0; lev<=finest_level; lev++) {
        max_norm_rho = amrex::max(max_norm_rho, rho[lev]->norm0());
    }
    amrex::ParallelDescriptor::ReduceRealMax(max_norm_rho);

    if (max_norm_rho == 0) {
        ablastr::warn_manager::WMRecordWarning(
                "ElectrostaticSolver",
                "Max norm of rho is 0",
                ablastr::warn_manager::WarnPriority::low
        );
    }
    return max_norm_rho;
}

/** Interpolate the potential `phi` from level `lev` to `lev+1`
 *
 * Needed to solve Poisson equation on `lev+1`
 * The coarser level `lev` provides both
 * the boundary values and initial guess for `phi`
 * on the finer level `lev+1`
 *
 * \param[in] phi_lev  The potential on a given mesh refinement level `lev`
 * \param[in] phi_lev_plus_one The potential on the next level `lev+1`
 * \param[in] geom_lev The geometry of level `lev`
 * \param[in] do_single_precision_comms perform communications in single precision
 * \param[in] refratio mesh refinement ratio between level `lev` and `lev+1`
 * \param[in] ncomp Number of components of the multifab (1)
 * \param[in] ng Number of ghost cells (1 if collocated, 0 otherwise)
 */
inline void interpolatePhiBetweenLevels (
    amrex::MultiFab const* phi_lev,
    amrex::MultiFab* phi_lev_plus_one,
    amrex::Geometry const & geom_lev,
    bool do_single_precision_comms,
    const amrex::IntVect& refratio,
    const int ncomp,
    const int ng)
{
    using namespace amrex::literals;

    // Allocate phi_cp for lev+1
    amrex::BoxArray ba = phi_lev_plus_one->boxArray();
    ba.coarsen(refratio);
    amrex::MultiFab phi_cp(ba, phi_lev_plus_one->DistributionMap(), ncomp, ng);
    if (ng > 0) {
        // Set all values outside the domain to zero
        phi_cp.setDomainBndry(0.0_rt, geom_lev);
    }

    // Copy from phi[lev] to phi_cp (in parallel)
    const amrex::Periodicity& crse_period = geom_lev.periodicity();

    ablastr::utils::communication::ParallelCopy(
        phi_cp,
        *phi_lev,
        0,
        0,
        1,
        amrex::IntVect(0),
        amrex::IntVect(ng),
        do_single_precision_comms,
        crse_period
    );

    // Local interpolation from phi_cp to phi[lev+1]
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*phi_lev_plus_one, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Array4<amrex::Real> const phi_fp_arr = phi_lev_plus_one->array(mfi);
        amrex::Array4<amrex::Real const> const phi_cp_arr = phi_cp.array(mfi);

        details::PoissonInterpCPtoFP const interp(phi_fp_arr, phi_cp_arr, refratio);

        amrex::Box const& b = mfi.growntilebox(ng);
        amrex::ParallelFor(b, interp);
    }
}

/** Compute the potential `phi` by solving the Poisson equation
 *
 * Uses `rho` as a source, assuming that the source moves at a
 * constant speed \f$\boldsymbol{\beta}\f$. This uses the AMReX solver.
 *
 * More specifically, this solves the equation
 * \f[
 *   \nabla^2 r \phi - (\boldsymbol{\beta}\cdot\boldsymbol{\nabla})^2 r \phi = -\frac{r \rho}{\epsilon_0}
 * \f]
 *
 * \tparam T_PostPhiCalculationFunctor a calculation per level directly after phi was calculated
 * \tparam T_BoundaryHandler handler for boundary conditions, for example @see ElectrostaticSolver::PoissonBoundaryHandler (EB ONLY)
 * \tparam T_FArrayBoxFactory usually nothing or an amrex::EBFArrayBoxFactory (EB ONLY)
 * \param[in] rho The charge density a given species
 * \param[out] phi The potential to be computed by this function
 * \param[in] beta Represents the velocity of the source of `phi`
 * \param[in] mlmg_options Settings for the MLMG solve, @see ablastr::fields::MLMGOptions
 * \param[in] geom the geometry per level (e.g., from AmrMesh)
 * \param[in] dmap the distribution mapping per level (e.g., from AmrMesh)
 * \param[in] grids the grids per level (e.g., from AmrMesh)
 * \param[in] grid_type Integer that corresponds to the type of grid used in the simulation (collocated, staggered, hybrid)
 * \param[in] is_solver_igf_on_lev0 boolean to select the Poisson solver: 1 for FFT on level 0 & Multigrid on other levels, 0 for Multigrid on all levels
 * \param[in] is_igf_2d boolean to select between fully 3D Poisson solver and quasi-3D, i.e. one 2D Poisson solve on every z slice (default: false)
 * \param[in] eb_enabled solve with embedded boundaries
 * \param[in] do_single_precision_comms perform communications in single precision
 * \param[in] rel_ref_ratio mesh refinement ratio between levels (default: 1)
 * \param[in] post_phi_calculation perform a calculation per level directly after phi was calculated; required for embedded boundaries (default: none)
 * \param[in] boundary_handler a handler for boundary conditions, for example @see ElectrostaticSolver::PoissonBoundaryHandler
 * \param[in] current_time the current time; required for embedded boundaries (default: none)
 * \param[in] eb_farray_box_factory a factory for field data, @see amrex::EBFArrayBoxFactory; required for embedded boundaries (default: none)
 */
template<
    typename T_PostPhiCalculationFunctor = std::nullopt_t,
    typename T_BoundaryHandler = std::nullopt_t,
    typename T_FArrayBoxFactory = void
>
void
computePhi (
    ablastr::fields::MultiLevelScalarField const& rho,
    ablastr::fields::MultiLevelScalarField const& phi,
    std::array<amrex::Real, 3> const beta,
    MLMGOptions const & mlmg_options,
    amrex::Vector<amrex::Geometry> const& geom,
    amrex::Vector<amrex::DistributionMapping> const& dmap,
    amrex::Vector<amrex::BoxArray> const& grids,
    utils::enums::GridType grid_type,
    bool is_solver_igf_on_lev0,
    [[maybe_unused]] bool const is_igf_2d,
    bool eb_enabled = false,
    bool do_single_precision_comms = false,
    std::optional<amrex::Vector<amrex::IntVect> > rel_ref_ratio = std::nullopt,
    [[maybe_unused]] T_PostPhiCalculationFunctor post_phi_calculation = std::nullopt,
    [[maybe_unused]] T_BoundaryHandler const& boundary_handler = std::nullopt,
    [[maybe_unused]] std::optional<amrex::Real const> current_time = std::nullopt, // only used for EB
    [[maybe_unused]] const std::optional<amrex::Vector<T_FArrayBoxFactory const *> >& eb_farray_box_factory = std::nullopt, // only used for EB
    std::optional<amrex::Vector<amrex::iMultiFab const *> > overset_masks = std::nullopt
#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
    , ablastr::fields::ConstMultiLevelScalarField const* relative_permittivity = nullptr,
    Insert::MaterialPoissonCache* material_solver_cache = nullptr,
    ablastr::fields::ConstMultiLevelScalarField const* fixed_potential = nullptr
#endif
)
{
    using namespace amrex::literals;

    ABLASTR_PROFILE("computePhi");

    auto const finest_level = static_cast<int>(rho.size() - 1);

#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
    if (relative_permittivity != nullptr) {
        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<int>(relative_permittivity->size()) == finest_level + 1,
            "Relative-permittivity data must cover every Poisson level.");
        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
            !is_solver_igf_on_lev0,
            "Variable relative permittivity is supported only by the MLMG Poisson solver.");
        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
            !eb_enabled,
            "Variable relative permittivity cannot be combined with embedded boundaries.");
        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
            std::all_of(beta.begin(), beta.end(), [](amrex::Real const value) {
                return value == 0.0;
            }),
            "Variable relative permittivity is supported only for a lab-frame solve.");
    }
    bool const use_material_solver_cache = material_solver_cache != nullptr;
    ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
        !use_material_solver_cache || relative_permittivity != nullptr,
        "The material Poisson cache requires relative permittivity.");
    if (use_material_solver_cache) {
        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
            !overset_masks || overset_masks->size() == rho.size(),
            "Material Poisson cache masks must match the number of levels.");
        if (material_solver_cache->empty()) {
            material_solver_cache->resize(finest_level + 1);
        }
        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<int>(material_solver_cache->size()) == finest_level + 1,
            "Material Poisson cache level count changed after initialization.");
    }
    ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
        fixed_potential == nullptr ||
            (overset_masks && fixed_potential->size() == rho.size() &&
             overset_masks->size() == rho.size()),
        "Fixed material potentials require a matching mask on every level.");
#endif
    if (!rel_ref_ratio.has_value()) {
        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(finest_level == 0u,
                                           "rel_ref_ratio must be set if mesh-refinement is used");
        rel_ref_ratio = amrex::Vector<amrex::IntVect>{{amrex::IntVect(AMREX_D_DECL(1, 1, 1))}};
    }

#if !defined(AMREX_USE_EB)
    ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(!eb_enabled,
                                       "Embedded boundary solve requested but not compiled in");
#endif

#if !defined(ABLASTR_USE_FFT)
    ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE( !is_solver_igf_on_lev0,
    "Must compile with FFT support to use the IGF solver!");
#endif

#if !defined(WARPX_DIM_3D)
    ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE( !is_solver_igf_on_lev0,
    "The FFT Poisson solver is currently only implemented for 3D!");
#endif

    // Set the value of beta
    amrex::Array<amrex::Real, AMREX_SPACEDIM> beta_solver =
#if defined(WARPX_DIM_1D_Z)
            {{ beta[2] }};  // beta_x and beta_z
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            {{ beta[0], beta[2] }};  // beta_x and beta_z
#elif defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
            {{ beta[0] }};  // beta_x
#else
            {{ beta[0], beta[1], beta[2] }};
#endif

    ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE (
        std::all_of(beta_solver.begin(), beta_solver.end(),
            [](const auto b){return (b<1.0_rt);}),
        "Components of beta_solver must be < 1.");

    amrex::LPInfo info;
    applyMLMGOptions(info, mlmg_options);

#if (defined(ABLASTR_USE_FFT) && defined(WARPX_DIM_3D))
    // determine if rho is zero everywhere
    amrex::Real const max_norm_b = getMaxNormRho(
        amrex::GetVecOfConstPtrs(rho), finest_level);

    // A charge density that is zero everywhere makes a purely relative
    // convergence criterion meaningless, so require a finite absolute
    // tolerance in that case.
    amrex::Real const absolute_tolerance =
        (max_norm_b == 0._rt && mlmg_options.absolute_tolerance == 0._rt) ?
            amrex::Real(1.e-6) : mlmg_options.absolute_tolerance;
#else
    amrex::Real const absolute_tolerance = mlmg_options.absolute_tolerance;
#endif

    for (int lev=0; lev<=finest_level; lev++) {
        amrex::Array<amrex::Real,AMREX_SPACEDIM> const dx_scaled
                {AMREX_D_DECL(geom[lev].CellSize(0)/std::sqrt(1._rt-beta_solver[0]*beta_solver[0]),
                              geom[lev].CellSize(1)/std::sqrt(1._rt-beta_solver[1]*beta_solver[1]),
                              geom[lev].CellSize(2)/std::sqrt(1._rt-beta_solver[2]*beta_solver[2]))};

#if (defined(ABLASTR_USE_FFT) && defined(WARPX_DIM_3D))
        // Use the Integrated Green Function solver (FFT) on the coarsest level if it was selected
        if(is_solver_igf_on_lev0 && lev==0){
            if ( max_norm_b == 0 ) {
                phi[lev]->setVal(0);
            } else {
                computePhiIGF( *rho[lev], *phi[lev], dx_scaled, is_igf_2d);
            }
            continue;
        }
#endif
        // Use the Multigrid (MLMG) solver if selected or on refined patches
        // but first scale rho appropriately
        rho[lev]->mult(-1._rt / ablastr::constant::SI::epsilon_0);

#ifdef WARPX_DIM_RZ
        constexpr bool is_rz = true;
#else
        constexpr bool is_rz = false;
#endif

        if (!eb_enabled && !is_rz) {
            // Determine whether to use semi-coarsening
            int max_semicoarsening_level = 0;
            int semicoarsening_direction = -1;
            const auto min_dir = static_cast<int>(std::distance(dx_scaled.begin(),
                                                                std::min_element(dx_scaled.begin(), dx_scaled.end())));
            const auto max_dir = static_cast<int>(std::distance(dx_scaled.begin(),
                                                                std::max_element(dx_scaled.begin(), dx_scaled.end())));
            if (dx_scaled[max_dir] > dx_scaled[min_dir]) {
                semicoarsening_direction = max_dir;
                max_semicoarsening_level = static_cast<int>(std::log2(dx_scaled[max_dir] / dx_scaled[min_dir]));
            }
            if (max_semicoarsening_level > 0) {
                info.setSemicoarsening(true);
                info.setMaxSemicoarseningLevel(max_semicoarsening_level);
                info.setSemicoarseningDirection(semicoarsening_direction);
            }
        }

#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
        std::unique_ptr<amrex::MLNodeLinOp> owned_linop;
        amrex::MLNodeLinOp* linop = nullptr;
        bool configure_linop = true;
        if (relative_permittivity != nullptr) {
            // The volume-anode overset mask remains valid for at most four
            // multigrid coarsening levels.  Set this before constructing the
            // operator: the material operator and its MLMG instance are cached.
            info.setMaxCoarseningLevel(4);

            amrex::MultiFab const* const epsilon_r = (*relative_permittivity)[lev];
            ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
                epsilon_r != nullptr,
                "Relative-permittivity data contains a null level pointer.");
            ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
                epsilon_r->ixType().cellCentered() && epsilon_r->nComp() == 1,
                "Relative permittivity must be a single-component cell-centered field.");
            ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
                epsilon_r->boxArray().CellEqual(grids[lev]) &&
                    epsilon_r->DistributionMap() == dmap[lev],
                "Relative permittivity must use the Poisson grid layout.");

            if (use_material_solver_cache) {
                auto& cached = (*material_solver_cache)[lev];
                auto const* mask = overset_masks ? overset_masks->at(lev) : nullptr;
                if (!cached) {
                    cached = std::make_unique<Insert::MaterialPoissonLevelCache>();
                    cached->m_geometry = geom[lev];
                    cached->m_grids = grids[lev];
                    cached->m_distribution = dmap[lev];
                    cached->m_epsilon = epsilon_r;
                    cached->m_mask = mask;
                    cached->m_grid_type = grid_type;
                    cached->m_linop = std::make_unique<amrex::MLNodeLaplacian>(
                        amrex::Vector<amrex::Geometry>{geom[lev]},
                        amrex::Vector<amrex::BoxArray>{grids[lev]},
                        amrex::Vector<amrex::DistributionMapping>{dmap[lev]}, info);
                    cached->m_linop->setSigma(0, *epsilon_r);
                } else {
                    ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
                        cached->m_linop && cached->m_mlmg &&
                        cached->m_grids == grids[lev] && cached->m_distribution == dmap[lev] &&
                        cached->m_epsilon == epsilon_r && cached->m_mask == mask &&
                        cached->m_grid_type == grid_type &&
                        cached->m_geometry.Domain() == geom[lev].Domain(),
                        "Material Poisson cache operator, material or layout changed.");
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
                            cached->m_geometry.ProbLo(d) == geom[lev].ProbLo(d) &&
                            cached->m_geometry.ProbHi(d) == geom[lev].ProbHi(d) &&
                            cached->m_geometry.isPeriodic(d) == geom[lev].isPeriodic(d),
                            "Material Poisson cache geometry changed.");
                    }
                    if constexpr (!std::is_same_v<T_BoundaryHandler, std::nullopt_t>) {
                        ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
                            cached->m_linop->m_lobc_orig[0] == boundary_handler.lobc &&
                            cached->m_linop->m_hibc_orig[0] == boundary_handler.hibc,
                            "Material Poisson cache boundary types changed.");
                    }
                    configure_linop = false;
                }
                linop = cached->m_linop.get();
            } else {
                auto linop_nodelap = std::make_unique<amrex::MLNodeLaplacian>(
                    amrex::Vector<amrex::Geometry>{geom[lev]},
                    amrex::Vector<amrex::BoxArray>{grids[lev]},
                    amrex::Vector<amrex::DistributionMapping>{dmap[lev]},
                    info);
                linop_nodelap->setSigma(0, *epsilon_r);
                owned_linop = std::move(linop_nodelap);
                linop = owned_linop.get();
            }
        } else {
#else
        std::unique_ptr<amrex::MLNodeLinOp> linop;
#endif
        if (eb_enabled || is_rz) {
            // In the presence of EB or RZ: the solver assumes that the beam is
            // propagating along  one of the axes of the grid, i.e. that only *one*
            // of the components of `beta` is non-negligible.
            auto linop_nodelap = std::make_unique<amrex::MLEBNodeFDLaplacian>();
            if (eb_enabled) {
#if defined(AMREX_USE_EB)
                if constexpr(std::is_same_v<void, T_FArrayBoxFactory>) {
                    throw std::runtime_error("EB requested by eb_farray_box_factory not provided!");
                } else {
                    linop_nodelap->define(
                        amrex::Vector<amrex::Geometry>{geom[lev]},
                        amrex::Vector<amrex::BoxArray>{grids[lev]},
                        amrex::Vector<amrex::DistributionMapping>{dmap[lev]},
                        info,
                        amrex::Vector<amrex::EBFArrayBoxFactory const*>{eb_farray_box_factory.value()[lev]}
                    );
                }
#endif
            }
            else {
                // TODO: rather use MLNodeTensorLaplacian (for RZ w/o EB) here? Semi-Coarsening would be nice here
                linop_nodelap->define(
                    amrex::Vector<amrex::Geometry>{geom[lev]},
                    amrex::Vector<amrex::BoxArray>{grids[lev]},
                    amrex::Vector<amrex::DistributionMapping>{dmap[lev]},
                    info
                );
            }

            // Note: this assumes that the beam is propagating along
            // one of the axes of the grid, i.e. that only *one* of the
            // components of `beta` is non-negligible. // we use this
#if defined(WARPX_DIM_RZ)
            linop_nodelap->setRZ(true);
            linop_nodelap->setSigma({0._rt, 1._rt-beta_solver[1]*beta_solver[1]});
#else
            linop_nodelap->setSigma({AMREX_D_DECL(
                1._rt-beta_solver[0]*beta_solver[0],
                1._rt-beta_solver[1]*beta_solver[1],
                1._rt-beta_solver[2]*beta_solver[2])});
#endif
#if defined(AMREX_USE_EB)
            if (eb_enabled) {
                if constexpr (!std::is_same_v<T_BoundaryHandler, std::nullopt_t>) {
                    // if the EB potential only depends on time, the potential can be passed
                    // as a float instead of a callable
                    if (boundary_handler.phi_EB_only_t) {
                        linop_nodelap->setEBDirichlet(boundary_handler.potential_eb_t(current_time.value()));
                    } else {
                        linop_nodelap->setEBDirichlet(boundary_handler.getPhiEB(current_time.value()));
                    }
                } else
                {
                    ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE( !is_solver_igf_on_lev0,
                        "EB Poisson solver enabled but no 'boundary_handler' passed!");
                }
            }
#endif
#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
            owned_linop = std::move(linop_nodelap);
            linop = owned_linop.get();
#else
            linop = std::move(linop_nodelap);
#endif
        } else {
            // In the absence of EB and RZ: use a more generic solver
            // that can handle beams propagating in any direction
            auto linop_tenslap = std::make_unique<amrex::MLNodeTensorLaplacian>(
                amrex::Vector<amrex::Geometry>{geom[lev]},
                amrex::Vector<amrex::BoxArray>{grids[lev]},
                amrex::Vector<amrex::DistributionMapping>{dmap[lev]},
                info
            );
            linop_tenslap->setBeta(beta_solver); // for the non-axis-aligned solver
#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
            owned_linop = std::move(linop_tenslap);
            linop = owned_linop.get();
#else
            linop = std::move(linop_tenslap);
#endif
        }
#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
        }
#endif

#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
        if (configure_linop) {
#endif
        // Level 0 domain boundary
        if constexpr (std::is_same_v<T_BoundaryHandler, std::nullopt_t>) {
            amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> const lobc = {AMREX_D_DECL(
                amrex::LinOpBCType::Dirichlet,
                amrex::LinOpBCType::Dirichlet,
                amrex::LinOpBCType::Dirichlet
            )};
            amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> const hibc = lobc;
            linop->setDomainBC(lobc, hibc);
        } else {
            linop->setDomainBC(boundary_handler.lobc, boundary_handler.hibc);
        }

        if (overset_masks.has_value() && lev < static_cast<int>(overset_masks->size()) &&
            overset_masks.value()[lev] != nullptr)
        {
            linop->setOversetMask(0, *overset_masks.value()[lev]);
        }
#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
        }
#endif

        // Solve the Poisson equation
#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
        std::unique_ptr<amrex::MLMG> owned_mlmg;
        amrex::MLMG* selected_mlmg = nullptr;
        if (use_material_solver_cache) {
            auto& cached = *(*material_solver_cache)[lev];
            if (!cached.m_mlmg) {
                cached.m_mlmg = std::make_unique<amrex::MLMG>(*linop);
            }
            selected_mlmg = cached.m_mlmg.get();
        } else
        {
            owned_mlmg = std::make_unique<amrex::MLMG>(*linop);
            selected_mlmg = owned_mlmg.get();
        }
        amrex::MLMG& mlmg = *selected_mlmg;
#else
        amrex::MLMG mlmg(*linop); // actual solver defined here
        applyMLMGOptions(mlmg, mlmg_options);
#endif

        mlmg.setMaxIter(max_iters);
#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
        if (relative_permittivity != nullptr) {
            mlmg.setPreSmooth(1);
            mlmg.setPostSmooth(2);
        }
#endif
        mlmg.setFinalSmooth(num_final_sweeps.value_or(8));
        mlmg.setConvergenceNormType(amrex::MLMGNormType::greater);
        mlmg.setNoGpuSync(true);

        const int ng = int(grid_type == utils::enums::GridType::Collocated); // ghost cells
        if (ng) {
            // In this case, computeE needs to use ghost nodes data. So we
            // ask MLMG to fill BC for us after it solves the problem.
            mlmg.setFinalFillBC(true);
        }

        rho[lev]->OverrideSync(geom[lev].periodicity());

#ifdef WARPX_USE_HALL_ELECTROSTATIC_MATERIALS
        // Coarse-to-fine interpolation overwrites the fine initial guess.
        // Restore the prescribed anode potential here, on every solve.
        if (fixed_potential != nullptr) {
            auto const* fixed = fixed_potential->at(lev);
            auto const* mask = overset_masks->at(lev);
            ABLASTR_ALWAYS_ASSERT_WITH_MESSAGE(
                fixed && mask && fixed->nComp() == 1 &&
                fixed->boxArray() == phi[lev]->boxArray() &&
                mask->boxArray() == phi[lev]->boxArray() &&
                fixed->DistributionMap() == phi[lev]->DistributionMap() &&
                mask->DistributionMap() == phi[lev]->DistributionMap(),
                "Fixed anode potentials and mask must share the phi layout.");
            for (amrex::MFIter mfi(*phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                auto const p = phi[lev]->array(mfi);
                auto const v = fixed->const_array(mfi);
                auto const m = mask->const_array(mfi);
                amrex::ParallelFor(mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        if (m(i,j,k) == 0) { p(i,j,k) = v(i,j,k); }
                    });
            }
        }
#endif
        // Solve Poisson equation at lev
        mlmg.solve( {phi[lev]}, {rho[lev]},
                     mlmg_options.relative_tolerance, absolute_tolerance );

        const amrex::IntVect& refratio = rel_ref_ratio.value()[lev];
        const int ncomp = linop->getNComp();

        // needed for solving the levels by levels:
        // - coarser level is initial guess for finer level
        // - coarser level provides boundary values for finer level patch
        // Interpolation from phi[lev] to phi[lev+1]
        // (This provides both the boundary conditions and initial guess for phi[lev+1])
        if (lev < finest_level) {
            interpolatePhiBetweenLevels(phi[lev],
                                        phi[lev+1],
                                        geom[lev],
                                        do_single_precision_comms,
                                        refratio,
                                        ncomp,
                                        ng);
        }

        // Run additional operations, such as calculation of the E field for embedded boundaries
        if constexpr (!std::is_same_v<T_PostPhiCalculationFunctor, std::nullopt_t>) {
            if (post_phi_calculation.has_value()) {
                post_phi_calculation.value()(mlmg, lev);
            }
        }
        rho[lev]->mult(-ablastr::constant::SI::epsilon_0);  // Multiply rho by epsilon again

    } // loop over lev(els)
} // computePhi
} // namespace ablastr::fields

#endif // ABLASTR_POISSON_SOLVER_H
