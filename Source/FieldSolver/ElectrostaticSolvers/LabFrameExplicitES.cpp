/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald, Arianna Formenti, Revathi Jambunathan
 *
 * License: BSD-3-Clause-LBNL
 */
#include "LabFrameExplicitES.H"
#include "EmbeddedBoundary/Enabled.H"
#include "Fields.H"
#include "Insert/Boundary/InsertBoundaryPhi.h"
#include "Insert/Config/WarpXFunctionConfig.h"
#include "Insert/Config/WarpXSimulationConfig.h"
#include "Insert/Core/WarpXInsert.h"
#include "Insert/Fields/ECDIChargeFilter.h"
#include "Python/callbacks.H"
#include "WarpX.H"
#ifdef WARPX_USE_ELECTROSTATIC_MATERIALS
#include "Insert/Fields/ElectrostaticMaterial.H"
#include "Insert/Fields/WallCharge.H"
#include "Insert/Fields/MaterialPoissonSolver.H"
#endif

#include <AMReX_ParmParse.H>

#include <algorithm>
#include <string>

using namespace amrex;

void LabFrameExplicitES::InitData() {
    auto & warpx = WarpX::GetInstance();
#ifdef WARPX_USE_ELECTROSTATIC_MATERIALS
    m_material_solver.readParameters(warpx.Geom(0));
    if (Insert::ElectrostaticMaterial::GetInstance().enabled()) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(warpx.Geom(0).Coord() == 0 && !warpx.Geom(0).isAnyPeriodic(),
            "The material solver requires nonperiodic Cartesian geometry.");
    }
    if (m_material_solver.isMonopole()) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrame &&
            warpx.evolve_scheme == EvolveScheme::Explicit,
            "Monopole Robin currently requires explicit lab-frame electrostatic evolution.");
        // The material backend owns the physical Robin boundary. The common
        // hierarchy and history still receive an explicitly initialized handler.
        m_poisson_boundary_handler->lobc.fill(amrex::LinOpBCType::Robin);
        m_poisson_boundary_handler->hibc.fill(amrex::LinOpBCType::Robin);
        m_poisson_boundary_handler->dirichlet_flag.fill(false);
        m_poisson_boundary_handler->has_non_periodic = true;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(Insert::ElectrostaticMaterial::GetInstance().enabled(),
            "Monopole Robin requires insert.use_electrostatic_materials = 1.");
        for (int d = 0; d < 3; ++d) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                WarpX::field_boundary_lo[d] == FieldBoundaryType::None &&
                WarpX::field_boundary_hi[d] == FieldBoundaryType::None,
                "Monopole Robin requires boundary.field_lo/hi = none none none; "
                "the Robin solver owns the electrostatic outer boundary.");
        }
        bool ecdi_enabled = false;
        amrex::ParmParse("insert.ecdi_control").query("enabled", ecdi_enabled);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(WarpX::nox == 1 && !WarpX::use_filter && !ecdi_enabled,
            "Monopole Robin currently requires algo.particle_shape = 1, "
            "warpx.use_filter = 0 and insert.ecdi_control.enabled = 0.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            WarpX::field_gathering_algo == GatheringAlgo::EnergyConserving,
            "Monopole Robin requires algo.field_gathering = energy-conserving; "
            "momentum-conserving gathering is not supported.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            WarpX::grid_type == ablastr::utils::enums::GridType::Collocated ||
            WarpX::galerkin_interpolation,
            "Staggered monopole Robin requires interpolation.galerkin_scheme = 1.");
        return;
    }
#else
    std::string boundary = "dirichlet";
    amrex::ParmParse("insert").query("electrostatic_boundary", boundary);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(boundary == "dirichlet",
        "insert.electrostatic_boundary = monopole requires 3D.");
#endif
    m_poisson_boundary_handler->DefinePhiBCs(warpx.Geom(0));
#ifdef WARPX_USE_ELECTROSTATIC_MATERIALS
    if (Insert::ElectrostaticMaterial::GetInstance().enabled()) {
        for (int d = 0; d < 3; ++d) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_poisson_boundary_handler->lobc[d] == amrex::LinOpBCType::Dirichlet &&
                m_poisson_boundary_handler->hibc[d] == amrex::LinOpBCType::Dirichlet,
                "Material Dirichlet mode requires Dirichlet conditions on all six faces.");
        }
    }
#endif
}

int LabFrameExplicitES::phiExtrapolationCacheDepth () const
{
    return self_fields_phi_extrapolation_order + 1;
}

void LabFrameExplicitES::definePhiExtrapolationCache (
    const ablastr::fields::MultiLevelScalarField& phi)
{
    bool reset_history = false;
    int const cache_depth = phiExtrapolationCacheDepth();

    for (int hist = 0; hist < cache_depth; ++hist) {
        if (m_phi_extrapolation_cache[hist].size() != phi.size()) {
            m_phi_extrapolation_cache[hist].resize(phi.size());
            reset_history = true;
        }
    }
    for (int hist = cache_depth;
         hist < static_cast<int>(m_phi_extrapolation_cache.size()); ++hist)
    {
        m_phi_extrapolation_cache[hist].clear();
    }

    for (int lev = 0; lev < static_cast<int>(phi.size()); ++lev) {
        auto const& phi_mf = *phi[lev];
        for (int hist = 0; hist < cache_depth; ++hist) {
            auto& cache_mf = m_phi_extrapolation_cache[hist][lev];
            bool const cache_needs_define = (
                !cache_mf ||
                cache_mf->nComp() != phi_mf.nComp() ||
                cache_mf->nGrowVect() != phi_mf.nGrowVect() ||
                !cache_mf->boxArray().CellEqual(phi_mf.boxArray()) ||
                !(cache_mf->DistributionMap() == phi_mf.DistributionMap())
            );

            if (!cache_needs_define) {
                continue;
            }

            cache_mf = std::make_unique<amrex::MultiFab>(
                phi_mf.boxArray(), phi_mf.DistributionMap(), phi_mf.nComp(),
                phi_mf.nGrowVect());
            cache_mf->setVal(0.0_rt);
            reset_history = true;
        }
    }

    if (reset_history) {
        m_phi_extrapolation_history_depth = 0;
    } else if (m_phi_extrapolation_history_depth > cache_depth) {
        m_phi_extrapolation_history_depth = cache_depth;
    }
}

void LabFrameExplicitES::preparePhiExtrapolatedInitialGuess (
    const ablastr::fields::MultiLevelScalarField& phi)
{
    definePhiExtrapolationCache(phi);

    if (m_phi_extrapolation_history_depth == 0) {
        return;
    }

    int const extrapolation_order =
        std::min(self_fields_phi_extrapolation_order,
                 m_phi_extrapolation_history_depth - 1);

    for (int lev = 0; lev < static_cast<int>(phi.size()); ++lev) {
        auto& phi_mf = *phi[lev];
        auto& phi_n = *m_phi_extrapolation_cache[0][lev];
        const int ncomp = phi_mf.nComp();
        const int ngrow = phi_mf.nGrow();

        amrex::MultiFab::Copy(phi_mf, phi_n, 0, 0, ncomp, ngrow);

        if (extrapolation_order == 2 &&
            (self_fields_phi_extrapolation_alpha != 0.0_rt ||
             self_fields_phi_extrapolation_beta != 0.0_rt))
        {
            auto& phi_nm1 = *m_phi_extrapolation_cache[1][lev];
            auto& phi_nm2 = *m_phi_extrapolation_cache[2][lev];
            phi_mf.mult(
                1.0_rt + self_fields_phi_extrapolation_alpha +
                self_fields_phi_extrapolation_beta,
                0, ncomp, ngrow);
            amrex::MultiFab::Saxpy(
                phi_mf,
                -(self_fields_phi_extrapolation_alpha +
                  2.0_rt * self_fields_phi_extrapolation_beta),
                phi_nm1, 0, 0, ncomp, ngrow);
            amrex::MultiFab::Saxpy(
                phi_mf, self_fields_phi_extrapolation_beta, phi_nm2, 0, 0,
                ncomp, ngrow);
        } else if (extrapolation_order == 1 &&
                   self_fields_phi_extrapolation_alpha != 0.0_rt)
        {
            auto& phi_nm1 = *m_phi_extrapolation_cache[1][lev];
            phi_mf.mult(
                1.0_rt + self_fields_phi_extrapolation_alpha, 0, ncomp, ngrow);
            amrex::MultiFab::Saxpy(
                phi_mf, -self_fields_phi_extrapolation_alpha, phi_nm1, 0, 0,
                ncomp, ngrow);
        }
    }
}

void LabFrameExplicitES::updatePhiExtrapolationHistory (
    const ablastr::fields::MultiLevelScalarField& phi)
{
    definePhiExtrapolationCache(phi);
    int const cache_depth = phiExtrapolationCacheDepth();
    int const shift_count = std::min(m_phi_extrapolation_history_depth,
                                     cache_depth - 1);

    for (int lev = 0; lev < static_cast<int>(phi.size()); ++lev) {
        auto const& phi_mf = *phi[lev];
        const int ncomp = phi_mf.nComp();
        const int ngrow = phi_mf.nGrow();

        for (int hist = shift_count; hist > 0; --hist) {
            amrex::MultiFab::Copy(
                *m_phi_extrapolation_cache[hist][lev],
                *m_phi_extrapolation_cache[hist - 1][lev], 0, 0, ncomp,
                ngrow);
        }
        amrex::MultiFab::Copy(
            *m_phi_extrapolation_cache[0][lev], phi_mf, 0, 0, ncomp, ngrow);
    }

    if (m_phi_extrapolation_history_depth < cache_depth) {
        ++m_phi_extrapolation_history_depth;
    }
}

void LabFrameExplicitES::ComputeSpaceChargeField (
    ablastr::fields::MultiFabRegister& fields,
    MultiParticleContainer& mpc,
    MultiFluidContainer* mfl,
    int max_level,
    bool verbose_step)
{
    using ablastr::fields::MultiLevelScalarField;
    using ablastr::fields::MultiLevelVectorField;
    using warpx::fields::FieldType;

    bool const skip_lev0_coarse_patch = true;
    bool const use_robin =
#ifdef WARPX_USE_ELECTROSTATIC_MATERIALS
        m_material_solver.isMonopole();
#else
        false;
#endif

    const MultiLevelScalarField rho_fp = fields.get_mr_levels(FieldType::rho_fp, max_level);
    const MultiLevelScalarField rho_cp = fields.get_mr_levels(FieldType::rho_cp, max_level, skip_lev0_coarse_patch);
    const MultiLevelScalarField phi_fp = fields.get_mr_levels(FieldType::phi_fp, max_level);
    const MultiLevelVectorField Efield_fp = fields.get_mr_levels_alldirs(FieldType::Efield_fp, max_level);

    ExecutePythonCallback("beforedeposition");
    mpc.DepositCharge(rho_fp, 0.0_rt);
    if (mfl) {
        const int lev = 0;
        mfl->DepositCharge(fields, *rho_fp[lev], lev);
    }
    ExecutePythonCallback("afterdeposition");

#ifdef WARPX_USE_ELECTROSTATIC_MATERIALS
    auto& electrostatic_material = Insert::ElectrostaticMaterial::GetInstance();
    if (electrostatic_material.enabled() && max_level > 0) {
        // Add raw persistent contributions once. Only this working rho is synchronized.
        Insert::WallCharge::GetInstance().addTo(rho_fp, max_level);
    }
#endif

    // Apply filter, perform MPI exchange, interpolate across levels
    const Vector<std::unique_ptr<MultiFab> > rho_buf(num_levels);
    auto & warpx = WarpX::GetInstance();
    warpx.SyncRho( rho_fp, rho_cp, amrex::GetVecOfPtrs(rho_buf) );

#ifndef WARPX_DIM_RZ
    for (int lev = 0; lev < num_levels && !use_robin; lev++) {
        // Reflect density over PEC boundaries, if needed.
        warpx.ApplyRhofieldBoundary(lev, rho_fp[lev], PatchType::fine);
    }
#endif

    Insert::FilterRhoForECDIControl(rho_fp, max_level);

    // beta is zero in lab frame
    // Todo: use simpler finite difference form with beta=0
    const std::array<Real, 3> beta = {0._rt};

    bool const use_level_history =
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_XZ)
        !EB::enabled() && !IsPythonCallbackInstalled("poissonsolver") &&
        WarpX::poisson_solver_id == PoissonSolverAlgo::Multigrid;
#else
        false;
#endif
    if (!use_level_history) { preparePhiExtrapolatedInitialGuess(phi_fp); }

    // set the boundary potentials appropriately
    if (!use_robin) { setPhiBC(phi_fp, warpx.gett_new(0)); }
#ifdef WARPX_USE_ELECTROSTATIC_MATERIALS
    std::optional<amrex::Vector<amrex::iMultiFab const*>> anode_masks = std::nullopt;
    ablastr::fields::ConstMultiLevelScalarField const* relative_permittivity = nullptr;
    if (electrostatic_material.enabled()) {
        if (max_level == 0) {
            // Preserve the single-level source ordering: wall charge bypasses
            // particle filtering and PEC image-charge processing.
            Insert::WallCharge::GetInstance().addTo(rho_fp, max_level);
        }
        electrostatic_material.prepare(rho_fp, phi_fp, max_level);
        anode_masks = *electrostatic_material.anodeMasks();
        relative_permittivity = electrostatic_material.relativePermittivity();
    }
#endif

    // Set the zmin-face ring anode potential, if requested (skipped when the
    // electrostatic material path owns the anode).
    Insert::SetBoundaryPhi();

    if (use_level_history) {
        ablastr::fields::NodalPoissonInitialGuess::BC lobc{}, hibc{};
        if (!use_robin) {
            lobc = m_poisson_boundary_handler->lobc;
            hibc = m_poisson_boundary_handler->hibc;
        }
        amrex::Vector<amrex::iMultiFab const*> fixed_masks;
#ifdef WARPX_USE_ELECTROSTATIC_MATERIALS
        if (anode_masks) { fixed_masks = *anode_masks; }
        if (use_robin) {
            // Monopole owns the physical boundary; its face nodes are unknowns.
            lobc.fill(amrex::LinOpBCType::Robin);
            hibc.fill(amrex::LinOpBCType::Robin);
        }
#endif
        m_nodal_initial_guess.beginSolve(
            phi_fp, warpx.Geom(), lobc, hibc, fixed_masks, warpx.gett_new(0),
            self_fields_phi_extrapolation_order,
            self_fields_phi_extrapolation_alpha, self_fields_phi_extrapolation_beta);
    }

    // Compute the potential phi, by solving the Poisson equation
    if (IsPythonCallbackInstalled("poissonsolver")) {

        // Use the Python level solver (user specified)
        ExecutePythonCallback("poissonsolver");

    } else {

#if defined(WARPX_DIM_1D_Z)
        // Use the tridiag solver with 1D
        amrex::ignore_unused(verbose_step);
        computePhiTriDiagonal(rho_fp, phi_fp);
#else
        // Use the AMREX MLMG or the FFT (IGF) solver otherwise
        ablastr::fields::MLMGOptions mlmg_options = m_mlmg_options;
        mlmg_options.verbosity = verbose_step ? m_mlmg_options.verbosity : 0;
        ablastr::fields::PoissonLevelPrepare prepare_level;
        if (use_level_history) {
            prepare_level = [this, verbosity = mlmg_options.verbosity]
                (int lev, amrex::MultiFab& phi) {
                    m_nodal_initial_guess.prepare(lev, phi, verbosity);
                };
        }
        ablastr::fields::PoissonLevelSolver level_solver;
#ifdef WARPX_USE_ELECTROSTATIC_MATERIALS
        if (relative_permittivity) {
            level_solver = [&](int lev, amrex::LPInfo const&, amrex::MultiFab& phi,
                               amrex::MultiFab& rhs,
                               ablastr::fields::MLMGOptions const& options,
                               amrex::Real absolute_tolerance) -> amrex::MLMG& {
                return m_material_solver.solveLevel(lev, rhs, phi,
                    *relative_permittivity->at(lev), *anode_masks->at(lev),
                    *electrostatic_material.anodePotentials()->at(lev), warpx.Geom(lev),
                    options, absolute_tolerance);
            };
        }
#endif
        computePhi(rho_fp, phi_fp, beta, mlmg_options, is_igf_2d_slices,
                   Efield_fp, std::nullopt, level_solver, prepare_level);
#endif
    }
    Insert::FillPhiGhosts(phi_fp
#ifdef WARPX_USE_ELECTROSTATIC_MATERIALS
                          , &m_material_solver
#endif
    );
    // Keep extrapolation history on the uncorrected Poisson potential.
    if (use_level_history) {
        m_nodal_initial_guess.commit(phi_fp);
    } else {
        updatePhiExtrapolationHistory(phi_fp);
    }
#if defined(WARPX_DIM_XZ) && defined(BENCHMARK_2D)
    Insert::VoltageAdjustment();
#endif

    // Compute the electric field. Note that if an EB is used the electric
    // field will be calculated in the computePhi call.
    if (!EB::enabled() || IsPythonCallbackInstalled("poissonsolver")) {
        computeE(Efield_fp, phi_fp, beta);
    }
#ifdef WARPX_USE_ELECTROSTATIC_MATERIALS
    if (use_robin) {
        for (auto const& level : Efield_fp) {
            for (auto* field : level) { field->FillBoundary(); }
        }
    } else {
        electrostatic_material.fillElectricFieldGhosts(Efield_fp);
    }
#endif
}

/* \brief Compute the potential by solving Poisson's equation with
          a 1D tridiagonal solve.

   \param[in] rho The charge density a given species
   \param[out] phi The potential to be computed by this function
*/
void LabFrameExplicitES::computePhiTriDiagonal (
    const ablastr::fields::MultiLevelScalarField& rho,
    const ablastr::fields::MultiLevelScalarField& phi)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(num_levels == 1,
    "The tridiagonal solver cannot be used with mesh refinement");

    const int lev = 0;
    auto & warpx = WarpX::GetInstance();

    const amrex::Real* dx = warpx.Geom(lev).CellSize();
    const amrex::Real xmin = warpx.Geom(lev).ProbLo(0);
    const amrex::Real xmax = warpx.Geom(lev).ProbHi(0);
    const int nx_full_domain = static_cast<int>( (xmax - xmin)/dx[0] + 0.5_rt );

    int nx_solve_min = 1;
    int nx_solve_max = nx_full_domain - 1;

    auto field_boundary_lo0 = WarpX::field_boundary_lo[0];
    auto field_boundary_hi0 = WarpX::field_boundary_hi[0];
    if (field_boundary_lo0 == FieldBoundaryType::Neumann || field_boundary_lo0 == FieldBoundaryType::Periodic) {
        // Neumann or periodic boundary condition
        // Solve for the point on the lower boundary
        nx_solve_min = 0;
    }
    if (field_boundary_hi0 == FieldBoundaryType::Neumann || field_boundary_hi0 == FieldBoundaryType::Periodic) {
        // Neumann or periodic boundary condition
        // Solve for the point on the upper boundary
        nx_solve_max = nx_full_domain;
    }

    // Create a 1-D MultiFab that covers all of x.
    // The tridiag solve will be done in this MultiFab and then copied out afterwards.
    const amrex::IntVect lo_full_domain(AMREX_D_DECL(0,0,0));
    const amrex::IntVect hi_full_domain(AMREX_D_DECL(nx_full_domain,0,0));
    const amrex::Box box_full_domain_node(lo_full_domain, hi_full_domain, amrex::IntVect::TheNodeVector());
    const BoxArray ba_full_domain_node(box_full_domain_node);
    const amrex::Vector<int> pmap = {0}; // The data will only be on processor 0
    const amrex::DistributionMapping dm_full_domain(pmap);

    // Put the data in the pinned arena since the tridiag solver will be done on the CPU, but have
    // the data readily accessible from the GPU.
    auto phi1d_mf = MultiFab(ba_full_domain_node, dm_full_domain, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    auto zwork1d_mf = MultiFab(ba_full_domain_node, dm_full_domain, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    auto rho1d_mf = MultiFab(ba_full_domain_node, dm_full_domain, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));

    if (field_boundary_lo0 == FieldBoundaryType::PEC || field_boundary_hi0 == FieldBoundaryType::PEC) {
        // Copy from phi to get the boundary values
        phi1d_mf.ParallelCopy(*phi[lev], 0, 0, 1);
    }
    rho1d_mf.ParallelCopy(*rho[lev], 0, 0, 1);

    // Multiplier on the charge density
    const amrex::Real norm = dx[0]*dx[0]/PhysConst::epsilon_0;
    rho1d_mf.mult(norm);

    // Use the MFIter loop since when parallel, only process zero has a FAB.
    // This skips the loop on all other processors.
    for (MFIter mfi(phi1d_mf); mfi.isValid(); ++mfi) {

        const auto& phi1d_arr = phi1d_mf[mfi].array();
        const auto& zwork1d_arr = zwork1d_mf[mfi].array();
        const auto& rho1d_arr = rho1d_mf[mfi].array();

        // The loops are always performed on the CPU

        amrex::Real diag = 2._rt;

        // The initial values depend on the boundary condition
        if (field_boundary_lo0 == FieldBoundaryType::PEC) {

            phi1d_arr(1,0,0) = (phi1d_arr(0,0,0) + rho1d_arr(1,0,0))/diag;

        } else if (field_boundary_lo0 == FieldBoundaryType::Neumann) {

            // Neumann boundary condition
            phi1d_arr(0,0,0) = rho1d_arr(0,0,0)/diag;

            zwork1d_arr(1,0,0) = 2._rt/diag;
            diag = 2._rt - zwork1d_arr(1,0,0);
            phi1d_arr(1,0,0) = (rho1d_arr(1,0,0) - (-1._rt)*phi1d_arr(1-1,0,0))/diag;

        } else if (field_boundary_lo0 == FieldBoundaryType::Periodic) {

            phi1d_arr(0,0,0) = rho1d_arr(0,0,0)/diag;

            zwork1d_arr(1,0,0) = 1._rt/diag;
            diag = 2._rt - zwork1d_arr(1,0,0);
            phi1d_arr(1,0,0) = (rho1d_arr(1,0,0) - (-1._rt)*phi1d_arr(1-1,0,0))/diag;

        }

        // Loop upward, calculating the Gaussian elimination multipliers and right hand sides
        for (int i_up = 2 ; i_up < nx_solve_max ; i_up++) {

            zwork1d_arr(i_up,0,0) = 1._rt/diag;
            diag = 2._rt - zwork1d_arr(i_up,0,0);
            phi1d_arr(i_up,0,0) = (rho1d_arr(i_up,0,0) - (-1._rt)*phi1d_arr(i_up-1,0,0))/diag;

        }

        // The last value depend on the boundary condition
        amrex::Real zwork_product = 1.; // Needed for parallel boundaries
        if (field_boundary_hi0 == FieldBoundaryType::PEC) {

            int const nxm1 = nx_full_domain - 1;
            zwork1d_arr(nxm1,0,0) = 1._rt/diag;
            diag = 2._rt - zwork1d_arr(nxm1,0,0);
            phi1d_arr(nxm1,0,0) = (phi1d_arr(nxm1+1,0,0) + rho1d_arr(nxm1,0,0) - (-1._rt)*phi1d_arr(nxm1-1,0,0))/diag;

        } else if (field_boundary_hi0 == FieldBoundaryType::Neumann) {

            // Neumann boundary condition
            zwork1d_arr(nx_full_domain,0,0) = 1._rt/diag;
            diag = 2._rt - 2._rt*zwork1d_arr(nx_full_domain,0,0);
            if (diag == 0._rt) {
                // This happens if the lower boundary is also Neumann.
                // It this case, the potential is relative to an arbitrary constant,
                // so set the upper boundary to zero to force a value.
                phi1d_arr(nx_full_domain,0,0) = 0.;
            } else {
                phi1d_arr(nx_full_domain,0,0) = (rho1d_arr(nx_full_domain,0,0) - (-2._rt)*phi1d_arr(nx_full_domain-1,0,0))/diag;
            }

        }


        for (int i_down = nx_solve_max-1 ; i_down >= nx_solve_min ; i_down--) {
            phi1d_arr(i_down,0,0) = phi1d_arr(i_down,0,0) + zwork1d_arr(i_down+1,0,0)*phi1d_arr(i_down+1,0,0);
        }

    }

    // Copy phi1d to phi
    phi[lev]->ParallelCopy(phi1d_mf, 0, 0, 1);
}

/* \brief Compute the potential by solving Poisson's equation
 *        with periodic boundaries using
          a 1D tridiagonal solve.
          This makes use of the Sherman–Morrison formula.
          The code is based on the code given in the Wikipedia page,
          https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm#Variants.

   \param[in] rho The charge density a given species
   \param[out] phi The potential to be computed by this function
*/
void LabFrameExplicitES::computePhiTriDiagonal_periodic (
    const ablastr::fields::MultiLevelScalarField& rho,
    const ablastr::fields::MultiLevelScalarField& phi)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(num_levels == 1,
    "The tridiagonal solver cannot be used with mesh refinement");

    const int lev = 0;
    auto & warpx = WarpX::GetInstance();

    const amrex::Real* dx = warpx.Geom(lev).CellSize();
    const amrex::Real xmin = warpx.Geom(lev).ProbLo(0);
    const amrex::Real xmax = warpx.Geom(lev).ProbHi(0);
    const int nx = static_cast<int>( (xmax - xmin)/dx[0] + 0.5_rt );


    // Create a 1-D MultiFab that covers all of x.
    // The tridiag solve will be done in this MultiFab and then copied out afterwards.
    const amrex::IntVect lo_full_domain(AMREX_D_DECL(0,0,0));
    const amrex::IntVect hi_full_domain(AMREX_D_DECL(nx,0,0));
    const amrex::Box box_full_domain_node(lo_full_domain, hi_full_domain, amrex::IntVect::TheNodeVector());
    const BoxArray ba_full_domain_node(box_full_domain_node);
    const amrex::Vector<int> pmap = {0}; // The data will only be on processor 0
    const amrex::DistributionMapping dm_full_domain(pmap);

    // Put the data in the pinned arena since the tridiag solver will be done on the CPU, but have
    // the data readily accessible from the GPU.
    auto phi1d_mf = MultiFab(ba_full_domain_node, dm_full_domain, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));

    // Work arrays
    auto cmod_mf = MultiFab(ba_full_domain_node, dm_full_domain, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    auto u_mf = MultiFab(ba_full_domain_node, dm_full_domain, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));

    // Copy rho into phi1d_mf to start
    phi1d_mf.ParallelCopy(*rho[lev], 0, 0, 1);

    // Multiplier on the charge density
    const amrex::Real norm = dx[0]*dx[0]/PhysConst::epsilon_0;
    phi1d_mf.mult(norm);

    // Use the MFIter loop since when parallel, only process zero has a FAB.
    // This skips the loop on all other processors.
    for (MFIter mfi(phi1d_mf); mfi.isValid(); ++mfi) {

        const auto& x = phi1d_mf[mfi].array();
        const auto& cmod = cmod_mf[mfi].array();
        const auto& u = u_mf[mfi].array();

        // The loops are always performed on the CPU

        {

        // This code is adapted from the Wikipedia page on the Tridiagonal matrix algorithm
        // https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm#Variants.
        // Licensed under CC BY-SA 4.0
        // Modifications: The a, b, and c inputs are replaced with the fixed values.

        const amrex::Real alpha = -1.0_rt;
        const amrex::Real beta = -1.0_rt;

        /* arbitrary, but chosen such that division by zero is avoided */
        const amrex::Real gamma = -2.0_rt;

        cmod(0,0,0) = -1.0_rt / (2.0_rt - gamma);
        u(0,0,0) = gamma / (2.0_rt - gamma);
        x(0,0,0) /= (2.0_rt - gamma);

        /* loop from 1 to nx - 2 inclusive */
        for (int ix = 1; ix + 1 < nx; ix++) {
            const amrex::Real m = 1.0_rt / (2.0_rt - -1.0_rt * cmod(ix - 1,0,0));
            cmod(ix,0,0) = -1.0_rt * m;
            u(ix,0,0) = (0.0f  - -1.0_rt * u(ix - 1,0,0)) * m;
            x(ix,0,0) = (x(ix,0,0) - -1.0_rt * x(ix - 1,0,0)) * m;
        }

        /* handle nx - 1 */
        const amrex::Real m = 1.0_rt / (2.0_rt - alpha * beta / gamma - -1.0_rt * cmod(nx - 2,0,0));
        u(nx - 1,0,0) = (alpha    - -1.0_rt * u(nx - 2,0,0)) * m;
        x(nx - 1,0,0) = (x(nx - 1,0,0) - -1.0_rt * x(nx - 2,0,0)) * m;

        /* loop from nx - 2 to 0 inclusive */
        for (int ix = nx - 2; ix >= 0; ix--) {
            u(ix,0,0) -= cmod(ix,0,0) * u(ix + 1,0,0);
            x(ix,0,0) -= cmod(ix,0,0) * x(ix + 1,0,0);
        }

        const amrex::Real fact = (x(0,0,0) + x(nx - 1,0,0) * alpha / gamma) / (1.0_rt + u(0,0,0) + u(nx - 1,0,0) * alpha / gamma);

        /* loop from 0 to nx - 1 inclusive */
        for (int ix = 0; ix < nx; ix++) {
            x(ix,0,0) -= fact * u(ix,0,0);
        }

        }

        x(nx,0,0) = x(0,0,0);

        // In a test case, this was giving an relative residual of around 1.e-10.
        // A dozen or so SOR iterations could improve that by a factor of 10.
        // Is it worth it?

    }

    // Copy phi1d to phi
    phi[lev]->ParallelCopy(phi1d_mf, 0, 0, 1);
}
