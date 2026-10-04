#include "HallFarFieldPoisson.H"

#ifdef WARPX_USE_HALL_FAR_FIELD

#include "EmbeddedBoundary/Enabled.H"
#include "Insert/Fields/HallElectrostaticMaterial.H"
#include "Python/callbacks.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <ablastr/coarsen/average.H>
#include <ablastr/fields/Interpolate.H>
#include <ablastr/fields/PoissonSolver.H>
#include <ablastr/utils/Communication.H>

#include <AMReX_CoordSys.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>
#include <AMReX_Reduce.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <utility>

namespace Insert {
namespace {

/** Convert a physical extension to whole coarse cells without silently rounding it. */
int ExtensionCells (amrex::Real extension, amrex::Real spacing)
{
    amrex::Real const cells = extension / spacing;
    amrex::Real const nearest = std::round(cells);
    amrex::Real const tolerance = 64 * std::numeric_limits<amrex::Real>::epsilon()
        * std::max(amrex::Real(1), std::abs(cells));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(cells) && cells >= -tolerance &&
            std::abs(cells - nearest) <= tolerance &&
            nearest < std::numeric_limits<int>::max() / 4,
        "Hall far-field domain must contain the fine domain with extensions of whole "
        "coarse cells in each direction.");
    return static_cast<int>(nearest);
}

void WriteBox (std::ostream& os, amrex::Box const& box)
{
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        os << ' ' << box.smallEnd(d);
    }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        os << ' ' << box.bigEnd(d);
    }
}

// NVCC extended device lambdas cannot be defined in private member functions.
void FillCoarseSupport (amrex::MultiFab& support, amrex::iMultiFab const& anode_mask)
{
    for (amrex::MFIter mfi(support, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const mask = anode_mask.const_array(mfi);
        auto const values = support.array(mfi);
        amrex::ParallelFor(mfi.tilebox(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                values(i, j, k) = amrex::Real(mask(i, j, k));
            });
    }
}

void ApplyFixedPotential (
    amrex::MultiFab& phi, amrex::iMultiFab const& mask,
    amrex::MultiFab const& anode_phi, amrex::Geometry const& geom)
{
    auto const lo = amrex::lbound(amrex::surroundingNodes(geom.Domain()));
    auto const hi = amrex::ubound(amrex::surroundingNodes(geom.Domain()));
    for (amrex::MFIter mfi(phi, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const values = phi.array(mfi);
        auto const active = mask.const_array(mfi);
        auto const fixed = anode_phi.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (i == lo.x || i == hi.x || j == lo.y || j == hi.y ||
                    k == lo.z || k == hi.z) {
                    values(i, j, k) = 0;
                } else if (active(i, j, k) == 0) {
                    values(i, j, k) = fixed(i, j, k);
                }
            });
    }
}

void CheckGroundedAnode (
    amrex::MultiFab const& anode_phi, amrex::Geometry const& geom)
{
    auto const lo = amrex::lbound(amrex::surroundingNodes(geom.Domain()));
    auto const hi = amrex::ubound(amrex::surroundingNodes(geom.Domain()));
    amrex::ReduceOps<amrex::ReduceOpMax> ops;
    amrex::ReduceData<amrex::Real> data(ops);
    using Tuple = decltype(data)::Type;
    for (amrex::MFIter mfi(anode_phi, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const fixed = anode_phi.const_array(mfi);
        ops.eval(mfi.tilebox(), data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept -> Tuple
            {
                bool const boundary = i == lo.x || i == hi.x || j == lo.y || j == hi.y ||
                    k == lo.z || k == hi.z;
                return {boundary ? std::abs(fixed(i, j, k)) : amrex::Real(0)};
            });
    }
    amrex::Real conflict = amrex::get<0>(data.value());
    amrex::ParallelDescriptor::ReduceRealMax(conflict);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(conflict == 0,
        "A nonzero Hall anode potential touches the grounded coarse-domain boundary. "
        "Extend insert.far_field.prob_lo/prob_hi or adjust the anode geometry.");
}

/** Check the coincident anode nodes retained by AMReX's nodal mask injection. */
int FirstUnresolvedAnodeLevel (
    amrex::MLNodeLaplacian const& linop, amrex::iMultiFab const& mask)
{
    for (int mg = 0; mg < linop.NMGLevels(0); ++mg) {
        auto const ratio = linop.Geom(0).Domain().length() / linop.Geom(0, mg).Domain().length();
        amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMax> ops;
        amrex::ReduceData<int, int> data(ops);
        using Tuple = decltype(data)::Type;
        constexpr int empty_min = std::numeric_limits<int>::max();
        constexpr int empty_max = std::numeric_limits<int>::lowest();
        for (amrex::MFIter mfi(mask, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto const values = mask.const_array(mfi);
            ops.eval(mfi.tilebox(), data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept -> Tuple
                {
                    bool const retained = values(i, j, k) == 0 &&
                        i % ratio[0] == 0 && j % ratio[1] == 0 && k % ratio[2] == 0;
                    return {retained ? k / ratio[2] : empty_min,
                            retained ? k / ratio[2] : empty_max};
                });
        }
        auto const result = data.value();
        int min_k = amrex::get<0>(result);
        int max_k = amrex::get<1>(result);
        amrex::ParallelDescriptor::ReduceIntMin(min_k);
        amrex::ParallelDescriptor::ReduceIntMax(max_k);
        if (min_k >= max_k) {
            return mg;
        }
    }
    return -1;
}

} // namespace

void HallFarFieldPoisson::initialize (
    WarpX const& warpx, amrex::MultiFab const& phi, amrex::MultiFab const& rho,
    SolverSettings const& defaults)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!m_defined,
        "Hall far-field grid has already been initialized.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(warpx.maxLevel() == 0 && warpx.finestLevel() == 0,
        "Hall far-field grid requires amr.max_level = 0.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrame ||
            WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrameElectroMagnetostatic,
        "Hall far-field grid requires a lab-frame electrostatic solver.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        WarpX::poisson_solver_id == PoissonSolverAlgo::Multigrid && !EB::enabled() &&
            !IsPythonCallbackInstalled("poissonsolver"),
        "Hall far-field grid requires MLMG without EB or a Python poissonsolver callback.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(HallElectrostaticMaterial::GetInstance().enabled(),
        "Hall far-field grid requires insert.use_electrostatic_materials = 1.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !warpx.getdo_moving_window() && !warpx.get_load_balance_intervals().isActivated(),
        "Hall far-field grid does not support moving windows or load balancing.");
    int regrid_interval = -1;
    utils::parser::queryWithParser(amrex::ParmParse("warpx"), "regrid_int", regrid_interval);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(regrid_interval <= 0,
        "Hall far-field grid does not support regridding.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        phi.ixType().nodeCentered() && rho.ixType().nodeCentered() &&
            phi.nComp() == 1 && rho.nComp() == 1 &&
            phi.boxArray() == rho.boxArray() &&
            phi.DistributionMap() == rho.DistributionMap(),
        "Hall far-field grid requires matching single-component nodal phi and rho.");

    m_fine_geom = warpx.Geom(0);
    m_fine_boxes = amrex::convert(phi.boxArray(), amrex::IntVect::TheZeroVector());
    m_fine_phi_ngrow = phi.nGrowVect();
    m_fine_rho_ngrow = rho.nGrowVect();
    readSolverSettings(defaults);
    defineGeometry(m_fine_geom);
    defineBoxes(m_fine_boxes, phi.DistributionMap());
    allocateFields(phi, rho, m_settings.m_extrapolation_order + 1);
    defineChargeTransfer(rho);
    defineSolver();
    m_defined = true;
    reportLayout();
}

void HallFarFieldPoisson::readSolverSettings (SolverSettings const& defaults)
{
    m_settings = defaults;
    amrex::ParmParse const pp("insert.far_field");
    utils::parser::queryWithParser(pp, "relative_tolerance", m_settings.m_relative_tolerance);
    utils::parser::queryWithParser(pp, "absolute_tolerance", m_settings.m_absolute_tolerance);
    utils::parser::queryWithParser(pp, "max_iters", m_settings.m_max_iters);
    utils::parser::queryWithParser(pp, "verbosity", m_settings.m_verbosity);
    utils::parser::queryWithParser(pp, "max_coarsening_level", m_settings.m_max_coarsening_level);
    utils::parser::queryWithParser(pp, "num_final_sweeps", m_settings.m_num_final_sweeps);
    utils::parser::queryWithParser(pp, "phi_extrapolation_order", m_settings.m_extrapolation_order);
    utils::parser::queryWithParser(pp, "phi_extrapolation_alpha", m_settings.m_extrapolation_alpha);
    utils::parser::queryWithParser(pp, "phi_extrapolation_beta", m_settings.m_extrapolation_beta);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(m_settings.m_relative_tolerance) &&
            std::isfinite(m_settings.m_absolute_tolerance) &&
            m_settings.m_relative_tolerance >= 0 && m_settings.m_absolute_tolerance >= 0 &&
            (m_settings.m_relative_tolerance > 0 || m_settings.m_absolute_tolerance > 0) &&
            m_settings.m_max_iters > 0 && m_settings.m_verbosity >= 0 &&
            m_settings.m_max_coarsening_level >= 0 && m_settings.m_num_final_sweeps > 0,
        "Hall far-field solver requires finite nonnegative tolerances (at least one positive), "
        "positive max_iters/num_final_sweeps and nonnegative verbosity/max_coarsening_level.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_settings.m_extrapolation_order >= 0 && m_settings.m_extrapolation_order <= 2 &&
            std::isfinite(m_settings.m_extrapolation_alpha) &&
            std::isfinite(m_settings.m_extrapolation_beta),
        "Hall far-field phi_extrapolation_order must be 0, 1 or 2, with finite alpha/beta.");
}

void HallFarFieldPoisson::defineGeometry (amrex::Geometry const& fine_geom)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        fine_geom.Coord() == amrex::CoordSys::cartesian && !fine_geom.isAnyPeriodic(),
        "Hall far-field grid requires nonperiodic Cartesian geometry.");
    amrex::ParmParse const pp("insert.far_field");
    amrex::Vector<amrex::Real> prob_lo;
    amrex::Vector<amrex::Real> prob_hi;
    amrex::Vector<int> ratio;
    utils::parser::getArrWithParser(pp, "prob_lo", prob_lo);
    utils::parser::getArrWithParser(pp, "prob_hi", prob_hi);
    utils::parser::getArrWithParser(pp, "coarsening_ratio", ratio);
    pp.query("layout_file", m_layout_file);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        prob_lo.size() == AMREX_SPACEDIM && prob_hi.size() == AMREX_SPACEDIM &&
            ratio.size() == AMREX_SPACEDIM,
        "insert.far_field.prob_lo, prob_hi and coarsening_ratio require three values each.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::all_of(ratio.begin(), ratio.end(), [](int r) { return r > 0; }) &&
            std::any_of(ratio.begin(), ratio.end(), [](int r) { return r > 1; }),
        "Hall far-field coarsening ratios must be positive, with at least one greater than one.");
    m_ratio = amrex::IntVect(ratio.data());
    amrex::Box const& fine_domain = fine_geom.Domain();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(fine_domain.coarsenable(m_ratio),
        "Hall fine-domain endpoints must align with the coarse mesh.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        fine_domain.length(0) == fine_domain.length(1) && ratio[0] == ratio[1],
        "Hall far-field grid requires equal x/y fine cell counts and coarsening ratios.");

    amrex::Box const embedded_domain = amrex::coarsen(fine_domain, m_ratio);
    amrex::Box coarse_domain = embedded_domain;
    amrex::IntVect high_extension;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(prob_lo[d]) && std::isfinite(prob_hi[d]) && prob_hi[d] > prob_lo[d],
            "Hall far-field physical bounds must be finite and strictly ordered.");
        amrex::Real const spacing = fine_geom.CellSize(d) * ratio[d];
        m_fine_offset[d] = ExtensionCells(fine_geom.ProbLo(d) - prob_lo[d], spacing);
        high_extension[d] = ExtensionCells(prob_hi[d] - fine_geom.ProbHi(d), spacing);
        auto const lo = static_cast<amrex::Long>(embedded_domain.smallEnd(d)) - m_fine_offset[d];
        auto const hi = static_cast<amrex::Long>(embedded_domain.bigEnd(d)) + high_extension[d];
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            lo > std::numeric_limits<int>::min() && hi < std::numeric_limits<int>::max() &&
                hi - lo + 1 < std::numeric_limits<int>::max(),
            "Hall far-field index extent exceeds the supported integer range.");
        coarse_domain.setSmall(d, static_cast<int>(lo));
        coarse_domain.setBig(d, static_cast<int>(hi));
        // Snap tolerated floating-point roundoff to the aligned physical bounds.
        prob_lo[d] = fine_geom.ProbLo(d) - m_fine_offset[d] * spacing;
        prob_hi[d] = fine_geom.ProbHi(d) + high_extension[d] * spacing;
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_fine_offset[0] == high_extension[0] && m_fine_offset[1] == high_extension[1] &&
            coarse_domain.length(0) == coarse_domain.length(1),
        "Hall fine grid must be centered in x/y with equal coarse x/y cell counts.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(coarse_domain != embedded_domain,
        "Hall far-field domain must extend beyond the fine domain.");

    // Retain natural coarse indices (possibly negative) so the existing WarpX
    // restriction/interpolation kernels still use fine_index = coarse_index * ratio.
    m_geom.define(coarse_domain, amrex::RealBox(prob_lo.data(), prob_hi.data()),
                  amrex::CoordSys::cartesian, {0, 0, 0});
}

void HallFarFieldPoisson::defineBoxes (
    amrex::BoxArray const& fine_boxes, amrex::DistributionMapping const& fine_dmap)
{
    amrex::Box const& domain = m_fine_geom.Domain();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        fine_boxes.size() > 0 && fine_boxes.size() == fine_dmap.size() &&
            fine_boxes.isDisjoint() && fine_boxes.contains(domain) &&
            fine_boxes.numPts() == domain.numPts(),
        "Hall fine Boxes must cover the domain exactly without overlaps.");

    bool const serial = amrex::ParallelDescriptor::NProcs() == 1;
    bool const full_xy = fine_boxes[0].length(0) == domain.length(0) &&
                         fine_boxes[0].length(1) == domain.length(1);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(serial || !full_xy,
        "Hall MPI layout requires x/y center-plane bisection.");
    amrex::Vector<amrex::Box> coarse_boxes(fine_boxes.size());
    for (int ibox = 0; ibox < fine_boxes.size(); ++ibox) {
        amrex::Box const fine_box = fine_boxes[ibox];
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            domain.contains(fine_box) && fine_box.coarsenable(m_ratio),
            "Hall fine Box boundaries must align with the coarse mesh.");
        for (int d = 0; d < 2; ++d) {
            int const half = domain.length(d) / 2;
            int const first = domain.smallEnd(d);
            bool const valid_full = serial && full_xy &&
                fine_box.smallEnd(d) == first && fine_box.length(d) == domain.length(d);
            bool const valid_half = !full_xy && domain.length(d) % 2 == 0 &&
                fine_box.length(d) == half &&
                (fine_box.smallEnd(d) == first || fine_box.smallEnd(d) == first + half);
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(valid_full || valid_half,
                "Hall fine Boxes must use x/y center-plane bisection; only serial runs "
                "may use full-width x/y Boxes. Other x/y subdivisions are unsupported.");
        }
        auto box = amrex::coarsen(fine_box, m_ratio);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (fine_box.smallEnd(d) == domain.smallEnd(d)) {
                box.setSmall(d, m_geom.Domain().smallEnd(d));
            }
            if (fine_box.bigEnd(d) == domain.bigEnd(d)) {
                box.setBig(d, m_geom.Domain().bigEnd(d));
            }
        }
        coarse_boxes[ibox] = box;
    }
    m_boxes = amrex::BoxArray(coarse_boxes.data(), static_cast<int>(coarse_boxes.size()));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_boxes.isDisjoint() && m_boxes.contains(m_geom.Domain()) &&
            m_boxes.numPts() == m_geom.Domain().numPts(),
        "Expanded Hall coarse Boxes must cover the coarse domain exactly without overlaps.");
    // Box numbering is unchanged: each overlapping coarse/fine pair has the same owner.
    m_dmap = fine_dmap;
}

void HallFarFieldPoisson::allocateFields (
    amrex::MultiFab const& phi, amrex::MultiFab const& rho, int const history_capacity)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(history_capacity >= 1 && history_capacity <= 3,
        "Hall far-field potential history requires one to three slots.");
    amrex::BoxArray const nodes = amrex::convert(m_boxes, amrex::IntVect::TheNodeVector());
    amrex::BoxArray const patch_nodes = amrex::coarsen(phi.boxArray(), m_ratio);
    amrex::IntVect const rho_ngrow = (rho.nGrowVect() + m_ratio - 1) / m_ratio;
    amrex::IntVect phi_ngrow = phi.nGrowVect();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        phi_ngrow[d] = std::max(phi_ngrow[d], 1);
    }
    m_rho.define(nodes, m_dmap, 1, rho_ngrow);
    m_phi.define(nodes, m_dmap, 1, phi_ngrow);
    m_phi_patch.define(patch_nodes, m_dmap, 1, phi_ngrow);
    m_rho.setVal(0);
    m_phi.setVal(0);
    m_phi_patch.setVal(0);
    for (int hist = 0; hist < history_capacity; ++hist) {
        m_phi_history[hist] = std::make_unique<amrex::MultiFab>(nodes, m_dmap, 1, phi_ngrow);
        m_phi_history[hist]->setVal(0);
    }
    m_history_depth = 0;
    m_history_times.fill(0);
}

void HallFarFieldPoisson::defineChargeTransfer (amrex::MultiFab const& rho)
{
    // Extend only physical faces, rounding outwards to complete coarse cells.
    // This makes formerly exterior deposition ghosts valid source nodes. Internal
    // Box overlaps are still only nodal overlaps and have a cached unique owner.
    auto const padding = ((rho.nGrowVect() + m_ratio - 1) / m_ratio) * m_ratio;
    amrex::Vector<amrex::Box> source_boxes(m_fine_boxes.size());
    for (int ibox = 0; ibox < m_fine_boxes.size(); ++ibox) {
        auto box = m_fine_boxes[ibox];
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (box.smallEnd(d) == m_fine_geom.Domain().smallEnd(d)) {
                box.growLo(d, padding[d]);
            }
            if (box.bigEnd(d) == m_fine_geom.Domain().bigEnd(d)) {
                box.growHi(d, padding[d]);
            }
        }
        source_boxes[ibox] = box;
    }
    amrex::BoxArray source_nodes(source_boxes.data(), static_cast<int>(source_boxes.size()));
    source_nodes.surroundingNodes();
    auto const patch_nodes = amrex::coarsen(source_nodes, m_ratio);
    auto const coarse_nodes = m_rho.boxArray();
    m_rho_source.define(source_nodes, m_dmap, 2, 0);
    m_rho_patch.define(patch_nodes, m_dmap, 2, 0);
    m_rho_raw.define(coarse_nodes, m_dmap, 2, 0);
    m_source_support.define(source_nodes, m_dmap, 1, 0);
    m_source_anode_mask.define(source_nodes, m_dmap, 1, 0);
    m_anode_mask.define(coarse_nodes, m_dmap, 1, 0);
    m_source_owner = m_rho_source.OwnerMask();
    m_coarse_owner = m_rho_raw.OwnerMask();
    m_rho_source.setVal(0);
    m_rho_patch.setVal(0);
    m_rho_raw.setVal(0);

    auto& material = HallElectrostaticMaterial::GetInstance();
    material.fillAnodeMask(m_source_anode_mask, m_fine_geom);
    material.fillAnodeMask(m_anode_mask, m_geom);

    // P(mask) > 0 iff at least one non-anode coarse recipient has nonzero
    // restriction weight. Use the existing interpolation kernel, not new weights.
    amrex::MultiFab coarse_support(coarse_nodes, m_dmap, 1, 0);
    amrex::MultiFab patch_support(patch_nodes, m_dmap, 1, 0);
    FillCoarseSupport(coarse_support, m_anode_mask);
    patch_support.setVal(0);
    ablastr::utils::communication::ParallelCopy(
        patch_support, coarse_support, 0, 0, 1,
        amrex::IntVect::TheZeroVector(), amrex::IntVect::TheZeroVector(), false);
    for (amrex::MFIter mfi(m_source_support, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        ablastr::fields::details::PoissonInterpCPtoFP const interp(
            m_source_support.array(mfi), patch_support.const_array(mfi), m_ratio);
        amrex::ParallelFor(mfi.tilebox(), interp);
    }
    amrex::Gpu::streamSynchronize();
}

void HallFarFieldPoisson::restrictCharge (amrex::MultiFab const& rho)
{
    m_charge_ready = false;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_defined && rho.boxArray().CellEqual(m_fine_boxes) &&
            rho.DistributionMap() == m_dmap && rho.nGrowVect() == m_fine_rho_ngrow &&
            rho.ixType().nodeCentered() && rho.nComp() == 1,
        "Hall charge transfer requires the initialized fine rho layout.");
    auto const zero = amrex::IntVect::TheZeroVector();
    // SyncRho has already summed deposition into valid AND exterior ghost nodes.
    // The second copy gives valid data priority over stale internal ghosts after
    // filtering, wall-charge addition and fine-anode clearing. Do not sum again.
    m_rho_source.setVal(0);
    ablastr::utils::communication::ParallelCopy(
        m_rho_source, rho, 0, 0, 1, rho.nGrowVect(), zero, false);
    ablastr::utils::communication::ParallelCopy(
        m_rho_source, rho, 0, 0, 1, zero, zero, false);
    m_rho_source.OverrideSync();

    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum> fine_ops;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real> fine_data(fine_ops);
    using FineTuple = decltype(fine_data)::Type;
    for (amrex::MFIter mfi(m_rho_source, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const source = m_rho_source.array(mfi);
        auto const owner = m_source_owner->const_array(mfi);
        auto const mask = m_source_anode_mask.const_array(mfi);
        auto const support = m_source_support.const_array(mfi);
        fine_ops.eval(mfi.tilebox(), fine_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept -> FineTuple
            {
                // Restrict disjoint contributions. Shared coarse nodes are summed later.
                amrex::Real const value = owner(i, j, k) && mask(i, j, k)
                    ? source(i, j, k, 0) : amrex::Real(0);
                amrex::Real const magnitude = std::abs(value);
                source(i, j, k, 0) = value;
                source(i, j, k, 1) = magnitude;
                amrex::Real const unsupported = magnitude != 0 && support(i, j, k) <= 0
                    ? amrex::Real(1) : amrex::Real(0);
                return {value, magnitude, unsupported};
            });
    }
    auto const fine_sums = fine_data.value();

    // With no source ghosts, Coarsen computes each Box's partial contribution.
    // Source Box endpoints are coarse-aligned, so every stencil recipient lies
    // in that Box's coarse patch. ParallelAdd accumulates those partial results.
    ablastr::coarsen::average::Coarsen(m_rho_patch, m_rho_source, m_ratio);
    m_rho_raw.setVal(0);
    ablastr::utils::communication::ParallelAdd(
        m_rho_raw, m_rho_patch, 0, 0, 2, zero, zero, false);
    m_rho_raw.OverrideSync();

    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> coarse_ops;
    amrex::ReduceData<amrex::Real, amrex::Real> coarse_data(coarse_ops);
    using CoarseTuple = decltype(coarse_data)::Type;
    m_rho.setVal(0);
    for (amrex::MFIter mfi(m_rho_raw, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const raw = m_rho_raw.const_array(mfi);
        auto const owner = m_coarse_owner->const_array(mfi);
        auto const mask = m_anode_mask.const_array(mfi);
        auto const rhs = m_rho.array(mfi);
        coarse_ops.eval(mfi.tilebox(), coarse_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept -> CoarseTuple
            {
                bool const anode = mask(i, j, k) == 0;
                rhs(i, j, k) = anode ? amrex::Real(0) : raw(i, j, k, 0);
                if (!owner(i, j, k)) {
                    return {0, 0};
                }
                amrex::Real const charge = raw(i, j, k, 0);
                amrex::Real const magnitude = raw(i, j, k, 1);
                return {charge, magnitude};
            });
    }
    auto const coarse_sums = coarse_data.value();
    // These are deposited nodal charges rho*dV, including exterior support, not
    // trapezoidal integration of a continuum density over the truncated fine box.
    amrex::Real const fine_volume = m_fine_geom.CellSize(0) * m_fine_geom.CellSize(1)
        * m_fine_geom.CellSize(2);
    amrex::Real const coarse_volume = m_geom.CellSize(0) * m_geom.CellSize(1)
        * m_geom.CellSize(2);
    std::array<amrex::Real, 5> sums{
        amrex::get<0>(fine_sums) * fine_volume,
        amrex::get<1>(fine_sums) * fine_volume,
        amrex::get<2>(fine_sums),
        amrex::get<0>(coarse_sums) * coarse_volume,
        amrex::get<1>(coarse_sums) * coarse_volume};
    amrex::ParallelDescriptor::ReduceRealSum(sums.data(), static_cast<int>(sums.size()));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::all_of(sums.begin(), sums.end(), [](amrex::Real v) { return std::isfinite(v); }),
        "Non-finite source encountered in Hall far-field charge transfer.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(sums[2] == 0,
        "A charged fine node has no non-anode coarse recipient. Reduce the coarse cell size "
        "or extend the coarse domain; source charge will not be silently discarded.");
    amrex::Real const tolerance = 512 * std::numeric_limits<amrex::Real>::epsilon()
        * (1 + std::log2(std::max(amrex::Real(1),
                                 amrex::Real(m_rho_source.boxArray().numPts()))))
        * std::max(sums[1], sums[4]);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::abs(sums[0] - sums[3]) <= tolerance && std::abs(sums[1] - sums[4]) <= tolerance,
        "Hall charge restriction lost source support or duplicated nodes. The coarse domain "
        "must contain the complete deposition/restriction support, including fine-domain ghosts.");
    // Fixed-potential nodes keep rho=0; the masked contribution is neither
    // redistributed nor recorded as physical wall absorption.
    m_rho.OverrideSync();
    ablastr::utils::communication::FillBoundary(m_rho, false);
    m_charge_ready = true;
}

void HallFarFieldPoisson::defineSolver ()
{
    m_rhs.define(m_rho.boxArray(), m_dmap, 1, m_rho.nGrowVect());
    m_epsilon_r.define(m_boxes, m_dmap, 1, 0);
    m_anode_phi.define(m_phi.boxArray(), m_dmap, 1, 0);
    m_rhs.setVal(0);
    HallElectrostaticMaterial::GetInstance().initializeMaterialFields(
        m_epsilon_r, m_anode_mask, m_anode_phi, m_geom);
    CheckGroundedAnode(m_anode_phi, m_geom);

    int max_coarsening = m_settings.m_max_coarsening_level;
    while (true) {
        amrex::LPInfo info;
        info.setMaxCoarseningLevel(max_coarsening);
        m_linop = std::make_unique<amrex::MLNodeLaplacian>(
            amrex::Vector<amrex::Geometry>{m_geom},
            amrex::Vector<amrex::BoxArray>{m_boxes},
            amrex::Vector<amrex::DistributionMapping>{m_dmap}, info);
        int const unresolved = FirstUnresolvedAnodeLevel(*m_linop, m_anode_mask);
        if (unresolved < 0) {
            break;
        }
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(unresolved > 0,
            "The Hall anode is unresolved on the coarse physical grid. "
            "Reduce insert.far_field.coarsening_ratio.");
        // An upper bound, not a demand to coarsen past the represented electrode.
        max_coarsening = unresolved - 1;
    }
    m_linop->setSigma(0, m_epsilon_r);
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> domain_bc;
    domain_bc.fill(amrex::LinOpBCType::Dirichlet);
    m_linop->setDomainBC(domain_bc, domain_bc);
    m_linop->setOversetMask(0, m_anode_mask);
    m_mlmg = std::make_unique<amrex::MLMG>(*m_linop);
}

void HallFarFieldPoisson::prepareInitialGuess (amrex::Real const time)
{
    if (m_history_depth > 0 && time < m_history_times[0]) {
        m_history_depth = 0;
    }
    if (m_history_depth == 0) {
        m_phi.setVal(0);
        return;
    }
    auto const ngrow = m_phi.nGrowVect();
    amrex::MultiFab::Copy(m_phi, *m_phi_history[0], 0, 0, 1, ngrow);
    if (time == m_history_times[0] || m_history_depth < 2) {
        return;
    }

    amrex::Real const interval = time - m_history_times[0];
    amrex::Real const previous_interval = m_history_times[0] - m_history_times[1];
    amrex::Real const time_scale = std::max({std::abs(time), std::abs(m_history_times[0]),
                                            std::abs(m_history_times[1])});
    amrex::Real const tolerance = 64 * std::numeric_limits<amrex::Real>::epsilon() * time_scale;
    if (previous_interval <= 0 || std::abs(interval - previous_interval) > tolerance) {
        // Fixed-step coefficients are not applied across a timestep change.
        m_history_depth = 1;
        return;
    }
    if (m_history_depth >= 3 &&
        std::abs(previous_interval - (m_history_times[1] - m_history_times[2])) > tolerance)
    {
        m_history_depth = 2;
    }

    int const order = std::min(m_settings.m_extrapolation_order, m_history_depth - 1);
    amrex::Real const alpha = m_settings.m_extrapolation_alpha;
    amrex::Real const beta = order == 2 ? m_settings.m_extrapolation_beta : amrex::Real(0);
    amrex::MultiFab::LinComb(
        m_phi, 1 + alpha + beta, *m_phi_history[0], 0,
        -(alpha + 2 * beta), *m_phi_history[1], 0, 0, 1, ngrow);
    if (order == 2) {
        amrex::MultiFab::Saxpy(m_phi, beta, *m_phi_history[2], 0, 0, 1, ngrow);
    }
}

void HallFarFieldPoisson::updateHistory (amrex::Real const time)
{
    int const capacity = m_settings.m_extrapolation_order + 1;
    if (m_history_depth == 0 || time != m_history_times[0]) {
        int const last = std::min(m_history_depth, capacity - 1);
        for (int hist = last; hist > 0; --hist) {
            std::swap(m_phi_history[hist], m_phi_history[hist - 1]);
            m_history_times[hist] = m_history_times[hist - 1];
        }
        m_history_depth = std::min(m_history_depth + 1, capacity);
    }
    // A repeated solve at the same physical time replaces the latest snapshot.
    amrex::MultiFab::Copy(*m_phi_history[0], m_phi, 0, 0, 1, m_phi.nGrowVect());
    m_history_times[0] = time;
}

void HallFarFieldPoisson::solve (amrex::Real const time, bool const verbose_step)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_defined && m_charge_ready && m_linop && m_mlmg && std::isfinite(time),
        "Hall coarse solve requires freshly restricted charge, a solver and a finite time.");
    m_charge_ready = false;
    prepareInitialGuess(time);
    ApplyFixedPotential(m_phi, m_anode_mask, m_anode_phi, m_geom);
    m_phi.OverrideSync();
    // computePhi scales its input rho in-place. A separate RHS preserves the
    // physical charge field even when the solver throws before restoring it.
    amrex::MultiFab::Copy(m_rhs, m_rho, 0, 0, 1, m_rho.nGrowVect());
    ablastr::fields::ConstMultiLevelScalarField const epsilon{&m_epsilon_r};
    // The preconfigured cache keeps the independently chosen geometry, mask and
    // MG hierarchy. No WarpX AMR level or boundary handler is used for this solve.
    ablastr::fields::computePhi<std::nullopt_t, std::nullopt_t, void>(
        {&m_rhs}, {&m_phi}, {0, 0, 0},
        m_settings.m_relative_tolerance, m_settings.m_absolute_tolerance,
        m_settings.m_max_iters, verbose_step ? m_settings.m_verbosity : 0,
        {m_geom}, {m_dmap}, {m_boxes}, ablastr::utils::enums::GridType::Collocated,
        false, false, false, false,
        std::nullopt, m_settings.m_num_final_sweeps,
        std::nullopt, std::nullopt, std::nullopt, std::nullopt,
        amrex::Vector<amrex::iMultiFab const*>{&m_anode_mask},
        &epsilon, &m_linop, &m_mlmg);
    bool const finite_phi = m_phi.is_finite(0, 1, 0);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(m_mlmg->getFinalResidual()) && finite_phi,
        "Hall coarse Poisson solve returned a non-finite residual or potential.");
    m_phi.OverrideSync();
    updateHistory(time);
}

void HallFarFieldPoisson::assertLayout (
    amrex::Geometry const& geom, amrex::MultiFab const& phi, amrex::MultiFab const& rho) const
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_defined && geom.Domain() == m_fine_geom.Domain() &&
            geom.Coord() == m_fine_geom.Coord() && !geom.isAnyPeriodic() &&
            phi.ixType().nodeCentered() && rho.ixType().nodeCentered() &&
            phi.nComp() == 1 && rho.nComp() == 1 &&
            phi.boxArray().CellEqual(m_fine_boxes) && rho.boxArray() == phi.boxArray() &&
            phi.DistributionMap() == m_dmap && rho.DistributionMap() == m_dmap &&
            phi.nGrowVect() == m_fine_phi_ngrow && rho.nGrowVect() == m_fine_rho_ngrow,
        "Hall far-field grid does not support changes to the fine field layout.");
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            geom.ProbLo(d) == m_fine_geom.ProbLo(d) &&
                geom.ProbHi(d) == m_fine_geom.ProbHi(d),
            "Hall far-field grid does not support moving fine-domain bounds.");
    }
}

void HallFarFieldPoisson::reportLayout () const
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }
    std::ostringstream os;
    auto const history_capacity = std::count_if(
        m_phi_history.begin(), m_phi_history.end(), [](auto const& field) { return bool(field); });
    os << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10);
    os << "Hall far-field layout initialized (stage three; fine boundary transfer not active).\n";
    os << "fine_domain";
    WriteBox(os, m_fine_geom.Domain());
    os << "\ncoarse_domain";
    WriteBox(os, m_geom.Domain());
    os << "\nratio " << m_ratio[0] << ' ' << m_ratio[1] << ' ' << m_ratio[2];
    os << "\nfine_offset " << m_fine_offset[0] << ' ' << m_fine_offset[1]
       << ' ' << m_fine_offset[2];
    os << "\nhistory_capacity " << history_capacity << "\nhistory_depth " << m_history_depth;
    os << "\nmg_levels " << m_linop->NMGLevels(0)
       << "\nmax_coarsening_level " << m_settings.m_max_coarsening_level;
    os << "\nranks " << amrex::ParallelDescriptor::NProcs();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        os << "\nphysical " << d << ' ' << m_fine_geom.ProbLo(d) << ' '
           << m_fine_geom.ProbHi(d) << ' ' << m_geom.ProbLo(d) << ' ' << m_geom.ProbHi(d);
    }
    os << "\n# box: id fine_rank coarse_rank fine_lo[3] fine_hi[3] coarse_lo[3] coarse_hi[3]";
    for (int ibox = 0; ibox < m_boxes.size(); ++ibox) {
        os << "\nbox " << ibox << ' ' << m_dmap[ibox] << ' ' << m_phi.DistributionMap()[ibox];
        WriteBox(os, m_fine_boxes[ibox]);
        WriteBox(os, m_boxes[ibox]);
    }
    os << '\n';
    amrex::Print() << os.str();
    if (!m_layout_file.empty()) {
        std::ofstream output(m_layout_file);
        output << os.str();
        output.close();
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(output.good(),
            "Could not write insert.far_field.layout_file: " + m_layout_file);
    }
}

} // namespace Insert

#endif // WARPX_USE_HALL_FAR_FIELD
