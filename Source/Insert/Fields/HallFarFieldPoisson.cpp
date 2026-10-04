#include "HallFarFieldPoisson.H"

#ifdef WARPX_USE_HALL_FAR_FIELD

#include "Insert/Fields/HallElectrostaticMaterial.H"
#include "Utils/Parser/ParserUtils.H"
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

/** Extend only physical ghosts from synchronized interior samples.
 * Quadratic phi extension makes centered boundary gradients second-order one-sided.
 * Linear E extension retains nonzero tangential fields at artificial boundaries.
 */
void ExtendPhysicalGhosts (
    amrex::MultiFab& field, amrex::MultiFab& work,
    amrex::Geometry const& geom, int const order)
{
    field.OverrideSync();
    work.setVal(0);
    ablastr::utils::communication::ParallelCopy(
        work, field, 0, 0, 1, amrex::IntVect::TheZeroVector(), work.nGrowVect(), false);
    auto const domain = amrex::convert(geom.Domain(), field.ixType());
    auto const lo = amrex::lbound(domain);
    auto const hi = amrex::ubound(domain);
    // Whole FABs: ghosts must not be updated repeatedly by neighboring tiles.
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
        auto const values = field.array(mfi);
        auto const samples = work.const_array(mfi);
        amrex::ParallelFor(mfi.fabbox(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (i >= lo.x && i <= hi.x && j >= lo.y && j <= hi.y &&
                    k >= lo.z && k <= hi.z) {
                    return;
                }
                int const point[3] = {i, j, k};
                int const lower[3] = {lo.x, lo.y, lo.z};
                int const upper[3] = {hi.x, hi.y, hi.z};
                int indices[3][3];
                amrex::Real weights[3][3];
                int count[3];
                for (int d = 0; d < 3; ++d) {
                    bool const outside = point[d] < lower[d] || point[d] > upper[d];
                    count[d] = outside ? order + 1 : 1;
                    indices[d][0] = point[d];
                    weights[d][0] = 1;
                    if (outside) {
                        int const side = point[d] < lower[d] ? 1 : -1;
                        int const edge = side == 1 ? lower[d] : upper[d];
                        amrex::Real const g = std::abs(point[d] - edge);
                        for (int n = 0; n <= order; ++n) {
                            indices[d][n] = edge + side * n;
                        }
                        if (order == 2) {
                            weights[d][0] = (g + 1) * (g + 2) / 2;
                            weights[d][1] = -g * (g + 2);
                            weights[d][2] = g * (g + 1) / 2;
                        } else {
                            weights[d][0] = g + 1;
                            weights[d][1] = -g;
                        }
                    }
                }
                amrex::Real value = 0;
                for (int a = 0; a < count[0]; ++a) {
                    for (int b = 0; b < count[1]; ++b) {
                        for (int c = 0; c < count[2]; ++c) {
                            value += weights[0][a] * weights[1][b] * weights[2][c] *
                                samples(indices[0][a], indices[1][b], indices[2][c]);
                        }
                    }
                }
                values(i, j, k) = value;
            });
    }
    ablastr::utils::communication::FillBoundary(field, false);
}

namespace {

/** Convert an aligned physical extension to whole coarse cells. */
int ExtensionCells (amrex::Real extension, amrex::Real spacing)
{
    amrex::Real const cells = extension / spacing;
    amrex::Real const nearest = std::round(cells);
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

    m_fine_geom = warpx.Geom(0);
    m_fine_boxes = amrex::convert(phi.boxArray(), amrex::IntVect::TheZeroVector());
    readSolverSettings(defaults);
    defineGeometry(m_fine_geom);
    defineBoxes(m_fine_boxes, phi.DistributionMap());
    allocateFields(phi, rho, m_settings.m_extrapolation_order + 1);
    defineChargeTransfer(rho);
    defineSolver();
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
}

void HallFarFieldPoisson::defineGeometry (amrex::Geometry const& fine_geom)
{
    amrex::ParmParse const pp("insert.far_field");
    amrex::Vector<amrex::Real> prob_lo;
    amrex::Vector<amrex::Real> prob_hi;
    amrex::Vector<int> ratio;
    utils::parser::getArrWithParser(pp, "prob_lo", prob_lo);
    utils::parser::getArrWithParser(pp, "prob_hi", prob_hi);
    utils::parser::getArrWithParser(pp, "coarsening_ratio", ratio);
    pp.query("layout_file", m_layout_file);
    m_ratio = amrex::IntVect(ratio.data());
    amrex::Box const& fine_domain = fine_geom.Domain();

    amrex::Box const embedded_domain = amrex::coarsen(fine_domain, m_ratio);
    amrex::Box coarse_domain = embedded_domain;
    amrex::IntVect high_extension;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        amrex::Real const spacing = fine_geom.CellSize(d) * ratio[d];
        m_fine_offset[d] = ExtensionCells(fine_geom.ProbLo(d) - prob_lo[d], spacing);
        high_extension[d] = ExtensionCells(prob_hi[d] - fine_geom.ProbHi(d), spacing);
        auto const lo = static_cast<amrex::Long>(embedded_domain.smallEnd(d)) - m_fine_offset[d];
        auto const hi = static_cast<amrex::Long>(embedded_domain.bigEnd(d)) + high_extension[d];
        coarse_domain.setSmall(d, static_cast<int>(lo));
        coarse_domain.setBig(d, static_cast<int>(hi));
        // Snap tolerated floating-point roundoff to the aligned physical bounds.
        prob_lo[d] = fine_geom.ProbLo(d) - m_fine_offset[d] * spacing;
        prob_hi[d] = fine_geom.ProbHi(d) + high_extension[d] * spacing;
    }

    // Retain natural coarse indices (possibly negative) so the existing WarpX
    // restriction/interpolation kernels still use fine_index = coarse_index * ratio.
    m_geom.define(coarse_domain, amrex::RealBox(prob_lo.data(), prob_hi.data()),
                  amrex::CoordSys::cartesian, {0, 0, 0});
}

void HallFarFieldPoisson::defineBoxes (
    amrex::BoxArray const& fine_boxes, amrex::DistributionMapping const& fine_dmap)
{
    amrex::Box const& domain = m_fine_geom.Domain();

    amrex::Vector<amrex::Box> coarse_boxes(fine_boxes.size());
    for (int ibox = 0; ibox < fine_boxes.size(); ++ibox) {
        amrex::Box const fine_box = fine_boxes[ibox];
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
    // Box numbering is unchanged: each overlapping coarse/fine pair has the same owner.
    m_dmap = fine_dmap;
}

void HallFarFieldPoisson::allocateFields (
    amrex::MultiFab const& phi, amrex::MultiFab const& rho, int const history_capacity)
{
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
    m_rho_source.define(source_nodes, m_dmap, 1, 0);
    m_rho_patch.define(patch_nodes, m_dmap, 1, 0);
    m_source_anode_mask.define(source_nodes, m_dmap, 1, 0);
    m_anode_mask.define(m_rho.boxArray(), m_dmap, 1, 0);
    m_source_owner = m_rho_source.OwnerMask();
    m_rho_source.setVal(0);
    m_rho_patch.setVal(0);

    auto& material = HallElectrostaticMaterial::GetInstance();
    material.fillAnodeMask(m_source_anode_mask, m_fine_geom);
    material.fillAnodeMask(m_anode_mask, m_geom);
}

void HallFarFieldPoisson::restrictCharge (amrex::MultiFab const& rho)
{
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

    for (amrex::MFIter mfi(m_rho_source, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const source = m_rho_source.array(mfi);
        auto const owner = m_source_owner->const_array(mfi);
        auto const mask = m_source_anode_mask.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Restrict disjoint contributions. Shared coarse nodes are summed later.
                if (!owner(i, j, k) || !mask(i, j, k)) {
                    source(i, j, k) = 0;
                }
            });
    }

    // Coarsen each Box's partial contribution, then sum shared recipients.
    ablastr::coarsen::average::Coarsen(m_rho_patch, m_rho_source, m_ratio);
    m_rho.setVal(0);
    ablastr::utils::communication::ParallelAdd(
        m_rho, m_rho_patch, 0, 0, 1, zero, zero, false);
    m_rho.OverrideSync();
    for (amrex::MFIter mfi(m_rho, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const mask = m_anode_mask.const_array(mfi);
        auto const rhs = m_rho.array(mfi);
        amrex::ParallelFor(mfi.tilebox(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (mask(i, j, k) == 0) {
                    rhs(i, j, k) = 0;
                }
            });
    }
    // Fixed-potential nodes keep rho=0; the masked contribution is neither
    // redistributed nor recorded as physical wall absorption.
    m_rho.OverrideSync();
    ablastr::utils::communication::FillBoundary(m_rho, false);
}

void HallFarFieldPoisson::defineSolver ()
{
    m_rhs.define(m_rho.boxArray(), m_dmap, 1, m_rho.nGrowVect());
    m_epsilon_r.define(m_boxes, m_dmap, 1, 0);
    m_anode_phi.define(m_phi.boxArray(), m_dmap, 1, 0);
    m_rhs.setVal(0);
    HallElectrostaticMaterial::GetInstance().initializeMaterialFields(
        m_epsilon_r, m_anode_mask, m_anode_phi, m_geom);

    int max_coarsening = m_settings.m_max_coarsening_level;
    while (true) {
        amrex::LPInfo info;
        info.setMaxCoarseningLevel(max_coarsening);
        m_linop = std::make_unique<amrex::MLNodeLaplacian>(
            amrex::Vector<amrex::Geometry>{m_geom},
            amrex::Vector<amrex::BoxArray>{m_boxes},
            amrex::Vector<amrex::DistributionMapping>{m_dmap}, info);
        int const unresolved = FirstUnresolvedAnodeLevel(*m_linop, m_anode_mask);
        if (unresolved <= 0) {
            break;
        }
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
    m_phi.OverrideSync();
    updateHistory(time);
}

void HallFarFieldPoisson::interpolateToFine (amrex::MultiFab& phi)
{
    m_phi.OverrideSync();
    ablastr::utils::communication::FillBoundary(m_phi, false);
    // Geometry uses natural coarse indices: physical offsets are already encoded
    // in the expanded domain's lower index, including asymmetric z extensions.
    m_phi_patch.setVal(0);
    ablastr::utils::communication::ParallelCopy(
        m_phi_patch, m_phi, 0, 0, 1, amrex::IntVect::TheZeroVector(),
        amrex::IntVect::TheZeroVector(), false);
    phi.setVal(0);
    for (amrex::MFIter mfi(phi, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        ablastr::fields::details::PoissonInterpCPtoFP const interp(
            phi.array(mfi), m_phi_patch.const_array(mfi), m_ratio);
        amrex::ParallelFor(mfi.tilebox(), interp);
    }
    phi.OverrideSync();
    // Fine material data take precedence where the anode meets a domain face.
    HallElectrostaticMaterial::GetInstance().restoreAnodePotential(phi, 0);
}

void HallFarFieldPoisson::fillElectricFieldGhosts (
    std::array<amrex::MultiFab*, 3> const& fields)
{
    for (int d = 0; d < 3; ++d) {
        auto& field = *fields[d];
        if (!m_fine_e_work[d]) {
            auto ngrow = field.nGrowVect();
            for (int axis = 0; axis < 3; ++axis) {
                ngrow[axis] = std::max(ngrow[axis], 1);
            }
            m_fine_e_work[d] = std::make_unique<amrex::MultiFab>(
                field.boxArray(), m_dmap, 1, ngrow);
        }
        ExtendPhysicalGhosts(field, *m_fine_e_work[d], m_fine_geom, 1);
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
    os << "Hall far-field layout initialized (coarse solve and fine Dirichlet coupling).\n";
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
    }
}

} // namespace Insert

#endif // WARPX_USE_HALL_FAR_FIELD
