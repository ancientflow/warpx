#include "HallFarFieldPoisson.H"

#ifdef WARPX_USE_HALL_FAR_FIELD

#include "EmbeddedBoundary/Enabled.H"
#include "Insert/Fields/HallElectrostaticMaterial.H"
#include "Python/callbacks.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <AMReX_CoordSys.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

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

} // namespace

void HallFarFieldPoisson::initialize (
    WarpX const& warpx, amrex::MultiFab const& phi, amrex::MultiFab const& rho,
    int const history_capacity)
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
    defineGeometry(m_fine_geom);
    defineBoxes(m_fine_boxes, phi.DistributionMap());
    allocateFields(phi, rho, history_capacity);
    m_defined = true;
    reportLayout();
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
    m_rho_patch.define(patch_nodes, m_dmap, 1, rho_ngrow);
    m_phi_patch.define(patch_nodes, m_dmap, 1, phi_ngrow);
    m_rho.setVal(0);
    m_phi.setVal(0);
    m_rho_patch.setVal(0);
    m_phi_patch.setVal(0);
    for (int hist = 0; hist < history_capacity; ++hist) {
        m_phi_history[hist] = std::make_unique<amrex::MultiFab>(nodes, m_dmap, 1, phi_ngrow);
        m_phi_history[hist]->setVal(0);
    }
    m_history_depth = 0;
    m_history_times.fill(0);
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
    os << "Hall far-field layout initialized (stage one; coarse solve not active).\n";
    os << "fine_domain";
    WriteBox(os, m_fine_geom.Domain());
    os << "\ncoarse_domain";
    WriteBox(os, m_geom.Domain());
    os << "\nratio " << m_ratio[0] << ' ' << m_ratio[1] << ' ' << m_ratio[2];
    os << "\nfine_offset " << m_fine_offset[0] << ' ' << m_fine_offset[1]
       << ' ' << m_fine_offset[2];
    os << "\nhistory_capacity " << history_capacity << "\nhistory_depth " << m_history_depth;
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
