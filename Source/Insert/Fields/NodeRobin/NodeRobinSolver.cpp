#include "NodeRobinSolver.H"

#include "CHT_MLNodeRobinLaplacian.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"

#include <ablastr/constant.H>

#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace Insert {
namespace {

amrex::Real uniformVoltage (amrex::iMultiFab const& mask, amrex::MultiFab const& voltage)
{
    amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMax> ops;
    amrex::ReduceData<amrex::Real, amrex::Real> data(ops);
    using Tuple = decltype(data)::Type;
    for (amrex::MFIter mfi(mask); mfi.isValid(); ++mfi) {
        auto const m = mask.const_array(mfi);
        auto const v = voltage.const_array(mfi);
        ops.eval(mfi.validbox(), data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> Tuple {
                return {m(i,j,k) == 0 ? v(i,j,k) : std::numeric_limits<amrex::Real>::max(),
                        m(i,j,k) == 0 ? v(i,j,k) : std::numeric_limits<amrex::Real>::lowest()};
            });
    }
    auto const result = data.value();
    auto lo = amrex::get<0>(result);
    auto hi = amrex::get<1>(result);
    amrex::ParallelDescriptor::ReduceRealMin(lo);
    amrex::ParallelDescriptor::ReduceRealMax(hi);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(lo) && std::isfinite(hi) && lo == hi,
        "Monopole Robin requires a nonempty anode with one uniform finite voltage.");
    return lo;
}

int coarseningLimit (amrex::iMultiFab const& mask, amrex::Geometry const& geom, int requested)
{
    int const limit = requested < 0 ? 30 : requested;
    int cap = 0;
    auto grids = amrex::convert(mask.boxArray(), amrex::IntVect(0));
    auto domain = geom.Domain();
    auto const* current = &mask;
    std::unique_ptr<amrex::iMultiFab> coarse;
    while (cap < limit && domain.coarsenable(2, 2) && grids.coarsenable(2, 2)) {
        grids.coarsen(2);
        domain.coarsen(2);
        auto next = std::make_unique<amrex::iMultiFab>(
            amrex::convert(grids, amrex::IntVect(1)), mask.DistributionMap(), 1, 0);
        // Match the operator's nodal mask injection, including thin electrodes.
        amrex::average_down_nodal(*current, *next, amrex::IntVect(2));
        if (next->min(0) != 0) { break; }
        coarse = std::move(next);
        current = coarse.get();
        ++cap;
    }
    return cap;
}

} // namespace

struct NodeRobinSolver::Cache
{
    amrex::Geometry geom;
    amrex::MultiFab const* epsilon = nullptr;
    amrex::iMultiFab const* mask = nullptr;
    ablastr::fields::MLMGOptions options;
    std::unique_ptr<amrex::MultiFab> rhs;
    std::unique_ptr<amrex::MultiFab> residual;
    // Reverse destruction order: MLMG must die before the operator it references.
    std::unique_ptr<cht::MLNodeRobinLaplacian> op;
    std::unique_ptr<amrex::MLMG> mlmg;
    int solves = 0;
};

NodeRobinSolver::NodeRobinSolver () = default;
NodeRobinSolver::~NodeRobinSolver () = default;

void NodeRobinSolver::readParameters (amrex::Geometry const& geom)
{
    amrex::ParmParse pp("insert");
    std::string boundary = "dirichlet";
    pp.query("electrostatic_boundary", boundary);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(boundary == "dirichlet" || boundary == "monopole",
        "insert.electrostatic_boundary must be dirichlet or monopole.");
    m_enabled = boundary == "monopole";
    if (!m_enabled) { return; }
    amrex::Vector<amrex::Real> center;
    utils::parser::queryArrWithParser(pp, "robin_center", center);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(center.empty() || center.size() == 3,
        "insert.robin_center requires three coordinates.");
    for (int d = 0; d < 3; ++d) {
        m_center[d] = center.empty() ? (geom.ProbLo(d) + geom.ProbHi(d))*0.5 : center[d];
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(m_center[d]) &&
            m_center[d] > geom.ProbLo(d) && m_center[d] < geom.ProbHi(d),
            "The Robin expansion center must lie strictly inside the domain.");
    }
    utils::parser::queryWithParser(pp, "robin_phi_infinity", m_phi_infinity);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(m_phi_infinity),
        "insert.robin_phi_infinity must be finite.");
    pp.query("robin_diagnostics", m_diagnostics);
}

void NodeRobinSolver::solve (
    amrex::MultiFab const& rho, amrex::MultiFab& phi,
    amrex::MultiFab const& epsilon, amrex::iMultiFab const& mask,
    amrex::MultiFab const& voltage, amrex::Geometry const& geom,
    ablastr::fields::MLMGOptions const& options)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_enabled && geom.Coord() == 0 && !geom.isAnyPeriodic(),
        "Monopole Robin requires nonperiodic 3D Cartesian geometry.");
    auto const nodal_grids = amrex::convert(epsilon.boxArray(), amrex::IntVect(1));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        rho.nComp() == 1 && phi.nComp() == 1 && voltage.nComp() == 1 && mask.nComp() == 1 &&
        rho.boxArray() == nodal_grids && phi.boxArray() == nodal_grids &&
        mask.boxArray() == nodal_grids && voltage.boxArray() == nodal_grids &&
        rho.DistributionMap() == epsilon.DistributionMap() &&
        phi.DistributionMap() == epsilon.DistributionMap() &&
        mask.DistributionMap() == epsilon.DistributionMap() &&
        voltage.DistributionMap() == epsilon.DistributionMap() && phi.nGrowVect().min() >= 1,
        "Robin rho, phi, mask, voltage and epsilon must share the solver layout.");
    auto const bottom = options.bottom_solver;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        bottom == amrex::BottomSolver::Default || bottom == amrex::BottomSolver::smoother ||
        bottom == amrex::BottomSolver::bicgstab || bottom == amrex::BottomSolver::cg ||
        bottom == amrex::BottomSolver::bicgcg || bottom == amrex::BottomSolver::cgbicg,
        "Robin supports only the native smoother and Krylov bottom solvers.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(voltage.is_finite(0, 1, 0),
        "Robin anode voltage must be finite.");
    auto const anode_voltage = uniformVoltage(mask, voltage);
    if (!m_cache) {
        m_cache = std::make_unique<Cache>();
        auto& c = *m_cache;
        c.geom = geom;
        c.epsilon = &epsilon;
        c.mask = &mask;
        c.options = options;
        amrex::LPInfo info;
        ablastr::fields::applyMLMGOptions(info, options);
        info.setSemicoarsening(false);
        int const cap = coarseningLimit(mask, geom, options.max_coarsening_level);
        info.setMaxCoarseningLevel(cap);
        c.op = std::make_unique<cht::MLNodeRobinLaplacian>(
            amrex::Vector<amrex::Geometry>{geom},
            amrex::Vector<amrex::BoxArray>{epsilon.boxArray()},
            amrex::Vector<amrex::DistributionMapping>{epsilon.DistributionMap()}, info);
        c.op->setSigma(0, epsilon);
        c.op->setMonopoleBC(m_center, m_phi_infinity);
        c.op->setOversetMask(0, mask, anode_voltage);
        c.mlmg = std::make_unique<amrex::MLMG>(*c.op);
        c.rhs = std::make_unique<amrex::MultiFab>(nodal_grids, phi.DistributionMap(), 1, 0);
        if (m_diagnostics) {
            c.residual = std::make_unique<amrex::MultiFab>(
                nodal_grids, phi.DistributionMap(), 1, 0);
            amrex::Print() << "Robin cache: max coarsening level = " << cap << '\n';
        }
    }
    auto& c = *m_cache;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(c.epsilon == &epsilon && c.mask == &mask &&
        c.rhs->boxArray() == phi.boxArray() && c.rhs->DistributionMap() == phi.DistributionMap() &&
        c.geom.Domain() == geom.Domain() &&
        c.options.max_coarsening_level == options.max_coarsening_level &&
        c.options.agglomeration == options.agglomeration &&
        c.options.consolidation == options.consolidation &&
        c.options.agglomeration_grid_size == options.agglomeration_grid_size &&
        c.options.consolidation_grid_size == options.consolidation_grid_size,
        "Robin cached layout, material, mask or hierarchy options changed.");
    for (int d = 0; d < 3; ++d) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(c.geom.ProbLo(d) == geom.ProbLo(d) &&
            c.geom.ProbHi(d) == geom.ProbHi(d), "Robin cached geometry changed.");
    }
    c.op->setInternalDirichletValue(anode_voltage);
    for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto const m = mask.const_array(mfi);
        auto const p = phi.array(mfi);
        amrex::ParallelFor(mfi.validbox(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (m(i,j,k) == 0) { p(i,j,k) = anode_voltage; }
            });
    }
    amrex::MultiFab::Copy(*c.rhs, rho, 0, 0, 1, 0);
    normalizeBoundarySource(*c.rhs, geom);
    c.rhs->mult(-1.0/ablastr::constant::SI::epsilon_0);
    c.rhs->OverrideSync(geom.periodicity());
    ablastr::fields::applyMLMGOptions(*c.mlmg, options);
    c.mlmg->setConvergenceNormType(amrex::MLMGNormType::greater);
    auto const abs_tol = options.absolute_tolerance == 0 && c.rhs->norm0() == 0
        ? amrex::Real(1.e-6) : options.absolute_tolerance;
    c.mlmg->solve({&phi}, {c.rhs.get()}, options.relative_tolerance, abs_tol);
    ++c.solves;
    if (m_diagnostics) {
        c.op->compPhysicalResidual(*c.residual, phi, *c.rhs);
        amrex::Print() << "Robin solve " << c.solves
                       << ": physical residual = " << c.residual->norm0()
                       << ", coefficient builds = " << c.op->coefficientBuilds() << '\n';
    }
}

void NodeRobinSolver::normalizeBoundarySource (
    amrex::MultiFab& source, amrex::Geometry const& geom)
{
    // Linear deposition stores q/dV at nodes. The nodal Robin operator uses
    // half/quarter/eighth dual volumes on physical faces/edges/corners.
    // Normalize this private source, never the persistent deposition or rho.
    auto const domain = amrex::surroundingNodes(geom.Domain());
    auto const lo = amrex::lbound(domain);
    auto const hi = amrex::ubound(domain);
    for (amrex::MFIter mfi(source); mfi.isValid(); ++mfi) {
        auto const rhs = source.array(mfi);
        amrex::ParallelFor(mfi.validbox(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                int const faces = int(i == lo.x || i == hi.x) +
                                  int(j == lo.y || j == hi.y) +
                                  int(k == lo.z || k == hi.z);
                rhs(i,j,k) *= amrex::Real(1 << faces);
            });
    }
}

void NodeRobinSolver::fillPhiGhosts (
    amrex::MultiFab& phi, amrex::Geometry const& geom) const
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(phi.nGrowVect().min() >= 1,
        "Collocated Robin fields require one phi ghost layer.");
    auto const domain = amrex::surroundingNodes(geom.Domain());
    auto const problo = geom.ProbLoArray();
    auto const dx = geom.CellSizeArray();
    auto const center = m_center;
    auto const far_phi = m_phi_infinity;
    for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto const p = phi.array(mfi);
        for (int d = 0; d < 3; ++d) {
            for (int sign : {-1, 1}) {
                auto face = domain;
                int const boundary = sign < 0 ? domain.smallEnd(d) : domain.bigEnd(d);
                face.setSmall(d, boundary + sign);
                face.setBig(d, boundary + sign);
                face &= mfi.fabbox();
                // Only one coordinate is outside: these are precisely the phi
                // nodes read by centered gradients on valid faces, edges and corners.
                amrex::ParallelFor(face,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        amrex::IntVect b(i,j,k);
                        b[d] = boundary;
                        auto inward = b;
                        inward[d] -= sign;
                        amrex::Real radius2 = 0;
                        amrex::Real rd = 0;
                        for (int axis = 0; axis < 3; ++axis) {
                            auto const r = problo[axis] +
                                (b[axis]-domain.smallEnd(axis))*dx[axis] - center[axis];
                            radius2 += r*r;
                            if (axis == d) { rd = r; }
                        }
                        // (phi_g-phi_in)/(2*h) = -q*(phi_b-phi_inf), q = n.r/r^2.
                        p(i,j,k) = p(inward) - 2*dx[d]*sign*rd/radius2*(p(b)-far_phi);
                    });
            }
        }
    }
}

} // namespace Insert
