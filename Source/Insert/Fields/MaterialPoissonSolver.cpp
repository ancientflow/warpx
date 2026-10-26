#include "MaterialPoissonSolver.H"

#include "NodeRobin/MLNodeRobinLaplacian.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"

#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace Insert {
namespace {

char const* bottomSolverName (amrex::BottomSolver solver)
{
    switch (solver) {
    case amrex::BottomSolver::Default: return "default";
    case amrex::BottomSolver::smoother: return "smoother";
    case amrex::BottomSolver::bicgstab: return "bicgstab";
    case amrex::BottomSolver::cg: return "cg";
    case amrex::BottomSolver::bicgcg: return "bicgcg";
    case amrex::BottomSolver::cgbicg: return "cgbicg";
    default: return "unsupported";
    }
}

// Conservative geometric guard: never average across a material transition,
// erase a mask feature inside a coarse cell, or bridge a resolved fine-grid gap.
// The candidate grids keep the fine layout's box order and distribution.
int coarseningLimit (amrex::iMultiFab const& mask, amrex::MultiFab const& epsilon,
                     amrex::Geometry const& geom, int requested)
{
    int const limit = requested < 0 ? 30 : requested;
    int cap = 0;
    auto grids = epsilon.boxArray();
    auto domain = geom.Domain();
    auto const* current = &mask;
    auto const* sigma = &epsilon;
    std::unique_ptr<amrex::iMultiFab> coarse;
    std::unique_ptr<amrex::MultiFab> coarse_sigma;
    while (cap < limit && domain.coarsenable(2, 2) && grids.coarsenable(2, 2)) {
        grids.coarsen(2);
        domain.coarsen(2);
        auto next = std::make_unique<amrex::iMultiFab>(
            amrex::convert(grids, amrex::IntVect(1)), mask.DistributionMap(), 1, 0);
        auto next_sigma = std::make_unique<amrex::MultiFab>(
            grids, epsilon.DistributionMap(), 1, 0);
        amrex::average_down_nodal(*current, *next, amrex::IntVect(2));
        amrex::average_down(*sigma, *next_sigma, 0, 1, amrex::IntVect(2));
        amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMin,
                         amrex::ReduceOpMax> ops;
        amrex::ReduceData<int, int, int> data(ops);
        using Tuple = decltype(data)::Type;
        for (amrex::MFIter mfi(*next); mfi.isValid(); ++mfi) {
            auto const m = next->const_array(mfi);
            ops.eval(mfi.validbox(), data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> Tuple {
                    return {0, m(i,j,k) == 0 ? k : std::numeric_limits<int>::max(),
                               m(i,j,k) == 0 ? k : std::numeric_limits<int>::lowest()};
                });
        }
        for (amrex::MFIter mfi(*next_sigma); mfi.isValid(); ++mfi) {
            auto const m = current->const_array(mfi);
            auto const e = sigma->const_array(mfi);
            ops.eval(mfi.validbox(), data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> Tuple {
                    int bad = 0;
                    int fine_min = 1, fine_max = 0, coarse_min = 1, coarse_max = 0;
                    for (int z = 0; z <= 2; ++z) {
                        for (int y = 0; y <= 2; ++y) {
                            for (int x = 0; x <= 2; ++x) {
                                int const v = m(2*i+x,2*j+y,2*k+z);
                                fine_min = amrex::min(fine_min, v);
                                fine_max = amrex::max(fine_max, v);
                                if (x != 1 && y != 1 && z != 1) {
                                    coarse_min = amrex::min(coarse_min, v);
                                    coarse_max = amrex::max(coarse_max, v);
                                }
                                // Preserve every two-edge connectivity/gap decision.
                                if (x == 1 && m(2*i,2*j+y,2*k+z) ==
                                    m(2*i+2,2*j+y,2*k+z) &&
                                    v != m(2*i,2*j+y,2*k+z)) { bad = 1; }
                                if (y == 1 && m(2*i+x,2*j,2*k+z) ==
                                    m(2*i+x,2*j+2,2*k+z) &&
                                    v != m(2*i+x,2*j,2*k+z)) { bad = 1; }
                                if (z == 1 && m(2*i+x,2*j+y,2*k) ==
                                    m(2*i+x,2*j+y,2*k+2) &&
                                    v != m(2*i+x,2*j+y,2*k)) { bad = 1; }
                                if (x < 2 && y < 2 && z < 2 &&
                                    e(2*i+x,2*j+y,2*k+z) != e(2*i,2*j,2*k)) {
                                    bad = 1;
                                }
                            }
                        }
                    }
                    bad |= fine_min != coarse_min || fine_max != coarse_max;
                    return {bad, std::numeric_limits<int>::max(),
                                 std::numeric_limits<int>::lowest()};
                });
        }
        auto const result = data.value();
        int bad = amrex::get<0>(result);
        int klo = amrex::get<1>(result);
        int khi = amrex::get<2>(result);
        amrex::ParallelDescriptor::ReduceIntMax(bad);
        amrex::ParallelDescriptor::ReduceIntMin(klo);
        amrex::ParallelDescriptor::ReduceIntMax(khi);
        bool const has_electrode = klo != std::numeric_limits<int>::max();
        if (bad || (has_electrode && khi <= klo)) { break; }
        coarse = std::move(next);
        coarse_sigma = std::move(next_sigma);
        current = coarse.get();
        sigma = coarse_sigma.get();
        ++cap;
    }
    return cap;
}

} // namespace

struct MaterialPoissonSolver::Cache
{
    amrex::Geometry geom;
    amrex::MultiFab const* epsilon = nullptr;
    amrex::iMultiFab const* mask = nullptr;
    ablastr::fields::MLMGOptions options;
    std::unique_ptr<amrex::MultiFab> rhs;
    std::unique_ptr<amrex::MultiFab> residual;
    // Reverse destruction order: MLMG must die before the operator it references.
    std::unique_ptr<Insert::Poisson::MLNodeRobinLaplacian> op;
    std::unique_ptr<amrex::MLMG> mlmg;
    int solves = 0;
    int coarsening_cap = 0;
};

MaterialPoissonSolver::MaterialPoissonSolver () = default;
MaterialPoissonSolver::~MaterialPoissonSolver () = default;

void MaterialPoissonSolver::readParameters (amrex::Geometry const& geom)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_caches.empty(),
        "Material boundary parameters cannot change after the first solve.");
    amrex::ParmParse pp("insert");
    std::string boundary = "dirichlet";
    pp.query("electrostatic_boundary", boundary);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(boundary == "dirichlet" || boundary == "monopole",
        "insert.electrostatic_boundary must be dirichlet or monopole.");
    m_monopole = boundary == "monopole";
    pp.query("robin_diagnostics", m_diagnostics);
    if (!m_monopole) { return; }
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
}

amrex::MLMG& MaterialPoissonSolver::solveLevel (
    int lev, amrex::MultiFab const& rho, amrex::MultiFab& phi,
    amrex::MultiFab const& epsilon, amrex::iMultiFab const& mask,
    amrex::MultiFab const& voltage, amrex::Geometry const& geom,
    ablastr::fields::MLMGOptions const& options, amrex::Real absolute_tolerance)
{
    auto const prepare_start = amrex::second();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(lev >= 0, "Invalid material AMR level.");
    bool const robin_level = m_monopole && lev == 0;
    if (static_cast<int>(m_caches.size()) <= lev) { m_caches.resize(lev+1); }
    auto& cache = m_caches[lev];
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(geom.Coord() == 0 && !geom.isAnyPeriodic(),
        "Material Poisson requires nonperiodic 3D Cartesian geometry.");
    auto const nodal_grids = amrex::convert(epsilon.boxArray(), amrex::IntVect(1));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        epsilon.nComp() == 1 && epsilon.ixType().cellCentered() &&
        rho.nComp() == 1 && phi.nComp() == 1 && voltage.nComp() == 1 && mask.nComp() == 1 &&
        rho.boxArray() == nodal_grids && phi.boxArray() == nodal_grids &&
        mask.boxArray() == nodal_grids && voltage.boxArray() == nodal_grids &&
        rho.DistributionMap() == epsilon.DistributionMap() &&
        phi.DistributionMap() == epsilon.DistributionMap() &&
        mask.DistributionMap() == epsilon.DistributionMap() &&
        voltage.DistributionMap() == epsilon.DistributionMap() && phi.nGrowVect().min() >= 1,
        "Material Poisson rho, phi, mask, voltage and epsilon must share the solver layout.");
    auto const bottom = options.bottom_solver;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        bottom == amrex::BottomSolver::Default || bottom == amrex::BottomSolver::smoother ||
        bottom == amrex::BottomSolver::bicgstab || bottom == amrex::BottomSolver::cg ||
        bottom == amrex::BottomSolver::bicgcg || bottom == amrex::BottomSolver::cgbicg,
        "Material Poisson supports only the native smoother and Krylov bottom solvers.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(voltage.is_finite(0, 1, 0),
        "Material Poisson anode voltage must be finite.");
    if (!cache) {
        cache = std::make_unique<Cache>();
        auto& c = *cache;
        c.geom = geom;
        c.epsilon = &epsilon;
        c.mask = &mask;
        c.options = options;
        amrex::LPInfo info;
        ablastr::fields::applyMLMGOptions(info, options);
        info.setSemicoarsening(false);
        int const cap = coarseningLimit(mask, epsilon, geom, options.max_coarsening_level);
        info.setMaxCoarseningLevel(cap);
        c.coarsening_cap = cap;
        c.op = std::make_unique<Insert::Poisson::MLNodeRobinLaplacian>(
            amrex::Vector<amrex::Geometry>{geom},
            amrex::Vector<amrex::BoxArray>{epsilon.boxArray()},
            amrex::Vector<amrex::DistributionMapping>{epsilon.DistributionMap()}, info);
        c.op->setSigma(0, epsilon);
        if (robin_level) {
            c.op->setMonopoleBC(m_center, m_phi_infinity);
        } else {
            Insert::Poisson::MLNodeRobinLaplacian::Faces faces;
            for (auto& face : faces) { face = {1, 0, 0}; }
            c.op->setRobinBC(faces);
        }
        c.op->setOversetMask(0, mask, 0);
        c.mlmg = std::make_unique<amrex::MLMG>(*c.op);
        c.rhs = std::make_unique<amrex::MultiFab>(nodal_grids, phi.DistributionMap(), 1, 0);
        if (m_diagnostics) {
            c.residual = std::make_unique<amrex::MultiFab>(
                nodal_grids, phi.DistributionMap(), 1, 0);
            amrex::Print() << "Material Poisson cache: max coarsening level = " << cap << '\n';
        }
    }
    auto& c = *cache;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(c.epsilon == &epsilon && c.mask == &mask &&
        c.rhs->boxArray() == phi.boxArray() && c.rhs->DistributionMap() == phi.DistributionMap() &&
        c.geom.Domain() == geom.Domain() &&
        c.options.max_coarsening_level == options.max_coarsening_level &&
        c.options.agglomeration == options.agglomeration &&
        c.options.consolidation == options.consolidation &&
        c.options.agglomeration_grid_size == options.agglomeration_grid_size &&
        c.options.consolidation_grid_size == options.consolidation_grid_size,
        "Material Poisson cached layout, material, mask or hierarchy options changed.");
    for (int d = 0; d < 3; ++d) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(c.geom.ProbLo(d) == geom.ProbLo(d) &&
            c.geom.ProbHi(d) == geom.ProbHi(d), "Material Poisson cached geometry changed.");
    }
    c.op->validateStaticMaterial(epsilon, mask);
    c.op->setInternalDirichletValues(voltage);
    if (!robin_level) {
        // setPhiBC has evaluated the current boundary values. Preserve its
        // edge/corner precedence and any spatially varying nodal face data.
        phi.OverrideSync(geom.periodicity());
        for (int face = 0; face < 6; ++face) { c.op->setRobinForcing(face, phi); }
    }
    amrex::MultiFab::Copy(*c.rhs, rho, 0, 0, 1, 0);
    // The common hierarchy has already scaled by -1/epsilon0. Dual-volume
    // normalization is dimensionless and applies only to physical Robin nodes.
    if (robin_level) { normalizeBoundarySource(*c.rhs, geom); }
    // Coverage, physical Dirichlet faces and electrodes all constrain only the
    // private solve RHS. The public deposition remains byte-for-byte unchanged.
    c.op->prepareForSolve();
    c.op->setDirichletNodesToZero(0, 0, *c.rhs);
    c.op->applyBC(0, 0, phi, amrex::MLLinOp::BCMode::Inhomogeneous,
                  amrex::MLLinOp::StateMode::Solution);
    if (m_diagnostics) { amrex::MultiFab::Copy(*c.residual, phi, 0, 0, 1, 0); }
    c.rhs->OverrideSync(geom.periodicity());
    auto effective_options = options;
    if (m_diagnostics) {
        // AMReX exposes bottom iteration counts but no stored return codes.
        // Verbosity >= 2 retains its bottom failure/fallback messages.
        effective_options.verbosity = std::max(2, options.verbosity);
        effective_options.bottom_verbosity = std::max(1, options.bottom_verbosity);
    }
    ablastr::fields::applyMLMGOptions(*c.mlmg, effective_options);
    // Keep coarse-interpolated exterior ghosts in the caller's phi. MLMG's
    // private nodal solution can have zero ghosts outside a local patch.
    c.mlmg->setFinalFillBC(false);
    c.mlmg->setConvergenceNormType(amrex::MLMGNormType::greater);
    auto const abs_tol = absolute_tolerance == 0 && c.rhs->norm0() == 0
        ? amrex::Real(1.e-6) : absolute_tolerance;
    amrex::Gpu::streamSynchronize();
    auto const solve_start = amrex::second();
    c.mlmg->solve({&phi}, {c.rhs.get()}, options.relative_tolerance, abs_tol);
    amrex::Gpu::streamSynchronize();
    auto const solve_end = amrex::second();
    ++c.solves;
    if (m_diagnostics || options.verbosity > 0) {
        amrex::Print() << "Material Poisson level " << lev << " ("
                       << (robin_level ? "monopole" : "dirichlet") << ") solve " << c.solves
                       << ": MG levels = " << c.op->NMGLevels(0)
                       << ", coarsening cap = " << c.coarsening_cap
                       << ", pre/post = " << options.num_pre_smooth
                       << '/' << options.num_post_smooth
                       << ", bottom solver = " << bottomSolverName(c.mlmg->getBottomSolver())
                       << ", iterations = " << c.mlmg->getNumIters()
                       << ", initial/final residual = " << c.mlmg->getInitResidual()
                       << '/' << c.mlmg->getFinalResidual()
                       << ", preparation seconds = " << solve_start - prepare_start
                       << ", solve seconds = " << solve_end - solve_start << '\n';
    }
    if (m_diagnostics) {
        auto const diagnostic_start = amrex::second();
        amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum> ops;
        amrex::ReduceData<amrex::Real, amrex::Long> data(ops);
        using Tuple = decltype(data)::Type;
        for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi) {
            auto const p = phi.const_array(mfi);
            auto const before = c.residual->const_array(mfi);
            auto const fixed = c.op->dirichletMask().const_array(mfi);
            auto const valid = mfi.validbox();
            ops.eval(valid, data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> Tuple {
                    bool const seam = !valid.strictly_contains(amrex::IntVect(i,j,k));
                    return {fixed(i,j,k) ? amrex::Math::abs(p(i,j,k)-before(i,j,k)) : 0,
                            amrex::Long(seam && !fixed(i,j,k))};
                });
        }
        auto const result = data.value();
        auto fixed_error = amrex::get<0>(result);
        auto seam_unknowns = amrex::get<1>(result);
        amrex::ParallelDescriptor::ReduceRealMax(fixed_error);
        amrex::ParallelDescriptor::ReduceLongSum(seam_unknowns);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(fixed_error == 0,
            "Material solve changed a prescribed physical, electrode or coarse/fine node.");
        amrex::Print() << "Material Poisson level " << lev
                       << ": fixed-node error = " << fixed_error
                       << ", unknown Box-boundary entries = " << seam_unknowns << '\n';
        c.op->compPhysicalResidual(*c.residual, phi, *c.rhs);
        auto const residual = c.residual->norm0();
        amrex::Gpu::streamSynchronize();
        amrex::Print() << "Material Poisson physical residual = " << residual
                       << ", diagnostic seconds = " << amrex::second() - diagnostic_start
                       << ", coefficient builds = " << c.op->coefficientBuilds()
                       << ", fixed nodes by MG level =";
        for (int mg = 0; mg < c.op->NMGLevels(0); ++mg) {
            amrex::Print() << ' ' << c.op->internalDirichletNodes(mg);
        }
        amrex::Print() << ", bottom Krylov iterations =";
        for (auto count : c.mlmg->getNumCGIters()) { amrex::Print() << ' ' << count; }
        amrex::Print() << '\n';
    }
    return *c.mlmg;
}

void MaterialPoissonSolver::normalizeBoundarySource (
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

void MaterialPoissonSolver::fillPhiGhosts (
    amrex::MultiFab& phi, amrex::Geometry const& geom) const
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_monopole && phi.nGrowVect().min() >= 1,
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
