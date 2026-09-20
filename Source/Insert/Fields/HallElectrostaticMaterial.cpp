#include "Insert/Fields/HallElectrostaticMaterial.H"

#include "EmbeddedBoundary/Enabled.H"
#include "Python/callbacks.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <AMReX_CoordSys.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <limits>
#include <string>

namespace Insert {

namespace {

/**
 * Initialize and validate the static electrostatic material fields on one level.
 *
 * The relative permittivity expression is evaluated at cell centers. The anode
 * implicit and potential expressions are evaluated at nodes, where a non-positive
 * implicit value produces an overset-mask value of zero. The routine also returns
 * the permittivity range, masked-node count, and masked z-index range used by the
 * caller to enforce the material invariants once when the cache is rebuilt.
 */
void
InitializeMaterialLevel (
    amrex::MultiFab& epsilon_r,
    amrex::iMultiFab& anode_mask,
    amrex::MultiFab& anode_phi,
    amrex::Geometry const& geom,
    amrex::ParserExecutor<3> const relative_permittivity,
    amrex::ParserExecutor<3> const anode_implicit,
    amrex::ParserExecutor<3> const anode_potential,
    amrex::Real& epsilon_min,
    amrex::Real& epsilon_max,
    int& min_masked_k,
    int& max_masked_k,
    amrex::Long& masked_nodes)
{
    auto const prob_lo = geom.ProbLoArray();
    auto const cell_size = geom.CellSizeArray();
    amrex::Dim3 const domain_lo = amrex::lbound(geom.Domain());

    for (amrex::MFIter mfi(epsilon_r, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& box = mfi.tilebox();
        auto const epsilon = epsilon_r.array(mfi);
        amrex::ParallelFor(
            box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real const x =
                    prob_lo[0] +
                    (i - domain_lo.x + static_cast<amrex::Real>(0.5)) * cell_size[0];
                amrex::Real const y =
                    prob_lo[1] +
                    (j - domain_lo.y + static_cast<amrex::Real>(0.5)) * cell_size[1];
                amrex::Real const z =
                    prob_lo[2] +
                    (k - domain_lo.z + static_cast<amrex::Real>(0.5)) * cell_size[2];
                epsilon(i, j, k) = relative_permittivity(x, y, z);
            });
    }

    for (amrex::MFIter mfi(anode_mask, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& box = mfi.tilebox();
        auto const mask = anode_mask.array(mfi);
        auto const fixed_phi = anode_phi.array(mfi);
        amrex::ParallelFor(
            box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real const x = prob_lo[0] + (i - domain_lo.x) * cell_size[0];
                amrex::Real const y = prob_lo[1] + (j - domain_lo.y) * cell_size[1];
                amrex::Real const z = prob_lo[2] + (k - domain_lo.z) * cell_size[2];
                bool const inside_anode = anode_implicit(x, y, z) <= 0.0;
                mask(i, j, k) = inside_anode ? 0 : 1;
                fixed_phi(i, j, k) = inside_anode ? anode_potential(x, y, z) : 0.0;
            });
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        epsilon_r.is_finite(0, 1, 0),
        "insert.relative_permittivity_function produced NaN or infinity.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        anode_phi.is_finite(0, 1, 0),
        "insert.anode_potential_function produced NaN or infinity inside the anode.");
    epsilon_min = epsilon_r.min(0);
    epsilon_max = epsilon_r.max(0);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        epsilon_min > 0.0,
        "insert.relative_permittivity_function must be strictly positive everywhere.");

    amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMax, amrex::ReduceOpSum>
        reduce_ops;
    amrex::ReduceData<int, int, amrex::Long> reduce_data(reduce_ops);
    using ReduceTuple = decltype(reduce_data)::Type;
    constexpr int no_masked_min = std::numeric_limits<int>::max();
    constexpr int no_masked_max = std::numeric_limits<int>::lowest();

    for (amrex::MFIter mfi(anode_mask, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& box = mfi.tilebox();
        auto const mask = anode_mask.const_array(mfi);
        reduce_ops.eval(
            box, reduce_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept -> ReduceTuple
            {
                bool const is_masked = mask(i, j, k) == 0;
                return {is_masked ? k : no_masked_min,
                        is_masked ? k : no_masked_max,
                        is_masked ? amrex::Long{1} : amrex::Long{0}};
            });
    }

    ReduceTuple const reduced = reduce_data.value();
    min_masked_k = amrex::get<0>(reduced);
    max_masked_k = amrex::get<1>(reduced);
    masked_nodes = amrex::get<2>(reduced);
    amrex::ParallelDescriptor::ReduceIntMin(min_masked_k);
    amrex::ParallelDescriptor::ReduceIntMax(max_masked_k);
    amrex::ParallelDescriptor::ReduceLongSum(masked_nodes);
}

} // namespace

HallElectrostaticMaterial&
HallElectrostaticMaterial::GetInstance ()
{
    static HallElectrostaticMaterial material;
    return material;
}

bool
HallElectrostaticMaterial::enabled ()
{
    readParametersOnce();
    return m_enabled;
}

void
HallElectrostaticMaterial::readParametersOnce ()
{
    if (m_parameters_read) {
        return;
    }

    amrex::ParmParse const pp("insert");
    int use_electrostatic_materials = 0;
    pp.query("use_electrostatic_materials", use_electrostatic_materials);
    m_enabled = use_electrostatic_materials != 0;

    if (!m_enabled) {
        m_parameters_read = true;
        return;
    }

#ifndef WARPX_DIM_3D
    WARPX_ABORT_WITH_MESSAGE(
        "insert.use_electrostatic_materials is currently supported only in 3D Cartesian geometry.");
#endif

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        WarpX::electrostatic_solver_id == ElectrostaticSolverAlgo::LabFrame ||
            WarpX::electrostatic_solver_id ==
                ElectrostaticSolverAlgo::LabFrameElectroMagnetostatic,
        "insert.use_electrostatic_materials requires a lab-frame electrostatic solver.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        WarpX::poisson_solver_id == PoissonSolverAlgo::Multigrid,
        "insert.use_electrostatic_materials requires warpx.poisson_solver = multigrid.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !EB::enabled(),
        "insert.use_electrostatic_materials cannot be combined with embedded boundaries.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !IsPythonCallbackInstalled("poissonsolver"),
        "insert.use_electrostatic_materials cannot be combined with a Python "
        "poissonsolver callback.");

    utils::parser::Store_parserString(
        pp, "relative_permittivity_function(x,y,z)", m_relative_permittivity_expression);
    utils::parser::Store_parserString(
        pp, "anode_implicit_function(x,y,z)", m_anode_implicit_expression);
    utils::parser::Store_parserString(
        pp, "anode_potential_function(x,y,z)", m_anode_potential_expression);

    m_relative_permittivity_parser = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_relative_permittivity_expression, {"x", "y", "z"}));
    m_anode_implicit_parser = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_anode_implicit_expression, {"x", "y", "z"}));
    m_anode_potential_parser = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(m_anode_potential_expression, {"x", "y", "z"}));

    m_relative_permittivity =
        utils::parser::compileParser<3>(m_relative_permittivity_parser.get());
    m_anode_implicit = utils::parser::compileParser<3>(m_anode_implicit_parser.get());
    m_anode_potential = utils::parser::compileParser<3>(m_anode_potential_parser.get());

    m_parameters_read = true;
}

void
HallElectrostaticMaterial::prepare (
    ablastr::fields::MultiLevelScalarField const& rho,
    ablastr::fields::MultiLevelScalarField const& phi,
    int const max_level)
{
    readParametersOnce();
    if (!m_enabled) {
        return;
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        max_level == 0,
        "insert.use_electrostatic_materials currently supports a single grid level only.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        static_cast<int>(rho.size()) == max_level + 1 &&
            static_cast<int>(phi.size()) == max_level + 1,
        "Electrostatic material preparation received an inconsistent number of grid levels.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        rho[0] != nullptr && phi[0] != nullptr,
        "Electrostatic material preparation requires valid rho and phi fields.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        phi[0]->ixType().nodeCentered(),
        "The electrostatic material anode mask requires a nodal phi field.");

    WarpX const& warpx = WarpX::GetInstance();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        warpx.Geom(0).Coord() == amrex::CoordSys::cartesian,
        "insert.use_electrostatic_materials requires Cartesian geometry.");

    if (!layoutMatches(phi, max_level)) {
        defineOrRebuild(phi, max_level);
    }

    for (int lev = 0; lev <= max_level; ++lev) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            rho[lev]->boxArray().CellEqual(m_anode_mask[lev]->boxArray()) &&
                rho[lev]->DistributionMap() == m_anode_mask[lev]->DistributionMap() &&
                rho[lev]->ixType().nodeCentered(),
            "The volume-anode mask requires rho and phi to share the nodal layout.");

        for (amrex::MFIter mfi(*m_anode_mask[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box const& box = mfi.tilebox();
            auto const mask = m_anode_mask[lev]->const_array(mfi);
            auto const fixed_phi = m_anode_phi[lev]->const_array(mfi);
            auto const phi_array = phi[lev]->array(mfi);
            auto const rho_array = rho[lev]->array(mfi);
            amrex::ParallelFor(
                box,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (mask(i, j, k) == 0) {
                        phi_array(i, j, k) = fixed_phi(i, j, k);
                        rho_array(i, j, k) = 0.0;
                    }
                });
        }
    }
}

bool
HallElectrostaticMaterial::layoutMatches (
    ablastr::fields::MultiLevelScalarField const& phi,
    int const max_level) const
{
    int const number_of_levels = max_level + 1;
    if (!m_defined || static_cast<int>(m_epsilon_r.size()) != number_of_levels ||
        static_cast<int>(m_anode_mask.size()) != number_of_levels ||
        static_cast<int>(m_anode_phi.size()) != number_of_levels ||
        static_cast<int>(m_domains.size()) != number_of_levels ||
        static_cast<int>(m_prob_lo.size()) != number_of_levels ||
        static_cast<int>(m_prob_hi.size()) != number_of_levels ||
        static_cast<int>(m_cell_size.size()) != number_of_levels)
    {
        return false;
    }

    WarpX const& warpx = WarpX::GetInstance();
    for (int lev = 0; lev < number_of_levels; ++lev) {
        amrex::Geometry const& geom = warpx.Geom(lev);
        auto const prob_lo = geom.ProbLoArray();
        auto const prob_hi = geom.ProbHiArray();
        auto const cell_size = geom.CellSizeArray();

        if (!m_epsilon_r[lev]->boxArray().CellEqual(warpx.boxArray(lev)) ||
            m_epsilon_r[lev]->DistributionMap() != warpx.DistributionMap(lev) ||
            !m_anode_mask[lev]->boxArray().CellEqual(phi[lev]->boxArray()) ||
            m_anode_mask[lev]->DistributionMap() != phi[lev]->DistributionMap() ||
            !m_anode_phi[lev]->boxArray().CellEqual(phi[lev]->boxArray()) ||
            m_anode_phi[lev]->DistributionMap() != phi[lev]->DistributionMap() ||
            m_domains[lev] != geom.Domain())
        {
            return false;
        }

        for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
            if (m_prob_lo[lev][dim] != prob_lo[dim] ||
                m_prob_hi[lev][dim] != prob_hi[dim] ||
                m_cell_size[lev][dim] != cell_size[dim])
            {
                return false;
            }
        }
    }
    return true;
}

void
HallElectrostaticMaterial::defineOrRebuild (
    ablastr::fields::MultiLevelScalarField const& phi,
    int const max_level)
{
    WarpX const& warpx = WarpX::GetInstance();
    int const number_of_levels = max_level + 1;

    m_defined = false;
    m_epsilon_r.clear();
    m_anode_mask.clear();
    m_anode_phi.clear();
    m_epsilon_r_ptrs.clear();
    m_anode_mask_ptrs.clear();
    m_epsilon_r.resize(number_of_levels);
    m_anode_mask.resize(number_of_levels);
    m_anode_phi.resize(number_of_levels);
    m_domains.resize(number_of_levels);
    m_prob_lo.resize(number_of_levels);
    m_prob_hi.resize(number_of_levels);
    m_cell_size.resize(number_of_levels);

    amrex::Real epsilon_min = std::numeric_limits<amrex::Real>::max();
    amrex::Real epsilon_max = std::numeric_limits<amrex::Real>::lowest();
    amrex::Long masked_nodes = 0;

    for (int lev = 0; lev < number_of_levels; ++lev) {
        amrex::Geometry const& geom = warpx.Geom(lev);
        auto const prob_lo = geom.ProbLoArray();
        auto const prob_hi = geom.ProbHiArray();
        auto const cell_size = geom.CellSizeArray();
        amrex::BoxArray const epsilon_ba = warpx.boxArray(lev);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            epsilon_ba.ixType().cellCentered(),
            "The relative permittivity cache requires a cell-centered WarpX grid.");

        m_epsilon_r[lev] = std::make_unique<amrex::MultiFab>(
            epsilon_ba, warpx.DistributionMap(lev), 1, 0);
        m_anode_mask[lev] = std::make_unique<amrex::iMultiFab>(
            phi[lev]->boxArray(), phi[lev]->DistributionMap(), 1, 0);
        m_anode_phi[lev] = std::make_unique<amrex::MultiFab>(
            phi[lev]->boxArray(), phi[lev]->DistributionMap(), 1, 0);

        amrex::Real level_epsilon_min;
        amrex::Real level_epsilon_max;
        int level_min_masked_k;
        int level_max_masked_k;
        amrex::Long level_masked_nodes;
        InitializeMaterialLevel(
            *m_epsilon_r[lev], *m_anode_mask[lev], *m_anode_phi[lev], geom,
            m_relative_permittivity, m_anode_implicit, m_anode_potential,
            level_epsilon_min, level_epsilon_max, level_min_masked_k,
            level_max_masked_k, level_masked_nodes);
        epsilon_min = std::min(epsilon_min, level_epsilon_min);
        epsilon_max = std::max(epsilon_max, level_epsilon_max);
        masked_nodes += level_masked_nodes;

        m_domains[lev] = geom.Domain();
        m_prob_lo[lev] = prob_lo;
        m_prob_hi[lev] = prob_hi;
        m_cell_size[lev] = cell_size;

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            level_masked_nodes > 0,
            "insert.anode_implicit_function did not select any nodal point.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            level_max_masked_k > level_min_masked_k,
            "The volume anode must cover at least two nodal layers in z.");
    }

    m_epsilon_r_ptrs.resize(number_of_levels);
    m_anode_mask_ptrs.resize(number_of_levels);
    for (int lev = 0; lev < number_of_levels; ++lev) {
        m_epsilon_r_ptrs[lev] = m_epsilon_r[lev].get();
        m_anode_mask_ptrs[lev] = m_anode_mask[lev].get();
    }
    m_defined = true;

    if (warpx.Verbose()) {
        amrex::Print() << "Hall electrostatic material initialized: epsilon_r in ["
                       << epsilon_min << ", " << epsilon_max << "], "
                       << masked_nodes << " masked nodal entries.\n"
                       << "  relative permittivity (dimensionless): "
                       << m_relative_permittivity_expression << "\n"
                       << "  anode implicit function: " << m_anode_implicit_expression << "\n"
                       << "  anode potential [V]: " << m_anode_potential_expression << "\n";
    }
}

ablastr::fields::ConstMultiLevelScalarField const*
HallElectrostaticMaterial::relativePermittivity () const noexcept
{
    return m_enabled && m_defined ? &m_epsilon_r_ptrs : nullptr;
}

amrex::Vector<amrex::iMultiFab const*> const*
HallElectrostaticMaterial::anodeMasks () const noexcept
{
    return m_enabled && m_defined ? &m_anode_mask_ptrs : nullptr;
}

} // namespace Insert
