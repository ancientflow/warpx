#include "Insert/Fields/HallElectrostaticMaterial.H"

#include "Utils/Parser/ParserUtils.H"
#include "WarpX.H"

#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>

namespace Insert {

namespace {

void RestoreAnodePotential (
    amrex::MultiFab& phi, amrex::iMultiFab const& mask, amrex::MultiFab const& anode_phi)
{
    for (amrex::MFIter mfi(phi, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const values = phi.array(mfi);
        auto const active = mask.const_array(mfi);
        auto const fixed = anode_phi.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (active(i, j, k) == 0) {
                    values(i, j, k) = fixed(i, j, k);
                }
            });
    }
    phi.OverrideSync();
}

void
FillAnodeMask (
    amrex::iMultiFab& mask, amrex::Geometry const& geom,
    amrex::ParserExecutor<3> const anode_implicit)
{
    auto const prob_lo = geom.ProbLoArray();
    auto const cell_size = geom.CellSizeArray();
    auto const domain_lo = amrex::lbound(geom.Domain());
    for (amrex::MFIter mfi(mask, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const values = mask.array(mfi);
        amrex::ParallelFor(mfi.tilebox(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real const x = prob_lo[0] + (i - domain_lo.x) * cell_size[0];
                amrex::Real const y = prob_lo[1] + (j - domain_lo.y) * cell_size[1];
                amrex::Real const z = prob_lo[2] + (k - domain_lo.z) * cell_size[2];
                values(i, j, k) = anode_implicit(x, y, z) <= 0.0 ? 0 : 1;
            });
    }
}

/**
 * Initialize the static electrostatic material fields on one level.
 *
 * The relative permittivity expression is evaluated at cell centers. The anode
 * implicit and potential expressions are evaluated at nodes, where a non-positive
 * implicit value produces an overset-mask value of zero.
 */
void
InitializeMaterialLevel (
    amrex::MultiFab& epsilon_r,
    amrex::iMultiFab& anode_mask,
    amrex::MultiFab& anode_phi,
    amrex::Geometry const& geom,
    amrex::ParserExecutor<3> const relative_permittivity,
    amrex::ParserExecutor<3> const anode_implicit,
    amrex::ParserExecutor<3> const anode_potential)
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

    FillAnodeMask(anode_mask, geom, anode_implicit);
    for (amrex::MFIter mfi(anode_mask, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& box = mfi.tilebox();
        auto const mask = anode_mask.const_array(mfi);
        auto const fixed_phi = anode_phi.array(mfi);
        amrex::ParallelFor(
            box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                amrex::Real const x = prob_lo[0] + (i - domain_lo.x) * cell_size[0];
                amrex::Real const y = prob_lo[1] + (j - domain_lo.y) * cell_size[1];
                amrex::Real const z = prob_lo[2] + (k - domain_lo.z) * cell_size[2];
                bool const inside_anode = mask(i, j, k) == 0;
                fixed_phi(i, j, k) = inside_anode ? anode_potential(x, y, z) : 0.0;
            });
    }
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
HallElectrostaticMaterial::fillAnodeMask (
    amrex::iMultiFab& mask, amrex::Geometry const& geom)
{
    readParametersOnce();
    FillAnodeMask(mask, geom, m_anode_implicit);
}

void
HallElectrostaticMaterial::initializeMaterialFields (
    amrex::MultiFab& epsilon_r, amrex::iMultiFab& anode_mask,
    amrex::MultiFab& anode_phi, amrex::Geometry const& geom)
{
    readParametersOnce();
    InitializeMaterialLevel(
        epsilon_r, anode_mask, anode_phi, geom,
        m_relative_permittivity, m_anode_implicit, m_anode_potential);
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

    if (!m_defined) {
        initialize(phi, max_level);
    }

    for (int lev = 0; lev <= max_level; ++lev) {
        for (amrex::MFIter mfi(*m_anode_mask[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box const& box = mfi.tilebox();
            auto const mask = m_anode_mask[lev]->const_array(mfi);
            //auto const fixed_phi = m_anode_phi[lev]->const_array(mfi);
            //auto const phi_array = phi[lev]->array(mfi);
            auto const rho_array = rho[lev]->array(mfi);
            amrex::ParallelFor(
                box,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (mask(i, j, k) == 0) {
                        // The overset mask keeps the initial anode potential fixed.
                        // phi_array(i, j, k) = fixed_phi(i, j, k);
                        rho_array(i, j, k) = 0.0;
                    }
                });
        }
    }
}

void
HallElectrostaticMaterial::restoreAnodePotential (amrex::MultiFab& phi, int const lev)
{
    RestoreAnodePotential(
        phi, *m_anode_mask[lev], *m_anode_phi[lev]);
}

void
HallElectrostaticMaterial::initialize (
    ablastr::fields::MultiLevelScalarField const& phi,
    int const max_level)
{
    WarpX const& warpx = WarpX::GetInstance();
    int const number_of_levels = max_level + 1;

    m_epsilon_r.resize(number_of_levels);
    m_anode_mask.resize(number_of_levels);
    m_anode_phi.resize(number_of_levels);

    for (int lev = 0; lev < number_of_levels; ++lev) {
        amrex::Geometry const& geom = warpx.Geom(lev);
        amrex::BoxArray const epsilon_ba = warpx.boxArray(lev);

        m_epsilon_r[lev] = std::make_unique<amrex::MultiFab>(
            epsilon_ba, warpx.DistributionMap(lev), 1, 0);
        m_anode_mask[lev] = std::make_unique<amrex::iMultiFab>(
            phi[lev]->boxArray(), phi[lev]->DistributionMap(), 1, 0);
        m_anode_phi[lev] = std::make_unique<amrex::MultiFab>(
            phi[lev]->boxArray(), phi[lev]->DistributionMap(), 1, 0);

        InitializeMaterialLevel(
            *m_epsilon_r[lev], *m_anode_mask[lev], *m_anode_phi[lev], geom,
            m_relative_permittivity, m_anode_implicit, m_anode_potential);

        for (amrex::MFIter mfi(*m_anode_mask[lev], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            amrex::Box const& box = mfi.tilebox();
            auto const mask = m_anode_mask[lev]->const_array(mfi);
            auto const fixed_phi = m_anode_phi[lev]->const_array(mfi);
            auto const phi_array = phi[lev]->array(mfi);
            amrex::ParallelFor(
                box,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (mask(i, j, k) == 0) {
                        phi_array(i, j, k) = fixed_phi(i, j, k);
                    }
                });
        }

    }

    m_epsilon_r_ptrs.resize(number_of_levels);
    m_anode_mask_ptrs.resize(number_of_levels);
    for (int lev = 0; lev < number_of_levels; ++lev) {
        m_epsilon_r_ptrs[lev] = m_epsilon_r[lev].get();
        m_anode_mask_ptrs[lev] = m_anode_mask[lev].get();
    }
    m_defined = true;

    if (warpx.Verbose()) {
        amrex::Print() << "Hall electrostatic material initialized.\n"
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
