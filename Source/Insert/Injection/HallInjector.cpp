#include "HallInjector.h"

#include "WarpX.H"

#include "Insert/Config/WarpXSimulationConfig.h"
#include "Insert/Injection/HallCoordinateDistribution.h"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Math.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Insert {
namespace {

constexpr amrex::ParticleReal L = 0.05;
constexpr amrex::ParticleReal XeMass = 2.179e-25;
constexpr amrex::ParticleReal Boltzmann = 1.38064852e-23;

amrex::ParticleReal
TwoPi ()
{
    return amrex::ParticleReal(2.0) * amrex::Math::pi<amrex::ParticleReal>();
}

amrex::Real
GetReal (
    amrex::ParmParse const& pp, std::string const& prefix,
    char const* name)
{
    amrex::Real value = 0.0;
    utils::parser::getWithParser(pp, prefix, name, value);
    return value;
}

amrex::Real
QueryReal (
    amrex::ParmParse const& pp, std::string const& prefix,
    char const* name, amrex::Real default_value)
{
    auto value = default_value;
    utils::parser::queryWithParser(pp, prefix, name, value);
    return value;
}

int
QueryInt (
    amrex::ParmParse const& pp, std::string const& prefix,
    char const* name, int default_value)
{
    auto value = default_value;
    utils::parser::queryWithParser(pp, prefix, name, value);
    return value;
}

bool
HasParameter (std::string const& name)
{
    amrex::ParmParse pp;
    return pp.contains(name);
}

std::unique_ptr<HallDistribution1D>
MakeUniformThetaDistribution ()
{
    return std::make_unique<HallUniformDistribution1D>(
        amrex::ParticleReal(0.0), TwoPi());
}

std::unique_ptr<HallDistribution1D>
MakePlasmaThetaDistribution ()
{
    amrex::ParmParse pp;
    constexpr char const* prefix = "channel_plasma.position.theta";
    if (HasParameter(std::string(prefix) + ".distribution")) {
        return MakeHallDistribution1D(pp, prefix);
    }
    return MakeUniformThetaDistribution();
}

std::unique_ptr<HallCoordinateDistribution>
MakeCylindricalPositionDistribution (
    amrex::ParticleReal rmin, amrex::ParticleReal rmax,
    std::unique_ptr<HallDistribution1D> theta_distribution,
    std::unique_ptr<HallDistribution1D> z_distribution)
{
    return std::make_unique<HallCoordinateDistribution>(
        HallCoordinateSystem::cylindrical,
        std::make_unique<HallAreaUniformDistribution1D>(rmin, rmax),
        std::move(theta_distribution), std::move(z_distribution));
}

std::unique_ptr<HallCoordinateDistribution>
MakeCartesianGaussianVelocityDistribution (
    amrex::ParticleReal sigma_x, amrex::ParticleReal sigma_y,
    amrex::ParticleReal sigma_z)
{
    return std::make_unique<HallCoordinateDistribution>(
        HallCoordinateSystem::cartesian,
        std::make_unique<HallGaussianDistribution1D>(
            amrex::ParticleReal(0.0), sigma_x),
        std::make_unique<HallGaussianDistribution1D>(
            amrex::ParticleReal(0.0), sigma_y),
        std::make_unique<HallGaussianDistribution1D>(
            amrex::ParticleReal(0.0), sigma_z));
}

std::unique_ptr<HallCoordinateDistribution>
MakeCartesianInflowVelocityDistribution (
    amrex::ParticleReal sigma_x, amrex::ParticleReal sigma_y,
    amrex::ParticleReal sigma_z, amrex::ParticleReal vz0)
{
    return std::make_unique<HallCoordinateDistribution>(
        HallCoordinateSystem::cartesian,
        std::make_unique<HallGaussianDistribution1D>(
            amrex::ParticleReal(0.0), sigma_x),
        std::make_unique<HallGaussianDistribution1D>(
            amrex::ParticleReal(0.0), sigma_y),
        std::make_unique<HallPositiveGaussianDistribution1D>(vz0, sigma_z));
}

std::unique_ptr<HallCoordinateDistribution>
MakeDummyPositionDistribution ()
{
    return std::make_unique<HallCoordinateDistribution>(
        HallCoordinateSystem::cartesian,
        std::make_unique<HallConstantDistribution1D>(amrex::ParticleReal(0.0)),
        std::make_unique<HallConstantDistribution1D>(amrex::ParticleReal(0.0)),
        std::make_unique<HallConstantDistribution1D>(amrex::ParticleReal(0.0)));
}

HallSpeciesVelocityConfig
MakeSpecies (
    std::string species_name, amrex::ParticleReal macro_weight,
    std::unique_ptr<HallCoordinateDistribution> velocity_space)
{
    HallSpeciesVelocityConfig species;
    species.species_name = std::move(species_name);
    species.macro_weight = macro_weight;
    species.velocity_space = std::move(velocity_space);
    return species;
}

std::vector<std::string>
ReadStringArray (amrex::ParmParse const& pp, char const* name)
{
    std::vector<std::string> values;
    pp.queryarr(name, values);
    return values;
}

std::vector<std::string>
ReadSourceSpecies (amrex::ParmParse const& pp, std::string const& prefix)
{
    std::vector<std::string> species;
    pp.queryarr(prefix + ".shared_species", species);
    if (species.empty()) {
        std::string single_species;
        utils::parser::query(pp, prefix, "species", single_species);
        if (!single_species.empty()) {
            species.push_back(single_species);
        }
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !species.empty(), "Hall injection source requires species.");
    return species;
}

std::unique_ptr<HallCoordinateDistribution>
MakeConfiguredPositionDistribution (
    amrex::ParmParse const& pp, std::string const& source_name)
{
    const auto prefix = source_name + ".position";
    std::string coupled_distribution;
    if (utils::parser::query(pp, prefix, "coupled_distribution",
                             coupled_distribution)) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            coupled_distribution == "hole_array_plane",
            "Unsupported Hall coupled position distribution: " +
                coupled_distribution);
        return MakeDummyPositionDistribution();
    }
    return MakeHallCoordinateDistribution(pp, prefix, HallCoordinateSpace::position);
}

HallInjectionSource
MakeConfiguredSource (amrex::ParmParse const& pp, std::string const& source_name)
{
    auto rate_model = MakeHallRateModel(pp, source_name);
    auto position = MakeConfiguredPositionDistribution(pp, source_name);

    std::vector<HallSpeciesVelocityConfig> species_configs;
    for (auto const& species_name : ReadSourceSpecies(pp, source_name)) {
        const auto species_prefix = source_name + "." + species_name;
        const auto weight = static_cast<amrex::ParticleReal>(
            QueryReal(pp, species_prefix, "weight",
                      QueryReal(pp, source_name, "macro_weight", 1.0)));
        const auto velocity_prefix =
            HasParameter(species_prefix + ".velocity.coordinate_system") ||
                    HasParameter(species_prefix + ".velocity.vx.distribution") ||
                    HasParameter(species_prefix + ".velocity.vnormal.distribution")
                ? species_prefix + ".velocity"
                : source_name + ".velocity";
        species_configs.push_back(MakeSpecies(
            species_name, weight,
            MakeHallCoordinateDistribution(
                pp, velocity_prefix, HallCoordinateSpace::velocity)));
    }

    const auto x_offset = static_cast<amrex::ParticleReal>(
        QueryReal(pp, source_name, "x_offset", 0.0));
    const auto y_offset = static_cast<amrex::ParticleReal>(
        QueryReal(pp, source_name, "y_offset", 0.0));
    HallInjectionSource source(
        source_name, std::move(rate_model), std::move(position),
        std::move(species_configs), x_offset, y_offset);

    const auto position_prefix = source_name + ".position";
    std::string coupled_distribution;
    if (utils::parser::query(pp, position_prefix, "coupled_distribution",
                             coupled_distribution)) {
        HallHoleArrayPlaneConfig hole_config;
        hole_config.hole_count = QueryInt(pp, source_name, "hole_count", 48);
        hole_config.hole_radius = static_cast<amrex::ParticleReal>(
            GetReal(pp, source_name, "hole_radius"));
        amrex::Real ring_radius = 0.0;
        if (!utils::parser::queryWithParser(
                pp, source_name, "hole_ring_radius", ring_radius)) {
            ring_radius = GetReal(pp, source_name, "ring_radius");
        }
        hole_config.ring_radius =
            static_cast<amrex::ParticleReal>(ring_radius);
        hole_config.z = static_cast<amrex::ParticleReal>(
            QueryReal(pp, source_name, "z", 0.0));
        source.setHoleArrayPlane(hole_config);
    }

    return source;
}

void
AppendConfiguredSources (
    std::vector<HallInjectionSource>& sources,
    std::vector<std::string> const& source_names)
{
    amrex::ParmParse pp;
    for (auto const& source_name : source_names) {
        sources.push_back(MakeConfiguredSource(pp, source_name));
    }
}

HallInjectionSource
MakeDefaultChannelPlasmaSource ()
{
    amrex::ParmParse pp_mc("my_constants");
    amrex::Real l_factor = 1.0;
    amrex::ParticleReal elec_weight = 0.0;
    pp_mc.get("l_factor", l_factor);
    pp_mc.getWithParser("elec_weight", elec_weight);

    const auto inv_l = amrex::ParticleReal(1.0) / l_factor;
    const auto r1 = amrex::ParticleReal(0.0105) * inv_l;
    const auto r2 = amrex::ParticleReal(0.0155) * inv_l;
    const auto z1 = amrex::ParticleReal(0.001) * inv_l;
    const auto z2 = amrex::ParticleReal(0.004) * inv_l;
    const auto dz = z2 - z1;
    const auto volume =
        (r2 * r2 - r1 * r1) * amrex::Math::pi<amrex::ParticleReal>() * dz;
    const auto half_l = L * inv_l / amrex::ParticleReal(2.0);

    std::vector<HallSpeciesVelocityConfig> species;
    species.push_back(MakeSpecies(
        "electrons", elec_weight,
        MakeCartesianGaussianVelocityDistribution(
            amrex::ParticleReal(592982.0), amrex::ParticleReal(592982.0),
            amrex::ParticleReal(592982.0))));
    species.push_back(MakeSpecies(
        "xe_ions", elec_weight,
        MakeCartesianGaussianVelocityDistribution(
            amrex::ParticleReal(1212.41), amrex::ParticleReal(1212.41),
            amrex::ParticleReal(1212.41))));

    return HallInjectionSource(
        "channel_plasma",
        std::make_unique<HallDensityVolumeRateModel>(1.0e18, volume, elec_weight),
        MakeCylindricalPositionDistribution(
            r1, r2, MakePlasmaThetaDistribution(),
            std::make_unique<HallUniformDistribution1D>(z1, z2)),
        std::move(species), half_l, half_l);
}

HallInjectionSource
MakeDefaultCathodeSource ()
{
    amrex::ParmParse pp_mc("my_constants");
    amrex::Real Ic = 0.0;
    amrex::Real l_factor = 1.0;
    amrex::ParticleReal elec_weight = 0.0;
    pp_mc.get("Ic", Ic);
    pp_mc.get("l_factor", l_factor);
    pp_mc.getWithParser("elec_weight", elec_weight);

    const auto inv_l = amrex::ParticleReal(1.0) / l_factor;
    const auto half_l = L * inv_l / amrex::ParticleReal(2.0);
    std::vector<HallSpeciesVelocityConfig> species;
    species.push_back(MakeSpecies(
        "electrons", elec_weight,
        MakeCartesianGaussianVelocityDistribution(
            amrex::ParticleReal(592982.0), amrex::ParticleReal(592982.0),
            amrex::ParticleReal(592982.0))));

    return HallInjectionSource(
        "cathode_electron",
        std::make_unique<HallCurrentRateModel>(Ic, elec_weight, l_factor),
        MakeCylindricalPositionDistribution(
            amrex::ParticleReal(0.016) * inv_l, amrex::ParticleReal(0.02) * inv_l,
            MakeUniformThetaDistribution(),
            std::make_unique<HallUniformDistribution1D>(
                amrex::ParticleReal(0.042) * inv_l,
                amrex::ParticleReal(0.046) * inv_l)),
        std::move(species), half_l, half_l);
}

HallInjectionSource
MakeDefaultXeNeutralSource ()
{
    amrex::ParmParse pp_mc("my_constants");
    amrex::Real l_factor = 1.0;
    amrex::Real m_dot = 0.0;
    amrex::ParticleReal atom_weight = 0.0;
    amrex::ParticleReal tx = 0.0;
    amrex::ParticleReal ty = 0.0;
    amrex::ParticleReal tz = 0.0;
    amrex::ParticleReal vz0 = 0.0;
    bool ifhole = false;
    pp_mc.get("l_factor", l_factor);
    pp_mc.getWithParser("xe_weight", atom_weight);
    pp_mc.get("m_dot", m_dot);
    pp_mc.query("ifhole", ifhole);
    pp_mc.get("Tx", tx);
    pp_mc.get("Ty", ty);
    pp_mc.get("Tz", tz);
    pp_mc.get("vz0", vz0);

    const auto inv_l = amrex::ParticleReal(1.0) / l_factor;
    const auto rmid =
        (amrex::ParticleReal(0.021) + amrex::ParticleReal(0.031)) * inv_l /
        amrex::ParticleReal(4.0);
    const auto rwidth = amrex::ParticleReal(0.001) * inv_l;
    const auto half_l = L * inv_l / amrex::ParticleReal(2.0);
    const auto sigmax = std::sqrt(Boltzmann * tx / XeMass);
    const auto sigmay = std::sqrt(Boltzmann * ty / XeMass);
    const auto sigmaz = std::sqrt(Boltzmann * tz / XeMass);

    std::vector<HallSpeciesVelocityConfig> species;
    species.push_back(MakeSpecies(
        "xe_netural", atom_weight,
        MakeCartesianInflowVelocityDistribution(sigmax, sigmay, sigmaz, vz0)));

    auto source = HallInjectionSource(
        "xe_neutral_inlet",
        std::make_unique<HallMassFlowRateModel>(
            m_dot, XeMass, atom_weight, l_factor),
        ifhole
            ? MakeDummyPositionDistribution()
            : MakeCylindricalPositionDistribution(
                  rmid - rwidth, rmid + rwidth, MakeUniformThetaDistribution(),
                  std::make_unique<HallConstantDistribution1D>(
                      amrex::ParticleReal(0.0))),
        std::move(species), half_l, half_l);

    if (ifhole) {
        HallHoleArrayPlaneConfig hole_config;
        hole_config.hole_count = 48;
        hole_config.ring_radius = rmid;
        hole_config.hole_radius = rwidth;
        hole_config.z = amrex::ParticleReal(0.0);
        source.setHoleArrayPlane(hole_config);
    }
    return source;
}

} // namespace

HallInjector&
HallInjector::GetInstance ()
{
    static HallInjector injector;
    return injector;
}

void
HallInjector::ReadParameters ()
{
    if (m_initialized) {
        return;
    }

    amrex::ParmParse pp_insert("insert");
    const auto initial_source_names = ReadStringArray(pp_insert, "initial_sources");
    const auto continuous_source_names =
        ReadStringArray(pp_insert, "continuous_sources");

    if (initial_source_names.empty()) {
#ifdef HALL3D
        m_initial_sources.push_back(MakeDefaultChannelPlasmaSource());
#endif
    } else {
        AppendConfiguredSources(m_initial_sources, initial_source_names);
    }

    if (continuous_source_names.empty()) {
#ifdef HALL3D
        m_continuous_sources.push_back(MakeDefaultCathodeSource());
#ifndef MCC_DENSITY_AVERAGE_USE
        m_continuous_sources.push_back(MakeDefaultXeNeutralSource());
#endif
#endif
#ifdef HALL3D_INIT
        auto fast_neutral_source = MakeDefaultXeNeutralSource();
        fast_neutral_source.setTimeStepOverride(5.6e-10);
        m_continuous_sources.push_back(std::move(fast_neutral_source));
#endif
    } else {
        AppendConfiguredSources(m_continuous_sources, continuous_source_names);
    }

    m_initialized = true;
}

void
HallInjector::InitializePlasma (WarpX& warpx)
{
    ReadParameters();
    for (auto& source : m_initial_sources) {
        source.Inject(warpx, amrex::Real(0.0));
    }
}

void
HallInjector::InjectParticles (WarpX& warpx, amrex::Real dt, int step)
{
    amrex::ignore_unused(step);
    ReadParameters();
    for (auto& source : m_continuous_sources) {
        source.Inject(warpx, dt);
    }
}

} // namespace Insert
