#include "HallInjector.h"

#include "WarpX.H"

#include "Insert/Config/WarpXFunctionConfig.h"
#include "Insert/Injection/HallCoordinateDistribution.h"
#include "Insert/Utils/InsertUtils.h"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"

#ifdef IONIZATION_SOURCE_INJECT
#include "Utils/WarpXConst.H"
#include "Insert/Injection/IonizationSourceSampler.h"
#endif

#include <AMReX_ParmParse.H>

#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Insert {
namespace {

#ifdef IONIZATION_SOURCE_INJECT
constexpr amrex::ParticleReal BirthElectronTemperatureEV = 1.0;

class IonizationSourceRateModel final : public HallRateModel
{
public:
    IonizationSourceRateModel (
        amrex::Real total_rate, amrex::ParticleReal macro_weight)
        : m_total_rate(total_rate), m_macro_weight(macro_weight)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_total_rate >= 0.0,
            "IonizationSourceRateModel requires non-negative A_tot.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_macro_weight > amrex::ParticleReal(0.0),
            "IonizationSourceRateModel requires positive macro_weight.");
    }

    [[nodiscard]] amrex::Real
    expectedMacroParticles (amrex::Real dt) const override
    {
        return m_total_rate * dt / static_cast<amrex::Real>(m_macro_weight);
    }

private:
    amrex::Real m_total_rate;
    amrex::ParticleReal m_macro_weight;
};
#endif

bool
HasParameter (std::string const& name)
{
    amrex::ParmParse pp;
    return pp.contains(name);
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
            QueryWithParser<amrex::Real>(
                pp, species_prefix, "weight",
                QueryWithParser<amrex::Real>(
                    pp, source_name, "macro_weight", 1.0)));
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
        QueryWithParser<amrex::Real>(pp, source_name, "x_offset", 0.0));
    const auto y_offset = static_cast<amrex::ParticleReal>(
        QueryWithParser<amrex::Real>(pp, source_name, "y_offset", 0.0));
    HallInjectionSource source(
        source_name, std::move(rate_model), std::move(position),
        std::move(species_configs), x_offset, y_offset);

    amrex::Real batch_multiplier = 0.0;
    if (utils::parser::queryWithParser(
            pp, source_name, "batch_multiplier", batch_multiplier)) {
        source.setBatchMultiplier(batch_multiplier);
    }

    const auto position_prefix = source_name + ".position";
    std::string coupled_distribution;
    if (utils::parser::query(pp, position_prefix, "coupled_distribution",
                             coupled_distribution)) {
        HallHoleArrayPlaneConfig hole_config;
        hole_config.hole_count =
            QueryWithParser<int>(pp, source_name, "hole_count", 48);
        hole_config.hole_radius = static_cast<amrex::ParticleReal>(
            GetWithParser<amrex::Real>(pp, source_name, "hole_radius"));
        amrex::Real ring_radius = 0.0;
        if (!utils::parser::queryWithParser(
                pp, source_name, "hole_ring_radius", ring_radius)) {
            ring_radius =
                GetWithParser<amrex::Real>(pp, source_name, "ring_radius");
        }
        hole_config.ring_radius =
            static_cast<amrex::ParticleReal>(ring_radius);
        hole_config.z = static_cast<amrex::ParticleReal>(
            QueryWithParser<amrex::Real>(pp, source_name, "z", 0.0));
        source.setHoleArrayPlane(hole_config);
    }

    return source;
}

#ifdef IONIZATION_SOURCE_INJECT
HallInjectionSource
MakeIonizationSource ();
#endif

void
AppendConfiguredSources (
    std::vector<HallInjectionSource>& sources,
    std::vector<std::string> const& source_names)
{
    amrex::ParmParse pp;
    for (auto const& source_name : source_names) {
#ifdef IONIZATION_SOURCE_INJECT
        if (source_name == "average_ionization_source") {
            sources.push_back(MakeIonizationSource());
            continue;
        }
#endif
        sources.push_back(MakeConfiguredSource(pp, source_name));
    }
}

#ifdef IONIZATION_SOURCE_INJECT
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

HallInjectionSource
MakeIonizationSource ()
{
    amrex::ParmParse pp_mc("my_constants");
    amrex::ParticleReal elec_weight = 0.0;
    pp_mc.getWithParser("elec_weight", elec_weight);

    auto position_sampler = std::make_unique<IonizationSourceSampler>();
    const auto source_weight = position_sampler->electronWeight();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::abs(static_cast<amrex::Real>(source_weight - elec_weight)) <=
            1.0e-12 *
                std::max(1.0, std::abs(static_cast<amrex::Real>(elec_weight))),
        "IONIZATION_SOURCE_INJECT requires my_constants.elec_weight to match "
        "ionization_source_fab/metadata.txt.");

    const auto electron_sigma = static_cast<amrex::ParticleReal>(
        std::sqrt(PhysConst::q_e * BirthElectronTemperatureEV / PhysConst::m_e));

    std::vector<HallSpeciesVelocityConfig> species;
    species.push_back(MakeSpecies(
        "electrons", elec_weight,
        MakeCartesianGaussianVelocityDistribution(
            electron_sigma, electron_sigma, electron_sigma)));
    species.push_back(MakeSpecies(
        "xe_ions", elec_weight,
        MakeCartesianGaussianVelocityDistribution(
            amrex::ParticleReal(1212.41), amrex::ParticleReal(1212.41),
            amrex::ParticleReal(1212.41))));

    const auto total_rate = position_sampler->totalRate();
    return HallInjectionSource(
        "average_ionization_source",
        std::make_unique<IonizationSourceRateModel>(total_rate, elec_weight),
        std::move(position_sampler), std::move(species));
}
#endif

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

    if (!initial_source_names.empty()) {
        AppendConfiguredSources(m_initial_sources, initial_source_names);
    }

    if (!continuous_source_names.empty()) {
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
