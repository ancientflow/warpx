#include "HallInjectionSource.h"

#include "WarpX.H"

#include "Particles/MultiParticleContainer.H"
#include "Utils/TextMsg.H"

#include <AMReX_Math.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <utility>

namespace Insert {
namespace {

amrex::RandomEngine
MakeRandomEngine ()
{
#ifdef AMREX_USE_GPU
    return amrex::RandomEngine(nullptr);
#else
    return amrex::RandomEngine{};
#endif
}

amrex::ParticleReal
TwoPi ()
{
    return amrex::ParticleReal(2.0) * amrex::Math::pi<amrex::ParticleReal>();
}

int
RandomHoleStart (int hole_count, amrex::RandomEngine const& engine)
{
    const int start = static_cast<int>(
        amrex::Random(engine) * static_cast<amrex::ParticleReal>(hole_count));
    return std::min(start, hole_count - 1);
}

} // namespace

HallInjectionSource::HallInjectionSource (
    std::string source_name,
    std::unique_ptr<HallRateModel> rate_model,
    std::unique_ptr<HallCoordinateDistribution> position_space,
    std::vector<HallSpeciesVelocityConfig> species,
    amrex::ParticleReal x_offset,
    amrex::ParticleReal y_offset)
    : m_source_name(std::move(source_name)),
      m_rate_model(std::move(rate_model)),
      m_position_space(std::move(position_space)),
      m_species(std::move(species)),
      m_x_offset(x_offset),
      m_y_offset(y_offset)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_rate_model != nullptr, "HallInjectionSource requires a rate model.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_position_space != nullptr,
        "HallInjectionSource requires a position distribution.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_species.empty(),
        "HallInjectionSource requires at least one species.");
    for (auto const& species_config : m_species) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            species_config.velocity_space != nullptr,
            "HallInjectionSource species requires a velocity distribution.");
    }
}

void
HallInjectionSource::setHoleArrayPlane (HallHoleArrayPlaneConfig config)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        config.hole_count > 0,
        "Hall hole_array_plane requires hole_count > 0.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        config.ring_radius >= amrex::ParticleReal(0.0),
        "Hall hole_array_plane requires ring_radius >= 0.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        config.hole_radius > amrex::ParticleReal(0.0),
        "Hall hole_array_plane requires hole_radius > 0.");
    m_hole_array_plane = config;
    m_use_hole_array_plane = true;
}

void
HallInjectionSource::setTimeStepOverride (amrex::Real dt)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(dt) && dt >= 0.0,
        "HallInjectionSource time step override must be finite and non-negative.");
    m_dt_override = dt;
    m_has_dt_override = true;
}

void
HallInjectionSource::setBatchMultiplier (amrex::Real batch_multiplier)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(batch_multiplier) && batch_multiplier >= 1.0,
        "HallInjectionSource batch multiplier must be larger than 1");
    m_batch_multiplier = batch_multiplier;
}

void
HallInjectionSource::Inject (WarpX& warpx, amrex::Real dt)
{
    const int count = consumeParticleCount(
        m_rate_model->expectedMacroParticles(m_has_dt_override ? m_dt_override : dt));
    if (count <= 0) {
        return;
    }

    amrex::RandomEngine uniform_engine(MakeRandomEngine());
    amrex::RandomEngine normal_engine(MakeRandomEngine());
    if (m_use_hole_array_plane) {
        m_hole_start = RandomHoleStart(m_hole_array_plane.hole_count, uniform_engine);
    }

    std::vector<EmissionSample> positions;
    positions.reserve(count);
    for (int i = 0; i < count; ++i) {
        positions.push_back(samplePosition(i, uniform_engine));
    }

    for (auto const& species : m_species) {
        addSpeciesParticles(warpx, species, positions, normal_engine);
    }
}

std::string const&
HallInjectionSource::name () const noexcept
{
    return m_source_name;
}

int
HallInjectionSource::consumeParticleCount (amrex::Real expected_macro_particles)
{
    if (m_batch_multiplier <= 0.0) {
        return m_accumulator.consume(expected_macro_particles);
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        expected_macro_particles >= 0.0,
        "HallInjectionSource requires non-negative particle counts.");
    const auto batch_count = static_cast<int>(
        expected_macro_particles * m_batch_multiplier + 1.0);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        batch_count > 0,
        "HallInjectionSource batch threshold must be positive.");

    m_batch_remainder += expected_macro_particles;
    if (m_batch_remainder > static_cast<amrex::Real>(batch_count)) {
        m_batch_remainder -= static_cast<amrex::Real>(batch_count);
        return batch_count;
    }
    return 0;
}

EmissionSample
HallInjectionSource::samplePosition (
    int particle_index, amrex::RandomEngine& engine)
{
    EmissionSample sample;
    if (m_use_hole_array_plane) {
        const auto radius =
            m_hole_array_plane.hole_radius * std::sqrt(amrex::Random(engine));
        const auto theta = TwoPi() * amrex::Random(engine);
        const int hole_index =
            (particle_index + m_hole_start) % m_hole_array_plane.hole_count;
        const auto hole_theta =
            static_cast<amrex::ParticleReal>(hole_index) * TwoPi() /
            static_cast<amrex::ParticleReal>(m_hole_array_plane.hole_count);

        sample.x =
            radius * std::cos(theta) +
            m_hole_array_plane.ring_radius * std::cos(hole_theta);
        sample.y =
            radius * std::sin(theta) +
            m_hole_array_plane.ring_radius * std::sin(hole_theta);
        sample.z = m_hole_array_plane.z;
    } else {
        sample = m_position_space->samplePosition(engine);
    }

    sample.x += m_x_offset;
    sample.y += m_y_offset;
    return sample;
}

void
HallInjectionSource::addSpeciesParticles (
    WarpX& warpx, HallSpeciesVelocityConfig const& species,
    std::vector<EmissionSample> const& positions,
    amrex::RandomEngine& normal_engine) const
{
    const auto count = static_cast<int>(positions.size());
    if (count <= 0) {
        return;
    }

    amrex::Vector<amrex::ParticleReal> px(count), py(count), pz(count);
    amrex::Vector<amrex::ParticleReal> vx(count), vy(count), vz(count);
    amrex::Vector<amrex::ParticleReal> pw(count);
    for (int i = 0; i < count; ++i) {
        auto const& position = positions[static_cast<std::size_t>(i)];
        const auto velocity =
            species.velocity_space->sampleVelocity(normal_engine, position);
        px[i] = position.x;
        py[i] = position.y;
        pz[i] = position.z;
        vx[i] = velocity.x;
        vy[i] = velocity.y;
        vz[i] = velocity.z;
        pw[i] = species.macro_weight * position.weight_factor;
    }

    static const amrex::Vector<amrex::Vector<int>> nattr;
    auto& mypc = warpx.GetPartContainer();
    auto& pc = mypc.GetParticleContainer(mypc.getSpeciesID(species.species_name));
    pc.AddNParticles(0, count, px, py, pz, vx, vy, vz, 1, {pw}, 0, nattr, 0);
}

} // namespace Insert
