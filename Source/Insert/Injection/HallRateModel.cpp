#include "HallRateModel.h"

#include "Insert/Utils/InsertUtils.h"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"

#include <AMReX.H>

#include <cmath>

namespace Insert {
namespace {

void
AssertPositive (amrex::Real value, std::string const& name)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(value > 0.0, name + " must be positive.");
}

} // namespace

HallFixedCountRateModel::HallFixedCountRateModel (amrex::Real count)
    : m_count(count)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_count >= 0.0, "HallFixedCountRateModel count must be non-negative.");
}

amrex::Real
HallFixedCountRateModel::expectedMacroParticles (amrex::Real dt) const
{
    amrex::ignore_unused(dt);
    return m_count;
}

HallDensityVolumeRateModel::HallDensityVolumeRateModel (
    amrex::Real density, amrex::Real volume, amrex::Real macro_weight)
    : m_density(density), m_volume(volume), m_macro_weight(macro_weight)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_density >= 0.0, "HallDensityVolumeRateModel density must be non-negative.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_volume >= 0.0, "HallDensityVolumeRateModel volume must be non-negative.");
    AssertPositive(m_macro_weight, "HallDensityVolumeRateModel macro_weight");
}

amrex::Real
HallDensityVolumeRateModel::expectedMacroParticles (amrex::Real dt) const
{
    amrex::ignore_unused(dt);
    return m_density * m_volume / m_macro_weight;
}

HallCurrentRateModel::HallCurrentRateModel (
    amrex::Real current, amrex::Real macro_weight, amrex::Real l_factor)
    : m_current(current), m_macro_weight(macro_weight), m_l_factor(l_factor)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_current >= 0.0, "HallCurrentRateModel current must be non-negative.");
    AssertPositive(m_macro_weight, "HallCurrentRateModel macro_weight");
    AssertPositive(m_l_factor, "HallCurrentRateModel l_factor");
}

amrex::Real
HallCurrentRateModel::expectedMacroParticles (amrex::Real dt) const
{
    return m_current * dt / (m_l_factor * m_l_factor * PhysConst::q_e * m_macro_weight);
}

HallMassFlowRateModel::HallMassFlowRateModel (
    amrex::Real mass_flow, amrex::Real particle_mass,
    amrex::Real macro_weight, amrex::Real l_factor)
    : m_mass_flow(mass_flow),
      m_particle_mass(particle_mass),
      m_macro_weight(macro_weight),
      m_l_factor(l_factor)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_mass_flow >= 0.0, "HallMassFlowRateModel mass_flow must be non-negative.");
    AssertPositive(m_particle_mass, "HallMassFlowRateModel particle_mass");
    AssertPositive(m_macro_weight, "HallMassFlowRateModel macro_weight");
    AssertPositive(m_l_factor, "HallMassFlowRateModel l_factor");
}

amrex::Real
HallMassFlowRateModel::expectedMacroParticles (amrex::Real dt) const
{
    return m_mass_flow * dt / (m_l_factor * m_l_factor * m_particle_mass * m_macro_weight);
}

int
HallFractionalParticleAccumulator::consume (
    amrex::Real expected_macro_particles)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        expected_macro_particles >= 0.0,
        "HallFractionalParticleAccumulator requires non-negative particle counts.");
    m_remainder += expected_macro_particles;
    const auto count = static_cast<int>(std::floor(m_remainder));
    m_remainder -= static_cast<amrex::Real>(count);
    return count;
}

amrex::Real
HallFractionalParticleAccumulator::remainder () const noexcept
{
    return m_remainder;
}

std::unique_ptr<HallRateModel>
MakeHallRateModel (
    amrex::ParmParse const& pp, std::string const& prefix)
{
    std::string rate;
    utils::parser::get(pp, prefix, "rate", rate);
    rate = ToLower(rate);

    if (rate == "fixed_count") {
        return std::make_unique<HallFixedCountRateModel>(
            GetWithParser<amrex::Real>(pp, prefix, "count"));
    }
    if (rate == "density_volume") {
        return std::make_unique<HallDensityVolumeRateModel>(
            GetWithParser<amrex::Real>(pp, prefix, "density"),
            GetWithParser<amrex::Real>(pp, prefix, "volume"),
            GetWithParser<amrex::Real>(pp, prefix, "macro_weight"));
    }
    if (rate == "current") {
        return std::make_unique<HallCurrentRateModel>(
            GetWithParser<amrex::Real>(pp, prefix, "current"),
            GetWithParser<amrex::Real>(pp, prefix, "macro_weight"),
            QueryWithParser<amrex::Real>(pp, prefix, "l_factor", 1.0));
    }
    if (rate == "mass_flow") {
        return std::make_unique<HallMassFlowRateModel>(
            GetWithParser<amrex::Real>(pp, prefix, "mass_flow"),
            GetWithParser<amrex::Real>(pp, prefix, "mass"),
            GetWithParser<amrex::Real>(pp, prefix, "macro_weight"),
            QueryWithParser<amrex::Real>(pp, prefix, "l_factor", 1.0));
    }

    WARPX_ABORT_WITH_MESSAGE("Unknown HallRateModel type: " + rate);
    return nullptr;
}

} // namespace Insert
