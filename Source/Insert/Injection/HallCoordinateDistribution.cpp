#include "HallCoordinateDistribution.h"

#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"

#include <AMReX_Dim3.H>

#include <string>
#include <utility>
#include <vector>

namespace Insert {
namespace {

std::string
DefaultCoordinateSystem (HallCoordinateSpace space)
{
    return (space == HallCoordinateSpace::position) ? "cylindrical" : "cartesian";
}

} // namespace

HallCoordinateDistribution::HallCoordinateDistribution (
    HallCoordinateSystem coordinate_system,
    std::unique_ptr<HallDistribution1D> q0,
    std::unique_ptr<HallDistribution1D> q1,
    std::unique_ptr<HallDistribution1D> q2,
    HallRotatingAxisFrame rotating_axis_frame)
    : m_coordinate_system(coordinate_system),
      m_q0(std::move(q0)),
      m_q1(std::move(q1)),
      m_q2(std::move(q2)),
      m_rotating_axis_frame(rotating_axis_frame)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_q0 && m_q1 && m_q2,
        "HallCoordinateDistribution requires three direction distributions.");
}

HallCoordinateSystem
HallCoordinateDistribution::coordinateSystem () const noexcept
{
    return m_coordinate_system;
}

amrex::XDim3
HallCoordinateDistribution::sampleCoordinates (
    amrex::RandomEngine const& engine) const
{
    return amrex::XDim3{
        m_q0->sample(engine),
        m_q1->sample(engine),
        m_q2->sample(engine)};
}

EmissionSample
HallCoordinateDistribution::samplePosition (
    amrex::RandomEngine const& engine) const
{
    return MakeEmissionSample(m_coordinate_system, sampleCoordinates(engine));
}

amrex::XDim3
HallCoordinateDistribution::sampleVelocity (
    amrex::RandomEngine const& engine,
    EmissionSample const& position) const
{
    return TransformVelocityToCartesian(
        m_coordinate_system, sampleCoordinates(engine), position,
        m_rotating_axis_frame);
}

std::unique_ptr<HallCoordinateDistribution>
MakeHallCoordinateDistribution (
    amrex::ParmParse const& pp, std::string const& prefix,
    HallCoordinateSpace space)
{
    std::string coupled_distribution;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !utils::parser::query(pp, prefix, "coupled_distribution", coupled_distribution),
        "Hall coupled coordinate distributions are not implemented yet: " +
            coupled_distribution);

    std::string coordinate_system = DefaultCoordinateSystem(space);
    utils::parser::query(pp, prefix, "coordinate_system", coordinate_system);
    const auto system = ParseHallCoordinateSystem(coordinate_system);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        system != HallCoordinateSystem::rotating_axis ||
            space == HallCoordinateSpace::velocity,
        "rotating_axis is only valid for Hall velocity distributions.");
    const auto axes = HallCoordinateAxisNames(system, space);

    HallRotatingAxisFrame rotating_axis_frame;
    if (system == HallCoordinateSystem::rotating_axis) {
        std::vector<amrex::Real> axis_at_theta0;
        utils::parser::getArrWithParser(
            pp, prefix, "axis_at_theta0", axis_at_theta0);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            axis_at_theta0.size() == 3,
            prefix + ".axis_at_theta0 must contain exactly three components.");
        rotating_axis_frame = MakeHallRotatingAxisFrame(amrex::XDim3{
            axis_at_theta0[0], axis_at_theta0[1], axis_at_theta0[2]});
    }

    auto q0 = MakeHallDistribution1D(pp, prefix + "." + axes[0]);
    auto q1 = MakeHallDistribution1D(pp, prefix + "." + axes[1]);
    auto q2 = MakeHallDistribution1D(pp, prefix + "." + axes[2]);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        system != HallCoordinateSystem::rotating_axis ||
            q2->min() >= amrex::ParticleReal(0.0),
        prefix + ".vz must sample non-negative velocities for rotating_axis.");

    return std::make_unique<HallCoordinateDistribution>(
        system, std::move(q0), std::move(q1), std::move(q2),
        rotating_axis_frame);
}

} // namespace Insert
