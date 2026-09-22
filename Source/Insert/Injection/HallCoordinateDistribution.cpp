#include "HallCoordinateDistribution.h"

#include "Insert/Utils/InsertUtils.h"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"

#include <AMReX_Dim3.H>
#include <AMReX_Random.H>

#include <cmath>
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

/** Sample one velocity axis as a (drift) Maxwellian whose mean/sigma are
 *  functions of the local radius r. Mirrors the constant-parameter
 *  semantics of HallGaussianDistribution1D/HallPositiveGaussianDistribution1D. */
amrex::ParticleReal
SampleRadialGaussianAxis (
    HallRadialGaussianOverride const& override_, amrex::ParticleReal r,
    amrex::RandomEngine const& engine)
{
    const amrex::ParticleReal sigma = (*override_.sigma_of_r)(r);
    amrex::ParticleReal mean = override_.mean_const;
    if (override_.mean_of_r.has_value()) {
        mean = (*override_.mean_of_r)(r);
    }
    if (sigma == amrex::ParticleReal(0.0)) {
        return mean;
    }
    amrex::ParticleReal value = amrex::RandomNormal(mean, sigma, engine);
    if (override_.positive) {
        while (value < amrex::ParticleReal(0.0)) {
            value = amrex::RandomNormal(mean, sigma, engine);
        }
    }
    return value;
}

amrex::ParticleReal
SampleVelocityAxis (
    HallDistribution1D const& base,
    std::optional<HallRadialGaussianOverride> const& override_,
    amrex::ParticleReal r, amrex::RandomEngine const& engine)
{
    if (!override_.has_value()) {
        return base.sample(engine);
    }
    return SampleRadialGaussianAxis(*override_, r, engine);
}

/** Parse optional <axis_prefix>.sigma(r) / .mean(r) parser functions.
 *  Returns std::nullopt when neither is present. */
std::optional<HallRadialGaussianOverride>
ParseRadialGaussianOverride (
    amrex::ParmParse const& pp, std::string const& axis_prefix)
{
    std::string sigma_expr;
    const bool has_sigma = utils::parser::Query_parserString(
        pp, axis_prefix + ".sigma(r)", sigma_expr);
    std::string mean_expr;
    const bool has_mean = utils::parser::Query_parserString(
        pp, axis_prefix + ".mean(r)", mean_expr);
    if (!has_sigma && !has_mean) {
        return std::nullopt;
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        has_sigma,
        axis_prefix + ": mean(r) requires sigma(r) for radially dependent "
                      "Gaussian velocity sampling.");

    std::string distribution;
    utils::parser::get(pp, axis_prefix, "distribution", distribution);
    distribution = ToLower(distribution);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        distribution == "gaussian" || distribution == "positive_gaussian",
        axis_prefix + ": sigma(r)/mean(r) require the gaussian or "
                      "positive_gaussian distribution, got '" + distribution +
            "'.");

    HallRadialGaussianOverride override_;
    override_.positive = (distribution == "positive_gaussian");
    override_.sigma_parser = utils::parser::makeParser(sigma_expr, {"r"});
    override_.sigma_of_r = override_.sigma_parser.compile<1>();
    if (has_mean) {
        override_.mean_parser = utils::parser::makeParser(mean_expr, {"r"});
        override_.mean_of_r = override_.mean_parser.compile<1>();
    } else {
        utils::parser::queryWithParser(
            pp, axis_prefix, "mean", override_.mean_const);
    }
    return override_;
}

} // namespace

HallCoordinateDistribution::HallCoordinateDistribution (
    HallCoordinateSystem coordinate_system,
    std::unique_ptr<HallDistribution1D> q0,
    std::unique_ptr<HallDistribution1D> q1,
    std::unique_ptr<HallDistribution1D> q2,
    HallRotatingAxisFrame rotating_axis_frame,
    std::array<std::optional<HallRadialGaussianOverride>, 3>
        radial_velocity)
    : m_coordinate_system(coordinate_system),
      m_q0(std::move(q0)),
      m_q1(std::move(q1)),
      m_q2(std::move(q2)),
      m_rotating_axis_frame(rotating_axis_frame),
      m_radial_velocity(std::move(radial_velocity))
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_q0 && m_q1 && m_q2,
        "HallCoordinateDistribution requires three direction distributions.");
    for (auto const& override_ : m_radial_velocity) {
        m_has_radial_velocity =
            m_has_radial_velocity || override_.has_value();
    }
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
    const amrex::XDim3 velocity =
        m_has_radial_velocity ? sampleVelocityRadial(engine, position)
                              : sampleCoordinates(engine);
    return TransformVelocityToCartesian(
        m_coordinate_system, velocity, position, m_rotating_axis_frame);
}

amrex::XDim3
HallCoordinateDistribution::sampleVelocityRadial (
    amrex::RandomEngine const& engine,
    EmissionSample const& position) const
{
    // Radius about the local injection axis: positions carry the source
    // offsets already applied, so subtract them back.
    const amrex::ParticleReal dx = position.x - position.x_offset;
    const amrex::ParticleReal dy = position.y - position.y_offset;
    const amrex::ParticleReal r = std::sqrt(dx * dx + dy * dy);
    return amrex::XDim3{
        SampleVelocityAxis(*m_q0, m_radial_velocity[0], r, engine),
        SampleVelocityAxis(*m_q1, m_radial_velocity[1], r, engine),
        SampleVelocityAxis(*m_q2, m_radial_velocity[2], r, engine)};
}

std::unique_ptr<HallCoordinateDistribution>
MakeHallCoordinateDistribution (
    amrex::ParmParse const& pp, std::string const& prefix,
    HallCoordinateSpace space, std::string const& allowed_coupled_distribution)
{
    std::string coupled_distribution;
    if (utils::parser::query(pp, prefix, "coupled_distribution",
                             coupled_distribution)) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            !allowed_coupled_distribution.empty() &&
                ToLower(coupled_distribution) == allowed_coupled_distribution,
            "Unsupported Hall coupled coordinate distribution: " +
                coupled_distribution);
    }

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

    // Optional radial (r) dependence of Gaussian velocity parameters,
    // configured via <prefix>.<axis>.sigma(r) / .mean(r).
    std::array<std::optional<HallRadialGaussianOverride>, 3>
        radial_velocity;
    if (space == HallCoordinateSpace::velocity) {
        for (int i = 0; i < 3; ++i) {
            radial_velocity[static_cast<std::size_t>(i)] =
                ParseRadialGaussianOverride(
                    pp, prefix + "." + axes[static_cast<std::size_t>(i)]);
        }
    }

    return std::make_unique<HallCoordinateDistribution>(
        system, std::move(q0), std::move(q1), std::move(q2),
        rotating_axis_frame, std::move(radial_velocity));
}

} // namespace Insert
