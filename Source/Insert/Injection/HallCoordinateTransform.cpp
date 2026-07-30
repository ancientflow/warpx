#include "HallCoordinateTransform.h"

#include "Insert/Utils/InsertUtils.h"
#include "Utils/TextMsg.H"

#include <AMReX_Math.H>

#include <cmath>
#include <limits>

namespace Insert {
namespace {

amrex::XDim3
Normalize (amrex::XDim3 value)
{
    const auto norm = std::sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(norm) && norm > amrex::Real(0.0),
        "Cannot normalize a zero or non-finite vector.");
    return amrex::XDim3{value.x / norm, value.y / norm, value.z / norm};
}

amrex::XDim3
Cross (amrex::XDim3 const& a, amrex::XDim3 const& b) noexcept
{
    return amrex::XDim3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x};
}

amrex::Real
LocalAzimuth (EmissionSample const& position) noexcept
{
    return std::atan2(
        static_cast<amrex::Real>(position.y - position.y_offset),
        static_cast<amrex::Real>(position.x - position.x_offset));
}

amrex::XDim3
RotateAroundZ (amrex::XDim3 const& value, amrex::Real theta) noexcept
{
    const auto cos_theta = std::cos(theta);
    const auto sin_theta = std::sin(theta);
    return amrex::XDim3{
        value.x * cos_theta - value.y * sin_theta,
        value.x * sin_theta + value.y * cos_theta,
        value.z};
}

} // namespace

HallRotatingAxisFrame
MakeHallRotatingAxisFrame (amrex::XDim3 const& axis_at_theta0)
{
    const auto z_axis = Normalize(axis_at_theta0);
    const auto denominator = amrex::Real(1.0) + z_axis.z;
    const auto threshold = amrex::Real(16.0) *
                           std::numeric_limits<amrex::Real>::epsilon();

    if (denominator <= threshold) {
        return HallRotatingAxisFrame{
            amrex::XDim3{1.0, 0.0, 0.0},
            amrex::XDim3{0.0, -1.0, 0.0},
            amrex::XDim3{0.0, 0.0, -1.0}};
    }

    const auto inverse_denominator = amrex::Real(1.0) / denominator;
    return HallRotatingAxisFrame{
        amrex::XDim3{
            amrex::Real(1.0) - z_axis.x * z_axis.x * inverse_denominator,
            -z_axis.x * z_axis.y * inverse_denominator,
            -z_axis.x},
        amrex::XDim3{
            -z_axis.x * z_axis.y * inverse_denominator,
            amrex::Real(1.0) - z_axis.y * z_axis.y * inverse_denominator,
            -z_axis.y},
        z_axis};
}

HallCoordinateSystem
ParseHallCoordinateSystem (std::string value)
{
    value = ToLower(value);
    if (value == "cartesian") {
        return HallCoordinateSystem::cartesian;
    }
    if (value == "cylindrical") {
        return HallCoordinateSystem::cylindrical;
    }
    if (value == "local_normal") {
        return HallCoordinateSystem::local_normal;
    }
    if (value == "rotating_axis") {
        return HallCoordinateSystem::rotating_axis;
    }
    WARPX_ABORT_WITH_MESSAGE("Unknown Hall coordinate system: " + value);
    return HallCoordinateSystem::cartesian;
}

std::array<std::string, 3>
HallCoordinateAxisNames (
    HallCoordinateSystem system, HallCoordinateSpace space)
{
    if (space == HallCoordinateSpace::position) {
        if (system == HallCoordinateSystem::cartesian) {
            return {"x", "y", "z"};
        }
        if (system == HallCoordinateSystem::cylindrical) {
            return {"r", "theta", "z"};
        }
        WARPX_ABORT_WITH_MESSAGE(
            "local_normal and rotating_axis are not valid position coordinate systems.");
    }

    if (system == HallCoordinateSystem::cartesian) {
        return {"vx", "vy", "vz"};
    }
    if (system == HallCoordinateSystem::cylindrical) {
        return {"vr", "vtheta", "vz"};
    }
    if (system == HallCoordinateSystem::rotating_axis) {
        return {"vx", "vy", "vz"};
    }
    return {"vnormal", "vt1", "vt2"};
}

EmissionSample
MakeEmissionSample (
    HallCoordinateSystem system, amrex::XDim3 const& coordinates)
{
    EmissionSample sample;
    if (system == HallCoordinateSystem::cartesian) {
        sample.x = static_cast<amrex::ParticleReal>(coordinates.x);
        sample.y = static_cast<amrex::ParticleReal>(coordinates.y);
        sample.z = static_cast<amrex::ParticleReal>(coordinates.z);
        return sample;
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        system == HallCoordinateSystem::cylindrical,
        "Only cartesian and cylindrical position coordinate systems are supported.");
    const auto radius = coordinates.x;
    const auto theta = coordinates.y;
    sample.x = static_cast<amrex::ParticleReal>(radius * std::cos(theta));
    sample.y = static_cast<amrex::ParticleReal>(radius * std::sin(theta));
    sample.z = static_cast<amrex::ParticleReal>(coordinates.z);
    return sample;
}

amrex::XDim3
TransformVelocityToCartesian (
    HallCoordinateSystem system, amrex::XDim3 const& velocity,
    EmissionSample const& position,
    HallRotatingAxisFrame const& rotating_axis_frame)
{
    if (system == HallCoordinateSystem::cartesian) {
        return velocity;
    }

    if (system == HallCoordinateSystem::cylindrical) {
        return RotateAroundZ(velocity, LocalAzimuth(position));
    }

    if (system == HallCoordinateSystem::rotating_axis) {
        const amrex::XDim3 velocity_at_theta0{
            velocity.x * rotating_axis_frame.x_axis.x +
                velocity.y * rotating_axis_frame.y_axis.x +
                velocity.z * rotating_axis_frame.z_axis.x,
            velocity.x * rotating_axis_frame.x_axis.y +
                velocity.y * rotating_axis_frame.y_axis.y +
                velocity.z * rotating_axis_frame.z_axis.y,
            velocity.x * rotating_axis_frame.x_axis.z +
                velocity.y * rotating_axis_frame.y_axis.z +
                velocity.z * rotating_axis_frame.z_axis.z};
        return RotateAroundZ(velocity_at_theta0, LocalAzimuth(position));
    }

    const auto normal = Normalize(amrex::XDim3{
        static_cast<amrex::Real>(position.nx),
        static_cast<amrex::Real>(position.ny),
        static_cast<amrex::Real>(position.nz)});
    const amrex::XDim3 reference =
        (std::abs(normal.z) < amrex::Real(0.9))
            ? amrex::XDim3{0.0, 0.0, 1.0}
            : amrex::XDim3{1.0, 0.0, 0.0};
    const auto tangent1 = Normalize(Cross(reference, normal));
    const auto tangent2 = Cross(normal, tangent1);

    return amrex::XDim3{
        velocity.x * normal.x + velocity.y * tangent1.x + velocity.z * tangent2.x,
        velocity.x * normal.y + velocity.y * tangent1.y + velocity.z * tangent2.y,
        velocity.x * normal.z + velocity.y * tangent1.z + velocity.z * tangent2.z};
}

} // namespace Insert
