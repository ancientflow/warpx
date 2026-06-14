#include "HallCoordinateTransform.h"

#include "Utils/TextMsg.H"

#include <AMReX_Math.H>

#include <algorithm>
#include <cctype>
#include <cmath>

namespace Insert {
namespace {

std::string
ToLower (std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

amrex::XDim3
Normalize (amrex::XDim3 value)
{
    const auto norm = std::sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        norm > amrex::Real(0.0), "Cannot normalize a zero vector.");
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

} // namespace

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
        WARPX_ABORT_WITH_MESSAGE("local_normal is not a valid position coordinate system.");
    }

    if (system == HallCoordinateSystem::cartesian) {
        return {"vx", "vy", "vz"};
    }
    if (system == HallCoordinateSystem::cylindrical) {
        return {"vr", "vtheta", "vz"};
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
    EmissionSample const& position)
{
    if (system == HallCoordinateSystem::cartesian) {
        return velocity;
    }

    if (system == HallCoordinateSystem::cylindrical) {
        const auto theta = std::atan2(
            static_cast<amrex::Real>(position.y),
            static_cast<amrex::Real>(position.x));
        const auto cos_theta = std::cos(theta);
        const auto sin_theta = std::sin(theta);
        return amrex::XDim3{
            velocity.x * cos_theta - velocity.y * sin_theta,
            velocity.x * sin_theta + velocity.y * cos_theta,
            velocity.z};
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
