#include "HallCoordinateTransform.h"

#include "Insert/Math/Sampling.h"
#include "Insert/Math/VectorOps.h"
#include "Insert/Utils/InsertUtils.h"
#include "Utils/TextMsg.H"

#include <AMReX_Math.H>

#include <cmath>
#include <limits>

namespace Insert {
namespace {

amrex::Real
LocalAzimuth (EmissionSample const& position) noexcept
{
    return std::atan2(
        static_cast<amrex::Real>(position.y - position.y_offset),
        static_cast<amrex::Real>(position.x - position.x_offset));
}

} // namespace

HallRotatingAxisFrame
MakeHallRotatingAxisFrame (amrex::XDim3 const& axis_at_theta0)
{
    const auto z_axis = Math::NormalizeChecked(axis_at_theta0);
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
    const auto planar = Math::PolarToCartesian(coordinates.x, coordinates.y);
    sample.x = static_cast<amrex::ParticleReal>(planar.x);
    sample.y = static_cast<amrex::ParticleReal>(planar.y);
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
        return Math::RotateAroundZ(velocity, LocalAzimuth(position));
    }

    if (system == HallCoordinateSystem::rotating_axis) {
        const amrex::XDim3 velocity_at_theta0 = Math::ExpandInBasis(
            velocity, rotating_axis_frame.x_axis, rotating_axis_frame.y_axis,
            rotating_axis_frame.z_axis);
        return Math::RotateAroundZ(velocity_at_theta0, LocalAzimuth(position));
    }

    const auto normal = Math::NormalizeChecked(amrex::XDim3{
        static_cast<amrex::Real>(position.nx),
        static_cast<amrex::Real>(position.ny),
        static_cast<amrex::Real>(position.nz)});
    // The reference-vector construction is kept here (rather than
    // Math::BuildOrthonormalBasis) because its tangent sign convention
    // differs: user-specified tangential velocity components would change
    // sign under the wall-operator convention.
    const amrex::XDim3 reference =
        (std::abs(normal.z) < amrex::Real(0.9))
            ? amrex::XDim3{0.0, 0.0, 1.0}
            : amrex::XDim3{1.0, 0.0, 0.0};
    const auto tangent1 = Math::NormalizeChecked(Math::Cross(reference, normal));
    const auto tangent2 = Math::Cross(normal, tangent1);

    // Local components are (vnormal, vt1, vt2); expand them in the
    // (normal, tangent1, tangent2) basis to get Cartesian components.
    return Math::ExpandInBasis(velocity, normal, tangent1, tangent2);
}

} // namespace Insert
