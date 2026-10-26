#include "CoordinateTransform.h"

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

RotatingAxisFrame
MakeRotatingAxisFrame (amrex::XDim3 const& axis_at_theta0)
{
    const auto z_axis = Math::NormalizeChecked(axis_at_theta0);
    const auto denominator = amrex::Real(1.0) + z_axis.z;
    const auto threshold = amrex::Real(16.0) *
                           std::numeric_limits<amrex::Real>::epsilon();

    if (denominator <= threshold) {
        return RotatingAxisFrame{
            amrex::XDim3{1.0, 0.0, 0.0},
            amrex::XDim3{0.0, -1.0, 0.0},
            amrex::XDim3{0.0, 0.0, -1.0}};
    }

    const auto inverse_denominator = amrex::Real(1.0) / denominator;
    return RotatingAxisFrame{
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

CoordinateSystem
ParseCoordinateSystem (std::string value)
{
    value = ToLower(value);
    if (value == "cartesian") {
        return CoordinateSystem::cartesian;
    }
    if (value == "cylindrical") {
        return CoordinateSystem::cylindrical;
    }
    if (value == "local_normal") {
        return CoordinateSystem::local_normal;
    }
    if (value == "rotating_axis") {
        return CoordinateSystem::rotating_axis;
    }
    WARPX_ABORT_WITH_MESSAGE("Unknown coordinate system: " + value);
    return CoordinateSystem::cartesian;
}

std::array<std::string, 3>
CoordinateAxisNames (
    CoordinateSystem system, CoordinateSpace space)
{
    if (space == CoordinateSpace::position) {
        if (system == CoordinateSystem::cartesian) {
            return {"x", "y", "z"};
        }
        if (system == CoordinateSystem::cylindrical) {
            return {"r", "theta", "z"};
        }
        WARPX_ABORT_WITH_MESSAGE(
            "local_normal and rotating_axis are not valid position coordinate systems.");
    }

    if (system == CoordinateSystem::cartesian) {
        return {"vx", "vy", "vz"};
    }
    if (system == CoordinateSystem::cylindrical) {
        return {"vr", "vtheta", "vz"};
    }
    if (system == CoordinateSystem::rotating_axis) {
        return {"vx", "vy", "vz"};
    }
    return {"vnormal", "vt1", "vt2"};
}

EmissionSample
MakeEmissionSample (
    CoordinateSystem system, amrex::XDim3 const& coordinates)
{
    EmissionSample sample;
    if (system == CoordinateSystem::cartesian) {
        sample.x = static_cast<amrex::ParticleReal>(coordinates.x);
        sample.y = static_cast<amrex::ParticleReal>(coordinates.y);
        sample.z = static_cast<amrex::ParticleReal>(coordinates.z);
        return sample;
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        system == CoordinateSystem::cylindrical,
        "Only cartesian and cylindrical position coordinate systems are supported.");
    const auto planar = Math::PolarToCartesian(coordinates.x, coordinates.y);
    sample.x = static_cast<amrex::ParticleReal>(planar.x);
    sample.y = static_cast<amrex::ParticleReal>(planar.y);
    sample.z = static_cast<amrex::ParticleReal>(coordinates.z);
    return sample;
}

amrex::XDim3
TransformVelocityToCartesian (
    CoordinateSystem system, amrex::XDim3 const& velocity,
    EmissionSample const& position,
    RotatingAxisFrame const& rotating_axis_frame)
{
    if (system == CoordinateSystem::cartesian) {
        return velocity;
    }

    if (system == CoordinateSystem::cylindrical) {
        return Math::RotateAroundZ(velocity, LocalAzimuth(position));
    }

    if (system == CoordinateSystem::rotating_axis) {
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
