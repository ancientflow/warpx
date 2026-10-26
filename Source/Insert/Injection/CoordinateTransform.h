#pragma once

#include <AMReX_Dim3.H>
#include <AMReX_REAL.H>

#include <array>
#include <string>

namespace Insert {

enum class CoordinateSystem
{
    cartesian,
    cylindrical,
    local_normal,
    rotating_axis
};

enum class CoordinateSpace
{
    position,
    velocity
};

struct EmissionSample
{
    amrex::ParticleReal x = amrex::ParticleReal(0.0);
    amrex::ParticleReal y = amrex::ParticleReal(0.0);
    amrex::ParticleReal z = amrex::ParticleReal(0.0);
    amrex::ParticleReal x_offset = amrex::ParticleReal(0.0);
    amrex::ParticleReal y_offset = amrex::ParticleReal(0.0);
    amrex::ParticleReal nx = amrex::ParticleReal(0.0);
    amrex::ParticleReal ny = amrex::ParticleReal(0.0);
    amrex::ParticleReal nz = amrex::ParticleReal(1.0);
    amrex::ParticleReal weight_factor = amrex::ParticleReal(1.0);
};

struct RotatingAxisFrame
{
    amrex::XDim3 x_axis{1.0, 0.0, 0.0};
    amrex::XDim3 y_axis{0.0, 1.0, 0.0};
    amrex::XDim3 z_axis{0.0, 0.0, 1.0};
};

[[nodiscard]] RotatingAxisFrame
MakeRotatingAxisFrame (amrex::XDim3 const& axis_at_theta0);

[[nodiscard]] CoordinateSystem
ParseCoordinateSystem (std::string value);

[[nodiscard]] std::array<std::string, 3>
CoordinateAxisNames (
    CoordinateSystem system, CoordinateSpace space);

[[nodiscard]] EmissionSample
MakeEmissionSample (
    CoordinateSystem system, amrex::XDim3 const& coordinates);

[[nodiscard]] amrex::XDim3
TransformVelocityToCartesian (
    CoordinateSystem system, amrex::XDim3 const& velocity,
    EmissionSample const& position,
    RotatingAxisFrame const& rotating_axis_frame = {});

} // namespace Insert
