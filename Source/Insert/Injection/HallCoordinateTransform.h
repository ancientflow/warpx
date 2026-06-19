#pragma once

#include <AMReX_Dim3.H>
#include <AMReX_REAL.H>

#include <array>
#include <string>

namespace Insert {

enum class HallCoordinateSystem
{
    cartesian,
    cylindrical,
    local_normal
};

enum class HallCoordinateSpace
{
    position,
    velocity
};

struct EmissionSample
{
    amrex::ParticleReal x = amrex::ParticleReal(0.0);
    amrex::ParticleReal y = amrex::ParticleReal(0.0);
    amrex::ParticleReal z = amrex::ParticleReal(0.0);
    amrex::ParticleReal nx = amrex::ParticleReal(0.0);
    amrex::ParticleReal ny = amrex::ParticleReal(0.0);
    amrex::ParticleReal nz = amrex::ParticleReal(1.0);
    amrex::ParticleReal weight_factor = amrex::ParticleReal(1.0);
};

[[nodiscard]] HallCoordinateSystem
ParseHallCoordinateSystem (std::string value);

[[nodiscard]] std::array<std::string, 3>
HallCoordinateAxisNames (
    HallCoordinateSystem system, HallCoordinateSpace space);

[[nodiscard]] EmissionSample
MakeEmissionSample (
    HallCoordinateSystem system, amrex::XDim3 const& coordinates);

[[nodiscard]] amrex::XDim3
TransformVelocityToCartesian (
    HallCoordinateSystem system, amrex::XDim3 const& velocity,
    EmissionSample const& position);

} // namespace Insert
