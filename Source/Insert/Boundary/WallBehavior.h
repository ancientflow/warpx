#pragma once

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

namespace Insert {

/** Common event vocabulary for parser-based and material collision models. */
enum class WallBehavior : int {
    none = -1, absorb, specular, diffuse, convert, secondary_electron_1,
    secondary_electron_2
};

constexpr int max_wall_behaviors = 6;
constexpr int max_wall_emissions = 2;

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int WallEmissionCount (WallBehavior const event) noexcept
{
    switch (event) {
    case WallBehavior::convert:
    case WallBehavior::secondary_electron_1:
        return 1;
    case WallBehavior::secondary_electron_2:
        return 2;
    default:
        return 0;
    }
}

} // namespace Insert
