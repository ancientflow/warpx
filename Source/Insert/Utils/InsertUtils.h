#pragma once

#include <AMReX_Extension.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <string>

namespace Insert {

struct HallAnodeRingConfig
{
    amrex::Real voltage = amrex::Real(0.0);
    amrex::ParticleReal center_x = amrex::ParticleReal(0.0);
    amrex::ParticleReal center_y = amrex::ParticleReal(0.0);
    amrex::ParticleReal r_min = amrex::ParticleReal(0.0);
    amrex::ParticleReal r_max = amrex::ParticleReal(0.0);
    amrex::ParticleReal r_min_sq = amrex::ParticleReal(0.0);
    amrex::ParticleReal r_max_sq = amrex::ParticleReal(0.0);
};

void CreateDirectoryTree (std::string const& dir);

HallAnodeRingConfig ReadHallAnodeRingConfig (amrex::Geometry const& geom);

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
IsHallAnodeRingHit (T const x, T const y,
                    HallAnodeRingConfig const config) noexcept
{
    amrex::ParticleReal const dx =
        static_cast<amrex::ParticleReal>(x) - config.center_x;
    amrex::ParticleReal const dy =
        static_cast<amrex::ParticleReal>(y) - config.center_y;
    amrex::ParticleReal const r_sq = dx * dx + dy * dy;
    return r_sq >= config.r_min_sq && r_sq <= config.r_max_sq;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
BacktraceParticleToZPlane (
    amrex::ParticleReal const x,
    amrex::ParticleReal const y,
    amrex::ParticleReal const z,
    amrex::ParticleReal const vx,
    amrex::ParticleReal const vy,
    amrex::ParticleReal const vz,
    amrex::ParticleReal const z_plane,
    amrex::ParticleReal& x_hit,
    amrex::ParticleReal& y_hit) noexcept
{
    x_hit = x;
    y_hit = y;

    if (vz == amrex::ParticleReal(0.0)) {
        return false;
    }

    amrex::ParticleReal const dt_back = (z - z_plane) / vz;
    if (dt_back < amrex::ParticleReal(0.0)) {
        return false;
    }

    x_hit = x - dt_back * vx;
    y_hit = y - dt_back * vy;
    return true;
}

} // namespace Insert
