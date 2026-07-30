#pragma once

#include "Utils/Parser/ParserUtils.H"

#include <AMReX_Extension.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Random.H>

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

int BoundaryParticleDiagInterval ();

bool DoBoundaryParticleDiag (int step);

std::string ParentPath (std::string const& path);

std::string PathJoin (std::string const& dir, std::string const& filename);

std::string ToLower (std::string value);

amrex::RandomEngine MakeRandomEngine ();

amrex::ParticleReal TwoPi ();

HallAnodeRingConfig ReadHallAnodeRingConfig (amrex::Geometry const& geom);

template <typename T>
T
GetWithParser (amrex::ParmParse const& pp, std::string const& prefix,
               char const* name)
{
    T value{};
    utils::parser::getWithParser(pp, prefix, name, value);
    return value;
}

template <typename T>
T
QueryWithParser (amrex::ParmParse const& pp, std::string const& prefix,
                 char const* name, T default_value)
{
    auto value = default_value;
    utils::parser::queryWithParser(pp, prefix, name, value);
    return value;
}

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
