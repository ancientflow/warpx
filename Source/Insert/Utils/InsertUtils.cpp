#include "InsertUtils.h"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Math.H>
#include <AMReX_Random.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cctype>

namespace Insert {

void
CreateDirectoryTree (std::string const& dir)
{
    if (dir.empty() || dir == ".") {
        return;
    }

    amrex::Vector<std::string> paths;
    std::string current;
    auto pos = std::string::size_type{0};
    if (dir[0] == '/') {
        current = "/";
        pos = 1;
    }
    while (pos < dir.size()) {
        auto const next = dir.find('/', pos);
        auto const part = dir.substr(pos, next - pos);
        if (!part.empty()) {
            if (!current.empty() && current != "/") {
                current += "/";
            }
            current += part;
            paths.push_back(current);
        }
        if (next == std::string::npos) {
            break;
        }
        pos = next + 1;
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        constexpr int permission_flag_rwxrxrx = 0755;
        for (auto const& path : paths) {
            if (!amrex::UtilCreateDirectory(path, permission_flag_rwxrxrx)) {
                amrex::CreateDirectoryFailed(path);
            }
        }
    }
    amrex::ParallelDescriptor::Barrier();
}

std::string
ParentPath (std::string const& path)
{
    auto const slash = path.find_last_of('/');
    if (slash == std::string::npos) {
        return ".";
    }
    if (slash == 0) {
        return "/";
    }
    return path.substr(0, slash);
}

std::string
PathJoin (std::string const& dir, std::string const& filename)
{
    if (!dir.empty() && dir.back() == '/') {
        return dir + filename;
    }
    return dir + "/" + filename;
}

std::string
ToLower (std::string value)
{
    std::transform(
        value.begin(), value.end(), value.begin(),
        [] (unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

amrex::RandomEngine
MakeRandomEngine ()
{
#ifdef AMREX_USE_GPU
    return amrex::RandomEngine(nullptr);
#else
    return amrex::RandomEngine{};
#endif
}

amrex::ParticleReal
TwoPi ()
{
    return amrex::ParticleReal(2.0) * amrex::Math::pi<amrex::ParticleReal>();
}

HallAnodeRingConfig
ReadHallAnodeRingConfig (amrex::Geometry const& geom)
{
    amrex::ParmParse pp_mc("my_constants");
    auto voltage = amrex::Real(0.0);
    auto l_factor = amrex::ParticleReal(1.0);
    auto length = static_cast<amrex::ParticleReal>(
        geom.ProbHi(0) - geom.ProbLo(0));
    pp_mc.query("voltage", voltage);
    pp_mc.query("l_factor", l_factor);
    pp_mc.query("L", length);

    amrex::ParticleReal const center_x =
        static_cast<amrex::ParticleReal>(geom.ProbLo(0)) +
        length / l_factor / amrex::ParticleReal(2.0);
    amrex::ParticleReal const center_y =
        static_cast<amrex::ParticleReal>(geom.ProbLo(1)) +
        length / l_factor / amrex::ParticleReal(2.0);
    amrex::ParticleReal const r_min =
        amrex::ParticleReal(0.021) / amrex::ParticleReal(2.0) / l_factor;
    amrex::ParticleReal const r_max =
        amrex::ParticleReal(0.031) / amrex::ParticleReal(2.0) / l_factor;

    return HallAnodeRingConfig{
        voltage,
        center_x,
        center_y,
        r_min,
        r_max,
        r_min * r_min,
        r_max * r_max};
}

} // namespace Insert
