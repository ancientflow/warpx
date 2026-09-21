#include "Insert/Boundary/AnalyticBoundaryInteraction.h"

#include "Insert/Boundary/AnalyticBoundaryGeometry.h"
#include "Insert/Boundary/BoundaryMathUtils.h"
#include "Insert/Core/WarpXInsert.h"
#include "Insert/Fields/HallWallCharge.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "WarpX.H"

#include <AMReX_BLassert.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelFor.H>

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Insert {
namespace {

#if defined(WARPX_DIM_3D)

struct AnalyticWall
{
    std::string name;
    std::unique_ptr<AnalyticBoundaryGeometry> geometry;
};

struct AnalyticWallConfiguration
{
    std::vector<AnalyticWall> walls;
    std::map<std::string, std::vector<bool>> absorbs;
};

[[nodiscard]] AnalyticWallConfiguration
ReadAnalyticWallConfiguration (MultiParticleContainer const& mpc)
{
    AnalyticWallConfiguration result;
    amrex::ParmParse const pp_insert("insert");
    std::vector<std::string> wall_names;
    pp_insert.queryarr("analytic_walls", wall_names);
    if (wall_names.empty()) {
        return result;
    }

    std::set<std::string> unique_names;
    for (std::string const& wall_name : wall_names) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            unique_names.insert(wall_name).second,
            "insert.analytic_walls must not contain duplicate wall names.");

        amrex::ParmParse const pp_wall("analytic_wall." + wall_name);
        std::string signed_value;
        std::string normal_x;
        std::string normal_y;
        std::string normal_z;
        utils::parser::Store_parserString(pp_wall, "signed_value", signed_value);
        utils::parser::Store_parserString(pp_wall, "normal_x", normal_x);
        utils::parser::Store_parserString(pp_wall, "normal_y", normal_y);
        utils::parser::Store_parserString(pp_wall, "normal_z", normal_z);
        result.walls.push_back({
            wall_name,
            std::make_unique<AnalyticBoundaryGeometry>(
                signed_value, normal_x, normal_y, normal_z)});
    }

    for (std::string const& species_name : mpc.GetSpeciesNames()) {
        std::vector<bool> absorbs(result.walls.size(), false);
        bool has_policy = false;
        amrex::ParmParse const pp_species(species_name);
        for (int wall_id = 0; wall_id < static_cast<int>(result.walls.size()); ++wall_id) {
            std::vector<std::string> behaviors;
            std::string const parameter =
                "analytic_wall." + result.walls[wall_id].name + ".behaviors";
            if (!pp_species.queryarr(parameter.c_str(), behaviors)) {
                continue;
            }
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                behaviors.size() == 1 && behaviors.front() == "absorb",
                species_name + "." + parameter +
                    " must be exactly 'absorb' in analytic-wall stage 2.");
            absorbs[wall_id] = true;
            has_policy = true;
        }
        if (has_policy) {
            result.absorbs.emplace(species_name, std::move(absorbs));
        }
    }
    return result;
}

template <typename Boundary>
void
ApplyAbsorbingWall (
    WarpXParticleContainer& pc, Boundary const& boundary,
    amrex::ParticleReal const species_charge, amrex::Real const dt,
    WallChargeGrid const& grid, amrex::MultiFab& wall_charge)
{
    for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {
            auto GetPosition = GetParticlePosition<PIdx>(pti);
            auto const charge_grid = wall_charge.array(pti);
            auto& soa = pti.GetParticleTile().GetStructOfArrays();
            uint64_t* const AMREX_RESTRICT idcpu = soa.GetIdCPUData().data();
            auto* const AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();
            auto* const AMREX_RESTRICT uy = soa.GetRealData(PIdx::uy).data();
            auto* const AMREX_RESTRICT uz = soa.GetRealData(PIdx::uz).data();
            auto* const AMREX_RESTRICT weight = soa.GetRealData(PIdx::w).data();

            // Deposition scatters several particles onto shared nodes. amrex::For
            // is required here; HostDevice atomics alone do not make ParallelFor
            // safe on CPU OpenMP execution.
            amrex::For(pti.numParticles(), [=] AMREX_GPU_DEVICE (long const i) noexcept {
                auto pidw = amrex::ParticleIDWrapper{idcpu[i]};
                if (!pidw.is_valid()) {
                    return;
                }

                amrex::ParticleReal x;
                amrex::ParticleReal y;
                amrex::ParticleReal z;
                GetPosition.AsStored(i, x, y, z);
                AnalyticBoundaryPosition const x_end{x, y, z};
                if (!BoundaryMath::IsOutsideDomain(boundary, x_end)) {
                    return;
                }

                AnalyticBoundaryPosition x_hit;
                amrex::Real dt_fraction_hit;
                BoundaryMath::BisectBoundaryIntersection(
                    boundary, x_end, ux[i], uy[i], uz[i], dt, x_hit,
                    dt_fraction_hit);
                DepositWallChargeToNodes(
                    charge_grid, grid, x_hit.x, x_hit.y, x_hit.z,
                    species_charge * weight[i]);
                pidw.make_invalid();
            });
        }
    }
}

#endif

} // namespace

void
AnalyticBoundaryInteraction ()
{
#if defined(WARPX_DIM_3D)
    WarpX& warpx_instance = WarpX::GetInstance();
    auto& mpc = warpx_instance.GetPartContainer();
    static AnalyticWallConfiguration const config = ReadAnalyticWallConfiguration(mpc);
    if (config.walls.empty()) {
        return;
    }

    HallWallCharge& persistent_wall_charge = HallWallCharge::GetInstance();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        persistent_wall_charge.isDefined(),
        "Analytic absorbing walls require the material Poisson path to initialize wall_charge.");
    WallChargeGrid const grid = MakeWallChargeGrid(warpx_instance.Geom(0));

    // Process walls in declaration order. Inputs must keep invalid regions
    // disjoint; this order is also the deterministic fallback at a shared edge.
    for (int wall_id = 0; wall_id < static_cast<int>(config.walls.size()); ++wall_id) {
        for (auto const& [species_name, policies] : config.absorbs) {
            if (!policies[wall_id]) {
                continue;
            }
            auto& pc = mpc.GetParticleContainerFromName(species_name);
            if (pc.getDoNotPush()) {
                continue;
            }
            amrex::Real const dt_effective =
                warpx_instance.getdt(0) * ParticleSubcyclingNdt(species_name);
            ApplyAbsorbingWall(
                pc, config.walls[wall_id].geometry->GetDeviceView(), pc.getCharge(),
                dt_effective, grid, persistent_wall_charge.get(0));
        }
    }
#endif
}

} // namespace Insert
