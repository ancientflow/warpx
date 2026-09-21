#include "Insert/Boundary/AnalyticBoundaryInteraction.h"

#include "Insert/Boundary/AnalyticBoundaryGeometry.h"
#include "Insert/Boundary/BoundaryMathUtils.h"
#include "Insert/Core/WarpXInsert.h"
#include "Insert/Fields/HallWallCharge.H"
#include "Insert/Math/ThermalVelocity.h"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleCreation/DefaultInitialization.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartUtils.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX_BLassert.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>
#include <AMReX_Scan.H>

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Insert {
namespace {

#if defined(WARPX_DIM_3D)

using namespace amrex::literals;

struct AnalyticWall
{
    std::string name;
    std::unique_ptr<AnalyticBoundaryGeometry> geometry;
};

constexpr int max_wall_behaviors = 6;
enum class WallBehavior : int {
    none = -1, absorb, specular, diffuse, neutralize, secondary_electron_1,
    secondary_electron_2
};

[[nodiscard]] WallBehavior
ParseBehavior (std::string const& name)
{
    if (name == "absorb") { return WallBehavior::absorb; }
    if (name == "specular") { return WallBehavior::specular; }
    if (name == "diffuse") { return WallBehavior::diffuse; }
    if (name == "neutralize") { return WallBehavior::neutralize; }
    if (name == "secondary_electron_1") { return WallBehavior::secondary_electron_1; }
    if (name == "secondary_electron_2") { return WallBehavior::secondary_electron_2; }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(false, "Unknown analytic-wall behavior: " + name);
    return WallBehavior::none;
}

struct WallPolicy
{
    bool active = false;
    int count = 0;
    std::array<WallBehavior, max_wall_behaviors> behaviors{};
    std::array<amrex::ParserExecutor<7>, max_wall_behaviors - 1> probabilities{};
    std::array<std::unique_ptr<amrex::Parser>, max_wall_behaviors - 1> parsers{};
    amrex::ParticleReal wall_temperature = 0.0_prt;
    amrex::ParticleReal secondary_temperature_eV = 0.0_prt;
    std::string neutral_species;
    std::string secondary_species;
};

/** Trivially copyable subset captured by particle kernels.  The parsers in
 * WallPolicy own the expression storage for these executors. */
struct WallPolicyView
{
    int count = 0;
    std::array<WallBehavior, max_wall_behaviors> behaviors{};
    std::array<amrex::ParserExecutor<7>, max_wall_behaviors - 1> probabilities{};
    amrex::ParticleReal wall_temperature = 0.0_prt;

    explicit WallPolicyView (WallPolicy const& policy)
        : count(policy.count), behaviors(policy.behaviors), probabilities(policy.probabilities),
          wall_temperature(policy.wall_temperature) {}
};

struct AnalyticWallConfiguration
{
    std::vector<AnalyticWall> walls;
    std::map<std::string, std::vector<WallPolicy>> policies;
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
        std::vector<WallPolicy> policies(result.walls.size());
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
                !behaviors.empty() && behaviors.size() <= max_wall_behaviors,
                species_name + "." + parameter + " must contain one to six behaviors.");
            WallPolicy& policy = policies[wall_id];
            policy.active = true;
            policy.count = static_cast<int>(behaviors.size());
            std::string const prefix = "analytic_wall." + result.walls[wall_id].name + ".";
            for (int i = 0; i < policy.count; ++i) {
                policy.behaviors[i] = ParseBehavior(behaviors[i]);
                if (i + 1 < policy.count) {
                    std::string expression;
                    utils::parser::Store_parserString(
                        pp_species, prefix + "p_" + behaviors[i], expression);
                    policy.parsers[i] = std::make_unique<amrex::Parser>(
                        utils::parser::makeParser(expression,
                            {"E_eV", "u_n", "u_t", "x", "y", "z", "t"}));
                    policy.probabilities[i] = policy.parsers[i]->compile<7>();
                }
            }
            bool need_temperature = false;
            bool need_neutral = false;
            bool need_secondary = false;
            for (int i = 0; i < policy.count; ++i) {
                need_temperature = need_temperature ||
                    policy.behaviors[i] == WallBehavior::diffuse ||
                    policy.behaviors[i] == WallBehavior::neutralize;
                need_neutral = need_neutral || policy.behaviors[i] == WallBehavior::neutralize;
                need_secondary = need_secondary ||
                    policy.behaviors[i] == WallBehavior::secondary_electron_1 ||
                    policy.behaviors[i] == WallBehavior::secondary_electron_2;
            }
            if (need_temperature) {
                pp_species.get((prefix + "wall_temperature").c_str(), policy.wall_temperature);
            }
            if (need_neutral) {
                pp_species.get((prefix + "neutral_species").c_str(), policy.neutral_species);
                amrex::ignore_unused(mpc.GetParticleContainerFromName(policy.neutral_species));
            }
            if (need_secondary) {
                pp_species.get(
                    (prefix + "secondary_electron_species").c_str(), policy.secondary_species);
                pp_species.get((prefix + "secondary_electron_temperature_eV").c_str(),
                    policy.secondary_temperature_eV);
                amrex::ignore_unused(mpc.GetParticleContainerFromName(policy.secondary_species));
            }
            has_policy = true;
        }
        if (has_policy) {
            result.policies.emplace(species_name, std::move(policies));
        }
    }
    return result;
}

template <typename Boundary>
void
CreateWallProducts (WarpXParticleContainer& source, WarpXParticleContainer& destination,
                    WarpXParIter const& pti, amrex::Gpu::DeviceVector<int> const& choices,
                    WallBehavior const requested, int const multiplicity, Boundary const& boundary,
                    amrex::Real const dt, amrex::ParticleReal const thermal_velocity, int const lev)
{
    int const np = pti.numParticles();
    auto* const choice = choices.dataPtr();
    amrex::Gpu::DeviceVector<int> counts(np), offsets(np);
    auto* const count = counts.dataPtr();
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept {
        count[i] = choice[i] == static_cast<int>(requested) ? multiplicity : 0;
    });
    amrex::Long const added = amrex::Scan::ExclusiveSum(np, count, offsets.data());
    if (added == 0) { return; }

    SmartCopyFactory const copy_factory(source, destination);
    auto const Copy = copy_factory.getSmartCopy();
    auto& dst_tile = destination.DefineAndReturnParticleTile(
        lev, pti.index(), pti.LocalTileIndex());
    amrex::Long const old_size = dst_tile.numParticles();
    dst_tile.resize(old_size + added);
    // Reacquire source data after resize: source and destination may be identical.
    auto const src_data = source.ParticlesAt(lev, pti).getParticleTileData();
    auto const dst_data = dst_tile.getParticleTileData();
    auto* const offset = offsets.dataPtr();
    amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (
        int i, amrex::RandomEngine const& engine) noexcept {
        if (count[i] == 0) { return; }
        AnalyticBoundaryPosition const end{
            src_data.m_rdata[PIdx::x][i], src_data.m_rdata[PIdx::y][i],
            src_data.m_rdata[PIdx::z][i]};
        AnalyticBoundaryPosition hit;
        amrex::Real fraction;
        BoundaryMath::BisectBoundaryIntersection(
            boundary, end, src_data.m_rdata[PIdx::ux][i], src_data.m_rdata[PIdx::uy][i],
            src_data.m_rdata[PIdx::uz][i], dt, hit, fraction);
        amrex::XDim3 const u = BoundaryMath::DiffuseVelocity(
            amrex::XDim3{boundary.Normal(hit).x, boundary.Normal(hit).y, boundary.Normal(hit).z},
            thermal_velocity, engine);
        for (int j = 0; j < count[i]; ++j) {
            int const idst = old_size + offset[i] + j;
            Copy(dst_data, src_data, i, idst, engine);
            dst_data.m_rdata[PIdx::ux][idst] = static_cast<amrex::ParticleReal>(u.x);
            dst_data.m_rdata[PIdx::uy][idst] = static_cast<amrex::ParticleReal>(u.y);
            dst_data.m_rdata[PIdx::uz][idst] = static_cast<amrex::ParticleReal>(u.z);
            dst_data.m_rdata[PIdx::x][idst] = hit.x +
                static_cast<amrex::ParticleReal>(u.x) * fraction * dt;
            dst_data.m_rdata[PIdx::y][idst] = hit.y +
                static_cast<amrex::ParticleReal>(u.y) * fraction * dt;
            dst_data.m_rdata[PIdx::z][idst] = hit.z +
                static_cast<amrex::ParticleReal>(u.z) * fraction * dt;
        }
    });
    ParticleCreation::DefaultInitializeRuntimeAttributes(
        dst_tile, destination, old_size, old_size + added);
    amrex::Gpu::synchronize();
    setNewParticleIDs(dst_tile, old_size, added);
}

template <typename Boundary>
void
ApplyWallPolicy (WarpXParticleContainer& pc, WallPolicy const& policy, Boundary const& boundary,
                 amrex::Real const time, amrex::Real const dt, WallChargeGrid const& grid,
                 amrex::MultiFab& wall_charge, MultiParticleContainer& mpc)
{
    amrex::ParticleReal const mass = pc.getMass();
    amrex::ParticleReal const charge = pc.getCharge();
    WallPolicyView const device_policy(policy);
    amrex::ParticleReal neutral_charge = 0.0_prt;
    amrex::ParticleReal secondary_charge = 0.0_prt;
    if (!policy.neutral_species.empty()) {
        neutral_charge = mpc.GetParticleContainerFromName(policy.neutral_species).getCharge();
    }
    if (!policy.secondary_species.empty()) {
        secondary_charge = mpc.GetParticleContainerFromName(policy.secondary_species).getCharge();
    }
    for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {
            auto& soa = pti.GetParticleTile().GetStructOfArrays();
            uint64_t* const AMREX_RESTRICT idcpu = soa.GetIdCPUData().data();
            auto* const AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();
            auto* const AMREX_RESTRICT uy = soa.GetRealData(PIdx::uy).data();
            auto* const AMREX_RESTRICT uz = soa.GetRealData(PIdx::uz).data();
            auto GetPosition = GetParticlePosition<PIdx>(pti);
            int const np = pti.numParticles();
            amrex::Gpu::DeviceVector<int> choices(np);
            auto* const choice = choices.dataPtr();
            amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (
                int i, amrex::RandomEngine const& engine) noexcept {
                auto pidw = amrex::ParticleIDWrapper{idcpu[i]};
                if (!pidw.is_valid()) {
                    choice[i] = static_cast<int>(WallBehavior::none);
                    return;
                }
                amrex::ParticleReal x;
                amrex::ParticleReal y;
                amrex::ParticleReal z;
                GetPosition.AsStored(i, x, y, z);
                AnalyticBoundaryPosition const x_end{x, y, z};
                if (!BoundaryMath::IsOutsideDomain(boundary, x_end)) {
                    choice[i] = static_cast<int>(WallBehavior::none);
                    return;
                }
                AnalyticBoundaryPosition x_hit;
                amrex::Real dt_fraction_hit;
                BoundaryMath::BisectBoundaryIntersection(
                    boundary, x_end, ux[i], uy[i], uz[i], dt, x_hit, dt_fraction_hit);
                auto const normal = boundary.Normal(x_hit);
                amrex::ParticleReal const n2 =
                    normal.x*normal.x + normal.y*normal.y + normal.z*normal.z;
                amrex::ParticleReal const inv_n = n2 > 0.0_prt ? 1.0_prt/std::sqrt(n2) : 0.0_prt;
                amrex::ParticleReal const un =
                    (ux[i]*normal.x + uy[i]*normal.y + uz[i]*normal.z)*inv_n;
                amrex::ParticleReal const u2 = ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i];
                amrex::ParticleReal const ut = std::sqrt(amrex::max(u2-un*un, 0.0_prt));
                amrex::ParticleReal const e_ev = mass*u2/(2.0_prt*PhysConst::q_e);
                amrex::ParticleReal p = 0.0_prt;
                amrex::ParticleReal const random = amrex::Random(engine);
                for (int j = 0; j + 1 < device_policy.count; ++j) {
                    p += static_cast<amrex::ParticleReal>(device_policy.probabilities[j](
                        e_ev, un, ut, x_hit.x, x_hit.y, x_hit.z, time));
                    if (random < p) {
                        choice[i] = static_cast<int>(device_policy.behaviors[j]);
                        return;
                    }
                }
                choice[i] = static_cast<int>(device_policy.behaviors[device_policy.count-1]);
            });

            if (!policy.neutral_species.empty()) {
                auto& target = mpc.GetParticleContainerFromName(policy.neutral_species);
                CreateWallProducts(
                    pc, target, pti, choices, WallBehavior::neutralize, 1, boundary, dt,
                    Math::ThermalVelocityFromTemperature(
                        policy.wall_temperature, target.getMass()), lev);
            }
            if (!policy.secondary_species.empty()) {
                auto& target = mpc.GetParticleContainerFromName(policy.secondary_species);
                auto const thermal = Math::ThermalVelocityFromEV(
                    policy.secondary_temperature_eV, target.getMass());
                CreateWallProducts(pc, target, pti, choices, WallBehavior::secondary_electron_1, 1,
                    boundary, dt, thermal, lev);
                CreateWallProducts(pc, target, pti, choices, WallBehavior::secondary_electron_2, 2,
                    boundary, dt, thermal, lev);
            }

            // Deposition scatters several particles onto shared nodes. amrex::For
            // is required here; HostDevice atomics alone do not make ParallelFor safe on CPU.
            auto const charge_grid = wall_charge.array(pti);
            // Product creation can resize this source tile when a product is emitted
            // into the same species. Reacquire all source pointers after that resize.
            auto& updated_soa = pc.ParticlesAt(lev, pti).GetStructOfArrays();
            uint64_t* const AMREX_RESTRICT updated_idcpu = updated_soa.GetIdCPUData().data();
            auto* const AMREX_RESTRICT updated_x = updated_soa.GetRealData(PIdx::x).data();
            auto* const AMREX_RESTRICT updated_y = updated_soa.GetRealData(PIdx::y).data();
            auto* const AMREX_RESTRICT updated_z = updated_soa.GetRealData(PIdx::z).data();
            auto* const AMREX_RESTRICT updated_ux = updated_soa.GetRealData(PIdx::ux).data();
            auto* const AMREX_RESTRICT updated_uy = updated_soa.GetRealData(PIdx::uy).data();
            auto* const AMREX_RESTRICT updated_uz = updated_soa.GetRealData(PIdx::uz).data();
            auto* const AMREX_RESTRICT updated_weight = updated_soa.GetRealData(PIdx::w).data();
            // Charge deposition is a scatter operation.  Keep it in amrex::For:
            // ParallelFor promises independent iterations even when atomics are used.
            amrex::For(np, [=] AMREX_GPU_DEVICE (long const i) noexcept {
                WallBehavior const selected = static_cast<WallBehavior>(choice[i]);
                if (selected == WallBehavior::none || selected == WallBehavior::specular ||
                    selected == WallBehavior::diffuse) { return; }
                AnalyticBoundaryPosition hit;
                amrex::Real fraction;
                BoundaryMath::BisectBoundaryIntersection(
                    boundary, {updated_x[i], updated_y[i], updated_z[i]}, updated_ux[i],
                    updated_uy[i], updated_uz[i], dt, hit, fraction);
                amrex::ParticleReal out_charge = 0.0_prt;
                if (selected == WallBehavior::neutralize) { out_charge = neutral_charge; }
                else if (selected == WallBehavior::secondary_electron_1) {
                    out_charge = secondary_charge;
                } else if (selected == WallBehavior::secondary_electron_2) {
                    out_charge = 2.0_prt * secondary_charge;
                }
                DepositWallChargeToNodes(
                    charge_grid, grid, hit.x, hit.y, hit.z,
                    (charge - out_charge) * updated_weight[i]);
            });
            amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (
                long const i, amrex::RandomEngine const& engine) noexcept {
                WallBehavior const selected = static_cast<WallBehavior>(choice[i]);
                auto pidw = amrex::ParticleIDWrapper{updated_idcpu[i]};
                if (!pidw.is_valid() || selected == WallBehavior::none) { return; }
                amrex::ParticleReal x = updated_x[i];
                amrex::ParticleReal y = updated_y[i];
                amrex::ParticleReal z = updated_z[i];
                AnalyticBoundaryPosition hit;
                amrex::Real fraction;
                BoundaryMath::BisectBoundaryIntersection(
                    boundary, {x,y,z}, updated_ux[i], updated_uy[i], updated_uz[i], dt, hit,
                    fraction);
                if (selected == WallBehavior::specular || selected == WallBehavior::diffuse) {
                    amrex::XDim3 const normal{
                        boundary.Normal(hit).x, boundary.Normal(hit).y, boundary.Normal(hit).z};
                    amrex::XDim3 u{updated_ux[i], updated_uy[i], updated_uz[i]};
                    if (selected == WallBehavior::specular) {
                        u = BoundaryMath::ReflectVelocity(normal, u, engine);
                    }
                    else { u = BoundaryMath::DiffuseVelocity(normal,
                        Math::ThermalVelocityFromTemperature(
                            device_policy.wall_temperature, mass), engine);
                    }
                    updated_ux[i] = static_cast<amrex::ParticleReal>(u.x);
                    updated_uy[i] = static_cast<amrex::ParticleReal>(u.y);
                    updated_uz[i] = static_cast<amrex::ParticleReal>(u.z);
                    updated_x[i] = hit.x + updated_ux[i] * fraction * dt;
                    updated_y[i] = hit.y + updated_uy[i] * fraction * dt;
                    updated_z[i] = hit.z + updated_uz[i] * fraction * dt;
                    return;
                }
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
        for (auto const& [species_name, policies] : config.policies) {
            WallPolicy const& policy = policies[wall_id];
            if (!policy.active) {
                continue;
            }
            auto& pc = mpc.GetParticleContainerFromName(species_name);
            if (pc.getDoNotPush()) {
                continue;
            }
            amrex::Real const dt_effective =
                warpx_instance.getdt(0) * ParticleSubcyclingNdt(species_name);
            ApplyWallPolicy(pc, policy, config.walls[wall_id].geometry->GetDeviceView(),
                warpx_instance.gett_new(0), dt_effective, grid,
                persistent_wall_charge.get(0), mpc);
        }
    }
#endif
}

} // namespace Insert
