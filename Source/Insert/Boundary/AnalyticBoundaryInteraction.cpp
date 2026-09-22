#include "Insert/Boundary/AnalyticBoundaryInteraction.h"

#include "Fields.H"
#include "Insert/Boundary/AnalyticBoundaryGeometry.h"
#include "Insert/Boundary/BoundaryMathUtils.h"
#include "Insert/Core/WarpXInsert.h"
#include "Insert/Fields/HallWallCharge.H"
#include "Insert/Math/ThermalVelocity.h"
#include "Particles/Gather/FieldGather.H"
#include "Particles/Gather/GetExternalFields.H"
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
#include <AMReX_FArrayBox.H>
#include <AMReX_IndexType.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <AMReX_Reduce.H>
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
    bool deposit_wall_charge = true;
    bool has_neutral_product_charge = false;
    amrex::ParticleReal neutral_product_charge = 0.0_prt;
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
            pp_species.query(
                (prefix + "deposit_wall_charge").c_str(), policy.deposit_wall_charge);
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
                policy.has_neutral_product_charge =
                    utils::parser::queryWithParser(pp_species,
                        (prefix + "neutral_product_charge").c_str(),
                        policy.neutral_product_charge) != 0;
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
                    amrex::Gpu::DeviceVector<amrex::ParticleReal> const& hit_x,
                    amrex::Gpu::DeviceVector<amrex::ParticleReal> const& hit_y,
                    amrex::Gpu::DeviceVector<amrex::ParticleReal> const& hit_z,
                    amrex::Gpu::DeviceVector<amrex::Real> const& hit_fraction,
                    WallBehavior const requested, int const multiplicity, Boundary const& boundary,
                    amrex::Real const dt, amrex::ParticleReal const thermal_velocity, int const lev)
{
    int const np = pti.numParticles();
    auto* const choice = choices.dataPtr();
    auto* const hx = hit_x.dataPtr();
    auto* const hy = hit_y.dataPtr();
    auto* const hz = hit_z.dataPtr();
    auto* const hfrac = hit_fraction.dataPtr();
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
        // Hit point and remaining time fraction were computed once by the
        // behavior-selection kernel; reuse them here.
        AnalyticBoundaryPosition const hit{hx[i], hy[i], hz[i]};
        amrex::Real const fraction = hfrac[i];
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
                 amrex::MultiFab& wall_charge, MultiParticleContainer& mpc,
                 std::array<amrex::Long, max_wall_behaviors>& counts)
{
    amrex::ParticleReal const mass = pc.getMass();
    amrex::ParticleReal const charge = pc.getCharge();
    WallPolicyView const device_policy(policy);
    amrex::ParticleReal neutral_charge = 0.0_prt;
    amrex::ParticleReal secondary_charge = 0.0_prt;
    if (!policy.neutral_species.empty()) {
        // The product charge defaults to the neutral species' charge, but can
        // be overridden: a species may carry a nonzero bookkeeping charge
        // (e.g. 1 C so that charge deposition yields a neutral density field)
        // while the physical neutralization product is uncharged.
        neutral_charge = policy.has_neutral_product_charge
            ? policy.neutral_product_charge
            : mpc.GetParticleContainerFromName(policy.neutral_species).getCharge();
    }
    if (!policy.secondary_species.empty()) {
        secondary_charge = mpc.GetParticleContainerFromName(policy.secondary_species).getCharge();
    }
    using ablastr::fields::Direction;
    using warpx::fields::FieldType;
    WarpX& warpx = WarpX::GetInstance();
    amrex::ParticleReal const q_over_m = mass != 0.0_prt ? charge/mass : 0.0_prt;
    for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
        // Incident-velocity correction: re-gather, at the back-traced
        // step-start position, the same fields that drove this step's
        // momentum push.
        bool const correct_impact_velocity =
            mass > 0.0_prt &&
            warpx.m_fields.has(FieldType::Efield_aux, Direction{0}, lev) &&
            warpx.m_fields.has(FieldType::Bfield_aux, Direction{0}, lev);
        amrex::XDim3 const dinv = WarpX::InvCellSize(std::max(lev, 0));
        int const nox = WarpX::nox;
        int const n_rz_azimuthal_modes = WarpX::n_rz_azimuthal_modes;
        bool const galerkin_interpolation = WarpX::galerkin_interpolation;
        amrex::GpuArray<amrex::Real, 3> const prob_lo = warpx.Geom(lev).ProbLoArray();
        amrex::GpuArray<amrex::Real, 3> const prob_hi = warpx.Geom(lev).ProbHiArray();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {
            amrex::Array4<amrex::Real const> ex_arr{};
            amrex::Array4<amrex::Real const> ey_arr{};
            amrex::Array4<amrex::Real const> ez_arr{};
            amrex::Array4<amrex::Real const> bx_arr{};
            amrex::Array4<amrex::Real const> by_arr{};
            amrex::Array4<amrex::Real const> bz_arr{};
            amrex::IndexType ex_type, ey_type, ez_type, bx_type, by_type, bz_type;
            amrex::XDim3 xyzmin{0.0, 0.0, 0.0};
            amrex::Dim3 lo{0, 0, 0};
            GetExternalEBField const getExternalEB(pti);
            if (correct_impact_velocity) {
                amrex::MultiFab const& Ex =
                    *warpx.m_fields.get(FieldType::Efield_aux, Direction{0}, lev);
                amrex::MultiFab const& Ey =
                    *warpx.m_fields.get(FieldType::Efield_aux, Direction{1}, lev);
                amrex::MultiFab const& Ez =
                    *warpx.m_fields.get(FieldType::Efield_aux, Direction{2}, lev);
                amrex::MultiFab const& Bx =
                    *warpx.m_fields.get(FieldType::Bfield_aux, Direction{0}, lev);
                amrex::MultiFab const& By =
                    *warpx.m_fields.get(FieldType::Bfield_aux, Direction{1}, lev);
                amrex::MultiFab const& Bz =
                    *warpx.m_fields.get(FieldType::Bfield_aux, Direction{2}, lev);
                amrex::Box const box = pti.tilebox();
                xyzmin = WarpX::LowerCorner(box, lev, 0.0_rt);
                lo = amrex::lbound(box);
                amrex::FArrayBox const& exfab = Ex[pti];
                amrex::FArrayBox const& eyfab = Ey[pti];
                amrex::FArrayBox const& ezfab = Ez[pti];
                amrex::FArrayBox const& bxfab = Bx[pti];
                amrex::FArrayBox const& byfab = By[pti];
                amrex::FArrayBox const& bzfab = Bz[pti];
                ex_arr = exfab.array();
                ey_arr = eyfab.array();
                ez_arr = ezfab.array();
                bx_arr = bxfab.array();
                by_arr = byfab.array();
                bz_arr = bzfab.array();
                ex_type = exfab.box().ixType();
                ey_type = eyfab.box().ixType();
                ez_type = ezfab.box().ixType();
                bx_type = bxfab.box().ixType();
                by_type = byfab.box().ixType();
                bz_type = bzfab.box().ixType();
            }
            auto& soa = pti.GetParticleTile().GetStructOfArrays();
            uint64_t* const AMREX_RESTRICT idcpu = soa.GetIdCPUData().data();
            auto* const AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();
            auto* const AMREX_RESTRICT uy = soa.GetRealData(PIdx::uy).data();
            auto* const AMREX_RESTRICT uz = soa.GetRealData(PIdx::uz).data();
            auto GetPosition = GetParticlePosition<PIdx>(pti);
            int const np = pti.numParticles();
            amrex::Gpu::DeviceVector<int> choices(np);
            auto* const choice = choices.dataPtr();
            // Hit point and remaining time fraction for each particle that
            // crosses the wall, computed once by the selection kernel and
            // reused by product creation, charge deposition and reflection.
            amrex::Gpu::DeviceVector<amrex::ParticleReal> hit_x(np), hit_y(np), hit_z(np);
            amrex::Gpu::DeviceVector<amrex::Real> hit_fraction(np);
            auto* const hx = hit_x.dataPtr();
            auto* const hy = hit_y.dataPtr();
            auto* const hz = hit_z.dataPtr();
            auto* const hfrac = hit_fraction.dataPtr();
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
                hx[i] = x_hit.x;
                hy[i] = x_hit.y;
                hz[i] = x_hit.z;
                hfrac[i] = dt_fraction_hit;
                // Incident velocity at the hit time. The stored velocity is
                // defined at mid-step (t^n + dt/2); this step's momentum push
                // was driven by the fields gathered at the step-start
                // position, so back-trace the straight-line trajectory to x^n
                // and re-gather the same fields there.
                amrex::ParticleReal ux_hit = ux[i];
                amrex::ParticleReal uy_hit = uy[i];
                amrex::ParticleReal uz_hit = uz[i];
                if (correct_impact_velocity) {
                    // Time from step start to the hit: dt_fraction_hit spans
                    // hit -> end, so start -> hit is (1 - dt_fraction_hit)*dt.
                    amrex::ParticleReal const dt_hit = static_cast<amrex::ParticleReal>(
                        (1.0 - dt_fraction_hit)*dt);
                    amrex::ParticleReal const xs = x_hit.x - ux[i]*dt_hit;
                    amrex::ParticleReal const ys = x_hit.y - uy[i]*dt_hit;
                    amrex::ParticleReal const zs = x_hit.z - uz[i]*dt_hit;
                    bool const start_inside =
                        xs >= static_cast<amrex::ParticleReal>(prob_lo[0]) &&
                        xs <= static_cast<amrex::ParticleReal>(prob_hi[0]) &&
                        ys >= static_cast<amrex::ParticleReal>(prob_lo[1]) &&
                        ys <= static_cast<amrex::ParticleReal>(prob_hi[1]) &&
                        zs >= static_cast<amrex::ParticleReal>(prob_lo[2]) &&
                        zs <= static_cast<amrex::ParticleReal>(prob_hi[2]);
                    if (start_inside) {
                        amrex::ParticleReal Exp = 0.0_prt;
                        amrex::ParticleReal Eyp = 0.0_prt;
                        amrex::ParticleReal Ezp = 0.0_prt;
                        amrex::ParticleReal Bxp = 0.0_prt;
                        amrex::ParticleReal Byp = 0.0_prt;
                        amrex::ParticleReal Bzp = 0.0_prt;
                        doGatherShapeN(xs, ys, zs, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                            ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                            ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                            dinv, xyzmin, lo, n_rz_azimuthal_modes, nox,
                            galerkin_interpolation);
                        getExternalEB(i, Exp, Eyp, Ezp, Bxp, Byp, Bzp);
                        // Non-relativistic Lorentz acceleration, matching the
                        // non-relativistic trajectory model of BoundaryMath.
                        amrex::ParticleReal const ax =
                            q_over_m*(Exp + uy[i]*Bzp - uz[i]*Byp);
                        amrex::ParticleReal const ay =
                            q_over_m*(Eyp + uz[i]*Bxp - ux[i]*Bzp);
                        amrex::ParticleReal const az =
                            q_over_m*(Ezp + ux[i]*Byp - uy[i]*Bxp);
                        // t_hit - t_midstep = ((1 - dt_fraction_hit) - 0.5) * dt
                        amrex::ParticleReal const dt_from_mid =
                            static_cast<amrex::ParticleReal>(
                                (0.5 - dt_fraction_hit)*dt);
                        ux_hit += ax*dt_from_mid;
                        uy_hit += ay*dt_from_mid;
                        uz_hit += az*dt_from_mid;
                    }
                }
                auto const normal = boundary.Normal(x_hit);
                amrex::ParticleReal const n2 =
                    normal.x*normal.x + normal.y*normal.y + normal.z*normal.z;
                amrex::ParticleReal const inv_n = n2 > 0.0_prt ? 1.0_prt/std::sqrt(n2) : 0.0_prt;
                amrex::ParticleReal const un =
                    (ux_hit*normal.x + uy_hit*normal.y + uz_hit*normal.z)*inv_n;
                amrex::ParticleReal const u2 = ux_hit*ux_hit + uy_hit*uy_hit + uz_hit*uz_hit;
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
                    pc, target, pti, choices, hit_x, hit_y, hit_z, hit_fraction,
                    WallBehavior::neutralize, 1, boundary, dt,
                    Math::ThermalVelocityFromTemperature(
                        policy.wall_temperature, target.getMass()), lev);
            }
            if (!policy.secondary_species.empty()) {
                auto& target = mpc.GetParticleContainerFromName(policy.secondary_species);
                auto const thermal = Math::ThermalVelocityFromEV(
                    policy.secondary_temperature_eV, target.getMass());
                CreateWallProducts(pc, target, pti, choices, hit_x, hit_y, hit_z, hit_fraction,
                    WallBehavior::secondary_electron_1, 1,
                    boundary, dt, thermal, lev);
                CreateWallProducts(pc, target, pti, choices, hit_x, hit_y, hit_z, hit_fraction,
                    WallBehavior::secondary_electron_2, 2,
                    boundary, dt, thermal, lev);
            }

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
            if (policy.deposit_wall_charge) {
                // Deposition scatters several particles onto shared nodes.
                // amrex::For is required here; HostDevice atomics alone do not
                // make ParallelFor safe on CPU.
                auto const charge_grid = wall_charge.array(pti);
                // Charge deposition is a scatter operation.  Keep it in
                // amrex::For: ParallelFor promises independent iterations even
                // when atomics are used.
                amrex::For(np, [=] AMREX_GPU_DEVICE (long const i) noexcept {
                    WallBehavior const selected = static_cast<WallBehavior>(choice[i]);
                    if (selected == WallBehavior::none || selected == WallBehavior::specular ||
                        selected == WallBehavior::diffuse) { return; }
                    // Reuse the hit point computed by the selection kernel.
                    AnalyticBoundaryPosition const hit{hx[i], hy[i], hz[i]};
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
            }
            amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (
                long const i, amrex::RandomEngine const& engine) noexcept {
                WallBehavior const selected = static_cast<WallBehavior>(choice[i]);
                auto pidw = amrex::ParticleIDWrapper{updated_idcpu[i]};
                if (!pidw.is_valid() || selected == WallBehavior::none) { return; }
                if (selected == WallBehavior::specular || selected == WallBehavior::diffuse) {
                    // Reuse the hit point computed by the selection kernel.
                    AnalyticBoundaryPosition const hit{hx[i], hy[i], hz[i]};
                    amrex::Real const fraction = hfrac[i];
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

            // Count the events selected for this tile. choice[i] still holds the
            // behavior even for particles that were invalidated above.
            amrex::Long tile_counts[max_wall_behaviors];
            for (int b = 0; b < max_wall_behaviors; ++b) {
                tile_counts[b] = amrex::Reduce::Sum<amrex::Long>(
                    np, [=] AMREX_GPU_DEVICE (long i) noexcept -> amrex::Long {
                        return choice[i] == b ? 1 : 0;
                    });
            }
#ifdef AMREX_USE_OMP
#pragma omp critical
#endif
            {
                for (int b = 0; b < max_wall_behaviors; ++b) {
                    counts[b] += tile_counts[b];
                }
            }
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

    std::array<amrex::Long, max_wall_behaviors> wall_event_counts{};

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
                persistent_wall_charge.get(0), mpc, wall_event_counts);
        }
    }

    amrex::Long counts[max_wall_behaviors];
    for (int b = 0; b < max_wall_behaviors; ++b) { counts[b] = wall_event_counts[b]; }
    amrex::ParallelDescriptor::ReduceLongSum(counts, max_wall_behaviors);

    amrex::Long const n_absorb = counts[static_cast<int>(WallBehavior::absorb)];
    amrex::Long const n_specular = counts[static_cast<int>(WallBehavior::specular)];
    amrex::Long const n_diffuse = counts[static_cast<int>(WallBehavior::diffuse)];
    amrex::Long const n_neutralize = counts[static_cast<int>(WallBehavior::neutralize)];
    amrex::Long const n_secondary_1 = counts[static_cast<int>(WallBehavior::secondary_electron_1)];
    amrex::Long const n_secondary_2 = counts[static_cast<int>(WallBehavior::secondary_electron_2)];
    amrex::Long const n_absorbed = n_absorb + n_neutralize + n_secondary_1 + n_secondary_2;
    amrex::Long const n_emitted = n_secondary_1 + 2*n_secondary_2;
    amrex::Print() << "AnalyticBoundaryInteraction: step " << warpx_instance.getistep(0)
        << '\n'
        << "  absorbed particles:  " << n_absorbed
        << " (absorb " << n_absorb << ", neutralize " << n_neutralize
        << ", secondary " << n_secondary_1 + n_secondary_2 << ")\n"
        << "  specular:            " << n_specular << '\n'
        << "  diffuse:             " << n_diffuse << '\n'
        << "  emitted secondaries: " << n_emitted
        << " (secondary_electron_1 " << n_secondary_1
        << ", secondary_electron_2 " << n_secondary_2 << ")\n";
#endif
}

} // namespace Insert
