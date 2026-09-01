#include "Insert/Boundary/AnalyticBoundaryInteraction.h"

#include "Insert/Boundary/AnalyticBoundaryGeometry.h"
#include "Insert/Boundary/BoundaryMathUtils.h"
#include "Insert/Core/WarpXInsert.h"
#include "Particles/MultiParticleContainer.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <utility>

namespace Insert {
namespace {

/** \brief Per-species analytic-boundary configuration, parsed once from
 *         <species>.analytic_boundary.* and cached for the whole run. */
struct SpeciesBoundaryConfig {
    amrex::ParticleReal absorb_fraction = 0.0;
    amrex::ParticleReal specular_fraction = 0.0;
    amrex::ParticleReal wall_temperature = 0.0;
    std::unique_ptr<AnalyticBoundaryGeometry> geometry;
};

/** \brief Apply the analytic-boundary interaction to all live particles of
 *         one species, in place.
 *
 *  For each valid particle on the solid side of the boundary, the wall
 *  behavior is drawn once from the configured probabilities:
 *    - absorb:   the particle is invalidated (removed by Redistribute);
 *    - specular: the proper velocity is mirrored about the boundary normal;
 *    - diffuse:  the velocity is re-drawn from the wall Maxwellian.
 *  Reflected particles are placed at the trajectory/boundary intersection
 *  and advanced with their new velocity for the remaining fraction of the
 *  timestep.
 *
 * \param pc                particle container of the interacting species
 * \param boundary          geometry operator (device view)
 * \param absorb_fraction   absorption probability in [0, 1]
 * \param specular_fraction specular reflection probability in [0, 1];
 *                          the remaining probability is diffuse re-emission
 * \param wall_vth          thermal velocity spread for diffuse re-emission
 * \param mass              species mass (for the u -> v conversion)
 * \param dt                timestep
 */
template <typename Boundary>
void
ApplyAnalyticBoundaryInteraction (
    WarpXParticleContainer& pc, Boundary const& boundary,
    amrex::ParticleReal const absorb_fraction,
    amrex::ParticleReal const specular_fraction,
    amrex::ParticleReal const wall_vth,
    amrex::ParticleReal const mass, amrex::Real const dt)
{
    for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {
            auto GetPosition = GetParticlePosition<PIdx>(pti);
            auto SetPosition = SetParticlePosition<PIdx>(pti);

            auto& ptile = pti.GetParticleTile();
            auto& soa = ptile.GetStructOfArrays();
            uint64_t* const AMREX_RESTRICT idcpu = soa.GetIdCPUData().data();
            auto* const AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();
            auto* const AMREX_RESTRICT uy = soa.GetRealData(PIdx::uy).data();
            auto* const AMREX_RESTRICT uz = soa.GetRealData(PIdx::uz).data();

            amrex::ParallelForRNG(
                pti.numParticles(),
                [=] AMREX_GPU_DEVICE(
                    long i, amrex::RandomEngine const& engine) noexcept {
#if defined(WARPX_DIM_3D)
                    auto pidw = amrex::ParticleIDWrapper{idcpu[i]};
                    if (!pidw.is_valid()) { return; }

                    amrex::ParticleReal x, y, z;
                    GetPosition.AsStored(i, x, y, z);
                    amrex::XDim3 const x_end{x, y, z};
                    amrex::XDim3 const u_in{ux[i], uy[i], uz[i]};

                    // Particles still on the domain side are untouched.
                    if (!BoundaryMath::IsOutsideDomain(boundary, x_end)) {
                        return;
                    }

                    // Policy: draw the wall behavior once per hit. Absorbed
                    // particles are invalidated in place and removed by the
                    // next Redistribute.
                    amrex::ParticleReal const selector =
                        static_cast<amrex::ParticleReal>(
                            amrex::Random(engine));
                    if (selector < absorb_fraction) {
                        pidw.make_invalid();
                        return;
                    }

                    // Locate the contact point on the boundary by bisecting
                    // along the trajectory, and evaluate the domain-pointing
                    // normal there.
                    amrex::XDim3 x_hit;
                    amrex::Real dt_fraction_hit;
                    BoundaryMath::BisectBoundaryIntersection(
                        boundary, x_end, u_in, mass, dt, x_hit,
                        dt_fraction_hit);
                    amrex::XDim3 const normal = boundary.Normal(x_hit);

                    amrex::XDim3 u_out;
                    if (selector < absorb_fraction + specular_fraction) {
                        u_out = BoundaryMath::ReflectVelocity(
                            normal, u_in, engine);
                    } else {
                        u_out = BoundaryMath::DiffuseVelocity(
                            normal, wall_vth, engine);
                    }

                    // Advance from the contact point with the new velocity
                    // for the remaining fraction of the timestep.
                    BoundaryMath::AdvancePosition(
                        x_hit, u_out, mass, dt_fraction_hit * dt);
                    SetPosition.AsStored(i, x_hit.x, x_hit.y, x_hit.z);
                    ux[i] = u_out.x;
                    uy[i] = u_out.y;
                    uz[i] = u_out.z;
#else
                    amrex::ignore_unused(
                        engine, boundary, mass, dt, absorb_fraction,
                        specular_fraction, wall_vth);
#endif
                });
        }
    }
}

} // namespace

void
AnalyticBoundaryInteraction ()
{
#if defined(WARPX_DIM_3D)
    WarpX& warpx_instance = WarpX::GetInstance();
    auto& mpc = warpx_instance.GetPartContainer();

    // Per-species configuration, built once on the first call: each species
    // may define its own analytic boundary and wall behavior under
    // <species>.analytic_boundary.*.
    static std::map<std::string, SpeciesBoundaryConfig> const configs =
        [&mpc] {
            std::map<std::string, SpeciesBoundaryConfig> result;
            for (auto const& name : mpc.GetSpeciesNames()) {
                amrex::ParmParse const pp(name);
                int enabled = 0;
                pp.query("analytic_boundary.enabled", enabled);
                if (!enabled) {
                    continue;
                }

                SpeciesBoundaryConfig cfg;
                pp.query("analytic_boundary.absorb_fraction",
                         cfg.absorb_fraction);
                pp.query("analytic_boundary.specular_fraction",
                         cfg.specular_fraction);
                pp.query("analytic_boundary.wall_temperature",
                         cfg.wall_temperature);
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    cfg.absorb_fraction >= amrex::ParticleReal(0.0) &&
                        cfg.absorb_fraction <= amrex::ParticleReal(1.0),
                    name + ".analytic_boundary.absorb_fraction must be in "
                    "[0, 1].");
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    cfg.specular_fraction >= amrex::ParticleReal(0.0) &&
                        cfg.absorb_fraction + cfg.specular_fraction <=
                            amrex::ParticleReal(1.0),
                    name + ".analytic_boundary.specular_fraction must be in "
                    "[0, 1] and absorb_fraction + specular_fraction must "
                    "not exceed 1; the remaining probability is diffuse "
                    "re-emission.");
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    cfg.absorb_fraction + cfg.specular_fraction >=
                            amrex::ParticleReal(1.0) ||
                        cfg.wall_temperature > amrex::ParticleReal(0.0),
                    name + ".analytic_boundary.wall_temperature must be "
                    "positive for a non-zero diffuse reflection fraction.");

                std::string signed_value;
                std::string normal_x;
                std::string normal_y;
                std::string normal_z;
                utils::parser::Store_parserString(
                    pp, "analytic_boundary.signed_value", signed_value);
                utils::parser::Store_parserString(
                    pp, "analytic_boundary.normal_x", normal_x);
                utils::parser::Store_parserString(
                    pp, "analytic_boundary.normal_y", normal_y);
                utils::parser::Store_parserString(
                    pp, "analytic_boundary.normal_z", normal_z);
                cfg.geometry = std::make_unique<AnalyticBoundaryGeometry>(
                    signed_value, normal_x, normal_y, normal_z);

                result.emplace(name, std::move(cfg));
            }
            return result;
        }();

    if (configs.empty()) {
        return;
    }

    for (auto const& [name, cfg] : configs) {
        auto& pc = mpc.GetParticleContainerFromName(name);

        // Species under particle subcycling are only pushed every ndt
        // steps. Mirror the official boundary treatment
        // (MultiParticleContainer::ApplyBoundaryConditions): skip the sweep
        // on steps where this species was not pushed, and use its effective
        // timestep ndt * dt for the backtrace and post-reflection advance.
        if (pc.getDoNotPush()) {
            continue;
        }
        amrex::Real const dt_effective =
            warpx_instance.getdt(0) * ParticleSubcyclingNdt(name);

        // The temperature to thermal-velocity conversion needs the species
        // mass, so it happens here in the driver layer.
        amrex::ParticleReal diffuse_vth = 0.0;
        if (cfg.absorb_fraction + cfg.specular_fraction <
            amrex::ParticleReal(1.0))
        {
            diffuse_vth = static_cast<amrex::ParticleReal>(
                std::sqrt(PhysConst::kb * cfg.wall_temperature /
                          pc.getMass()));
        }

        ApplyAnalyticBoundaryInteraction(
            pc, cfg.geometry->GetDeviceView(), cfg.absorb_fraction,
            cfg.specular_fraction, diffuse_vth,
            static_cast<amrex::ParticleReal>(pc.getMass()), dt_effective);
    }
#endif
}

} // namespace Insert
