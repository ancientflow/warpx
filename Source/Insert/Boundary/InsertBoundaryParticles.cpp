#include "InsertBoundaryParticles.h"

#include "WarpX.H"

#include "EmbeddedBoundary/Enabled.H"
#include "Initialization/SampleGaussianFluxDistribution.H"
#include "Insert/Boundary/NeutralAtomEBGeometry.h"
#include "Insert/Boundary/ZMinWallCharge.h"
#include "Insert/Config/WarpXSimulationConfig.h"
#include "Insert/Diagnostics/InsertRuntimeDiagnostics.h"
#include "Insert/Utils/InsertUtils.h"
#include "Particles/Algorithms/KineticEnergy.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Particles/ParticleCreation/FilterCopyTransform.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartUtils.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Array4.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <algorithm>
#include <cmath>
#include <string>

namespace Insert {
namespace {

using NeutralAtomEBGeometry::GetTruncatedConeNormal;
using NeutralAtomEBGeometry::TruncatedConeGeometry;

template <int N = 1, typename SrcPC, typename DstPC, typename FilterFunc,
          typename TransformFunc>
amrex::Long
FilterCopyTransformParticleTiles (SrcPC& src_pc, DstPC& dst_pc,
                                  FilterFunc const& filter,
                                  TransformFunc const& transform,
                                  int lev_min = 0, int lev_max = -1) {
#ifdef AMREX_USE_GPU
    auto const* src_arena = src_pc.arena();
    auto const* dst_arena = dst_pc.arena();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_arena->isDevice() || src_arena->isManaged() ||
            src_arena == amrex::The_Pinned_Arena(),
        "FilterCopyTransformParticleTiles source particles must be in a "
        "device, managed, or AMReX pinned arena for GPU execution.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        dst_arena->isDevice() || dst_arena->isManaged() ||
            dst_arena == amrex::The_Pinned_Arena(),
        "FilterCopyTransformParticleTiles destination particles must be in a "
        "device, managed, or AMReX pinned arena for GPU execution.");
#endif

    const SmartCopyFactory copy_factory(src_pc, dst_pc);
    const auto Copy = copy_factory.getSmartCopy();

    if (lev_max < 0) {
        lev_max = src_pc.numLevels() - 1;
    }
    lev_max = std::min(lev_max, src_pc.numLevels() - 1);

    amrex::Long total_added = 0;

    for (int lev = lev_min; lev <= lev_max; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())                      \
    reduction(+ : total_added)
#endif
        for (WarpXParIter pti(src_pc, lev); pti.isValid(); ++pti) {
            auto& src_tile = pti.GetParticleTile();
            if (src_tile.numParticles() == 0) {
                continue;
            }

            auto& dst_tile = dst_pc.DefineAndReturnParticleTile(
                lev, pti.index(), pti.LocalTileIndex());
            const auto np_dst = dst_tile.numParticles();

            const auto num_added = filterCopyTransformParticles<N>(
                dst_pc, dst_tile, src_tile, np_dst, filter, Copy, transform);
            setNewParticleIDs(dst_tile, np_dst, num_added);

            total_added += static_cast<amrex::Long>(num_added);
        }
    }

    return total_added;
}

template <int N = 1, typename DstPC, typename FilterFunc,
          typename TransformFunc>
amrex::Long
FilterCopyTransformBoundaryBuffer (ParticleBoundaryBuffer& boundary_buffer,
                                   std::string const& src_species,
                                   int const boundary, DstPC& dst_pc,
                                   FilterFunc const& filter,
                                   TransformFunc const& transform,
                                   int lev_min = 0, int lev_max = -1) {
    auto* src_pc =
        boundary_buffer.getParticleBufferPointer(src_species, boundary);

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_pc != nullptr && src_pc->isDefined(),
        "The requested absorbing-boundary particle buffer is not defined. "
        "Enable the corresponding save_particles_at_* option and call this "
        "after ParticleBoundaryBuffer has gathered absorbed particles.");

    return FilterCopyTransformParticleTiles<N>(*src_pc, dst_pc, filter,
                                               transform, lev_min, lev_max);
}

template <typename SrcPC, typename DstPC, typename DecisionFunc,
          typename TransformFunc>
amrex::Long
VariableCountCopyTransformParticleTiles (SrcPC& src_pc, DstPC& dst_pc,
                                         DecisionFunc const& decision_func,
                                         TransformFunc const& transform,
                                         int lev_min = 0, int lev_max = -1) {
    // Some wall interactions create a variable number of particles from one
    // absorbed particle, so the fixed-N helper is not expressive enough.
#ifdef AMREX_USE_GPU
    auto const* src_arena = src_pc.arena();
    auto const* dst_arena = dst_pc.arena();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_arena->isDevice() || src_arena->isManaged() ||
            src_arena == amrex::The_Pinned_Arena(),
        "VariableCountCopyTransformParticleTiles source particles must be in "
        "a device, managed, or AMReX pinned arena for GPU execution.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        dst_arena->isDevice() || dst_arena->isManaged() ||
            dst_arena == amrex::The_Pinned_Arena(),
        "VariableCountCopyTransformParticleTiles destination particles must be "
        "in a device, managed, or AMReX pinned arena for GPU execution.");
#endif

    const SmartCopyFactory copy_factory(src_pc, dst_pc);
    const auto Copy = copy_factory.getSmartCopy();

    if (lev_max < 0) {
        lev_max = src_pc.numLevels() - 1;
    }
    lev_max = std::min(lev_max, src_pc.numLevels() - 1);

    amrex::Long total_added = 0;

    for (int lev = lev_min; lev <= lev_max; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())                      \
    reduction(+ : total_added)
#endif
        for (WarpXParIter pti(src_pc, lev); pti.isValid(); ++pti) {
            auto& src_tile = pti.GetParticleTile();
            const auto np_src = src_tile.numParticles();
            if (np_src == 0) {
                continue;
            }

            auto& dst_tile = dst_pc.DefineAndReturnParticleTile(
                lev, pti.index(), pti.LocalTileIndex());
            const auto np_dst = dst_tile.numParticles();

            amrex::Gpu::DeviceVector<int> emit_count(np_src);
            amrex::Gpu::DeviceVector<int> behavior(np_src);

            auto* p_emit_count = emit_count.dataPtr();
            auto* p_behavior = behavior.dataPtr();
            const auto src_data = src_tile.getParticleTileData();

            // Draw the physical decision once before the prefix sum. Reusing
            // this decision in transform keeps random behavior and emit count
            // consistent for the same incident particle.
            amrex::ParallelForRNG(
                np_src, [=] AMREX_GPU_DEVICE(
                            int i, amrex::RandomEngine const& engine) noexcept {
                    const auto decision = decision_func(src_data, i, engine);
                    p_emit_count[i] = decision.m_emit_count;
                    p_behavior[i] = decision.m_behavior;
                });

            amrex::Gpu::DeviceVector<int> offsets(np_src);
            // emit_count is the scan input, not a boolean mask: its value
            // directly determines how many destination slots each source owns.
            const int num_added =
                amrex::Scan::ExclusiveSum(np_src, p_emit_count, offsets.data());
            if (num_added == 0) {
                continue;
            }

            const auto old_np = dst_tile.size();
            const auto new_np =
                std::max(np_dst + num_added, dst_tile.numParticles());
            dst_tile.resize(new_np);

            auto* p_offsets = offsets.dataPtr();
            const auto dst_data = dst_tile.getParticleTileData();

            amrex::ParallelForRNG(
                np_src, [=] AMREX_GPU_DEVICE(
                            int i, amrex::RandomEngine const& engine) noexcept {
                    const int count = p_emit_count[i];
                    for (int j = 0; j < count; ++j) {
                        const int i_dst = np_dst + p_offsets[i] + j;
                        Copy(dst_data, src_data, i, i_dst, engine);
                        transform(dst_data, src_data, i, i_dst, p_behavior[i],
                                  j, count, engine);
                    }
                });

            ParticleCreation::DefaultInitializeRuntimeAttributes(
                dst_tile, dst_pc, old_np, new_np);
            amrex::Gpu::synchronize();

            setNewParticleIDs(dst_tile, np_dst, num_added);
            total_added += static_cast<amrex::Long>(num_added);
        }
    }

    return total_added;
}

template <typename DstPC, typename DecisionFunc, typename TransformFunc>
amrex::Long
VariableCountCopyTransformBoundaryBuffer (
    ParticleBoundaryBuffer& boundary_buffer, std::string const& src_species,
    int const boundary, DstPC& dst_pc, DecisionFunc const& decision_func,
    TransformFunc const& transform, int lev_min = 0, int lev_max = -1) {
    auto* src_pc =
        boundary_buffer.getParticleBufferPointer(src_species, boundary);

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_pc != nullptr && src_pc->isDefined(),
        "The requested absorbing-boundary particle buffer is not defined. "
        "Enable the corresponding save_particles_at_* option and call this "
        "after ParticleBoundaryBuffer has gathered absorbed particles.");

    return VariableCountCopyTransformParticleTiles(
        *src_pc, dst_pc, decision_func, transform, lev_min, lev_max);
}

enum class NeutralAtomReflectionBehavior : int {
    Diffuse = 0,
    Specular = 1
};

struct NeutralAtomReflectionDecision {
    int m_emit_count = 0;
    int m_behavior =
        static_cast<int>(NeutralAtomReflectionBehavior::Diffuse);
};

struct NeutralAtomReflectionDecisionFunc {
    int m_behavior =
        static_cast<int>(NeutralAtomReflectionBehavior::Diffuse);
    TruncatedConeGeometry m_cone;

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    NeutralAtomReflectionDecision
    operator()(PData const& ptd, int const i,
               amrex::RandomEngine const& engine) const noexcept {
#if defined(WARPX_DIM_3D)
        amrex::XDim3 const x_hit{
            ptd.m_rdata[PIdx::x][i],
            ptd.m_rdata[PIdx::y][i],
            ptd.m_rdata[PIdx::z][i]};
        amrex::XDim3 normal_to_domain;
        if (!GetTruncatedConeNormal(m_cone, x_hit, normal_to_domain))
        {
            // Plane-section impacts are deliberately ignored for this case.
            return {};
        }
#else
        amrex::ignore_unused(ptd, i, m_cone);
        return {};
#endif

        // TODO: Add per-particle model selection here if the reflection
        // behavior later depends on probability, incident state or wall data.
        amrex::ignore_unused(engine);
        return {1, m_behavior};
    }
};

struct SpecularReflectionOperator {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator()(
        amrex::XDim3 const& normal_to_domain,
        amrex::XDim3 const& u_in, amrex::XDim3& u_out,
        amrex::RandomEngine const& engine) const noexcept {
        amrex::ParticleReal const normal_norm_sq =
            normal_to_domain.x * normal_to_domain.x
            + normal_to_domain.y * normal_to_domain.y
            + normal_to_domain.z * normal_to_domain.z;
        if (normal_norm_sq <= amrex::ParticleReal(0.0)) {
            u_out = u_in;
            amrex::ignore_unused(engine);
            return;
        }

        using std::sqrt;
        amrex::ParticleReal const inv_normal_norm =
            amrex::ParticleReal(1.0) / sqrt(normal_norm_sq);
        amrex::XDim3 const normal{
            normal_to_domain.x * inv_normal_norm,
            normal_to_domain.y * inv_normal_norm,
            normal_to_domain.z * inv_normal_norm};
        amrex::ParticleReal const u_dot_normal =
            u_in.x * normal.x + u_in.y * normal.y + u_in.z * normal.z;

        u_out = {
            u_in.x - amrex::ParticleReal(2.0) * u_dot_normal * normal.x,
            u_in.y - amrex::ParticleReal(2.0) * u_dot_normal * normal.y,
            u_in.z - amrex::ParticleReal(2.0) * u_dot_normal * normal.z};

        amrex::ignore_unused(engine);
    }
};

struct DiffuseReemissionOperator {
    amrex::ParticleReal m_vth = 0.0;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator()(
        amrex::XDim3 const& normal_to_domain,
        amrex::XDim3 const& u_in, amrex::XDim3& u_out,
        amrex::RandomEngine const& engine) const noexcept {
        amrex::ParticleReal const normal_norm_sq =
            normal_to_domain.x * normal_to_domain.x
            + normal_to_domain.y * normal_to_domain.y
            + normal_to_domain.z * normal_to_domain.z;
        if (normal_norm_sq <= amrex::ParticleReal(0.0)) {
            u_out = u_in;
            return;
        }

        using std::abs;
        using std::sqrt;
        amrex::ParticleReal const inv_normal_norm =
            amrex::ParticleReal(1.0) / sqrt(normal_norm_sq);
        amrex::XDim3 const normal{
            normal_to_domain.x * inv_normal_norm,
            normal_to_domain.y * inv_normal_norm,
            normal_to_domain.z * inv_normal_norm};

        amrex::XDim3 tangent_one;
        if (abs(normal.z) < amrex::ParticleReal(0.9)) {
            amrex::ParticleReal const inv_tangent_norm =
                amrex::ParticleReal(1.0)
                / sqrt(normal.x * normal.x + normal.y * normal.y);
            tangent_one = {
                -normal.y * inv_tangent_norm,
                normal.x * inv_tangent_norm,
                amrex::ParticleReal(0.0)};
        } else {
            amrex::ParticleReal const inv_tangent_norm =
                amrex::ParticleReal(1.0)
                / sqrt(normal.y * normal.y + normal.z * normal.z);
            tangent_one = {
                amrex::ParticleReal(0.0),
                normal.z * inv_tangent_norm,
                -normal.y * inv_tangent_norm};
        }
        amrex::XDim3 const tangent_two{
            normal.y * tangent_one.z - normal.z * tangent_one.y,
            normal.z * tangent_one.x - normal.x * tangent_one.z,
            normal.x * tangent_one.y - normal.y * tangent_one.x};

        amrex::ParticleReal const u_tangent_one =
            amrex::RandomNormal(amrex::ParticleReal(0.0), m_vth, engine);
        amrex::ParticleReal const u_tangent_two =
            amrex::RandomNormal(amrex::ParticleReal(0.0), m_vth, engine);
        amrex::ParticleReal const u_normal =
            generateGaussianFluxDist(
                amrex::ParticleReal(0.0), m_vth, engine);

        u_out = {
            u_tangent_one * tangent_one.x
                + u_tangent_two * tangent_two.x + u_normal * normal.x,
            u_tangent_one * tangent_one.y
                + u_tangent_two * tangent_two.y + u_normal * normal.y,
            u_tangent_one * tangent_one.z
                + u_tangent_two * tangent_two.z + u_normal * normal.z};

        amrex::ignore_unused(u_in);
    }
};

struct NeutralAtomReflectionTransform {
    TruncatedConeGeometry m_cone;
    amrex::ParticleReal m_position_epsilon = 0.0;
    SpecularReflectionOperator m_specular;
    DiffuseReemissionOperator m_diffuse;

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator()(DstData& dst, SrcData const& src, int const i_src,
                    int const i_dst, int const behavior,
                    int const emission_index, int const emit_count,
                    amrex::RandomEngine const& engine) const noexcept {
#if defined(WARPX_DIM_3D)
        amrex::XDim3 const u_in{
            src.m_rdata[PIdx::ux][i_src],
            src.m_rdata[PIdx::uy][i_src],
            src.m_rdata[PIdx::uz][i_src]};
        amrex::XDim3 const x_hit{
            src.m_rdata[PIdx::x][i_src],
            src.m_rdata[PIdx::y][i_src],
            src.m_rdata[PIdx::z][i_src]};
        amrex::XDim3 normal_to_domain;
        bool const on_cone =
            GetTruncatedConeNormal(m_cone, x_hit, normal_to_domain);
        amrex::XDim3 x_out = x_hit;
        amrex::XDim3 u_out = u_in;

        if (behavior ==
            static_cast<int>(NeutralAtomReflectionBehavior::Specular))
        {
            m_specular(normal_to_domain, u_in, u_out, engine);
        } else {
            m_diffuse(normal_to_domain, u_in, u_out, engine);
        }

        // Move a fixed distance along the reflected direction. Proper velocity
        // and physical velocity have the same direction, so no gamma
        // conversion is required for this displacement.
        using std::sqrt;
        amrex::ParticleReal const u_out_norm =
            sqrt(u_out.x * u_out.x + u_out.y * u_out.y +
                 u_out.z * u_out.z);
        if (u_out_norm > amrex::ParticleReal(0.0)) {
            amrex::ParticleReal const displacement_scale =
                m_position_epsilon / u_out_norm;
            x_out.x += displacement_scale * u_out.x;
            x_out.y += displacement_scale * u_out.y;
            x_out.z += displacement_scale * u_out.z;
        } else if (on_cone) {
            x_out.x += m_position_epsilon * normal_to_domain.x;
            x_out.y += m_position_epsilon * normal_to_domain.y;
            x_out.z += m_position_epsilon * normal_to_domain.z;
        }

        dst.m_rdata[PIdx::x][i_dst] = x_out.x;
        dst.m_rdata[PIdx::y][i_dst] = x_out.y;
        dst.m_rdata[PIdx::z][i_dst] = x_out.z;
        dst.m_rdata[PIdx::ux][i_dst] = u_out.x;
        dst.m_rdata[PIdx::uy][i_dst] = u_out.y;
        dst.m_rdata[PIdx::uz][i_dst] = u_out.z;

        amrex::ignore_unused(emission_index, emit_count);
#else
        amrex::ignore_unused(
            dst, src, i_src, i_dst, behavior, emission_index, emit_count,
            engine, m_cone, m_position_epsilon, m_specular, m_diffuse);
#endif
    }
};

enum class SecondaryEmissionBehavior : int {
    Absorb = 0,
    Reflect = 1,
    EmitOne = 2,
    EmitTwo = 3
};

struct SecondaryEmissionDecision {
    // m_behavior selects the transform; m_emit_count selects allocation size.
    int m_emit_count = 0;
    int m_behavior = static_cast<int>(SecondaryEmissionBehavior::Absorb);
};

struct SecondaryEmissionDecisionFunc {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;
    HallAnodeRingConfig m_anode_ring;
    amrex::ParticleReal m_mass = PhysConst::m_e;
    amrex::ParticleReal m_absorb_a = 0.0;
    amrex::ParticleReal m_absorb_x0_eV = 1.0;
    amrex::ParticleReal m_reflect_a = 0.0;
    amrex::ParticleReal m_reflect_x0_eV = 1.0;
    amrex::ParticleReal m_emit_one_a = 0.0;
    amrex::ParticleReal m_emit_one_x0_eV = 1.0;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::ParticleReal
    probabilityTerm (amrex::ParticleReal a, amrex::ParticleReal x0_eV,
                     amrex::ParticleReal energy_eV) noexcept {
        using std::exp;

        const amrex::ParticleReal e_over_x0 = energy_eV / x0_eV;
        return a * exp(-e_over_x0 * e_over_x0);
    }

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE SecondaryEmissionDecision
    operator()(PData const& ptd, int const i,
               amrex::RandomEngine const& engine) const noexcept {
        if (!amrex::ParticleIDWrapper{ptd.m_idcpu[i]}.is_valid()) {
            return {};
        }
#if defined(WARPX_DIM_3D)
        const auto& p = ptd.getSuperParticle(i);
        amrex::ParticleReal x, y, z;
        get_particle_position(p, x, y, z);

        const amrex::ParticleReal ux_inc = ptd.m_rdata[PIdx::ux][i];
        const amrex::ParticleReal uy_inc = ptd.m_rdata[PIdx::uy][i];
        const amrex::ParticleReal uz_inc = ptd.m_rdata[PIdx::uz][i];
        const amrex::ParticleReal z_boundary = m_plo[2];

        // The boundary buffer stores the scraped particle after it crossed the
        // domain boundary. Backtrace to the actual wall-hit position before
        // applying ring geometry.
        amrex::ParticleReal x_emit, y_emit;
        BacktraceParticleToZPlane(x, y, z, ux_inc, uy_inc, uz_inc, z_boundary,
                                  x_emit, y_emit);
        // The anode ring is a hard absorber: classify as behavior 1 directly.
        if (IsHallAnodeRingHit(x_emit, y_emit, m_anode_ring)) {
            return {};
        }

        constexpr auto eV = PhysConst::q_e;
        // SEE probabilities are parameterized in terms of incident electron
        // energy in eV.
        const auto energy_eV = static_cast<amrex::ParticleReal>(
            Algorithms::KineticEnergy<double>(ux_inc, uy_inc, uz_inc, m_mass) /
            eV);

        // Sample the four behaviors in order. The selected constants are
        // assumed to satisfy p_absorb + p_reflect + p_emit_one < 1, leaving the
        // final else branch for two-electron emission.
        const amrex::ParticleReal p_absorb =
            probabilityTerm(m_absorb_a, m_absorb_x0_eV, energy_eV);
        const amrex::ParticleReal p_reflect =
            probabilityTerm(m_reflect_a, m_reflect_x0_eV, energy_eV);
        const amrex::ParticleReal p_emit_one =
            amrex::ParticleReal(1.0) -
            probabilityTerm(m_emit_one_a, m_emit_one_x0_eV, energy_eV);

        const amrex::ParticleReal selector = amrex::Random(engine);
        if (selector < p_absorb) {
            return {};
        }
        if (selector < p_absorb + p_reflect) {
            return {1, static_cast<int>(SecondaryEmissionBehavior::Reflect)};
        }
        if (selector < p_absorb + p_reflect + p_emit_one) {
            return {1, static_cast<int>(SecondaryEmissionBehavior::EmitOne)};
        }
        return {2, static_cast<int>(SecondaryEmissionBehavior::EmitTwo)};
#else
        amrex::ignore_unused(ptd, i, engine);
        return {};
#endif
    }
};

struct SecondaryEmissionTransform {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_phi;
    amrex::ParticleReal m_eps = 0.0;
    amrex::ParticleReal m_vth = 0.0;
    ZMinWallChargeGrid m_wall_charge_grid;
    amrex::Array4<amrex::Real> m_wall_charge_density;
    bool m_deposit_wall_charge = false;
    amrex::Real m_charge_density_per_weight = amrex::Real(0.0);

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    operator()(DstData& dst, SrcData const& src, int const i_src,
               int const i_dst, int const behavior, int const emission_index,
               int const emit_count,
               amrex::RandomEngine const& engine) const noexcept {
#if defined(WARPX_DIM_3D)
        using std::abs;

        const auto& p = src.getSuperParticle(i_src);
        amrex::ParticleReal x, y, z;
        get_particle_position(p, x, y, z);

        const amrex::ParticleReal ux_inc = src.m_rdata[PIdx::ux][i_src];
        const amrex::ParticleReal uy_inc = src.m_rdata[PIdx::uy][i_src];
        const amrex::ParticleReal uz_inc = src.m_rdata[PIdx::uz][i_src];
        const amrex::ParticleReal z_boundary = m_plo[2];

        amrex::ParticleReal x_hit, y_hit;
        BacktraceParticleToZPlane(x, y, z, ux_inc, uy_inc, uz_inc, z_boundary,
                                  x_hit, y_hit);

        amrex::ParticleReal x_emit = x_hit;
        amrex::ParticleReal y_emit = y_hit;
        const amrex::ParticleReal xlo = m_plo[0] + m_eps;
        const amrex::ParticleReal xhi = m_phi[0] - m_eps;
        const amrex::ParticleReal ylo = m_plo[1] + m_eps;
        const amrex::ParticleReal yhi = m_phi[1] - m_eps;
        if (x_emit < xlo) {
            x_emit = xlo;
        }
        if (x_emit > xhi) {
            x_emit = xhi;
        }
        if (y_emit < ylo) {
            y_emit = ylo;
        }
        if (y_emit > yhi) {
            y_emit = yhi;
        }

        dst.m_rdata[PIdx::x][i_dst] = x_emit;
        dst.m_rdata[PIdx::y][i_dst] = y_emit;
        dst.m_rdata[PIdx::z][i_dst] = z_boundary + m_eps;

        if (m_deposit_wall_charge) {
            amrex::Real const emitted_charge_density =
                static_cast<amrex::Real>(src.m_rdata[PIdx::w][i_src]) *
                m_charge_density_per_weight;
            // Electrons have negative charge, so subtracting emitted charge
            // raises the wall charge density.
            DepositZMinWallChargeToNodes(m_wall_charge_density,
                                         m_wall_charge_grid, x_hit, y_hit,
                                         -emitted_charge_density);
        }

        // Reflection keeps the incident tangential velocity and flips the
        // normal component back into the domain.
        if (behavior == static_cast<int>(SecondaryEmissionBehavior::Reflect)) {
            dst.m_rdata[PIdx::ux][i_dst] = ux_inc;
            dst.m_rdata[PIdx::uy][i_dst] = uy_inc;
            dst.m_rdata[PIdx::uz][i_dst] = abs(uz_inc);
            return;
        }

        amrex::ignore_unused(emission_index, emit_count);
        // True secondary electrons use a half-Maxwellian normal distribution
        // and thermal tangential components.
        dst.m_rdata[PIdx::ux][i_dst] =
            amrex::RandomNormal(amrex::ParticleReal(0.0), m_vth, engine);
        dst.m_rdata[PIdx::uy][i_dst] =
            amrex::RandomNormal(amrex::ParticleReal(0.0), m_vth, engine);
        dst.m_rdata[PIdx::uz][i_dst] =
            abs(amrex::RandomNormal(amrex::ParticleReal(0.0), m_vth, engine));
#else
        amrex::ignore_unused(dst, src, i_src, i_dst, behavior, emission_index,
                             emit_count, engine, m_wall_charge_grid,
                             m_wall_charge_density,
                             m_charge_density_per_weight);
#endif
    }
};

struct AnodeRingIonFilter {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;
    HallAnodeRingConfig m_anode_ring;
    bool m_select_ring = true;

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
    operator()(PData const& ptd, int const i,
               amrex::RandomEngine const& /*engine*/) const noexcept {
        if (!amrex::ParticleIDWrapper{ptd.m_idcpu[i]}.is_valid()) {
            return false;
        }
#if defined(WARPX_DIM_3D)
        const auto& p = ptd.getSuperParticle(i);
        amrex::ParticleReal x, y, z;
        get_particle_position(p, x, y, z);

        const amrex::ParticleReal ux_inc = ptd.m_rdata[PIdx::ux][i];
        const amrex::ParticleReal uy_inc = ptd.m_rdata[PIdx::uy][i];
        const amrex::ParticleReal uz_inc = ptd.m_rdata[PIdx::uz][i];
        const amrex::ParticleReal z_boundary = m_plo[2];

        amrex::ParticleReal x_emit, y_emit;
        BacktraceParticleToZPlane(x, y, z, ux_inc, uy_inc, uz_inc, z_boundary,
                                  x_emit, y_emit);
        const bool is_ring_hit =
            IsHallAnodeRingHit(x_emit, y_emit, m_anode_ring);
        return m_select_ring ? is_ring_hit : !is_ring_hit;
#else
        amrex::ignore_unused(ptd, i);
        return false;
#endif
    }
};

struct AnodeIonEmissionTransform {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_phi;
    amrex::ParticleReal m_eps = 0.0;
    amrex::ParticleReal m_vth = 0.0;

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    operator()(DstData& dst, SrcData const& src, int const i_src,
               int const i_dst,
               amrex::RandomEngine const& engine) const noexcept {
#if defined(WARPX_DIM_3D)
        const auto& p = src.getSuperParticle(i_src);
        amrex::ParticleReal x, y, z;
        get_particle_position(p, x, y, z);

        const amrex::ParticleReal ux_inc = src.m_rdata[PIdx::ux][i_src];
        const amrex::ParticleReal uy_inc = src.m_rdata[PIdx::uy][i_src];
        const amrex::ParticleReal uz_inc = src.m_rdata[PIdx::uz][i_src];
        const amrex::ParticleReal z_boundary = m_plo[2];

        amrex::ParticleReal x_emit, y_emit;
        BacktraceParticleToZPlane(x, y, z, ux_inc, uy_inc, uz_inc, z_boundary,
                                  x_emit, y_emit);

        const amrex::ParticleReal xlo = m_plo[0] + m_eps;
        const amrex::ParticleReal xhi = m_phi[0] - m_eps;
        const amrex::ParticleReal ylo = m_plo[1] + m_eps;
        const amrex::ParticleReal yhi = m_phi[1] - m_eps;
        if (x_emit < xlo) {
            x_emit = xlo;
        }
        if (x_emit > xhi) {
            x_emit = xhi;
        }
        if (y_emit < ylo) {
            y_emit = ylo;
        }
        if (y_emit > yhi) {
            y_emit = yhi;
        }

        dst.m_rdata[PIdx::x][i_dst] = x_emit;
        dst.m_rdata[PIdx::y][i_dst] = y_emit;
        dst.m_rdata[PIdx::z][i_dst] = z_boundary + m_eps;

        dst.m_rdata[PIdx::ux][i_dst] =
            amrex::RandomNormal(amrex::ParticleReal(0.0), m_vth, engine);
        dst.m_rdata[PIdx::uy][i_dst] =
            amrex::RandomNormal(amrex::ParticleReal(0.0), m_vth, engine);
        dst.m_rdata[PIdx::uz][i_dst] =
            generateGaussianFluxDist(amrex::ParticleReal(0.0), m_vth, engine);
#else
        amrex::ignore_unused(dst, src, i_src, i_dst, engine);
#endif
    }
};

} // namespace

void
NeutralAtomEBInteraction () {
#if defined(HALL3D) && defined(AMREX_USE_EB)
    static bool const enabled = [] {
        amrex::ParmParse const pp("insert.neutral_atom_eb");
        int value = 0;
        pp.query("enabled", value);
        return value != 0;
    }();

    if (!enabled) {
        return;
    }

    static std::string const species = [] {
        amrex::ParmParse const pp("insert.neutral_atom_eb");
        std::string value;
        pp.query("species", value);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            !value.empty(),
            "insert.neutral_atom_eb.species must be specified when the neutral "
            "atom EB interaction is enabled.");
        return value;
    }();

    static int const reflection_behavior = [] {
        amrex::ParmParse const pp("insert.neutral_atom_eb");
        std::string model = "diffuse";
        pp.query("model", model);
        model = ToLower(model);

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            model == "diffuse" || model == "specular",
            "insert.neutral_atom_eb.model must be either diffuse or specular.");

        return model == "specular"
                   ? static_cast<int>(
                         NeutralAtomReflectionBehavior::Specular)
                   : static_cast<int>(
                         NeutralAtomReflectionBehavior::Diffuse);
    }();

    static TruncatedConeGeometry const cone = [] {
        amrex::ParmParse const pp("insert.neutral_atom_eb");
        TruncatedConeGeometry value;
        pp.get("k", value.m_k);
        pp.get("a1", value.m_a1);
        pp.get("b1", value.m_b1);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            value.m_k != amrex::ParticleReal(0.0),
            "insert.neutral_atom_eb.k must be non-zero.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            value.m_a1 > amrex::ParticleReal(0.0),
            "insert.neutral_atom_eb.a1 must be positive.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            value.m_b1 > value.m_a1,
            "insert.neutral_atom_eb.b1 must be greater than a1.");
        return value;
    }();

    static amrex::ParticleReal const wall_temperature = [] {
        amrex::ParmParse const pp("insert.neutral_atom_eb");
        amrex::ParticleReal value = 0.0;
        pp.query("wall_temperature", value);
        return value;
    }();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        reflection_behavior !=
                static_cast<int>(NeutralAtomReflectionBehavior::Diffuse) ||
            wall_temperature > amrex::ParticleReal(0.0),
        "insert.neutral_atom_eb.wall_temperature must be positive for "
        "diffuse reflection.");

    static amrex::ParticleReal const configured_position_epsilon = [] {
        amrex::ParmParse const pp("insert.neutral_atom_eb");
        amrex::ParticleReal value = -1.0;
        pp.query("position_epsilon", value);
        return value;
    }();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        configured_position_epsilon == amrex::ParticleReal(-1.0) ||
            configured_position_epsilon > amrex::ParticleReal(0.0),
        "insert.neutral_atom_eb.position_epsilon must be positive when "
        "specified.");

    WarpX& warpx_instance = WarpX::GetInstance();
    if (!DoBoundaryParticleDiag(warpx_instance.getistep(0))) {
        return;
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        EB::enabled(),
        "insert.neutral_atom_eb.enabled requires an embedded boundary.");

    auto& boundary_buffer = warpx_instance.GetParticleBoundaryBuffer();
    auto& particle_container = warpx_instance.GetPartContainer();
    auto& neutral_atoms =
        particle_container.GetParticleContainerFromName(species);
    constexpr int eb_boundary = AMREX_SPACEDIM * 2;
    auto* neutral_buffer =
        boundary_buffer.getParticleBufferPointer(species, eb_boundary);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        neutral_buffer != nullptr && neutral_buffer->isDefined(),
        "The neutral atom EB buffer is not defined. Set "
        "<species>.save_particles_at_eb = 1.");

    auto const dx = warpx_instance.Geom(0).CellSizeArray();
    amrex::ParticleReal min_dx =
        static_cast<amrex::ParticleReal>(dx[0]);
    for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
        min_dx = std::min(
            min_dx, static_cast<amrex::ParticleReal>(dx[idim]));
    }
    amrex::ParticleReal position_epsilon =
        amrex::ParticleReal(0.1) * min_dx;
    if (configured_position_epsilon > amrex::ParticleReal(0.0)) {
        position_epsilon = configured_position_epsilon;
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        position_epsilon > amrex::ParticleReal(0.0),
        "insert.neutral_atom_eb.position_epsilon must be positive.");

    amrex::ParticleReal diffuse_vth = 0.0;
    if (reflection_behavior ==
        static_cast<int>(NeutralAtomReflectionBehavior::Diffuse))
    {
        diffuse_vth = static_cast<amrex::ParticleReal>(
            std::sqrt(PhysConst::kb * wall_temperature /
                      neutral_atoms.getMass()));
    }

    NeutralAtomReflectionDecisionFunc const decision_func{
        reflection_behavior, cone};
    NeutralAtomReflectionTransform const transform{
        cone,
        position_epsilon,
        SpecularReflectionOperator{},
        DiffuseReemissionOperator{diffuse_vth}};

    amrex::Long const num_reflected =
        VariableCountCopyTransformBoundaryBuffer(
            boundary_buffer, species, eb_boundary, neutral_atoms,
            decision_func, transform);
    amrex::Print() << "Reflected neutral atoms: " << num_reflected << "\n";

    // The neutral buffer is intentionally left untouched here. It is cleared
    // together with the electron and ion boundary buffers by the existing
    // boundary-particle cleanup path.
#endif
}

void
SecondaryEmission () {
#ifdef HALL3D
    WarpX& warpx_instance = WarpX::GetInstance();
    if (DoBoundaryParticleDiag(warpx_instance.getistep(0))) {
        auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();

        auto& mypc = warpx_instance.GetPartContainer();
        auto& elec_pc = mypc.GetParticleContainerFromName("electrons");

        const auto plo = warpx_instance.Geom(0).ProbLoArray();
        const auto phi = warpx_instance.Geom(0).ProbHiArray();
        const auto dx = warpx_instance.Geom(0).CellSizeArray();
        amrex::Real min_dx = dx[0];
        for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
            min_dx = std::min(min_dx, dx[idim]);
        }

        HallAnodeRingConfig const anode_ring =
            ReadHallAnodeRingConfig(warpx_instance.Geom(0));

        constexpr amrex::ParticleReal absorb_a = amrex::ParticleReal(0.5);
        constexpr amrex::ParticleReal absorb_x0_eV = amrex::ParticleReal(43.5);
        constexpr amrex::ParticleReal reflect_a = amrex::ParticleReal(0.5);
        constexpr amrex::ParticleReal reflect_x0_eV = amrex::ParticleReal(30.0);
        constexpr amrex::ParticleReal emit_one_a = amrex::ParticleReal(1);
        constexpr amrex::ParticleReal emit_one_x0_eV =
            amrex::ParticleReal(127.9);

        constexpr amrex::ParticleReal emission_temperature_eV =
            amrex::ParticleReal(3.0);
        const amrex::ParticleReal emission_vth =
            static_cast<amrex::ParticleReal>(std::sqrt(
                emission_temperature_eV * PhysConst::q_e / elec_pc.getMass()));

        ZMinWallChargeGrid wall_charge_grid{};
        amrex::Array4<amrex::Real> wall_charge_density;
        bool deposit_wall_charge = false;
        amrex::Real charge_density_per_weight = amrex::Real(0.0);
        if (g_accumulated_wall_charge_density != nullptr) {
            wall_charge_grid = MakeZMinWallChargeGrid(warpx_instance.Geom(0));
            for (amrex::MFIter mfi(*g_accumulated_wall_charge_density);
                 mfi.isValid(); ++mfi)
            {
                wall_charge_density =
                    g_accumulated_wall_charge_density->array(mfi);
            }
            deposit_wall_charge = true;
            charge_density_per_weight =
                static_cast<amrex::Real>(elec_pc.getCharge()) /
                (wall_charge_grid.dx * wall_charge_grid.dy);
        }

        const SecondaryEmissionDecisionFunc decision_func{
            plo,
            anode_ring,
            elec_pc.getMass(),
            absorb_a,
            absorb_x0_eV,
            reflect_a,
            reflect_x0_eV,
            emit_one_a,
            emit_one_x0_eV};
        const SecondaryEmissionTransform transform{
            plo,
            phi,
            amrex::ParticleReal(0.1 * min_dx),
            emission_vth,
            wall_charge_grid,
            wall_charge_density,
            deposit_wall_charge,
            charge_density_per_weight};

        amrex::Long num_added = VariableCountCopyTransformBoundaryBuffer(
            mybpc, "electrons", 4, elec_pc, decision_func, transform);
        amrex::Print() << "Emission electron: " << num_added << "\n";
    }
#endif
}

void
AnodeIonNeutralization () {
#ifdef HALL3D
    WarpX& warpx_instance = WarpX::GetInstance();
    if (DoBoundaryParticleDiag(warpx_instance.getistep(0))) {
        auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();

        auto& mypc = warpx_instance.GetPartContainer();
        auto& atom_pc = mypc.GetParticleContainerFromName("xe_netural");

        const auto plo = warpx_instance.Geom(0).ProbLoArray();
        const auto phi = warpx_instance.Geom(0).ProbHiArray();
        const auto dx = warpx_instance.Geom(0).CellSizeArray();
        amrex::Real min_dx = dx[0];
        for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
            min_dx = std::min(min_dx, dx[idim]);
        }

        constexpr amrex::ParticleReal atom_temperature_K =
            amrex::ParticleReal(400.0);
        const amrex::ParticleReal atom_vth = static_cast<amrex::ParticleReal>(
            std::sqrt(PhysConst::kb * atom_temperature_K / atom_pc.getMass()));

        auto& ion_pc = mypc.GetParticleContainerFromName("xe_ions");

        HallAnodeRingConfig const anode_ring =
            ReadHallAnodeRingConfig(warpx_instance.Geom(0));

        constexpr amrex::ParticleReal ion_temperature_eV =
            amrex::ParticleReal(3.0);
        const amrex::ParticleReal ion_vth = static_cast<amrex::ParticleReal>(
            std::sqrt(ion_temperature_eV * PhysConst::q_e / ion_pc.getMass()));

        const AnodeRingIonFilter ring_filter{plo, anode_ring, true};
        const AnodeRingIonFilter non_ring_filter{plo, anode_ring, false};
        const AnodeIonEmissionTransform neutral_transform{
            plo, phi, amrex::ParticleReal(0.1 * min_dx), atom_vth};
        const AnodeIonEmissionTransform ion_transform{
            plo, phi, amrex::ParticleReal(0.1 * min_dx), ion_vth};

        const amrex::Long num_neutralized = FilterCopyTransformBoundaryBuffer(
            mybpc, "xe_ions", 4, atom_pc, ring_filter, neutral_transform);
        const amrex::Long num_ion_emitted = FilterCopyTransformBoundaryBuffer(
            mybpc, "xe_ions", 4, ion_pc, non_ring_filter, ion_transform);
        amrex::Print() << "Neutralized anode-ring Xe ions: " << num_neutralized
                       << "\n";
        amrex::Print() << "Re-emitted non-ring Xe ions: " << num_ion_emitted
                       << "\n";
    }
#endif
}

} // namespace Insert
