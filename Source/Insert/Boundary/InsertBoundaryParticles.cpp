#include "InsertBoundaryParticles.h"

#include "WarpX.H"

#include "Initialization/SampleGaussianFluxDistribution.H"
#include "Insert/Config/WarpXSimulationConfig.h"
#include "Particles/Algorithms/KineticEnergy.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Particles/ParticleCreation/FilterCopyTransform.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartUtils.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Utils/WarpXConst.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <algorithm>
#include <cmath>
#include <string>

namespace Insert {
namespace {

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
    // SEE can create 0, 1, or 2 particles from one absorbed particle, so the
    // fixed-N FilterCopyTransform helper is not expressive enough here.
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
            // emit_count is the scan input, not a boolean mask: values 0/1/2
            // directly determine how many destination slots each source owns.
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
    amrex::ParticleReal m_anode_ring_x = 0.0;
    amrex::ParticleReal m_anode_ring_y = 0.0;
    amrex::ParticleReal m_anode_ring_rmin = 0.0;
    amrex::ParticleReal m_anode_ring_rmax = 0.0;
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
        using std::sqrt;

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
        amrex::ParticleReal x_emit = x;
        amrex::ParticleReal y_emit = y;
        if (uz_inc != amrex::ParticleReal(0.0)) {
            const amrex::ParticleReal dt_back = (z - z_boundary) / uz_inc;
            if (dt_back >= amrex::ParticleReal(0.0)) {
                x_emit -= dt_back * ux_inc;
                y_emit -= dt_back * uy_inc;
            }
        }

        const amrex::ParticleReal dx = x_emit - m_anode_ring_x;
        const amrex::ParticleReal dy = y_emit - m_anode_ring_y;
        const amrex::ParticleReal r = sqrt(dx * dx + dy * dy);
        // The anode ring is a hard absorber: classify as behavior 1 directly.
        if (r >= m_anode_ring_rmin && r <= m_anode_ring_rmax) {
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

        amrex::ParticleReal x_emit = x;
        amrex::ParticleReal y_emit = y;
        if (uz_inc != amrex::ParticleReal(0.0)) {
            const amrex::ParticleReal dt_back = (z - z_boundary) / uz_inc;
            if (dt_back >= amrex::ParticleReal(0.0)) {
                x_emit -= dt_back * ux_inc;
                y_emit -= dt_back * uy_inc;
            }
        }

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
                             emit_count, engine);
#endif
    }
};

struct AnodeRingIonFilter {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;
    amrex::ParticleReal m_anode_ring_x = 0.0;
    amrex::ParticleReal m_anode_ring_y = 0.0;
    amrex::ParticleReal m_anode_ring_rmin = 0.0;
    amrex::ParticleReal m_anode_ring_rmax = 0.0;
    bool m_select_ring = true;

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
    operator()(PData const& ptd, int const i,
               amrex::RandomEngine const& /*engine*/) const noexcept {
        if (!amrex::ParticleIDWrapper{ptd.m_idcpu[i]}.is_valid()) {
            return false;
        }
#if defined(WARPX_DIM_3D)
        using std::sqrt;

        const auto& p = ptd.getSuperParticle(i);
        amrex::ParticleReal x, y, z;
        get_particle_position(p, x, y, z);

        const amrex::ParticleReal ux_inc = ptd.m_rdata[PIdx::ux][i];
        const amrex::ParticleReal uy_inc = ptd.m_rdata[PIdx::uy][i];
        const amrex::ParticleReal uz_inc = ptd.m_rdata[PIdx::uz][i];
        const amrex::ParticleReal z_boundary = m_plo[2];

        amrex::ParticleReal x_emit = x;
        amrex::ParticleReal y_emit = y;
        if (uz_inc != amrex::ParticleReal(0.0)) {
            const amrex::ParticleReal dt_back = (z - z_boundary) / uz_inc;
            if (dt_back >= amrex::ParticleReal(0.0)) {
                x_emit -= dt_back * ux_inc;
                y_emit -= dt_back * uy_inc;
            }
        }

        const amrex::ParticleReal dx = x_emit - m_anode_ring_x;
        const amrex::ParticleReal dy = y_emit - m_anode_ring_y;
        const amrex::ParticleReal r = sqrt(dx * dx + dy * dy);
        const bool is_ring_hit =
            r >= m_anode_ring_rmin && r <= m_anode_ring_rmax;
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

        amrex::ParticleReal x_emit = x;
        amrex::ParticleReal y_emit = y;
        if (uz_inc != amrex::ParticleReal(0.0)) {
            const amrex::ParticleReal dt_back = (z - z_boundary) / uz_inc;
            if (dt_back >= amrex::ParticleReal(0.0)) {
                x_emit -= dt_back * ux_inc;
                y_emit -= dt_back * uy_inc;
            }
        }

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
SecondaryEmission () {
#ifdef HALL3D
    static int times = 0;
    static const int gap = 10;
    times++;

    if (times == gap) {
        times = 0;

        WarpX& warpx_instance = WarpX::GetInstance();

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

        amrex::ParticleReal l_factor, anode_length;
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("l_factor", l_factor);
        pp_mc.get("L", anode_length);

        // Use the same scaled anode-ring geometry as anode current collection.
        const amrex::ParticleReal anode_ring_center =
            anode_length / l_factor / amrex::ParticleReal(2.0);
        const amrex::ParticleReal anode_ring_inner_radius =
            amrex::ParticleReal(0.021) / amrex::ParticleReal(2.0) / l_factor;
        const amrex::ParticleReal anode_ring_outer_radius =
            amrex::ParticleReal(0.031) / amrex::ParticleReal(2.0) / l_factor;

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

        const SecondaryEmissionDecisionFunc decision_func{
            plo,
            anode_ring_center,
            anode_ring_center,
            anode_ring_inner_radius,
            anode_ring_outer_radius,
            elec_pc.getMass(),
            absorb_a,
            absorb_x0_eV,
            reflect_a,
            reflect_x0_eV,
            emit_one_a,
            emit_one_x0_eV};
        const SecondaryEmissionTransform transform{
            plo, phi, amrex::ParticleReal(0.1 * min_dx), emission_vth};

        amrex::Long num_added = VariableCountCopyTransformBoundaryBuffer(
            mybpc, "electrons", 4, elec_pc, decision_func, transform);
        amrex::Print() << "Emission electron: " << num_added << "\n";
    }
#endif
}

void
AnodeIonNeutralization () {
#ifdef HALL3D
    static int times = 0;
    static const int gap = 10;
    times++;

    if (times == gap) {
        times = 0;

        WarpX& warpx_instance = WarpX::GetInstance();

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

        amrex::ParticleReal l_factor, anode_length;
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("l_factor", l_factor);
        pp_mc.get("L", anode_length);

        const amrex::ParticleReal anode_ring_center =
            anode_length / l_factor / amrex::ParticleReal(2.0);
        const amrex::ParticleReal anode_ring_inner_radius =
            amrex::ParticleReal(0.021) / amrex::ParticleReal(2.0) / l_factor;
        const amrex::ParticleReal anode_ring_outer_radius =
            amrex::ParticleReal(0.031) / amrex::ParticleReal(2.0) / l_factor;

        constexpr amrex::ParticleReal ion_temperature_eV =
            amrex::ParticleReal(3.0);
        const amrex::ParticleReal ion_vth = static_cast<amrex::ParticleReal>(
            std::sqrt(ion_temperature_eV * PhysConst::q_e / ion_pc.getMass()));

        const AnodeRingIonFilter ring_filter{
            plo, anode_ring_center, anode_ring_center,
            anode_ring_inner_radius, anode_ring_outer_radius, true};
        const AnodeRingIonFilter non_ring_filter{
            plo, anode_ring_center, anode_ring_center,
            anode_ring_inner_radius, anode_ring_outer_radius, false};
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
