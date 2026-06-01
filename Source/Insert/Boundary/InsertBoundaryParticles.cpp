#include "InsertBoundaryParticles.h"

#include "WarpX.H"

#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Particles/ParticleCreation/FilterCopyTransform.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartUtils.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Utils/WarpXConst.H"
#include "Insert/Config/WarpXSimulationConfig.h"

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
                                  int lev_min = 0, int lev_max = -1)
{
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
                                   int lev_min = 0, int lev_max = -1)
{
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

struct SecondaryEmissionFilter {
    amrex::ParticleReal m_probability = 0.1;

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
    operator() (PData const& ptd, int const i,
                amrex::RandomEngine const& engine) const noexcept
    {
        if (!amrex::ParticleIDWrapper{ptd.m_idcpu[i]}.is_valid()) {
            return false;
        }
        return amrex::Random(engine) < m_probability;
    }
};

struct SecondaryEmissionTransform {
    int m_normal_index = 0;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_phi;
    amrex::ParticleReal m_eps = 0.0;
    amrex::ParticleReal m_vz_std = 0.0;

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    operator() (DstData& dst, SrcData const& src, int const i_src,
                int const i_dst,
                amrex::RandomEngine const& engine) const noexcept
    {
#if defined(WARPX_DIM_3D)
        const auto& p = src.getSuperParticle(i_src);
        amrex::ParticleReal x, y, z;
        get_particle_position(p, x, y, z);

        const amrex::ParticleReal ux_inc = src.m_rdata[PIdx::ux][i_src];
        const amrex::ParticleReal uy_inc = src.m_rdata[PIdx::uy][i_src];
        const amrex::ParticleReal uz_inc = src.m_rdata[PIdx::uz][i_src];
        const amrex::ParticleReal nz =
            src.m_runtime_rdata[m_normal_index + 2][i_src];

        const amrex::ParticleReal normal_sign = (nz >= amrex::ParticleReal(0.0))
                                                    ? amrex::ParticleReal(1.0)
                                                    : amrex::ParticleReal(-1.0);
        const amrex::ParticleReal z_boundary =
            (normal_sign > amrex::ParticleReal(0.0)) ? m_plo[2] : m_phi[2];

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
        dst.m_rdata[PIdx::z][i_dst] = z_boundary + normal_sign * m_eps;

        amrex::ParticleReal vz_emit =
            amrex::RandomNormal(amrex::ParticleReal(0.0), m_vz_std, engine);
        if (vz_emit < amrex::ParticleReal(0.0)) {
            vz_emit = -vz_emit;
        }

        dst.m_rdata[PIdx::ux][i_dst] = amrex::ParticleReal(0.0);
        dst.m_rdata[PIdx::uy][i_dst] = amrex::ParticleReal(0.0);
        dst.m_rdata[PIdx::uz][i_dst] = normal_sign * vz_emit;
#else
        amrex::ignore_unused(dst, src, i_src, i_dst, engine);
#endif
    }
};

} // namespace

void
SecondaryEmission ()
{
#ifdef HALL3D
    static int times = 0;
    static const int gap = 10;
    times++;

    if (times == gap) {
        times = 0;

        WarpX& warpx_instance = WarpX::GetInstance();

        auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();
        auto& elec_zmin = mybpc.getParticleBuffer("electrons", 4);
        const int normal_index = elec_zmin.GetRealCompIndex("nx") -
                                 WarpXParticleContainer::NArrayReal;

        auto& mypc = warpx_instance.GetPartContainer();
        auto& elec_pc = mypc.GetParticleContainerFromName("electrons");

        const auto plo = warpx_instance.Geom(0).ProbLoArray();
        const auto phi = warpx_instance.Geom(0).ProbHiArray();
        const auto dx = warpx_instance.Geom(0).CellSizeArray();
        amrex::Real min_dx = dx[0];
        for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
            min_dx = std::min(min_dx, dx[idim]);
        }

        constexpr amrex::ParticleReal emission_energy_eV =
            amrex::ParticleReal(2.0);
        const amrex::ParticleReal vz_std =
            static_cast<amrex::ParticleReal>(std::sqrt(
                2.0 * emission_energy_eV * PhysConst::q_e / PhysConst::m_e));

        const SecondaryEmissionFilter filter{amrex::ParticleReal(0.5)};
        const SecondaryEmissionTransform transform{
            normal_index, plo, phi, amrex::ParticleReal(0.1 * min_dx), vz_std};

        int num_added = FilterCopyTransformBoundaryBuffer(
            mybpc, "electrons", 4, elec_pc, filter, transform);
        amrex::Print() << "Emission electron: " << num_added << "\n";
    }
#endif
}

} // namespace Insert
