#ifndef WARPX_INSERT_FILTERCOPYTRANSFORMCOUPLED_H_
#define WARPX_INSERT_FILTERCOPYTRANSFORMCOUPLED_H_

#include "Particles/Algorithms/KineticEnergy.H"
#include "Particles/ParticleCreation/FilterCopyTransform.H"
#include "Particles/ShapeFactors.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Dim3.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_Random.H>

#include <utility>

namespace Insert::CollisionDetail {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
DenseBinIndex (amrex::Dim3 const cell, amrex::Dim3 const ncell) noexcept {
    return (cell.z * ncell.y + cell.y) * ncell.x + cell.x;
}

/**
 * \brief Transform functor that applies only the incident-electron energy loss
 * and momentum change of an impact ionization event, without creating the
 * secondary electron and ion.
 *
 * This is used together with filterCopyTransformParticles when
 * IONIZATION_SOURCE_INJECT is enabled: the average ionization source is
 * injected elsewhere, so the coupled MCC path only needs to account for the
 * energy loss of the colliding electron.
 */
class IonizationSourceTransformFunc {
  public:
    IonizationSourceTransformFunc (amrex::ParticleReal energy_cost,
                                   double mass1)
        : m_energy_cost(energy_cost), m_mass1(mass1) {}

    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (SrcData& src, int const i_src,
                     amrex::RandomEngine const& engine) const noexcept {
        using namespace amrex;
        using std::sqrt;

        auto& ux = src.m_rdata[PIdx::ux][i_src];
        auto& uy = src.m_rdata[PIdx::uy][i_src];
        auto& uz = src.m_rdata[PIdx::uz][i_src];

        constexpr auto eV = PhysConst::q_e;
        const double E_coll =
            Algorithms::KineticEnergy<double>(ux, uy, uz, m_mass1) / eV;

        const auto E_out = static_cast<amrex::ParticleReal>(
            (E_coll - m_energy_cost) / 2.0_prt * eV);

        constexpr auto c2 = PhysConst::c2;
        const auto mc2 = m_mass1 * c2;

        const amrex::ParticleReal up =
            sqrt(E_out * (E_out + 2.0_prt * mc2) / c2) / m_mass1;

        ParticleUtils::RandomizeVelocity(ux, uy, uz, up, engine);
    }

  private:
    amrex::ParticleReal m_energy_cost;
    double m_mass1;
};

} // namespace Insert::CollisionDetail

/**
 * \brief Apply a filter, copy, and transform operation using a local
 * background-density field and record how many background particles should be
 * deleted from each cell.
 */
template <int N, int depos_order, typename DstPC, typename DstTile,
          typename SrcTile, typename Index, typename PredFunc,
          typename TransFunc, typename CopyFunc1, typename CopyFunc2>
Index
filterCopyTransformParticles (DstPC& pc1, DstPC& pc2, DstTile& dst1,
                              DstTile& dst2, SrcTile& src, Index dst1_index,
                              Index dst2_index, PredFunc&& filter,
                              CopyFunc1&& copy1, CopyFunc2&& copy2,
                              TransFunc&& transform, Index* p_delete,
                              const amrex::Array4<amrex::Real>& rho_array,
                              Index* p_particle_num, const amrex::XDim3 xyzmin,
                              const amrex::Box& box, const amrex::IntVect ngrow,
                              const amrex::XDim3& inv_cell_size,
                              const bool consume_background = true) noexcept {
    using namespace amrex;

    auto np = src.numParticles();
    if (np == 0) {
        return 0;
    }

    Gpu::DeviceVector<Index> mask(np);

    auto* p_mask = mask.dataPtr();
    const auto src_data = src.getParticleTileData();

    const Dim3 lo = lbound(box), hi = ubound(box), grow_dim = ngrow.dim3();
    Dim3 const ncell{hi.x - lo.x - (2 * grow_dim.x),
                     hi.y - lo.y - (2 * grow_dim.y),
                     hi.z - lo.z - (2 * grow_dim.z)};
    amrex::ParallelForRNG(
        np, [=] AMREX_GPU_DEVICE(int i, amrex::RandomEngine const& engine) {
            const auto& p = src_data.getSuperParticle(i);

            ParticleReal x, y, z;
            get_particle_position(p, x, y, z);
            const amrex::ParticleReal rpx = (x - xyzmin.x) * inv_cell_size.x,
                                      rpy = (y - xyzmin.y) * inv_cell_size.y,
                                      rpz = (z - xyzmin.z) * inv_cell_size.z;
            amrex::ParticleReal n_a = 0;
            const int pi = static_cast<int>(rpx - grow_dim.x),
                      pj = static_cast<int>(rpy - grow_dim.y),
                      pk = static_cast<int>(rpz - grow_dim.z);
            Compute_shape_factor<depos_order> const compute_shape_factor;
            amrex::Real sx[depos_order + 1] = {0._rt},
                                         sy[depos_order + 1] = {0._rt},
                                         sz[depos_order + 1] = {0._rt};

            int px = compute_shape_factor(sx, rpx),
                py = compute_shape_factor(sy, rpy),
                pz = compute_shape_factor(sz, rpz);

#if defined(WARPX_DIM_1D_Z)
            for (int ix = 0; ix <= depos_order; ix++) {
                n_a += sx[ix] * rho_array(lo.x + px + ix, 0, 0);
            }
#elif defined(WARPX_DIM_XZ)
            for (int iy = 0; iy <= depos_order; iy++) {
                for (int ix = 0; ix <= depos_order; ix++) {
                    n_a += sx[ix] * sy[iy] *
                           rho_array(lo.x + px + ix, lo.y + py + iy, 0);
                }
            }
#elif defined(WARPX_DIM_3D)
            for (int iz = 0; iz <= depos_order; iz++) {
                for (int iy = 0; iy <= depos_order; iy++) {
                    for (int ix = 0; ix <= depos_order; ix++) {
                        n_a += sx[ix] * sy[iy] * sz[iz] *
                               rho_array(lo.x + px + ix, lo.y + py + iy,
                                         lo.z + pz + iz);
                    }
                }
            }
#endif
            p_mask[i] = filter(src_data, i, engine, n_a);

            if (p_mask[i] && consume_background) {
                int pos = Insert::CollisionDetail::DenseBinIndex(
                    amrex::Dim3{pi, pj, pk}, ncell);
                Index np_old = p_particle_num[pos];
                while (np_old > 0) {
                    auto returned = amrex::Gpu::Atomic::CAS(
                        &p_particle_num[pos], np_old, np_old - 1);
                    if (returned == np_old) {
                        amrex::Gpu::Atomic::Add(&p_delete[pos], 1);
                        break;
                    }
                    np_old = p_particle_num[pos];
                }
                if (np_old <= 0) {
                    p_mask[i] = false;
                }
            }
        });

#ifdef IONIZATION_SOURCE_INJECT
    // When the average ionization source is injected externally, the MCC path
    // only applies the energy loss / momentum change to the incident electron
    // and does not create new electron-ion pairs.
    amrex::ignore_unused(pc1, pc2, dst1, dst2, dst1_index, dst2_index, copy1,
                         copy2);
    amrex::ParallelForRNG(
        np, [=] AMREX_GPU_HOST_DEVICE(int i,
                                      amrex::RandomEngine const& engine) {
            if (p_mask[i]) {
                transform(src_data, i, engine);
            }
        });
    return 0;
#else
    return filterCopyTransformParticles<N>(
        pc1, pc2, dst1, dst2, src, mask.dataPtr(), dst1_index, dst2_index,
        std::forward<CopyFunc1>(copy1), std::forward<CopyFunc2>(copy2),
        std::forward<TransFunc>(transform));
#endif
}

#endif // WARPX_INSERT_FILTERCOPYTRANSFORMCOUPLED_H_
