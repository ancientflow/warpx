#ifndef WARPX_INSERT_HALLWALLCHARGE_H_
#define WARPX_INSERT_HALLWALLCHARGE_H_

#include "Insert/Math/InterpUtils.h"

#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Dim3.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <memory>

namespace Insert {

struct WallChargeGrid {
    amrex::Dim3 domain_lo;
    amrex::Dim3 domain_hi;
    amrex::GpuArray<amrex::Real, 3> prob_lo;
    amrex::GpuArray<amrex::Real, 3> inv_cell_size;
    amrex::Real inv_cell_volume = 0.0;
};

[[nodiscard]] WallChargeGrid MakeWallChargeGrid (amrex::Geometry const& geom);

/** Deposit a charge in C as a trilinear nodal source in C/m^3. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DepositWallChargeToNodes (amrex::Array4<amrex::Real> const& wall_charge,
                          WallChargeGrid const& grid,
                          amrex::ParticleReal const x,
                          amrex::ParticleReal const y,
                          amrex::ParticleReal const z,
                          amrex::ParticleReal const charge) noexcept {
#if defined(WARPX_DIM_3D)
    if (grid.domain_hi.x <= grid.domain_lo.x ||
        grid.domain_hi.y <= grid.domain_lo.y ||
        grid.domain_hi.z <= grid.domain_lo.z) {
        return;
    }

    amrex::ParticleReal const x_node =
        (x - static_cast<amrex::ParticleReal>(grid.prob_lo[0])) *
        static_cast<amrex::ParticleReal>(grid.inv_cell_size[0]);
    amrex::ParticleReal const y_node =
        (y - static_cast<amrex::ParticleReal>(grid.prob_lo[1])) *
        static_cast<amrex::ParticleReal>(grid.inv_cell_size[1]);
    amrex::ParticleReal const z_node =
        (z - static_cast<amrex::ParticleReal>(grid.prob_lo[2])) *
        static_cast<amrex::ParticleReal>(grid.inv_cell_size[2]);
    int i_left;
    int j_left;
    int k_left;
    amrex::ParticleReal wx;
    amrex::ParticleReal wy;
    amrex::ParticleReal wz;
    if (!Math::LinearHatWeightsDiscardOutside(
            x_node, grid.domain_hi.x - grid.domain_lo.x + 1, i_left, wx) ||
        !Math::LinearHatWeightsDiscardOutside(
            y_node, grid.domain_hi.y - grid.domain_lo.y + 1, j_left, wy) ||
        !Math::LinearHatWeightsDiscardOutside(
            z_node, grid.domain_hi.z - grid.domain_lo.z + 1, k_left, wz)) {
        return;
    }

    int const i = grid.domain_lo.x + i_left;
    int const j = grid.domain_lo.y + j_left;
    int const k = grid.domain_lo.z + k_left;

    amrex::ParticleReal const charge_density =
        charge * static_cast<amrex::ParticleReal>(grid.inv_cell_volume);
    amrex::ParticleReal const wx0 = amrex::ParticleReal(1.0) - wx;
    amrex::ParticleReal const wy0 = amrex::ParticleReal(1.0) - wy;
    amrex::ParticleReal const wz0 = amrex::ParticleReal(1.0) - wz;
    // consider changing to multi loop
    amrex::HostDevice::Atomic::Add(
        &wall_charge(i, j, k),
        static_cast<amrex::Real>(charge_density * wx0 * wy0 * wz0));
    amrex::HostDevice::Atomic::Add(
        &wall_charge(i + 1, j, k),
        static_cast<amrex::Real>(charge_density * wx * wy0 * wz0));
    amrex::HostDevice::Atomic::Add(
        &wall_charge(i, j + 1, k),
        static_cast<amrex::Real>(charge_density * wx0 * wy * wz0));
    amrex::HostDevice::Atomic::Add(
        &wall_charge(i + 1, j + 1, k),
        static_cast<amrex::Real>(charge_density * wx * wy * wz0));
    amrex::HostDevice::Atomic::Add(
        &wall_charge(i, j, k + 1),
        static_cast<amrex::Real>(charge_density * wx0 * wy0 * wz));
    amrex::HostDevice::Atomic::Add(
        &wall_charge(i + 1, j, k + 1),
        static_cast<amrex::Real>(charge_density * wx * wy0 * wz));
    amrex::HostDevice::Atomic::Add(
        &wall_charge(i, j + 1, k + 1),
        static_cast<amrex::Real>(charge_density * wx0 * wy * wz));
    amrex::HostDevice::Atomic::Add(
        &wall_charge(i + 1, j + 1, k + 1),
        static_cast<amrex::Real>(charge_density * wx * wy * wz));
#else
    amrex::ignore_unused(wall_charge, grid, x, y, z, charge);
#endif
}

/** Persistent, nodal grid source for charge retained by analytic walls. */
class HallWallCharge {
  public:
    static HallWallCharge& GetInstance ();

    HallWallCharge (HallWallCharge const&) = delete;
    HallWallCharge (HallWallCharge&&) = delete;
    HallWallCharge& operator=(HallWallCharge const&) = delete;
    HallWallCharge& operator=(HallWallCharge&&) = delete;

    /** Define the persistent source on the electrostatic rho layout once. */
    void prepare (ablastr::fields::MultiLevelScalarField const& rho,
                  int max_level);

    /** Add the synchronized persistent source to the Poisson right-hand side.
     */
    void addTo (ablastr::fields::MultiLevelScalarField const& rho,
                int max_level);

    /** Return the per-level storage used by wall-interaction deposition
     * kernels. */
    [[nodiscard]] amrex::MultiFab& get (int lev);

    [[nodiscard]] bool isDefined () const noexcept { return m_defined; }

    /** Record that a deposition kernel has written into the persistent source.
     */
    void
    markDirty () noexcept {
        m_dirty = true;
    }

  private:
    HallWallCharge () = default;
    ~HallWallCharge () = default;

    void synchronize ();
    void assertLayout (ablastr::fields::MultiLevelScalarField const& rho,
                       int max_level) const;

    bool m_defined = false;
    bool m_dirty = false;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_wall_charge;
};

} // namespace Insert

#endif // WARPX_INSERT_HALLWALLCHARGE_H_
