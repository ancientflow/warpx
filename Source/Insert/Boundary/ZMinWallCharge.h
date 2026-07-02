#pragma once

#include <AMReX.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

class WarpX;

namespace amrex
{
class Geometry;
}

namespace Insert {

struct ZMinWallChargeGrid
{
    int nx = 0;
    int ny = 0;
    amrex::Real problo_x = amrex::Real(0.0);
    amrex::Real problo_y = amrex::Real(0.0);
    amrex::Real inv_dx = amrex::Real(0.0);
    amrex::Real inv_dy = amrex::Real(0.0);
    amrex::Real dx = amrex::Real(0.0);
    amrex::Real dy = amrex::Real(0.0);
    amrex::Real zmin = amrex::Real(0.0);
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
DepositZMinWallChargeToNodes (amrex::Real* const AMREX_RESTRICT wall_charge,
                              ZMinWallChargeGrid const grid,
                              amrex::ParticleReal const x,
                              amrex::ParticleReal const y,
                              amrex::ParticleReal const charge) noexcept
{
    if (wall_charge == nullptr || grid.nx < 2 || grid.ny < 2) {
        return;
    }

    amrex::Real const x_node = (x - grid.problo_x) * grid.inv_dx;
    amrex::Real const y_node = (y - grid.problo_y) * grid.inv_dy;
    if (x_node < amrex::Real(0.0) ||
        x_node > static_cast<amrex::Real>(grid.nx - 1) ||
        y_node < amrex::Real(0.0) ||
        y_node > static_cast<amrex::Real>(grid.ny - 1))
    {
        return;
    }

    int i_left = static_cast<int>(amrex::Math::floor(x_node));
    int j_left = static_cast<int>(amrex::Math::floor(y_node));
    amrex::Real wx_right = x_node - static_cast<amrex::Real>(i_left);
    amrex::Real wy_right = y_node - static_cast<amrex::Real>(j_left);

    if (i_left >= grid.nx - 1) {
        i_left = grid.nx - 2;
        wx_right = amrex::Real(1.0);
    }
    if (j_left >= grid.ny - 1) {
        j_left = grid.ny - 2;
        wy_right = amrex::Real(1.0);
    }

    amrex::Real const wx_left = amrex::Real(1.0) - wx_right;
    amrex::Real const wy_left = amrex::Real(1.0) - wy_right;
    long const offset_ll =
        static_cast<long>(j_left) * static_cast<long>(grid.nx) +
        static_cast<long>(i_left);
    long const offset_lr = offset_ll + 1;
    long const offset_ul = offset_ll + static_cast<long>(grid.nx);
    long const offset_ur = offset_ul + 1;

    amrex::HostDevice::Atomic::Add(&wall_charge[offset_ll],
                                   charge * wx_left * wy_left);
    amrex::HostDevice::Atomic::Add(&wall_charge[offset_lr],
                                   charge * wx_right * wy_left);
    amrex::HostDevice::Atomic::Add(&wall_charge[offset_ul],
                                   charge * wx_left * wy_right);
    amrex::HostDevice::Atomic::Add(&wall_charge[offset_ur],
                                   charge * wx_right * wy_right);
}

ZMinWallChargeGrid MakeZMinWallChargeGrid (amrex::Geometry const& geom);

long ZMinWallChargeSize (ZMinWallChargeGrid const& grid);

amrex::Gpu::DeviceVector<amrex::Real>
DepositZMinWallCharge (WarpX& warpx_instance, ZMinWallChargeGrid const& grid);

} // namespace Insert
