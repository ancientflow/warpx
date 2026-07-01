#pragma once

#include <AMReX_GpuContainers.H>
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

ZMinWallChargeGrid MakeZMinWallChargeGrid (amrex::Geometry const& geom);

long ZMinWallChargeSize (ZMinWallChargeGrid const& grid);

amrex::Gpu::DeviceVector<amrex::Real>
DepositZMinWallCharge (WarpX& warpx_instance, ZMinWallChargeGrid const& grid);

} // namespace Insert
