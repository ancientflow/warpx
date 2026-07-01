#include "ZMinWallCharge.h"

#include "WarpX.H"

#include "Insert/Config/WarpXSimulationConfig.h"
#include "Insert/Utils/InsertUtils.h"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"

#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_Math.H>

#include <limits>

namespace {

#ifdef HALL3D

using Insert::BacktraceParticleToZPlane;
using Insert::HallAnodeRingConfig;
using Insert::IsHallAnodeRingHit;
using Insert::ReadHallAnodeRingConfig;

constexpr int zlo_boundary = 4;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
IsZMinAnodeParticle (amrex::ParticleReal const x,
                     amrex::ParticleReal const y,
                     HallAnodeRingConfig const anode_ring) noexcept
{
    return IsHallAnodeRingHit(x, y, anode_ring);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
DepositZMinWallChargeToNodes (amrex::Real* const AMREX_RESTRICT wall_charge,
                               Insert::ZMinWallChargeGrid const grid,
                               amrex::ParticleReal const x,
                               amrex::ParticleReal const y,
                               amrex::ParticleReal const charge) noexcept
{
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

struct ZMinWallChargeDepositFunctor
{
    amrex::ParticleReal const* AMREX_RESTRICT px;
    amrex::ParticleReal const* AMREX_RESTRICT py;
    amrex::ParticleReal const* AMREX_RESTRICT pw;
    amrex::ParticleReal const* AMREX_RESTRICT pz;
    amrex::ParticleReal const* AMREX_RESTRICT pvx;
    amrex::ParticleReal const* AMREX_RESTRICT pvy;
    amrex::ParticleReal const* AMREX_RESTRICT pvz;
    amrex::ParticleReal species_charge;
    Insert::ZMinWallChargeGrid grid;
    HallAnodeRingConfig anode_ring;
    amrex::Real* AMREX_RESTRICT wall_charge;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void operator() (long const ip) const noexcept
    {
        amrex::ParticleReal x_hit, y_hit;
        if (!BacktraceParticleToZPlane(px[ip], py[ip], pz[ip], pvx[ip],
                                       pvy[ip], pvz[ip], grid.zmin, x_hit,
                                       y_hit))
        {
            return;
        }

        if (IsZMinAnodeParticle(x_hit, y_hit, anode_ring)) {
            return;
        }

        DepositZMinWallChargeToNodes(wall_charge, grid, x_hit, y_hit,
                                     species_charge * pw[ip]);
    }
};

void
DepositSpeciesZMinWallCharge (WarpXParticleContainer::Base* const particles,
                               amrex::ParticleReal const species_charge,
                               Insert::ZMinWallChargeGrid const grid,
                               HallAnodeRingConfig const anode_ring,
                               amrex::Real* const AMREX_RESTRICT wall_charge)
{
    if (particles == nullptr || !particles->isDefined()) {
        return;
    }

    for (auto pti = WarpXParIter(*particles, 0); pti.isValid(); ++pti) {
        auto& arr = pti.GetStructOfArrays().GetRealData();
        auto px = arr[PIdx::x].dataPtr();
        auto py = arr[PIdx::y].dataPtr();
        auto pw = arr[PIdx::w].dataPtr();
        auto pz = arr[PIdx::z].dataPtr();
        auto pvx = arr[PIdx::ux].dataPtr();
        auto pvy = arr[PIdx::uy].dataPtr();
        auto pvz = arr[PIdx::uz].dataPtr();
        int const np = pti.numParticles();

        ZMinWallChargeDepositFunctor const deposit{
            px,   py,         pw,         pz, pvx, pvy, pvz, species_charge,
            grid, anode_ring, wall_charge};
        amrex::ParallelFor(np, deposit);
    }
}

#endif

} // namespace

namespace Insert {

ZMinWallChargeGrid
MakeZMinWallChargeGrid (amrex::Geometry const& geom)
{
#ifdef HALL3D
    amrex::Box domain = geom.Domain();
    domain.surroundingNodes();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(domain.length(0) >= 2 &&
                                         domain.length(1) >= 2,
                                     "zmin wall charge deposition requires at "
                                     "least two zmin face nodes in x and y.");

    return ZMinWallChargeGrid{domain.length(0),
                              domain.length(1),
                              geom.ProbLo(0),
                              geom.ProbLo(1),
                              amrex::Real(1.0) / geom.CellSize(0),
                              amrex::Real(1.0) / geom.CellSize(1),
                              geom.CellSize(0),
                              geom.CellSize(1),
                              geom.ProbLo(2)};
#else
    (void)geom;
    return ZMinWallChargeGrid{};
#endif
}

long
ZMinWallChargeSize (ZMinWallChargeGrid const& grid)
{
    return static_cast<long>(grid.nx) * static_cast<long>(grid.ny);
}

amrex::Gpu::DeviceVector<amrex::Real>
DepositZMinWallCharge (WarpX& warpx_instance, ZMinWallChargeGrid const& grid)
{
#ifdef HALL3D
    long const data_size = ZMinWallChargeSize(grid);

    amrex::Gpu::DeviceVector<amrex::Real> device_wall_charge(data_size,
                                                             amrex::Real(0.0));
    amrex::Real* const wall_charge_ptr = device_wall_charge.dataPtr();

    auto& mypc = warpx_instance.GetPartContainer();
    auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();
    HallAnodeRingConfig const anode_ring =
        ReadHallAnodeRingConfig(warpx_instance.Geom(0));

    auto* elec_zmin = mybpc.getParticleBufferPointer("electrons", zlo_boundary);
    auto& elec_pc = mypc.GetParticleContainerFromName("electrons");
    DepositSpeciesZMinWallCharge(elec_zmin, elec_pc.getCharge(), grid,
                                 anode_ring, wall_charge_ptr);

    auto* xe_ion_zmin = mybpc.getParticleBufferPointer("xe_ions", zlo_boundary);
    auto& xe_ion_pc = mypc.GetParticleContainerFromName("xe_ions");
    DepositSpeciesZMinWallCharge(xe_ion_zmin, xe_ion_pc.getCharge(), grid,
                                 anode_ring, wall_charge_ptr);

    return device_wall_charge;
#else
    (void)warpx_instance;
    (void)grid;
    return amrex::Gpu::DeviceVector<amrex::Real>{};
#endif
}

} // namespace Insert
