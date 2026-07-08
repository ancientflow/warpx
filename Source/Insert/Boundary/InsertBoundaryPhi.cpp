#include "InsertBoundaryPhi.h"

#include "WarpX.H"

#include "Fields.H"
#include "Insert/Boundary/ZMinWallCharge.h"
#include "Insert/Config/WarpXFunctionConfig.h"
#include "Insert/Config/WarpXSimulationConfig.h"
#include "Insert/Diagnostics/InsertRuntimeDiagnostics.h"
#include "Insert/Fields/SpectralBoundarySchur.h"
#include "Insert/Utils/InsertUtils.h"
#include "Utils/WarpXConst.H"

#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <iostream>

namespace {

#ifdef HALL3D
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
IsHallAnodeRingNode (
    int i, int j, int k,
    int zlo, int xlo, int ylo,
    amrex::Real problo_x, amrex::Real problo_y,
    amrex::Real dx, amrex::Real dy,
    Insert::HallAnodeRingConfig const config)
{
    amrex::Real const x = problo_x + (i - xlo) * dx;
    amrex::Real const y = problo_y + (j - ylo) * dy;
    return k == zlo && Insert::IsHallAnodeRingHit(x, y, config);
}

#endif

} // namespace

namespace Insert {

void
VoltageAdjustment ()
{
#ifdef BENCHMARK_2D
    WarpX& warpx_instance = WarpX::GetInstance();
    auto phi_field = warpx_instance.m_fields.get(FieldType::phi_fp, 0);
    amrex::Real phisum = 0;

    static amrex::MultiFab mf(phi_field->boxArray(), phi_field->DistributionMap(),
                              phi_field->nComp(), phi_field->nGrow());
    static bool DotInit = false;
    if (!DotInit) {
        mf.setVal(0);
        for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            const amrex::Box& box = mfi.tilebox();
            amrex::Array4<amrex::Real> const& arr = mf.array(mfi);
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j) {
                if (i == 491) {
                    arr(i, j, 0) = 0.48;
                } else if (i == 492) {
                    arr(i, j, 0) = 0.52;
                }
            });
        }
        DotInit = true;
    }
    phisum = amrex::MultiFab::Dot(mf, 0, *phi_field, 0, 1, 0);

    phisum /= 256;

    std::cout << "voltage adjustment: " << phisum << std::endl;

    for (amrex::MFIter mfi(*phi_field, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& box = mfi.tilebox();

        amrex::Array4<amrex::Real> const& phi = phi_field->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j) {
            phi(i, j, 0) -= phisum * i / 491.52;
        });
    }
#endif
}

void
AnodeVoltage ()
{
#ifdef HALL3D
    WarpX& warpx_instance = WarpX::GetInstance();
    auto phi_field =
        warpx_instance.m_fields.get(warpx::fields::FieldType::phi_fp, 0);
    auto const config = ReadHallAnodeRingConfig(warpx_instance.Geom(0));
    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();
    amrex::Real const problo_x = warpx_instance.Geom(0).ProbLo(0);
    amrex::Real const problo_y = warpx_instance.Geom(0).ProbLo(1);
    amrex::Real const dx = warpx_instance.Geom(0).CellSize(0);
    amrex::Real const dy = warpx_instance.Geom(0).CellSize(1);
    for (amrex::MFIter mfi(*phi_field, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        amrex::Array4<amrex::Real> const& phi = phi_field->array(mfi);
        if (!domain.strictly_contains(box)) {
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (IsHallAnodeRingNode(i, j, k,
                                         domain.smallEnd(2),
                                         domain.smallEnd(0),
                                         domain.smallEnd(1),
                                         problo_x, problo_y, dx, dy, config))
                {
                    phi(i, j, k) = config.voltage;
                }
            });
        }
    }
    phi_field->FillBoundary(warpx_instance.Geom(0).periodicity());
#endif
}

void
DirichletPhiGuardSet ()
{
#ifndef WAVE1D
    using namespace amrex::literals;

    WarpX& warpx_instance = WarpX::GetInstance();
    auto rho = warpx_instance.m_fields.get(warpx::fields::FieldType::rho_fp, 0);
    auto phi = warpx_instance.m_fields.get(warpx::fields::FieldType::phi_fp, 0);
    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();
    const amrex::Real* dx_host = warpx_instance.Geom(0).CellSize();
    amrex::Gpu::DeviceVector<amrex::Real> dx_device(3);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, dx_host, dx_host + 3,
                     dx_device.begin());
    amrex::Real* dx = dx_device.dataPtr();
#ifdef HALL3D
    auto const anode_config = ReadHallAnodeRingConfig(warpx_instance.Geom(0));
    amrex::Real const problo_x = warpx_instance.Geom(0).ProbLo(0);
    amrex::Real const problo_y = warpx_instance.Geom(0).ProbLo(1);
#endif

    for (amrex::MFIter mfi(*phi, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& box = mfi.validbox();
        const auto& phi_arr = phi->array(mfi);
        const auto& rho_arr = rho->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (i == domain.smallEnd(0)) {
                phi_arr(i - 1, j, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i + 1, j, k) -
                    rho_arr(i, j, k) * dx[0] * dx[0] / PhysConst::epsilon_0;
            }
            if (i == domain.bigEnd(0)) {
                phi_arr(i + 1, j, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i - 1, j, k) -
                    rho_arr(i, j, k) * dx[0] * dx[0] / PhysConst::epsilon_0;
            }
#if defined(WARPX_DIM_XZ)
            if (j == domain.smallEnd(1)) {
                phi_arr(i, j - 1, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j + 1, k) -
                    rho_arr(i, j, k) * dx[1] * dx[1] / PhysConst::epsilon_0;
            }
            if (j == domain.bigEnd(1)) {
                phi_arr(i, j + 1, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j - 1, k) -
                    rho_arr(i, j, k) * dx[1] * dx[1] / PhysConst::epsilon_0;
            }
#endif
#if defined(WARPX_DIM_3D)
#if defined(HALL3D)
            if (IsHallAnodeRingNode(i, j, k,
                                    domain.smallEnd(2),
                                    domain.smallEnd(0),
                                    domain.smallEnd(1),
                                    problo_x, problo_y, dx[0], dx[1], anode_config))
#else
            if (k == domain.smallEnd(2))
#endif
            {
                phi_arr(i, j, k - 1) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j, k + 1) -
                    rho_arr(i, j, k) * dx[2] * dx[2] / PhysConst::epsilon_0;
            }
            if (k == domain.bigEnd(2)) {
                phi_arr(i, j, k + 1) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j, k - 1) -
                    rho_arr(i, j, k) * dx[2] * dx[2] / PhysConst::epsilon_0;
            }
#endif
        });
    }
#endif
}

void
HallThrusterPhiGuardSet ()
{
#ifdef HALL3D
    WarpX& warpx_instance = WarpX::GetInstance();
    auto phi = warpx_instance.m_fields.get(warpx::fields::FieldType::phi_fp, 0);

    auto const zmin_bc = WarpX::field_boundary_lo[WARPX_ZINDEX];
    bool const is_dirichlet = zmin_bc == FieldBoundaryType::PEC;
    bool const is_neumann = zmin_bc == FieldBoundaryType::Neumann;
    bool const use_schur = SpectralBoundarySchur::Enabled();
    if (!is_dirichlet && !is_neumann && !use_schur) {
        return;
    }

    amrex::Geometry const& geom = warpx_instance.Geom(0);
    amrex::Box domain = geom.Domain();
    domain.surroundingNodes();
    int const xlo = domain.smallEnd(0);
    int const ylo = domain.smallEnd(1);
    int const zlo = domain.smallEnd(2);

    amrex::Real const problo_x = geom.ProbLo(0);
    amrex::Real const problo_y = geom.ProbLo(1);
    amrex::Real const dx = geom.CellSize(0);
    amrex::Real const dy = geom.CellSize(1);
    amrex::Real const dz = geom.CellSize(2);
    auto const anode_config = ReadHallAnodeRingConfig(geom);

    amrex::Array4<amrex::Real const> wall_charge_density;
    bool const use_wall_charge_guard = is_neumann || use_schur;
    if (use_wall_charge_guard) {
        if (amrex::ParallelDescriptor::NProcs() != 1) {
            amrex::Abort("Hall zmin phi guard correction requires one MPI rank");
        }

        ZMinWallChargeGrid const grid = MakeZMinWallChargeGrid(geom);
        InitializeAccumulatedZMinWallChargeDensity(grid);
        for (amrex::MFIter mfi(*g_accumulated_wall_charge_density);
             mfi.isValid(); ++mfi)
        {
            wall_charge_density =
                g_accumulated_wall_charge_density->const_array(mfi);
        }
    }

    amrex::Real const inv_epsilon0 =
        amrex::Real(1.0) / PhysConst::epsilon_0;

    for (amrex::MFIter mfi(*phi, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi)
    {
        amrex::Box const& box = mfi.validbox();
        if (box.smallEnd(2) > zlo || box.bigEnd(2) < zlo) {
            continue;
        }

        amrex::Array4<amrex::Real> const& phi_arr = phi->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (k != zlo) {
                return;
            }

            if (IsHallAnodeRingNode(i, j, k, zlo, xlo, ylo, problo_x,
                                    problo_y, dx, dy, anode_config))
            {
                phi_arr(i, j, k - 1) =
                    use_schur ? amrex::Real(0.0) : anode_config.voltage;
                return;
            }

            if (use_wall_charge_guard) {
                amrex::Real const sigma_s =
                    wall_charge_density(i - xlo, j - ylo, 0);
                // zmin outward normal n=-z, so dphi/dz = -sigma_s/epsilon0.
                phi_arr(i, j, k - 1) =
                    phi_arr(i, j, k + 1) +
                    amrex::Real(2.0) * dz * sigma_s * inv_epsilon0;
            }
        });
    }
#endif
}

} // namespace Insert
