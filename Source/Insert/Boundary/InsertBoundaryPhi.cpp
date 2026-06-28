#include "InsertBoundaryPhi.h"

#include "WarpX.H"

#include "Fields.H"
#include "Insert/Config/WarpXFunctionConfig.h"
#include "Insert/Config/WarpXSimulationConfig.h"
#include "Insert/Utils/InsertUtils.h"
#include "Utils/WarpXConst.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>

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

#if defined(WARPX_DIM_3D)
struct PhiOversetMaskCache
{
    amrex::Vector<std::unique_ptr<amrex::iMultiFab> > masks;
    bool initialized = false;
};
#endif

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

amrex::Vector<std::unique_ptr<amrex::iMultiFab> > const&
BuildPhiOversetMasks (ablastr::fields::MultiLevelScalarField const& phi)
{
    static amrex::Vector<std::unique_ptr<amrex::iMultiFab> > empty_masks;

#if defined(HALL3D) && defined(WARPX_DIM_3D)
    static PhiOversetMaskCache cache;
    WarpX& warpx_instance = WarpX::GetInstance();

    if (!cache.initialized) {
        auto const config = ReadHallAnodeRingConfig(warpx_instance.Geom(0));
        cache.masks.reserve(phi.size());

        for (int lev = 0; lev < static_cast<int>(phi.size()); ++lev) {
            auto* phi_field = phi[lev];
            amrex::BoxArray mask_ba(phi_field->boxArray());
            mask_ba.convert(amrex::IntVect::TheNodeVector());
            auto mask = std::make_unique<amrex::iMultiFab>(
                mask_ba, phi_field->DistributionMap(), 1, 0);
            mask->setVal(1);

            amrex::Box domain = warpx_instance.Geom(lev).Domain();
            domain.surroundingNodes();
            amrex::Real const problo_x = warpx_instance.Geom(lev).ProbLo(0);
            amrex::Real const problo_y = warpx_instance.Geom(lev).ProbLo(1);
            amrex::Real const dx = warpx_instance.Geom(lev).CellSize(0);
            amrex::Real const dy = warpx_instance.Geom(lev).CellSize(1);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(*mask, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Array4<int> const& mask_arr = mask->array(mfi);
                amrex::Array4<amrex::Real> const& phi_arr = phi_field->array(mfi);
                amrex::Box const& box = mfi.tilebox();
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    if (IsHallAnodeRingNode(i, j, k,
                                             domain.smallEnd(2),
                                             domain.smallEnd(0),
                                             domain.smallEnd(1),
                                             problo_x, problo_y, dx, dy, config))
                    {
                        mask_arr(i, j, k) = 0;
                    }
                });
            }

            phi_field->FillBoundary(warpx_instance.Geom(lev).periodicity());
            cache.masks.push_back(std::move(mask));
        }

        cache.initialized = true;
    }

    return cache.masks;
#else
    return empty_masks;
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

} // namespace Insert
