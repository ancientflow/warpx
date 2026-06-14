#pragma once

#include "Insert/Config/WarpXFunctionConfig.h"

#include "Particles/Pusher/GetAndSetPosition.H"

#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

#if defined(IONIZATION_SOURCE_RECORD)

#ifndef IONIZATION_SOURCE_NR
#define IONIZATION_SOURCE_NR 256
#endif

#ifndef IONIZATION_SOURCE_NZ
#define IONIZATION_SOURCE_NZ 256
#endif

#endif

#if defined(IONIZATION_SOURCE_RECORD)

namespace Insert {

struct IonizationSourceRecordView {
    amrex::Real* node_count = nullptr;
    int nr = 0;
    int nz = 0;
    amrex::Real x_center = 0.0;
    amrex::Real y_center = 0.0;
    amrex::Real r_min = 0.0;
    amrex::Real z_min = 0.0;
    amrex::Real inv_dr = 0.0;
    amrex::Real inv_dz = 0.0;
};

void IonizationSourceRecordCollisionSample (amrex::Geometry const& geom,
                                            amrex::Real dt,
                                            amrex::ParticleReal elec_weight,
                                            int max_level);

IonizationSourceRecordView IonizationSourceGetRecordView ();

void IonizationSourceFinalize ();

template <typename ParticleTile>
void
IonizationSourceDepositNewElectrons (ParticleTile const& particle_tile,
                                     long first_new, long num_added,
                                     amrex::ParticleReal event_weight) {
    if (num_added <= 0) {
        return;
    }

    auto const view = IonizationSourceGetRecordView();
    auto const get_position =
        GetParticlePosition<PIdx>(particle_tile, first_new);
    amrex::Real* const node_count = view.node_count;
    int const nr = view.nr;
    int const nz = view.nz;
    amrex::Real const x_center = view.x_center;
    amrex::Real const y_center = view.y_center;
    amrex::Real const r_min = view.r_min;
    amrex::Real const z_min = view.z_min;
    amrex::Real const inv_dr = view.inv_dr;
    amrex::Real const inv_dz = view.inv_dz;
    amrex::Real const weight = static_cast<amrex::Real>(event_weight);

    amrex::ParallelFor(num_added, [=] AMREX_GPU_DEVICE(long ip) noexcept {
        amrex::ParticleReal x, y, z;
        get_position.AsStored(ip, x, y, z);

        amrex::Real const dx = static_cast<amrex::Real>(x) - x_center;
        amrex::Real const dy = static_cast<amrex::Real>(y) - y_center;
        amrex::Real const r = std::sqrt(dx * dx + dy * dy);
        amrex::Real const z_axis = static_cast<amrex::Real>(z);

        amrex::Real const rz = (z_axis - z_min) * inv_dz;
        amrex::Real const rr = (r - r_min) * inv_dr;

        if (rz < 0.0 || rr < 0.0 || rz > static_cast<amrex::Real>(nz) ||
            rr > static_cast<amrex::Real>(nr)) {
            return;
        }

        int iz = static_cast<int>(std::floor(rz));
        int ir = static_cast<int>(std::floor(rr));
        amrex::Real xi = rz - static_cast<amrex::Real>(iz);
        amrex::Real eta = rr - static_cast<amrex::Real>(ir);

        if (iz == nz) {
            iz = nz - 1;
            xi = 1.0;
        }
        if (ir == nr) {
            ir = nr - 1;
            eta = 1.0;
        }

        auto const node_index = [=] AMREX_GPU_DEVICE(int izn,
                                                     int irn) noexcept {
            return izn * (nr + 1) + irn;
        };

        amrex::Gpu::Atomic::AddNoRet(&node_count[node_index(iz, ir)],
                                     weight * (1.0 - xi) * (1.0 - eta));
        amrex::Gpu::Atomic::AddNoRet(&node_count[node_index(iz + 1, ir)],
                                     weight * xi * (1.0 - eta));
        amrex::Gpu::Atomic::AddNoRet(&node_count[node_index(iz, ir + 1)],
                                     weight * (1.0 - xi) * eta);
        amrex::Gpu::Atomic::AddNoRet(&node_count[node_index(iz + 1, ir + 1)],
                                     weight * xi * eta);
    });
}

} // namespace Insert

#endif
