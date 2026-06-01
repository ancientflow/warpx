#include "InsertRuntimeDiagnostics.h"

#include "WarpX.H"

#include "Insert/Config/WarpXFunctionConfig.h"
#include "Insert/Config/WarpXSimulationConfig.h"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <fstream>

namespace Insert {

void
ParticleNumber () {
#ifdef NUMP
    WarpX& warpx_instance = WarpX::GetInstance();
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();

    auto species_names = mypc.GetSpeciesNames();

    for (const auto& name : species_names) {
        auto& pc = mypc.GetParticleContainerFromName(name);
        amrex::Print() << name << " number: " << pc.TotalNumberOfParticles()
                       << "\n";
    }
#endif
}

void
ShowAndWriteIonzationNum (amrex::Vector<int> num) {
#ifdef COLLISION_RECORD
    amrex::Print() << "produce macro electron: " << num[0]
                   << "\nproduce macro xe ions: " << num[1] << "\n";

    static bool ifinit = false;
    std::fstream fileout;
    if (!ifinit) {
        fileout.open("collision_record.dat", std::ios::out);
        ifinit = true;
    } else {
        fileout.open("collision_record.dat", std::ios::app);
    }
    fileout << num[0] << "\t" << num[1] << "\n";
    fileout.close();
#else
    amrex::ignore_unused(num);
#endif
}

void
AnodeCurrentCalc () {
#ifdef HALL3D
    static int times = 0;
    static const int gap = 10;
    times++;
    static bool ifinit = false;
    if (times == gap) {
        std::fstream fileout;
        if (!ifinit) {
            fileout.open("anode_current.dat", std::ios::out);
            fileout << "time\t"
                       "zmin_electron\tzmin_ion\tanode_electron\t"
                       "anode_electron_cut\n";
            ifinit = true;
        } else {
            fileout.open("anode_current.dat", std::ios::app);
        }
        times = 0;

        WarpX& warpx_instance = WarpX::GetInstance();
        auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();

        auto& elec_zmin = mybpc.getParticleBuffer("electrons", 4);
        auto& xe_ion_zmin = mybpc.getParticleBuffer("xe_ions", 4);

        const int data_size = 4;
        amrex::Gpu::DeviceVector<amrex::ParticleReal> device_charge(data_size,
                                                                    0);
        amrex::Vector<amrex::ParticleReal> host_charge(data_size, 0);

        amrex::ParticleReal* device_ptr = device_charge.dataPtr();

        amrex::ParticleReal l_factor, half_L, L, rn1, rn2;
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("l_factor", l_factor);
        pp_mc.get("L", L);

        half_L = L / l_factor / 2;
        rn1 = 0.021 / 2 / l_factor;
        rn2 = 0.031 / 2 / l_factor;

        for (auto pti = WarpXParIter(elec_zmin, 0); pti.isValid(); ++pti) {
            auto& arr = pti.GetStructOfArrays().GetRealData();
            auto px = arr[PIdx::x].dataPtr();
            auto py = arr[PIdx::y].dataPtr();
            auto pw = arr[PIdx::w].dataPtr();
            auto pz = arr[PIdx::z].dataPtr();
            auto pvx = arr[PIdx::ux].dataPtr();
            auto pvy = arr[PIdx::uy].dataPtr();
            auto pvz = arr[PIdx::uz].dataPtr();
            int np = pti.numParticles();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) {
                amrex::ParticleReal x = px[ip], y = py[ip], w = pw[ip],
                                    z = pz[ip], vx = pvx[ip], vy = pvy[ip],
                                    vz = pvz[ip];
                amrex::ParticleReal dt = z / vz;
                x -= half_L;
                y -= half_L;
                amrex::ParticleReal r = std::sqrt(x * x + y * y);
                if (r >= rn1 && r <= rn2) {
                    amrex::Gpu::Atomic::Add(&device_ptr[data_size - 2], w);
                }
                x -= dt * vx;
                y -= dt * vy;
                r = std::sqrt(x * x + y * y);
                if (r >= rn1 && r <= rn2) {
                    amrex::Gpu::Atomic::Add(&device_ptr[data_size - 1], w);
                }
                amrex::Gpu::Atomic::Add(&device_ptr[0], w);
            });
        }

        for (auto pti = WarpXParIter(xe_ion_zmin, 0); pti.isValid(); ++pti) {
            auto& arr = pti.GetStructOfArrays().GetRealData();
            auto px = arr[PIdx::x].dataPtr();
            auto py = arr[PIdx::y].dataPtr();
            auto pw = arr[PIdx::w].dataPtr();
            auto pz = arr[PIdx::z].dataPtr();
            auto pvx = arr[PIdx::ux].dataPtr();
            auto pvy = arr[PIdx::uy].dataPtr();
            auto pvz = arr[PIdx::uz].dataPtr();
            int np = pti.numParticles();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) {
                amrex::ParticleReal x = px[ip], y = py[ip], w = pw[ip],
                                    z = pz[ip], vx = pvx[ip], vy = pvy[ip],
                                    vz = pvz[ip];
                amrex::ParticleReal dt = z / vz;
                x -= half_L;
                y -= half_L;
                x -= dt * vx;
                y -= dt * vy;
                amrex::ParticleReal r = std::sqrt(x * x + y * y);
                if (r >= rn1 && r <= rn2) {
                    amrex::Gpu::Atomic::Add(&device_ptr[1], w);
                }
            });
        }

        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_charge.begin(),
                         device_charge.end(), host_charge.begin());
        fileout << warpx_instance.gett_new(0);
        for (int i = 0; i < host_charge.size(); i++) {
            fileout << "\t" << host_charge[i];
        }
        fileout << "\n";
        fileout.close();
        mybpc.clearParticles();
    }
#endif
}

} // namespace Insert
