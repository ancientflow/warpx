#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#include "Fields.H"
#ifdef WARPX_USE_FFT
#ifdef WARPX_DIM_RZ
#include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#else
#include "FieldSolver/SpectralSolver/SpectralSolver.H"
#endif
#endif
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include "BackgroundCoupledDensity.h"
#include "WarpXFunctionConfig.h"
#include "WarpXSimulationConfig.h"
#include <algorithm>
#include <array>
#include <ostream>
#include <vector>



#ifdef NUMP
// 粒子数监测
void
ParticleNumber () {
    WarpX& warpx_instance = WarpX::GetInstance();
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();

    auto species_names = mypc.GetSpeciesNames();

    for (const auto& name : species_names) {
        auto& pc = mypc.GetParticleContainerFromName(name);
        amrex::Print() << name << " number: " << pc.TotalNumberOfParticles()
                       << "\n";
    }
}
#endif

#ifdef EXAMINE
void
PhiExamine () {
    WarpX& warpx_instance = WarpX::GetInstance();
    auto phi = warpx_instance.m_fields.get(FieldType::phi_fp, 0);
    auto step = warpx_instance.getistep(0);
    if (step != 1)
        return;
    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();

    auto lo = domain.smallEnd();
    auto hi = domain.bigEnd();
    const int sizex = hi[0] - lo[0] + 1, sizey = hi[1] - lo[1] + 1,
              sizez = hi[2] - lo[2];
    const int datasize = sizex * sizey * sizez;

    Vector<Real> phi3d(datasize);
    amrex::Gpu::DeviceVector<Real> phi3d_device(datasize);
    Real* pphi = phi3d_device.dataPtr();

    for (MFIter mfi(*phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();

        Array4<Real> const& phiarr = phi->array(mfi);
        if (!domain.strictly_contains(box)) {
            ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k == domain.smallEnd(2)) {
                    pphi[i * sizey * sizez + j * sizez + k] = phiarr(i, j, k);
                }
            });
        }
    }

    amrex::Gpu::copy(amrex::Gpu::deviceToHost, phi3d_device.begin(),
                     phi3d_device.end(), phi3d.begin());
    std::fstream fileout("phiexam.dat", std::ios::out);

    for (int i = 0; i < sizex; i++) {
        for (int j = 0; j < sizey; j++) {
            for (int k = 0; k < sizez; k++) {
                fileout << i << "\t" << j << "\t" << k << "\t"
                        << phi3d[i * sizey * sizez + j * sizez + k] << "\n";
            }
        }
    }
    fileout.close();
}

void
XeRhoExamine () {
    WarpX& warpx_instance = WarpX::GetInstance();
    auto rho = warpx_instance.m_fields.get(FieldType::rho_fp, 0);
    MultiFab xe_rho(rho->boxArray(), rho->DistributionMap(), rho->nComp(),
                    rho->nGrow());
    xe_rho.setVal(0.0);
    auto& xe_pc = mypc->GetParticleContainer(mypc->getSpeciesID("xe_netural"));
    AtomDepostiAPI(xe_pc, xe_rho);

    auto step = warpx_instance.getistep(0);
    if (step != 1)
        return;
    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();

    auto lo = domain.smallEnd();
    auto hi = domain.bigEnd();
    const int sizex = hi[0] - lo[0] + 1, sizey = hi[1] - lo[1] + 1,
              sizez = hi[2] - lo[2];
    const int datasize = sizex * sizey * sizez;

    Vector<Real> phi3d(datasize);
    amrex::Gpu::DeviceVector<Real> phi3d_device(datasize);
    Real* pphi = phi3d_device.dataPtr();

    for (MFIter mfi(*phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();

        Array4<Real> const& phiarr = phi->array(mfi);
        if (!domain.strictly_contains(box)) {
            ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k == domain.smallEnd(2)) {
                    pphi[i * sizey * sizez + j * sizez + k] = phiarr(i, j, k);
                }
            });
        }
    }

    amrex::Gpu::copy(amrex::Gpu::deviceToHost, phi3d_device.begin(),
                     phi3d_device.end(), phi3d.begin());
    std::fstream fileout("phiexam.dat", std::ios::out);

    for (int i = 0; i < sizex; i++) {
        for (int j = 0; j < sizey; j++) {
            for (int k = 0; k < sizez; k++) {
                fileout << i << "\t" << j << "\t" << k << "\t"
                        << phi3d[i * sizey * sizez + j * sizez + k] << "\n";
            }
        }
    }
    fileout.close();
}
#endif

#ifdef COLLISION_RECORD
void
ShowAndWriteIonzationNum (amrex::Vector<int> num) {
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
}
#endif

#ifndef PUSH_GAP
/**
 * 初始化各粒子推进步长
 */
void
PushGapInit () {
    const ParmParse pp_particles("particles");
    std::vector<std::string> species_names;
    pp_particles.queryarr("species_names", species_names);

    int gap;
    for (auto i : species_names) {
        ParmParse pp_species(i);
        pp_species.get("ndt", gap);
        AMREX_ASSERT_WITH_MESSAGE(gap > 0, "The gap should larger than 1.");
        push_control.push_back(std::make_pair(gap, 0));
    }
}
#endif

void
inline DirichletPhiGuardSet () {
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

    for (amrex::MFIter mfi(*phi, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& box = mfi.validbox();
        const auto& phi_arr = phi->array(mfi);
        const auto& rho_arr = rho->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // x direction
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
            // y direction
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
            // z direction
            if (k == domain.smallEnd(2)) {
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
}

#ifdef MCC_DENSITY
amrex::Vector<BackgroundCoupledDensity> global_background_density;
void
inline GlobalBackgroundDensityInit () {
    amrex::ParmParse const pp_coll("collisions");
    amrex::Vector<std::string> species_names;
    pp_coll.queryarr("background_species", species_names);

    global_background_density.resize(species_names.size());
    for (int i = 0; i < species_names.size(); i++) {
        global_background_density[i].m_ground_species = species_names[i];
        global_background_density[i].backgroundDensityInit();
    }
}

inline void
GlobalBackgroundDensityUpdate (const int step, const bool& if_split) {
    if (global_background_density.empty()) {
        return;
    }
    amrex::ParmParse pp_mc("my_constants");
    amrex::ParticleReal elec_weight;
    pp_mc.getWithParser("elec_weight", elec_weight);

    WarpX& warpx_instance = WarpX::GetInstance();
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();
    const bool if_update_sort = warpx_instance.sort_intervals.contains(step);
    for (auto& density : global_background_density) {
        const int ndt =
            mypc.GetParticleContainerFromName(density.m_ground_species)
                .get_ndt();
        const bool if_update_push =
            ((if_split && ((step % ndt == 0) || (step % ndt == 1))) ||
             (!if_split && (step % ndt == 1)));
        if (if_update_push || if_update_sort) {
            density.backgroundDensityUpdate(mypc, elec_weight);
        }
    }
}

inline void
GlobalBackgroundDensityClean (const int step, const bool& if_split) {
    if (global_background_density.empty()) {
        return;
    }
    amrex::ParmParse const pp_mc("my_constants");
    amrex::ParticleReal elec_weight;
    pp_mc.getWithParser("elec_weight", elec_weight);

    WarpX& warpx_instance = WarpX::GetInstance();
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();
    /**
     * deleteInvalidParticles在处理边界条件之后被调用，push会在当前
     * 时步的碰撞之后进行，下一时步将进行背景密度的重新计算，因此这一
     * 时步的所有碰撞完成之后，进行粒子清理
     */
    for (auto& density : global_background_density) {
        int const ndt =
            mypc.GetParticleContainerFromName(density.m_ground_species)
                .get_ndt();
        const bool if_clean = ((if_split && ((step + 1) % ndt == 0)) ||
                               (!if_split && (step % ndt == 0)));
        if (if_clean) {
            density.backgroundSpeciesClean(mypc);
        }
    }
}
#endif

#define WAVE1D
#ifdef WAVE1D
inline void
InitDisturbance () {
    WarpX& warpx_instance = WarpX::GetInstance();
    MultiParticleContainer const& mypc = warpx_instance.GetPartContainer();

    amrex::ParmParse const pp_mc("my_constants");
    double n0, l, mode = 1, Ti, Te;
    int ncell, nppc;
    pp_mc.getWithParser("n0", n0);
    pp_mc.getWithParser("l", l);
    pp_mc.query("mode", mode);
    pp_mc.getWithParser("ncell", ncell);
    pp_mc.get("nppc", nppc);
    pp_mc.get("Te", Te);
    pp_mc.get("Ti", Ti);

    double ratio = 0.01, nd = n0 * ratio, weight = nd * l / nppc, eVToK = 11605;

    int np = ncell * nppc;

    amrex::Vector<amrex::ParticleReal> pz(np, 0), px(np), py(np, 0), vx(np),
        vy(np), vz(np, 0), pw(np, weight);
    for (int i = 0; i < np; i++) {
    }
}
#endif

void
EmbededBoundaryParticleReflect () {
    WarpX& warpx_instance = WarpX::GetInstance();
}