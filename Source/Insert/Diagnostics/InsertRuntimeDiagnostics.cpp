#include "InsertRuntimeDiagnostics.h"

#include "WarpX.H"

#include "Insert/Config/WarpXFunctionConfig.h"
#include "Insert/Config/WarpXSimulationConfig.h"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Utils/WarpXConst.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

using namespace amrex::literals;

namespace {

constexpr int xlo_boundary = 0;
constexpr int xhi_boundary = 1;
constexpr int ylo_boundary = 2;
constexpr int yhi_boundary = 3;
constexpr int zlo_boundary = 4;
constexpr int zhi_boundary = 5;
constexpr int outlet_boundaries[] = {xlo_boundary, xhi_boundary, ylo_boundary,
                                     yhi_boundary, zhi_boundary};

int
HallDiagInterval () {
    int gap = 10;
    amrex::ParmParse pp_mc("my_constants");
    pp_mc.query("hall_diag_interval", gap);
    return std::max(gap, 1);
}

amrex::Real
SampleDt (amrex::Real const time, amrex::Real& last_sample_time,
          int const gap) {
    amrex::Real sample_dt = time - last_sample_time;
    if (sample_dt <= 0.0_rt) {
        sample_dt =
            static_cast<amrex::Real>(gap) * WarpX::GetInstance().getdt(0);
    }
    last_sample_time = time;
    return sample_dt;
}

bool
DiagEnabled (const char* const name) {
    bool enabled = false;
    amrex::ParmParse pp_mc("my_constants");
    pp_mc.query(name, enabled);
    return enabled;
}

} // namespace

namespace Insert {

void
ParticleNumber () {
#ifdef NUMP
    static bool const diag_enabled = DiagEnabled("particle_number_diag");
    if (!diag_enabled) { return; }

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
    static bool const diag_enabled = DiagEnabled("collision_record_diag");
    if (!diag_enabled) { return; }

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
    static bool const diag_enabled = DiagEnabled("anode_current_diag");
    if (!diag_enabled) { return; }

    static int times = 0;
    static int gap = 10;
    times++;

    static bool ifinit = false;
    static std::string anode_current_path = "anode_current.dat";

    if (!ifinit) {
        gap = HallDiagInterval();
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.query("anode_current_path", anode_current_path);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::fstream anode_file(anode_current_path, std::ios::out);
            anode_file << "time\t"
                          "zmin_electron\tzmin_ion\tanode_electron\t"
                          "anode_electron_cut\n";
        }

        ifinit = true;
    }

    if (times == gap) {
        times = 0;
        WarpX& warpx_instance = WarpX::GetInstance();
        auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();

        auto* elec_zmin =
            mybpc.getParticleBufferPointer("electrons", zlo_boundary);
        auto* xe_ion_zmin =
            mybpc.getParticleBufferPointer("xe_ions", zlo_boundary);

        constexpr int data_size = 4;
        amrex::Gpu::DeviceVector<amrex::Real> device_charge(data_size, 0.0_rt);
        amrex::Vector<amrex::Real> host_charge(data_size, 0.0_rt);
        amrex::Real* device_ptr = device_charge.dataPtr();

        amrex::ParticleReal l_factor, half_L, L, rn1, rn2;
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("l_factor", l_factor);
        pp_mc.get("L", L);

        half_L = L / l_factor / 2;
        rn1 = 0.021 / 2 / l_factor;
        rn2 = 0.031 / 2 / l_factor;

        if (elec_zmin != nullptr && elec_zmin->isDefined()) {
            for (auto pti = WarpXParIter(*elec_zmin, 0); pti.isValid(); ++pti) {
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
        }

        if (xe_ion_zmin != nullptr && xe_ion_zmin->isDefined()) {
            for (auto pti = WarpXParIter(*xe_ion_zmin, 0); pti.isValid();
                 ++pti) {
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
        }

        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_charge.begin(),
                         device_charge.end(), host_charge.begin());
        amrex::ParallelDescriptor::ReduceRealSum(
            host_charge.data(), static_cast<int>(host_charge.size()),
            amrex::ParallelDescriptor::IOProcessorNumber());

        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::fstream anode_file(anode_current_path, std::ios::app);
            anode_file << warpx_instance.gett_new(0);
            for (auto const value : host_charge) {
                anode_file << "\t" << value;
            }
            anode_file << "\n";
        }
    }
#endif
}

void
ThrustCalc () {
#ifdef HALL3D
    static bool const diag_enabled = DiagEnabled("thrust_diag");
    if (!diag_enabled) { return; }

    static int times = 0;
    static int gap = 10;
    times++;

    static bool ifinit = false;
    static std::string thrust_diag_path = "thrust.dat";
    static amrex::Real last_sample_time = 0.0_rt;

    if (!ifinit) {
        gap = HallDiagInterval();
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.query("thrust_diag_path", thrust_diag_path);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::fstream thrust_file(thrust_diag_path, std::ios::out);
            thrust_file << "time\tdt\tion_weight\tthrust_N\t"
                           "axial_momentum_kg_m_per_s\n";
        }

        ifinit = true;
    }

    if (times == gap) {
        times = 0;
        WarpX& warpx_instance = WarpX::GetInstance();
        auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();
        constexpr int data_size = 2;
        amrex::Gpu::DeviceVector<amrex::Real> device_data(data_size, 0.0_rt);
        amrex::Vector<amrex::Real> host_data(data_size, 0.0_rt);
        amrex::Real* data_ptr = device_data.dataPtr();

        auto& ion_pc =
            warpx_instance.GetPartContainer().GetParticleContainerFromName(
                "xe_ions");
        const amrex::Real ion_mass = ion_pc.getMass();

        for (int const outlet_boundary : outlet_boundaries) {
            auto* xe_ions =
                mybpc.getParticleBufferPointer("xe_ions", outlet_boundary);
            if (xe_ions != nullptr && xe_ions->isDefined()) {
                for (auto pti = WarpXParIter(*xe_ions, 0); pti.isValid();
                     ++pti) {
                    auto& arr = pti.GetStructOfArrays().GetRealData();
                    auto pw = arr[PIdx::w].dataPtr();
                    auto pvz = arr[PIdx::uz].dataPtr();
                    int np = pti.numParticles();

                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) {
                        amrex::Real const w = pw[ip];
                        amrex::Real const vz = pvz[ip];
                        amrex::Gpu::Atomic::Add(&data_ptr[0], w);
                        amrex::Gpu::Atomic::Add(&data_ptr[1],
                                                w * ion_mass * vz);
                    });
                }
            }
        }

        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_data.begin(),
                         device_data.end(), host_data.begin());
        amrex::ParallelDescriptor::ReduceRealSum(
            host_data.data(), static_cast<int>(host_data.size()),
            amrex::ParallelDescriptor::IOProcessorNumber());

        const amrex::Real time = warpx_instance.gett_new(0);
        const amrex::Real sample_dt = SampleDt(time, last_sample_time, gap);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            const amrex::Real thrust = host_data[1] / sample_dt;
            std::fstream thrust_file(thrust_diag_path, std::ios::app);
            thrust_file << time << "\t" << sample_dt << "\t" << host_data[0]
                        << "\t" << thrust << "\t" << host_data[1] << "\n";
        }
    }
#endif
}

void
BeamDivergenceCalc () {
#ifdef HALL3D
    static bool const diag_enabled = DiagEnabled("beam_divergence_diag");
    if (!diag_enabled) { return; }

    static int times = 0;
    static int gap = 10;
    times++;

    static bool ifinit = false;
    static std::string beam_divergence_path = "beam_divergence.dat";
    static amrex::Real last_sample_time = 0.0_rt;

    if (!ifinit) {
        gap = HallDiagInterval();
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.query("beam_divergence_path", beam_divergence_path);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::fstream divergence_file(beam_divergence_path, std::ios::out);
            divergence_file << "time\tdt\tion_weight\t"
                               "divergence_angle_rad\tdivergence_angle_deg\t"
                               "axial_momentum_kg_m_per_s\t"
                               "transverse_momentum_kg_m_per_s\n";
        }

        ifinit = true;
    }

    if (times == gap) {
        times = 0;
        WarpX& warpx_instance = WarpX::GetInstance();
        auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();
        constexpr int data_size = 3;
        amrex::Gpu::DeviceVector<amrex::Real> device_data(data_size, 0.0_rt);
        amrex::Vector<amrex::Real> host_data(data_size, 0.0_rt);
        amrex::Real* data_ptr = device_data.dataPtr();

        auto& ion_pc =
            warpx_instance.GetPartContainer().GetParticleContainerFromName(
                "xe_ions");
        const amrex::Real ion_mass = ion_pc.getMass();

        for (int const outlet_boundary : outlet_boundaries) {
            auto* xe_ions =
                mybpc.getParticleBufferPointer("xe_ions", outlet_boundary);
            if (xe_ions != nullptr && xe_ions->isDefined()) {
                for (auto pti = WarpXParIter(*xe_ions, 0); pti.isValid();
                     ++pti) {
                    auto& arr = pti.GetStructOfArrays().GetRealData();
                    auto pw = arr[PIdx::w].dataPtr();
                    auto pvx = arr[PIdx::ux].dataPtr();
                    auto pvy = arr[PIdx::uy].dataPtr();
                    auto pvz = arr[PIdx::uz].dataPtr();
                    int np = pti.numParticles();

                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) {
                        amrex::Real const w = pw[ip];
                        amrex::Real const vx = pvx[ip];
                        amrex::Real const vy = pvy[ip];
                        amrex::Real const vz = pvz[ip];
                        amrex::Real const v_perp = std::sqrt(vx * vx + vy * vy);
                        amrex::Gpu::Atomic::Add(&data_ptr[0], w);
                        amrex::Gpu::Atomic::Add(&data_ptr[1],
                                                w * ion_mass * vz);
                        amrex::Gpu::Atomic::Add(&data_ptr[2],
                                                w * ion_mass * v_perp);
                    });
                }
            }
        }

        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_data.begin(),
                         device_data.end(), host_data.begin());
        amrex::ParallelDescriptor::ReduceRealSum(
            host_data.data(), static_cast<int>(host_data.size()),
            amrex::ParallelDescriptor::IOProcessorNumber());

        const amrex::Real time = warpx_instance.gett_new(0);
        const amrex::Real sample_dt = SampleDt(time, last_sample_time, gap);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            const amrex::Real divergence_angle =
                std::atan2(host_data[2], host_data[1]);
            const amrex::Real divergence_angle_deg =
                divergence_angle * amrex::Real(180.0) / MathConst::pi;

            std::fstream divergence_file(beam_divergence_path, std::ios::app);
            divergence_file << time << "\t" << sample_dt << "\t" << host_data[0]
                            << "\t" << divergence_angle << "\t"
                            << divergence_angle_deg << "\t" << host_data[1]
                            << "\t" << host_data[2] << "\n";
        }
    }
#endif
}

void
IEDFCalc () {
#ifdef HALL3D
    static bool const diag_enabled = DiagEnabled("iedf_diag");
    if (!diag_enabled) { return; }

    static int times = 0;
    static int gap = 10;
    times++;

    static bool ifinit = false;
    static std::string iedf_path = "iedf.dat";
    static int iedf_bins = 200;
    static amrex::Real iedf_min_eV = 0.0_rt;
    static amrex::Real iedf_max_eV = 500.0_rt;
    static amrex::Real last_sample_time = 0.0_rt;

    if (!ifinit) {
        gap = HallDiagInterval();
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.query("iedf_path", iedf_path);
        pp_mc.query("iedf_bins", iedf_bins);
        pp_mc.query("iedf_min_eV", iedf_min_eV);
        pp_mc.query("iedf_max_eV", iedf_max_eV);

        iedf_bins = std::max(iedf_bins, 1);
        if (iedf_max_eV <= iedf_min_eV) {
            iedf_max_eV = iedf_min_eV + 1.0_rt;
        }

        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::fstream iedf_file(iedf_path, std::ios::out);
            iedf_file << "time\tbin_center_eV\tbin_lo_eV\tbin_hi_eV\t"
                         "ion_weight\tcurrent_A\tpdf_per_eV\n";
        }

        ifinit = true;
    }

    if (times == gap) {
        times = 0;
        WarpX& warpx_instance = WarpX::GetInstance();
        auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();
        amrex::Gpu::DeviceVector<amrex::Real> device_iedf(iedf_bins, 0.0_rt);
        amrex::Gpu::DeviceVector<amrex::Real> device_total(1, 0.0_rt);
        std::vector<amrex::Real> host_iedf(iedf_bins, 0.0_rt);
        amrex::Vector<amrex::Real> host_total(1, 0.0_rt);

        amrex::Real* iedf_ptr = device_iedf.dataPtr();
        amrex::Real* total_ptr = device_total.dataPtr();

        auto& ion_pc =
            warpx_instance.GetPartContainer().GetParticleContainerFromName(
                "xe_ions");
        const amrex::Real ion_mass = ion_pc.getMass();
        const amrex::Real e_charge = PhysConst::q_e;
        const amrex::Real inv_bin_width =
            static_cast<amrex::Real>(iedf_bins) / (iedf_max_eV - iedf_min_eV);
        const amrex::Real e_min = iedf_min_eV;
        const amrex::Real e_max = iedf_max_eV;
        const int num_bins = iedf_bins;

        for (int const outlet_boundary : outlet_boundaries) {
            auto* xe_ions =
                mybpc.getParticleBufferPointer("xe_ions", outlet_boundary);
            if (xe_ions != nullptr && xe_ions->isDefined()) {
                for (auto pti = WarpXParIter(*xe_ions, 0); pti.isValid();
                     ++pti) {
                    auto& arr = pti.GetStructOfArrays().GetRealData();
                    auto pw = arr[PIdx::w].dataPtr();
                    auto pvx = arr[PIdx::ux].dataPtr();
                    auto pvy = arr[PIdx::uy].dataPtr();
                    auto pvz = arr[PIdx::uz].dataPtr();
                    int np = pti.numParticles();

                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) {
                        amrex::Real const w = pw[ip];
                        amrex::Real const vx = pvx[ip];
                        amrex::Real const vy = pvy[ip];
                        amrex::Real const vz = pvz[ip];
                        amrex::Real const v2 = vx * vx + vy * vy + vz * vz;
                        amrex::Real const energy_eV =
                            amrex::Real(0.5) * ion_mass * v2 / e_charge;

                        amrex::Gpu::Atomic::Add(total_ptr, w);
                        if (energy_eV >= e_min && energy_eV < e_max) {
                            int const bin = static_cast<int>(
                                (energy_eV - e_min) * inv_bin_width);
                            if (bin >= 0 && bin < num_bins) {
                                amrex::Gpu::Atomic::Add(&iedf_ptr[bin], w);
                            }
                        }
                    });
                }
            }
        }

        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_iedf.begin(),
                         device_iedf.end(), host_iedf.begin());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_total.begin(),
                         device_total.end(), host_total.begin());
        amrex::ParallelDescriptor::ReduceRealSum(
            host_iedf.data(), static_cast<int>(host_iedf.size()),
            amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::ReduceRealSum(
            host_total.data(), static_cast<int>(host_total.size()),
            amrex::ParallelDescriptor::IOProcessorNumber());

        const amrex::Real time = warpx_instance.gett_new(0);
        const amrex::Real sample_dt = SampleDt(time, last_sample_time, gap);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            const amrex::Real bin_width = (iedf_max_eV - iedf_min_eV) /
                                          static_cast<amrex::Real>(iedf_bins);
            const amrex::Real total_weight = host_total[0];
            std::fstream iedf_file(iedf_path, std::ios::app);
            for (int i = 0; i < iedf_bins; ++i) {
                const amrex::Real bin_lo =
                    iedf_min_eV + static_cast<amrex::Real>(i) * bin_width;
                const amrex::Real bin_hi = bin_lo + bin_width;
                const amrex::Real bin_center =
                    amrex::Real(0.5) * (bin_lo + bin_hi);
                const amrex::Real current = e_charge * host_iedf[i] / sample_dt;
                const amrex::Real pdf =
                    total_weight > 0.0_rt
                        ? host_iedf[i] / (total_weight * bin_width)
                        : 0.0_rt;
                iedf_file << time << "\t" << bin_center << "\t" << bin_lo
                          << "\t" << bin_hi << "\t" << host_iedf[i] << "\t"
                          << current << "\t" << pdf << "\n";
            }
            iedf_file << "\n";
        }
    }
#endif
}

void
ClearHallBoundaryParticleCache () {
#ifdef HALL3D
    static bool const diag_enabled = DiagEnabled("clear_hall_boundary_particle_cache_diag");
    if (!diag_enabled) { return; }

    static int times = 0;
    static int gap = 10;
    times++;

    static bool ifinit = false;
    if (!ifinit) {
        gap = HallDiagInterval();
        ifinit = true;
    }

    if (times == gap) {
        times = 0;
        WarpX::GetInstance().GetParticleBoundaryBuffer().clearParticles();
    }
#endif
}

} // namespace Insert
