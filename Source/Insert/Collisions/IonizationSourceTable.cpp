#include "IonizationSourceTable.h"

#if defined(IONIZATION_SOURCE_RECORD)

#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include "Insert/Utils/InsertUtils.h"

#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <string>
#include <vector>

namespace {

using namespace amrex::literals;

constexpr amrex::Real pi = 3.141592653589793238462643383279502884_rt;

struct IonizationSourceTable {
    static constexpr int nr = IONIZATION_SOURCE_NR;
    static constexpr int nz = IONIZATION_SOURCE_NZ;

    void
    initialize (amrex::Geometry const& geom, amrex::ParticleReal elec_weight,
                int max_level) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            amrex::ParallelDescriptor::NProcs() == 1,
            "Average ionization source only supports one MPI rank in the first "
            "version.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            max_level == 0, "Average ionization source only supports one mesh "
                            "level in the first version.");

        auto const& prob_lo = geom.ProbLoArray();
        auto const& prob_hi = geom.ProbHiArray();

        x_center = 0.5_rt * (prob_lo[0] + prob_hi[0]);
        y_center = 0.5_rt * (prob_lo[1] + prob_hi[1]);
        r_min = 0.0_rt;
        r_max =
            std::min(std::min(x_center - prob_lo[0], prob_hi[0] - x_center),
                     std::min(y_center - prob_lo[1], prob_hi[1] - y_center));
        z_min = prob_lo[2];
        z_max = prob_hi[2];

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(r_max > r_min && z_max > z_min,
                                         "Average ionization source requires a "
                                         "non-empty r-z source domain.");

        dr = (r_max - r_min) / static_cast<amrex::Real>(nr);
        dz = (z_max - z_min) / static_cast<amrex::Real>(nz);
        inv_dr = 1.0_rt / dr;
        inv_dz = 1.0_rt / dz;
        recorded_elec_weight = elec_weight;
        initialized = true;

        node_count.resize((nr + 1) * (nz + 1), 0.0_rt);
    }

    void
    recordSample (amrex::Geometry const& geom, amrex::Real dt,
                  amrex::ParticleReal elec_weight, int max_level) {
        if (!initialized) {
            initialize(geom, elec_weight, max_level);
        } else {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::abs(static_cast<amrex::Real>(elec_weight -
                                                  recorded_elec_weight)) <=
                    1.0e-12_rt *
                        std::max(1.0_rt, std::abs(static_cast<amrex::Real>(
                                             recorded_elec_weight))),
                "Average ionization source recording requires a constant "
                "my_constants.elec_weight.");
        }

        sample_count += 1;
        total_time += dt;
    }

    Insert::IonizationSourceRecordView
    view () {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            initialized && !node_count.empty(),
            "Average ionization source recorder must be initialized before deposition.");

        return Insert::IonizationSourceRecordView{node_count.dataPtr(),
                                                  nr,
                                                  nz,
                                                  x_center,
                                                  y_center,
                                                  r_min,
                                                  z_min,
                                                  inv_dr,
                                                  inv_dz};
    }

    amrex::Real
    nodeVolume (int iz, int ir) const {
        // This is the integral of the nodal bilinear basis over 2*pi*r dr dz,
        // not a geometric half-cell control volume.
        amrex::Real volume = 0.0_rt;

        for (int cell_z = std::max(0, iz - 1); cell_z <= std::min(nz - 1, iz);
             ++cell_z) {
            for (int cell_r = std::max(0, ir - 1);
                 cell_r <= std::min(nr - 1, ir); ++cell_r) {
                amrex::Real const r0 =
                    r_min + static_cast<amrex::Real>(cell_r) * dr;
                bool const upper_r_node = (ir == cell_r + 1);
                amrex::Real const radial_basis_integral =
                    upper_r_node ? (0.5_rt * r0 + (1.0_rt / 3.0_rt) * dr)
                                 : (0.5_rt * r0 + (1.0_rt / 6.0_rt) * dr);

                volume += pi * dz * dr * radial_basis_integral;
            }
        }

        return volume;
    }

    amrex::Real
    cellRate (std::vector<amrex::Real> const& node_rate, int iz, int ir) const {
        auto const node_index = [] (int izn, int irn) {
            return izn * (nr + 1) + irn;
        };

        amrex::Real const s00 = node_rate[node_index(iz, ir)];
        amrex::Real const s10 = node_rate[node_index(iz + 1, ir)];
        amrex::Real const s01 = node_rate[node_index(iz, ir + 1)];
        amrex::Real const s11 = node_rate[node_index(iz + 1, ir + 1)];

        amrex::Real const a = s00;
        amrex::Real const b = s10 - s00;
        amrex::Real const c = s01 - s00;
        amrex::Real const d = s11 - s10 - s01 + s00;
        amrex::Real const alpha = a + 0.5_rt * b;
        amrex::Real const beta = c + 0.5_rt * d;
        amrex::Real const r0 = r_min + static_cast<amrex::Real>(ir) * dr;
        amrex::Real const q0 = alpha * r0;
        amrex::Real const q1 = alpha * dr + beta * r0;
        amrex::Real const q2 = beta * dr;
        amrex::Real const integral = q0 + 0.5_rt * q1 + (1.0_rt / 3.0_rt) * q2;

        return std::max(0.0_rt, 2.0_rt * pi * dz * dr * integral);
    }

    void
    writeArray (std::string const& path, std::vector<amrex::Real> const& data,
                int nrow, int ncol) const {
        if (!amrex::ParallelDescriptor::IOProcessor()) {
            return;
        }

        std::ofstream os(path);
        os << std::setprecision(17);
        os << nrow << " " << ncol << "\n";
        for (int i = 0; i < nrow; ++i) {
            for (int j = 0; j < ncol; ++j) {
                os << data[i * ncol + j];
                os << (j + 1 == ncol ? '\n' : ' ');
            }
        }
    }

    void
    writeMetadata (std::string const& path, amrex::Real a_tot) const {
        if (!amrex::ParallelDescriptor::IOProcessor()) {
            return;
        }

        std::ofstream os(path);
        os << std::setprecision(17);
        os << "version = 1\n";
        os << "axis = z\n";
        os << "center = " << x_center << " " << y_center << "\n";
        os << "r_min = " << r_min << "\n";
        os << "r_max = " << r_max << "\n";
        os << "z_min = " << z_min << "\n";
        os << "z_max = " << z_max << "\n";
        os << "nr = " << nr << "\n";
        os << "nz = " << nz << "\n";
        os << "T_avg = " << total_time << "\n";
        os << "N_sample = " << sample_count << "\n";
        os << "A_tot = " << a_tot << "\n";
        os << "elec_weight = " << recorded_elec_weight << "\n";
        os << "coordinate_units = simulation_units\n";
        os << "node_source_units = real_pairs_per_m3_per_second\n";
        os << "cell_rate_units = real_pairs_per_second_per_rz_cell\n";
        os << "plotfile = ionization_source_plt\n";
    }

    void
    writePlotfile (std::vector<amrex::Real> const& cell_rate,
                   std::vector<amrex::Real> const& cell_count) const {
        amrex::Gpu::DeviceVector<amrex::Real> d_cell_rate(cell_rate.size());
        amrex::Gpu::DeviceVector<amrex::Real> d_cell_count(cell_count.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, cell_rate.begin(),
                         cell_rate.end(), d_cell_rate.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, cell_count.begin(),
                         cell_count.end(), d_cell_count.begin());

        amrex::Box const domain(amrex::IntVect(0, 0, 0),
                                amrex::IntVect(nr - 1, nz - 1, 0));
        amrex::BoxArray ba(domain);
        amrex::DistributionMapping dm(ba);
        amrex::MultiFab plot_mf(ba, dm, 4, 0);

        amrex::Real const* const cell_rate_ptr = d_cell_rate.dataPtr();
        amrex::Real const* const cell_count_ptr = d_cell_count.dataPtr();
        amrex::Real const rmin = r_min;
        amrex::Real const dr_plot = dr;
        amrex::Real const dz_plot = dz;

        for (amrex::MFIter mfi(plot_mf); mfi.isValid(); ++mfi) {
            auto const& arr = plot_mf.array(mfi);
            amrex::ParallelFor(
                mfi.validbox(),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    amrex::ignore_unused(k);
                    int const ir = i;
                    int const iz = j;
                    int const idx = iz * nr + ir;
                    amrex::Real const r0 =
                        rmin + static_cast<amrex::Real>(ir) * dr_plot;
                    amrex::Real const r1 = r0 + dr_plot;
                    amrex::Real const volume =
                        pi * (r1 * r1 - r0 * r0) * dz_plot;
                    arr(i, j, k, 0) = cell_rate_ptr[idx];
                    arr(i, j, k, 1) = cell_rate_ptr[idx] / volume;
                    arr(i, j, k, 2) = cell_count_ptr[idx];
                    arr(i, j, k, 3) = volume;
                });
        }
        amrex::Gpu::streamSynchronize();

        amrex::RealBox const real_box({AMREX_D_DECL(r_min, z_min, 0.0_rt)},
                                      {AMREX_D_DECL(r_max, z_max, 1.0_rt)});
        amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
        amrex::Geometry source_geom(domain, &real_box, 0, is_periodic.data());

        amrex::Vector<amrex::MultiFab const*> mfs{&plot_mf};
        amrex::Vector<std::string> varnames{"source_rate", "source_density",
                                            "source_count", "cell_volume"};
        amrex::Vector<amrex::Geometry> geoms{source_geom};
        amrex::Vector<int> level_steps{0};
        amrex::Vector<amrex::IntVect> ref_ratio{amrex::IntVect(2)};

        amrex::WriteMultiLevelPlotfile("ionization_source_plt", 1, mfs,
                                       varnames, geoms, 0.0_rt, level_steps,
                                       ref_ratio);
    }

    void
    finalize () {
        if (!initialized) {
            return;
        }

        if (total_time <= 0.0_rt) {
            ablastr::warn_manager::WMRecordWarning(
                "Average ionization source",
                "Skipping ionization source output because no collision "
                "samples were recorded.");
            return;
        }

        amrex::Gpu::streamSynchronize();
        std::vector<amrex::Real> h_node_count(node_count.size());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, node_count.begin(),
                         node_count.end(), h_node_count.begin());

        std::vector<amrex::Real> node_rate(h_node_count.size(), 0.0_rt);
        for (int iz = 0; iz <= nz; ++iz) {
            for (int ir = 0; ir <= nr; ++ir) {
                int const idx = iz * (nr + 1) + ir;
                amrex::Real const volume = nodeVolume(iz, ir);
                node_rate[idx] = (volume > 0.0_rt)
                                     ? h_node_count[idx] / (total_time * volume)
                                     : 0.0_rt;
            }
        }

        std::vector<amrex::Real> cell_rate(nr * nz, 0.0_rt);
        std::vector<amrex::Real> cell_cdf(nr * nz, 0.0_rt);
        std::vector<amrex::Real> cell_count(nr * nz, 0.0_rt);
        amrex::Real a_tot = 0.0_rt;
        for (int iz = 0; iz < nz; ++iz) {
            for (int ir = 0; ir < nr; ++ir) {
                int const idx = iz * nr + ir;
                cell_rate[idx] = cellRate(node_rate, iz, ir);
                cell_count[idx] = cell_rate[idx] * total_time;
                a_tot += cell_rate[idx];
                cell_cdf[idx] = a_tot;
            }
        }

        if (a_tot > 0.0_rt) {
            for (auto& value : cell_cdf) {
                value /= a_tot;
            }
            cell_cdf.back() = 1.0_rt;
        }

        Insert::CreateDirectoryTree("ionization_source_fab");
        writeMetadata("ionization_source_fab/metadata.txt", a_tot);
        writeArray("ionization_source_fab/node_source_count", h_node_count,
                   nz + 1, nr + 1);
        writeArray("ionization_source_fab/node_source_rate", node_rate, nz + 1,
                   nr + 1);
        writeArray("ionization_source_fab/cell_rate", cell_rate, nz, nr);
        writeArray("ionization_source_fab/cell_cdf", cell_cdf, nz, nr);
        writePlotfile(cell_rate, cell_count);

        amrex::Print() << Utils::TextMsg::Info(
            "Average ionization source output: T_avg = " +
            std::to_string(total_time) + ", A_tot = " + std::to_string(a_tot));
    }

    bool initialized = false;
    long sample_count = 0;
    amrex::Real total_time = 0.0_rt;
    amrex::ParticleReal recorded_elec_weight = 0.0_prt;
    amrex::Real x_center = 0.0_rt;
    amrex::Real y_center = 0.0_rt;
    amrex::Real r_min = 0.0_rt;
    amrex::Real r_max = 0.0_rt;
    amrex::Real z_min = 0.0_rt;
    amrex::Real z_max = 0.0_rt;
    amrex::Real dr = 0.0_rt;
    amrex::Real dz = 0.0_rt;
    amrex::Real inv_dr = 0.0_rt;
    amrex::Real inv_dz = 0.0_rt;
    amrex::Gpu::DeviceVector<amrex::Real> node_count;
};

IonizationSourceTable&
GlobalIonizationSourceTable () {
    static IonizationSourceTable table;
    return table;
}

} // namespace

namespace Insert {

void
IonizationSourceRecordCollisionSample (amrex::Geometry const& geom,
                                       amrex::Real dt,
                                       amrex::ParticleReal elec_weight,
                                       int max_level) {
    GlobalIonizationSourceTable().recordSample(geom, dt, elec_weight,
                                               max_level);
}

IonizationSourceRecordView
IonizationSourceGetRecordView () {
    return GlobalIonizationSourceTable().view();
}

void
IonizationSourceFinalize () {
    GlobalIonizationSourceTable().finalize();
}

} // namespace Insert

#endif
