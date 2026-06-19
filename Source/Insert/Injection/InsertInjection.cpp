#include "InsertInjection.h"

#include "WarpX.H"

#include "Insert/Config/WarpXSimulationConfig.h"
#include "Insert/Injection/DistributionSampler1D.H"
#include "Insert/Injection/HallInjector.h"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace Insert {

namespace {

amrex::RandomEngine
MakeRandomEngine () {
#ifdef AMREX_USE_GPU
    return amrex::RandomEngine(nullptr);
#else
    return amrex::RandomEngine{};
#endif
}

} // namespace

#ifdef BENCHMARK_2D
namespace {

int
ElectronsCollection (WarpX& warpx_instance) {
    auto& boundary_buffer = warpx_instance.GetParticleBoundaryBuffer();
    int e_xlo =
            boundary_buffer.getNumParticlesInContainer("electrons", 0, false),
        xe_xlo =
            boundary_buffer.getNumParticlesInContainer("xe_ions", 0, false);
    return std::max(e_xlo - xe_xlo, 0);
}

void
ParticleInjection (WarpX& warpx_instance) {
    static const int nx = 512, ny = 256, nppc = 75;
    static const double lx = 0.025, ly = 0.0128, NPlasma = 5e16, s0 = 5.23e23,
                        const_dt = 5e-12;

    static const int once_injection = 1000;
    static double rest_macro_particles = 0;
    static const double one_step_real_particles =
        s0 * const_dt * 0.0128 * 2 / 3.1415926 * 0.0075;
    static const amrex::Real global_weight = NPlasma / nx / ny * lx * ly / nppc;
    static const double one_step_macro_particles =
        one_step_real_particles / global_weight;

    int this_step_pair = 0;
    int const this_step_cathode_electron = ElectronsCollection(warpx_instance);
    rest_macro_particles += one_step_macro_particles;
    if (rest_macro_particles > once_injection) {
        rest_macro_particles -= once_injection;
        this_step_pair = once_injection;
    }

    auto& mypc = warpx_instance.GetPartContainer();
    auto& electons_pc =
        mypc.GetParticleContainer(mypc.getSpeciesID("electrons"));
    auto& xe_pc = mypc.GetParticleContainer(mypc.getSpeciesID("xe_ions"));

    const amrex::Vector<amrex::Vector<int>> nattr;
    amrex::RandomEngine uniform_engine(MakeRandomEngine());
    amrex::RandomEngine normal_engine(MakeRandomEngine());

    int const electron_size = this_step_pair + this_step_cathode_electron;
    std::cout << electron_size << " " << this_step_pair << " " << std::endl;
    amrex::Vector<amrex::Real> pxe(electron_size), pye(electron_size, 0),
        pze(electron_size), vxe(electron_size), vye(electron_size),
        vze(electron_size), we(electron_size, global_weight),
        pxxe(this_step_pair), pyxe(this_step_pair, 0), pzxe(this_step_pair),
        vxxe(this_step_pair), vyxe(this_step_pair), vzxe(this_step_pair),
        wxe(this_step_pair, global_weight);

    if (this_step_pair > 0) {
        for (int i = 0; i < this_step_pair; i++) {
            double const r1 = amrex::Random(uniform_engine);
            double const r2 = amrex::Random(uniform_engine);
            pxe[i] = 0.00625 + std::asin(2 * r1 - 1) / 3.14159 * 0.0075;
            pze[i] = 0.0128 * r2;
            pxxe[i] = pxe[i];
            pzxe[i] = pze[i];

            vxe[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vye[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vze[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vxxe[i] = amrex::RandomNormal(0, 606.34, normal_engine);
            vyxe[i] = amrex::RandomNormal(0, 606.34, normal_engine);
            vzxe[i] = amrex::RandomNormal(0, 606.34, normal_engine);
        }
    }

    if (this_step_cathode_electron > 0) {
        for (int i = this_step_pair; i < electron_size; i++) {
            pxe[i] = 0.024;
            pze[i] = 0.0128 * amrex::Random(uniform_engine);
            vxe[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vye[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vze[i] = amrex::RandomNormal(0, 1326232, normal_engine);
        }
    }
    if (this_step_pair > 0) {
        std::cout << "add pairs number: " << this_step_pair << std::endl;
        xe_pc.AddNParticles(0, this_step_pair, pxxe, pyxe, pzxe, vxxe, vyxe,
                            vzxe, 1, {wxe}, 0, nattr, false);
    }

    if (electron_size > 0) {
        electons_pc.AddNParticles(0, electron_size, pxe, pye, pze, vxe, vye,
                                  vze, 1, {we}, 0, nattr, false);
        if (this_step_cathode_electron > 0) {
            auto& boundary_buffer = warpx_instance.GetParticleBoundaryBuffer();
            boundary_buffer.clearParticles();
        }
    }
}

} // namespace
#endif

#ifdef HALL3D
namespace {

[[maybe_unused]] amrex::ParticleReal
single_spoke_distribution (amrex::ParticleReal theta) {
    constexpr amrex::ParticleReal pi = amrex::Math::pi<amrex::ParticleReal>();
    constexpr amrex::ParticleReal center = pi;
    constexpr amrex::ParticleReal sigma = pi / amrex::ParticleReal(4.0);
    constexpr amrex::ParticleReal normalization =
        amrex::ParticleReal(0.5079812642688625);

    const amrex::ParticleReal normalized_distance = (theta - center) / sigma;
    return normalization * std::exp(-amrex::ParticleReal(0.5) *
                                    normalized_distance * normalized_distance);
}

template <int spoke_count>
amrex::ParticleReal
multi_spoke_distribution (amrex::ParticleReal theta) {
    constexpr amrex::ParticleReal pi = amrex::Math::pi<amrex::ParticleReal>();
    constexpr amrex::ParticleReal period = amrex::ParticleReal(2.0) * pi;
    constexpr amrex::ParticleReal spoke_interval =
        period / static_cast<amrex::ParticleReal>(spoke_count);

    amrex::ParticleReal local_theta = std::fmod(theta, spoke_interval);
    if (local_theta < amrex::ParticleReal(0.0)) {
        local_theta += spoke_interval;
    }

    // Map each spoke interval back to [0, 2*pi] and reuse the normalized
    // single-spoke distribution. No extra factor is needed: each compressed
    // interval contributes 1/spoke_count of the single-spoke integral.
    const amrex::ParticleReal equivalent_single_spoke_theta =
        local_theta * static_cast<amrex::ParticleReal>(spoke_count);
    return single_spoke_distribution(equivalent_single_spoke_theta);
}

amrex::ParticleReal
uniform_distribution (amrex::ParticleReal theta) {
    amrex::ignore_unused(theta);
    return 0.5 / amrex::Math::pi<amrex::ParticleReal>();
}

DistributionSampler1D
MakePlasmaThetaSampler () {
    constexpr int theta_num_bins = 1024;
    constexpr amrex::ParticleReal theta_min = 0.0;
    constexpr amrex::ParticleReal theta_max =
        amrex::ParticleReal(2.0) *
        amrex::ParticleReal(3.1415926535897932384626433832795);

    // Hard-coded single-spoke probability density. The spoke peaks at the
    // interval midpoint, and both interval boundaries are exactly 4 sigma away.
    return DistributionSampler1D(theta_min, theta_max, theta_num_bins,
                                 uniform_distribution);
}

} // namespace
#endif

void
LegacyCathodeInjection3D () {
#ifdef HALL3D
    static const double L = 0.05;
    static double Ic, dt, l_factor, elec_weight;

    static bool ifinit = false;

    if (!ifinit) {
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("dt", dt);
        pp_mc.get("Ic", Ic);
        pp_mc.get("l_factor", l_factor);
        pp_mc.getWithParser("elec_weight", elec_weight);
        ifinit = true;
    }

    const double num_of_real_electrons =
                     Ic * dt / l_factor / l_factor / 1.60217e-19,
                 num_of_marco_particle = num_of_real_electrons / elec_weight;

    WarpX& warpx_instance = WarpX::GetInstance();

    static double electron_rest = 0;
    const int one_step_injection = static_cast<int>(num_of_marco_particle) + 1;

    electron_rest += num_of_marco_particle;

    static const double r1 = 0.016 / l_factor, r2 = 0.02 / l_factor,
                        dr = r2 - r1, sumr = r1 + r2, multr = r1 * r2,
                        z1 = 0.042 / l_factor, z2 = 0.046 / l_factor,
                        dz = z2 - z1, Pi2 = 3.1415926 * 2, sigma = 592982,
                        half_l = L / 2 / l_factor;

    double r, theta;
    static amrex::Vector<amrex::Vector<int>> nattr;
    if (electron_rest > one_step_injection) {

        electron_rest -= one_step_injection;

        amrex::Vector<amrex::ParticleReal> pz(one_step_injection),
            px(one_step_injection), py(one_step_injection),
            vx(one_step_injection), vy(one_step_injection),
            vz(one_step_injection), pw(one_step_injection, elec_weight);

        amrex::RandomEngine uniform_engine(MakeRandomEngine()),
            normal_engine(MakeRandomEngine());

        for (int i = 0; i < one_step_injection; i++) {
            r = amrex::Random(uniform_engine) * dr + r1;
            theta = amrex::Random(uniform_engine) * Pi2;
            pz[i] = amrex::Random(uniform_engine) * dz + z1;

            r = std::sqrt(sumr * r - multr);
            px[i] = r * std::cos(theta) + half_l;
            py[i] = r * std::sin(theta) + half_l;

            vx[i] = amrex::RandomNormal(0, sigma, normal_engine);
            vy[i] = amrex::RandomNormal(0, sigma, normal_engine);
            vz[i] = amrex::RandomNormal(0, sigma, normal_engine);
        }

        auto& mypc = warpx_instance.GetPartContainer();
        auto& e_pc = mypc.GetParticleContainer(mypc.getSpeciesID("electrons"));
        e_pc.AddNParticles(0, one_step_injection, px, py, pz, vx, vy, vz, 1,
                           {pw}, 0, nattr, 0);
    }
#endif
}

void
LegacyPlasmaInit () {
#ifdef HALL3D
    const double L = 0.05;
    double l_factor, elec_weight;

    amrex::ParmParse pp_mc("my_constants");
    pp_mc.get("l_factor", l_factor);
    pp_mc.getWithParser("elec_weight", elec_weight);

    const double r1 = 0.0105 / l_factor, r2 = 0.0155 / l_factor, dr = r2 - r1,
                 sumr = r1 + r2, multr = r1 * r2, z1 = 0.001 / l_factor,
                 z2 = 0.004 / l_factor, dz = z2 - z1, sigma_e = 592982,
                 sigma_xe = 1212.41, half_l = L / 2 / l_factor;

    const double volume = (r2 * r2 - r1 * r1) * 3.14159 * dz;
    const auto theta_sampler = MakePlasmaThetaSampler();
    const auto num_elec = static_cast<int>(volume * 1e18 / elec_weight);

    WarpX& warpx_instance = WarpX::GetInstance();

    double r, theta;
    static amrex::Vector<amrex::Vector<int>> nattr;

    amrex::Vector<amrex::ParticleReal> pz(num_elec), px(num_elec), py(num_elec),
        vx(num_elec), vy(num_elec), vz(num_elec), pw(num_elec, elec_weight);

    amrex::RandomEngine uniform_engine(MakeRandomEngine()),
        normal_engine(MakeRandomEngine());

    for (int i = 0; i < num_elec; i++) {
        r = amrex::Random(uniform_engine) * dr + r1;
        theta = theta_sampler.sample(uniform_engine);
        pz[i] = amrex::Random(uniform_engine) * dz + z1;

        r = std::sqrt(sumr * r - multr);
        px[i] = r * std::cos(theta) + half_l;
        py[i] = r * std::sin(theta) + half_l;

        vx[i] = amrex::RandomNormal(0, sigma_e, normal_engine);
        vy[i] = amrex::RandomNormal(0, sigma_e, normal_engine);
        vz[i] = amrex::RandomNormal(0, sigma_e, normal_engine);
    }
    auto& mypc = warpx_instance.GetPartContainer();
    auto& e_pc = mypc.GetParticleContainer(mypc.getSpeciesID("electrons"));
    e_pc.AddNParticles(0, num_elec, px, py, pz, vx, vy, vz, 1, {pw}, 0, nattr,
                       0);
    for (int i = 0; i < num_elec; i++) {
        vx[i] = amrex::RandomNormal(0, sigma_xe, normal_engine);
        vy[i] = amrex::RandomNormal(0, sigma_xe, normal_engine);
        vz[i] = amrex::RandomNormal(0, sigma_xe, normal_engine);
    }

    auto& xe_pc = mypc.GetParticleContainer(mypc.getSpeciesID("xe_ions"));
    xe_pc.AddNParticles(0, num_elec, px, py, pz, vx, vy, vz, 1, {pw}, 0, nattr,
                        0);
#endif
}

void
LegacyXeInjection () {
#ifdef HALL3D
    static const double L = 0.05, NA = 6.03e23;
    static bool ifinit = false, ifhole = false;
    static int hole_num = 48;
    static double dt, l_factor, atom_weight, m_dot, Tx, Ty, Tz, vz0;
    static amrex::Vector<double> hole_x, hole_y;

    if (!ifinit) {
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("dt", dt);
        pp_mc.get("l_factor", l_factor);
        pp_mc.getWithParser("xe_weight", atom_weight);
        pp_mc.get("m_dot", m_dot);
        pp_mc.query("ifhole", ifhole);
        pp_mc.get("Tx", Tx);
        pp_mc.get("Ty", Ty);
        pp_mc.get("Tz", Tz);
        pp_mc.get("vz0", vz0);
        ifinit = true;
    }
    static const double V_per_sec = 0.06 / 60. / 1e6,
                        n_per_sec = V_per_sec * 101325 / 8.314 / 273.15 * NA,
                        n_per_step = n_per_sec * dt, mxe = 2.179e-25,
                        n_marco_per_step = n_per_step / atom_weight,
                        Pi2 = 3.1415926 * 2;

    const double n_marco_per_step_m =
        m_dot * dt / l_factor / l_factor / mxe / atom_weight;

    const unsigned one_times_inject_particle =
        static_cast<int>((n_marco_per_step_m * 100) + 1);
    static double xe_rest = 0;

    const double rmid = (0.021 + 0.031) / 4 / l_factor,
                 rwidth = 0.001 / l_factor, r1 = rmid - rwidth,
                 r2 = rmid + rwidth, sqdr = (r2 * r2) - (r1 * r1),
                 sqr1 = r1 * r1, kb = 1.38064852e-23,
                 sigmax = std::sqrt(kb * Tx / mxe),
                 sigmay = std::sqrt(kb * Ty / mxe),
                 sigmaz = std::sqrt(kb * Tz / mxe), half_L = L / 2 / l_factor,
                 sqdr_hole = rwidth * rwidth;
    amrex::ignore_unused(n_marco_per_step, sqdr_hole);
    double r, theta;
    static const amrex::Vector<amrex::Vector<int>> nattr;

    static bool if_hole_init = false;
    if (!if_hole_init) {
        const double per_angle = Pi2 / hole_num;
        double angle = 0;
        for (int i = 0; i < hole_num; i++) {
            hole_x.push_back(rmid * std::cos(angle));
            hole_y.push_back(rmid * std::sin(angle));
            angle += per_angle;
        }
        if_hole_init = true;
    }

    xe_rest += n_marco_per_step_m;

    if (xe_rest > one_times_inject_particle) {
        WarpX& warpx_instance = WarpX::GetInstance();

        xe_rest -= one_times_inject_particle;

        amrex::Vector<amrex::ParticleReal> pz(one_times_inject_particle, 0),
            px(one_times_inject_particle), py(one_times_inject_particle),
            vx(one_times_inject_particle), vy(one_times_inject_particle),
            vz(one_times_inject_particle),
            pw(one_times_inject_particle, atom_weight);

        amrex::RandomEngine uniform_engine(MakeRandomEngine()),
            normal_engine(MakeRandomEngine());
        int hole_start = std::rand() % hole_num;
        for (unsigned i = 0; i < one_times_inject_particle; i++) {
            if (!ifhole) {
                theta = amrex::Random(uniform_engine) * Pi2;
                r = std::sqrt(sqdr * amrex::Random(uniform_engine) + sqr1);
                px[i] = r * std::cos(theta) + half_L;
                py[i] = r * std::sin(theta) + half_L;
            } else {
                theta = amrex::Random(uniform_engine) * Pi2;
                r = rwidth * std::sqrt(amrex::Random(uniform_engine));
                px[i] = r * std::cos(theta) + half_L +
                        hole_x[(static_cast<int>(i) + hole_start) % hole_num];
                py[i] = r * std::sin(theta) + half_L +
                        hole_y[(static_cast<int>(i) + hole_start) % hole_num];
            }

            vx[i] = amrex::RandomNormal(0, sigmax, normal_engine);
            vy[i] = amrex::RandomNormal(0, sigmay, normal_engine);
            do {
                vz[i] = amrex::RandomNormal(0, sigmaz, normal_engine) + vz0;
            } while (vz[i] < 0);
        }

        auto& mypc = warpx_instance.GetPartContainer();
        auto& xe_pc =
            mypc.GetParticleContainer(mypc.getSpeciesID("xe_netural"));
        xe_pc.AddNParticles(0, one_times_inject_particle, px, py, pz, vx, vy,
                            vz, 1, {pw}, 0, nattr, 0);
        amrex::Print() << "Injection Xe Atom\n";
    }
#endif
}

void
LegacyXeFastInjection () {
#ifdef HALL3D_INIT
    double atom_weight, l_factor, m_dot, Tx, Ty, Tz, vz0;
    bool ifhole = false;
    static amrex::Vector<double> hole_x, hole_y;
    static int hole_num = 48;

    amrex::ParmParse pp_mc("my_constants");
    pp_mc.get("l_factor", l_factor);
    pp_mc.getWithParser("xe_weight", atom_weight);
    pp_mc.get("m_dot", m_dot);
    pp_mc.query("ifhole", ifhole);
    pp_mc.get("Tx", Tx);
    pp_mc.get("Ty", Ty);
    pp_mc.get("Tz", Tz);
    pp_mc.get("vz0", vz0);

    static const double L = 0.05, NA = 6.03e23;
    const double V_per_sec = 0.06 / 60. / 1e6, dt = 5.6e-10,
                 n_per_sec = V_per_sec * 101325 / 8.314 / 273.15 * NA,
                 n_per_step = n_per_sec * dt,
                 n_marco_per_step = n_per_step / atom_weight, mxe = 2.179e-25,
                 Pi2 = 3.1415926 * 2;
    amrex::ignore_unused(n_marco_per_step);
    const double n_marco_per_step_m =
        m_dot * dt / mxe / l_factor / l_factor / atom_weight;

    const double one_times_inject_particle =
        static_cast<int>(n_marco_per_step_m + 1);
    static double xe_rest = 0;
    static const double rmid = (0.021 + 0.031) / 4 / l_factor,
                        rwidth = 0.001 / l_factor, r1 = rmid - rwidth,
                        r2 = rmid + rwidth, sqdr = r2 * r2 - r1 * r1,
                        sqr1 = r1 * r1, kb = 1.38064852e-23,
                        sigmax = std::sqrt(kb * Tx / mxe),
                        sigmay = std::sqrt(kb * Ty / mxe),
                        sigmaz = std::sqrt(kb * Tz / mxe),
                        half_L = L / 2 / l_factor;
    double r, theta;
    static const amrex::Vector<amrex::Vector<int>> nattr;

    static bool if_hole_init = false;
    if (!if_hole_init) {
        double per_angle = Pi2 / hole_num, angle = 0;
        for (int i = 0; i < hole_num; i++) {
            hole_x.push_back(rmid * std::cos(angle));
            hole_y.push_back(rmid * std::sin(angle));
            angle += per_angle;
        }
        if_hole_init = true;
    }

    xe_rest += n_marco_per_step_m;
    if (xe_rest > one_times_inject_particle) {
        WarpX& warpx_instance = WarpX::GetInstance();

        xe_rest -= one_times_inject_particle;

        amrex::Vector<amrex::ParticleReal> pz(one_times_inject_particle, 0),
            px(one_times_inject_particle), py(one_times_inject_particle),
            vx(one_times_inject_particle), vy(one_times_inject_particle),
            vz(one_times_inject_particle, vz0),
            pw(one_times_inject_particle, atom_weight);

        amrex::RandomEngine uniform_engine(MakeRandomEngine()),
            normal_engine(MakeRandomEngine());

        int hole_start = std::rand() % hole_num;
        for (int i = 0; i < one_times_inject_particle; i++) {
            theta = amrex::Random(uniform_engine) * Pi2;
            r = std::sqrt(sqdr * amrex::Random(uniform_engine) + sqr1);

            if (!ifhole) {
                theta = amrex::Random(uniform_engine) * Pi2;
                r = std::sqrt(sqdr * amrex::Random(uniform_engine) + sqr1);
                px[i] = r * std::cos(theta) + half_L;
                py[i] = r * std::sin(theta) + half_L;
            } else {
                theta = amrex::Random(uniform_engine) * Pi2;
                r = rwidth * std::sqrt(amrex::Random(uniform_engine));
                px[i] = r * std::cos(theta) + half_L +
                        hole_x[(i + hole_start) % hole_num];
                py[i] = r * std::sin(theta) + half_L +
                        hole_y[(i + hole_start) % hole_num];
            }

            vx[i] = amrex::RandomNormal(0, sigmax, normal_engine);
            vy[i] = amrex::RandomNormal(0, sigmay, normal_engine);
            do {
                vz[i] = amrex::RandomNormal(0, sigmaz, normal_engine) + vz0;
            } while (vz[i] < 0);
        }

        auto& mypc = warpx_instance.GetPartContainer();
        auto& xe_pc =
            mypc.GetParticleContainer(mypc.getSpeciesID("xe_netural"));
        xe_pc.AddNParticles(0, one_times_inject_particle, px, py, pz, vx, vy,
                            vz, 1, {pw}, 0, nattr, 0);
        amrex::Print() << "Injection Xe Atom\n";
    }
#endif
}

void
CathodeInjection3D ()
{
#ifdef HALL3D
    InjectHallParticles();
#endif
}

void
PlasmaInit ()
{
#ifdef HALL3D
    InitializeHallInjection();
#endif
}

void
XeInjection ()
{
#ifdef HALL3D
    InjectHallParticles();
#endif
}

void
XeFastInjection ()
{
#ifdef HALL3D_INIT
    InjectHallParticles();
#endif
}

void
InitializeHallInjection ()
{
#if defined(HALL3D) || defined(HALL3D_INIT)
    HallInjector::GetInstance().InitializePlasma(WarpX::GetInstance());
#endif
}

void
InjectHallParticles ()
{
#if defined(HALL3D) || defined(HALL3D_INIT)
    amrex::Real dt = 0.0;
    amrex::ParmParse pp_mc("my_constants");
    pp_mc.query("dt", dt);
    HallInjector::GetInstance().InjectParticles(WarpX::GetInstance(), dt, 0);
#endif
}

} // namespace Insert
