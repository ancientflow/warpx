/* Copyright 2024-2025 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unistd.h>

/**
 * Generate a Xe-Xe neutral elastic-collision cross-section file for the
 * WarpX DSMC module, using the Variable Hard Sphere (VHS) model.
 *
 * Usage:
 *     XeXeVHSGenerator [output_file] [num_points] [E_min_eV] [E_max_eV]
 *
 * Arguments (all optional):
 *     output_file  Output path (default: Xe_Xe_VHS_elastic.dat)
 *     num_points   Number of log-spaced energy points (default: 1000)
 *     E_min_eV     Minimum center-of-mass energy in eV (default: 1e-6)
 *     E_max_eV     Maximum center-of-mass energy in eV (default: 10)
 *
 * VHS parameters for Xe-Xe (see xe_coll.md):
 *     m_Xe  = 131.293 u,  mu = m_Xe/2
 *     T_ref = 273 K
 *     d_ref = 5.74e-10 m
 *     omega = 0.85
 *
 * The cross-section as a function of the center-of-mass energy E [eV]:
 *
 *     sigma(E) = pi d_ref^2 / Gamma(2.5 - omega)
 *                * (E_ref / E)^(omega - 1/2),
 *     E_ref = k_B T_ref / e  [eV]
 *
 * The output file matches WarpX's ScatteringProcess::readCrossSectionFile:
 * two whitespace-separated columns (energy[eV] sigma[m^2]), strictly
 * increasing energy, no comments or unit strings.
 */

std::string
trim (const std::string& str) {
    auto start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return "";
    }
    auto end = str.find_last_not_of(" \t\r\n");
    return str.substr(start, end - start + 1);
}

[[nodiscard]]
std::string
readPathConfig () {
    char exe_path[4096];
    ssize_t len = readlink("/proc/self/exe", exe_path, sizeof(exe_path) - 1);
    if (len == -1) {
        return "";
    }
    exe_path[len] = '\0';

    std::string exe(exe_path);
    auto pos = exe.find_last_of('/');
    if (pos == std::string::npos) {
        return "";
    }

    std::string exe_dir = exe.substr(0, pos + 1);
    std::ifstream conf(exe_dir + "paths.conf");
    if (!conf.is_open()) {
        return "";
    }

    std::string line;
    if (std::getline(conf, line)) {
        std::string base = trim(line);
        if (!base.empty() && base[0] != '/') {
            base = exe_dir + base;
        }
        return base;
    }
    return "";
}

[[nodiscard]]
std::string
resolvePath (const std::string& base_dir, const std::string& filepath) {
    if (filepath.empty() || filepath[0] == '/') {
        return filepath;
    }
    if (base_dir.empty()) {
        return filepath;
    }
    return base_dir + "/" + filepath;
}

// Physical constants (SI)
constexpr double k_B = 1.380649e-23;        // Boltzmann constant [J/K]
constexpr double e_charge = 1.602176634e-19; // elementary charge [J/eV]
constexpr double u_mass = 1.66053906660e-27; // atomic mass unit [kg]

// Xe-Xe VHS parameters
constexpr double m_Xe = 131.293 * u_mass; // Xe mass [kg]
constexpr double T_ref = 273.0;           // reference temperature [K]
constexpr double d_ref = 5.74e-10;        // reference diameter [m]
constexpr double omega = 0.85;            // viscosity index

/**
 * VHS elastic cross-section [m^2] at center-of-mass energy E [eV].
 * E must be > 0 (sigma diverges as E -> 0).
 */
[[nodiscard]]
double
sigmaVHS (double energy_eV) {
    const double E_ref_eV = k_B * T_ref / e_charge;
    const double prefactor =
        M_PI * d_ref * d_ref / std::tgamma(2.5 - omega);
    return prefactor * std::pow(E_ref_eV / energy_eV, omega - 0.5);
}

int
main (int argc, char* argv[]) {
    if (argc > 5) {
        std::cerr << "Usage: " << argv[0]
                  << " [output_file] [num_points] [E_min_eV] [E_max_eV]"
                  << std::endl;
        return 1;
    }

    std::string base_dir = readPathConfig();

    std::string output_file = "Xe_Xe_VHS_elastic.dat";
    int num_points = 1000;
    double e_min = 1e-6;
    double e_max = 10.0;

    if (argc >= 2) {
        output_file = argv[1];
    }
    if (argc >= 3) {
        num_points = std::stoi(argv[2]);
    }
    if (argc >= 4) {
        e_min = std::stod(argv[3]);
    }
    if (argc >= 5) {
        e_max = std::stod(argv[4]);
    }
    output_file = resolvePath(base_dir, output_file);

    if (num_points < 2) {
        std::cerr << "Error: num_points must be >= 2." << std::endl;
        return 1;
    }
    if (e_min <= 0.0) {
        std::cerr << "Error: E_min_eV must be positive (VHS sigma diverges "
                     "at E = 0)."
                  << std::endl;
        return 1;
    }
    if (e_max <= e_min) {
        std::cerr << "Error: E_max_eV must be larger than E_min_eV."
                  << std::endl;
        return 1;
    }

    std::ofstream outfile(output_file);
    if (!outfile.is_open()) {
        std::cerr << "Error: failed to open output file: " << output_file
                  << std::endl;
        return 1;
    }

    // Log-spaced energy grid, strictly increasing
    const double log_ratio =
        std::log(e_max / e_min) / static_cast<double>(num_points - 1);

    outfile << std::scientific << std::setprecision(6);
    double prev_e = 0.0;
    double prev_sigma = 0.0;
    double first_e = 0.0, first_sigma = 0.0, last_e = 0.0, last_sigma = 0.0;
    for (int i = 0; i < num_points; ++i) {
        const double e = e_min * std::exp(static_cast<double>(i) * log_ratio);
        const double sigma = sigmaVHS(e);

        if (sigma <= 0.0) {
            std::cerr << "Error: non-positive cross-section at E = " << e
                      << " eV." << std::endl;
            return 1;
        }
        if (i > 0 && (e <= prev_e || sigma >= prev_sigma)) {
            std::cerr << "Error: grid not strictly increasing or "
                         "cross-section not monotonically decreasing at E = "
                      << e << " eV." << std::endl;
            return 1;
        }

        outfile << e << " " << sigma << "\n";

        if (i == 0) {
            first_e = e;
            first_sigma = sigma;
        }
        prev_e = e;
        prev_sigma = sigma;
        last_e = e;
        last_sigma = sigma;
    }
    outfile.close();

    const double mu = m_Xe / 2.0;
    std::cout << "Xe-Xe VHS elastic cross-section" << std::endl;
    std::cout << "  m_Xe  = " << m_Xe << " kg, mu = " << mu << " kg"
              << std::endl;
    std::cout << "  T_ref = " << T_ref << " K, d_ref = " << d_ref
              << " m, omega = " << omega << std::endl;
    std::cout << "  E_ref = " << k_B * T_ref / e_charge << " eV" << std::endl;
    std::cout << "Wrote " << num_points << " points, range [" << e_min << ", "
              << e_max << "] eV -> " << output_file << std::endl;

    // Numerical checks (xe_coll.md section 6)
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Test values:" << std::endl;
    for (const double e_test : {1e-4, 1e-3, 1e-2, 1e-1, 1.0}) {
        std::cout << "  E = " << e_test << " eV, sigma = " << sigmaVHS(e_test)
                  << " m^2" << std::endl;
    }
    const double slope = std::log(last_sigma / first_sigma) /
                         std::log(last_e / first_e);
    std::cout << "log-log slope = " << slope << " (expected "
              << -(omega - 0.5) << ")" << std::endl;

    return 0;
}
