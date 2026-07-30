/* Copyright 2024-2025 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

/**
 * Interpolate a cross-section data file onto an evenly-spaced energy grid
 * required by WarpX's ScatteringProcess.
 *
 * Usage:
 *     CrossSectionInterpolator <input_file> <energy_step> [output_file]
 *
 * Arguments:
 *     input_file   Path to the original cross-section file (two columns:
 * energy[eV] sigma[m^2]) energy_step  Desired uniform energy spacing in eV
 *     output_file  Optional output file path (default:
 * <input_file>_uniform.dat)
 *
 * The input file format matches WarpX's
 * ScatteringProcess::readCrossSectionFile:
 *     - Each line contains two whitespace-separated values: energy and
 * cross-section
 *     - Energy values do NOT need to be evenly spaced
 *     - Lines starting with '#' are treated as comments and ignored
 *
 * The output file will have:
 *     - Evenly-spaced energy grid from min(energy) to max(energy)
 *     - Energy step equal to the specified value (or very close to it)
 *     - Linear interpolation for sigma values; extrapolation uses the nearest
 * endpoint
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

struct CrossSectionData {
    std::vector<double> energies;
    std::vector<double> sigmas;
};

[[nodiscard]]
CrossSectionData
readCrossSectionFile (const std::string& filepath) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        std::cerr << "Error: failed to open input file: " << filepath
                  << std::endl;
        std::exit(1);
    }

    CrossSectionData data;
    std::string line;
    int line_number = 0;
    while (std::getline(infile, line)) {
        ++line_number;
        // Trim leading whitespace
        auto start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) {
            continue;
        } // empty line
        if (line[start] == '#') {
            continue;
        } // comment line

        std::istringstream iss(line);
        double energy = 0.0, sigma = 0.0;
        if (!(iss >> energy >> sigma)) {
            std::cerr << "Error: invalid line format at line " << line_number
                      << " (expected two columns): " << line << std::endl;
            std::exit(1);
        }
        data.energies.push_back(energy);
        data.sigmas.push_back(sigma);
    }

    if (data.energies.size() < 2) {
        std::cerr << "Error: at least two data points are required."
                  << std::endl;
        std::exit(1);
    }

    return data;
}

[[nodiscard]]
double
interpolateLinear (const std::vector<double>& energies,
                   const std::vector<double>& sigmas, double target_energy) {
    if (target_energy <= energies.front()) {
        return sigmas.front();
    }
    if (target_energy >= energies.back()) {
        return sigmas.back();
    }

    // Binary search for the bracketing interval
    auto it = std::upper_bound(energies.begin(), energies.end(), target_energy);
    std::size_t hi =
        static_cast<std::size_t>(std::distance(energies.begin(), it));
    std::size_t lo = hi - 1;

    const double e0 = energies[lo];
    const double e1 = energies[hi];
    const double s0 = sigmas[lo];
    const double s1 = sigmas[hi];

    const double t = (target_energy - e0) / (e1 - e0);
    return s0 + t * (s1 - s0);
}

[[nodiscard]]
CrossSectionData
generateUniformGrid (const CrossSectionData& data, double dE) {
    const double e_min = data.energies.front();
    const double e_max = data.energies.back();

    int n_points = static_cast<int>(std::round((e_max - e_min) / dE)) + 1;
    // Ensure e_max is included
    if (std::abs((e_min + (n_points - 1) * dE) - e_max) > 1e-12) {
        ++n_points;
    }

    CrossSectionData result;
    result.energies.reserve(static_cast<std::size_t>(n_points));
    result.sigmas.reserve(static_cast<std::size_t>(n_points));

    for (int i = 0; i < n_points; ++i) {
        double e = e_min + static_cast<double>(i) * dE;
        if (e > e_max) {
            e = e_max;
        }
        result.energies.push_back(e);
        result.sigmas.push_back(
            interpolateLinear(data.energies, data.sigmas, e));
    }

    return result;
}

void
writeCrossSectionFile (const std::string& filepath,
                       const CrossSectionData& data) {
    std::ofstream outfile(filepath);
    if (!outfile.is_open()) {
        std::cerr << "Error: failed to open output file: " << filepath
                  << std::endl;
        std::exit(1);
    }

    outfile << std::setprecision(12);
    for (std::size_t i = 0; i < data.energies.size(); ++i) {
        outfile << data.energies[i] << "    " << data.sigmas[i] << "\n";
    }
}

int
main (int argc, char* argv[]) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_file> <energy_step> [output_file]" << std::endl;
        return 1;
    }

    std::string base_dir = readPathConfig();

    const std::string input_arg = argv[1];
    const std::string input_file = resolvePath(base_dir, input_arg);
    const double energy_step = std::stod(argv[2]);

    if (energy_step <= 0.0) {
        std::cerr << "Error: energy_step must be positive." << std::endl;
        return 1;
    }

    std::string output_file;
    if (argc == 4) {
        output_file = resolvePath(base_dir, argv[3]);
    } else {
        // Default: append _uniform before the extension
        const std::size_t last_dot = input_file.find_last_of('.');
        if (last_dot != std::string::npos) {
            output_file = input_file.substr(0, last_dot) + "_uniform" +
                          input_file.substr(last_dot);
        } else {
            output_file = input_file + "_uniform";
        }
    }

    CrossSectionData data = readCrossSectionFile(input_file);
    CrossSectionData uniform = generateUniformGrid(data, energy_step);

    // Sanity check: verify uniform spacing
    const double actual_dE = uniform.energies[1] - uniform.energies[0];
    for (std::size_t i = 2; i < uniform.energies.size(); ++i) {
        const double diff = uniform.energies[i] - uniform.energies[i - 1];
        if (std::abs(diff - actual_dE) > actual_dE * 1e-6) {
            std::cerr
                << "Warning: energy grid is not perfectly uniform at index "
                << i << ". Expected dE=" << actual_dE << ", got " << diff << "."
                << std::endl;
        }
    }

    writeCrossSectionFile(output_file, uniform);

    std::cout << "Input:  " << data.energies.size() << " points, range ["
              << data.energies.front() << ", " << data.energies.back() << "] eV"
              << std::endl;
    std::cout << "Output: " << uniform.energies.size()
              << " points, dE = " << actual_dE << " eV -> " << output_file
              << std::endl;

    return 0;
}
