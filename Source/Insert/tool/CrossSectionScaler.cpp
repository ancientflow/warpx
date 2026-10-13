/* Copyright 2024-2025 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

/**
 * Scale cross-section values in a data file by a constant factor.
 *
 * Usage:
 *     CrossSectionScaler <input_file> <output_file> <scale_factor>
 *
 * Arguments:
 *     input_file   Path to the cross-section file (two columns: energy[eV]
 * sigma[m^2]) output_file  Path to the scaled output file scale_factor
 * Multiplicative factor applied to every sigma value
 *
 * The input file format matches WarpX's
 * ScatteringProcess::readCrossSectionFile:
 *     - Each line contains two whitespace-separated values: energy and
 * cross-section
 *     - Lines starting with '#' are treated as comments and ignored
 *
 * Energy values are preserved; only cross-sections are scaled.
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

    if (data.energies.size() < 1) {
        std::cerr << "Error: at least one data point is required." << std::endl;
        std::exit(1);
    }

    return data;
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
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_file> <output_file> <scale_factor>" << std::endl;
        return 1;
    }

    std::string base_dir = readPathConfig();

    const std::string input_file = resolvePath(base_dir, argv[1]);
    const std::string output_file = resolvePath(base_dir, argv[2]);
    const double scale_factor = std::stod(argv[3]);

    CrossSectionData data = readCrossSectionFile(input_file);

    for (auto& sigma : data.sigmas) {
        sigma *= scale_factor;
    }

    writeCrossSectionFile(output_file, data);

    std::cout << "Input:  " << data.energies.size() << " points from "
              << input_file << std::endl;
    std::cout << "Output: " << data.energies.size() << " points to "
              << output_file << std::endl;
    std::cout << "Scale factor: " << scale_factor << std::endl;

    return 0;
}
