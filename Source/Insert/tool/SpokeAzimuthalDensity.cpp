/* Copyright 2025 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

/**
 * Export normalized azimuthal density profiles of the spoke distributions
 * used by the Hall injection framework, for 1-4 spokes, as an ASCII table
 * suitable for plotting (e.g. in Origin).
 *
 * Usage:
 *     SpokeAzimuthalDensity [output_file] [num_points]
 *
 * Arguments:
 *     output_file  Output path (default: spoke_azimuthal_density.dat)
 *     num_points   Number of theta samples over [0, 360) deg (default: 3601)
 *
 * The profiles replicate the exact formulas in
 * Source/Insert/Injection/HallDistribution1D.cpp:
 *     - neutral_spoke: azimuthal depletion profile of the neutral gas
 *     - multi_spoke:   superposition of spoke_count Gaussian peaks (plasma)
 *
 * Default parameters match the current spoke setup in Script/3d_hall_spoke:
 *     spoke_phase        = 20 deg   (neutral density peak / ionization front)
 *     spoke_ion_width    = 0.30*pi  (depletion width of the neutral profile)
 *     spoke_min_ratio    = 1/6      (density ratio at the depletion minimum)
 *     spoke_drop_exponent= 4.0
 *     spoke_reverse      = 1        (clockwise recovery, opposite to CCW
 *                                    propagation)
 *     spoke_sigma        = pi/8     (Gaussian width of each plasma peak)
 *     spoke_plasma_shift = 15 deg   (plasma peak shifted from the neutral
 *                                    minimum toward the neutral peak)
 *     plasma phase       = spoke_phase - spoke_ion_width + spoke_plasma_shift
 *
 * For plotting, the phase is re-centered so that the neutral peak and the
 * plasma peak straddle the middle of the [0, 360] deg domain: the plasma
 * peak sits 39 deg behind the neutral peak, so the neutral peak is placed
 * at 199.5 deg and the plasma peak at 160.5 deg. This is a pure phase
 * shift and does not change the profile shapes.
 *
 * Each curve is normalized by its own maximum so that the peak value is 1.
 */

namespace {

constexpr double pi = 3.14159265358979323846264338327950288;
constexpr double two_pi = 2.0 * pi;

// Current spoke parameters (Script/3d_hall_spoke, Script/3d_hall_spoke_init)
constexpr double spoke_ion_width = 0.30 * pi;
constexpr double spoke_min_ratio = 1.0 / 6.0;
constexpr double spoke_drop_exponent = 4.0;
constexpr bool spoke_reverse = true;
constexpr double spoke_sigma = pi / 8.0;
constexpr double spoke_plasma_shift = pi / 12.0; // 15 deg
constexpr int max_spoke_count = 4;

// The physical phase is 20 deg (pi/9). For plotting, the phase is re-centered
// so that the neutral peak and the plasma peak straddle the middle of the
// plot domain (180 deg): the plasma peak sits (ion_width - shift) = 39 deg
// behind the neutral peak, so the neutral peak is placed 19.5 deg past 180.
constexpr double spoke_phase =
    pi + 0.5 * (spoke_ion_width - spoke_plasma_shift);
constexpr double plasma_phase =
    spoke_phase - spoke_ion_width + spoke_plasma_shift;

[[nodiscard]]
double
periodicDistance (double x, double center) noexcept
{
    double distance = std::fmod(x - center + pi, two_pi);
    if (distance < 0.0) {
        distance += two_pi;
    }
    return distance - pi;
}

[[nodiscard]]
double
periodicPhase (double theta, double phase, bool reverse) noexcept
{
    double local_phase = reverse ? phase - theta : theta - phase;
    local_phase = std::fmod(local_phase, two_pi);
    if (local_phase < 0.0) {
        local_phase += two_pi;
    }
    return local_phase;
}

[[nodiscard]]
double
gaussianDensity (double x, double center, double sigma) noexcept
{
    const double normalized = periodicDistance(x, center) / sigma;
    return std::exp(-0.5 * normalized * normalized);
}

// Mirrors NeutralSpokeDensity() in HallDistribution1D.cpp
[[nodiscard]]
double
neutralSpokeDensity (double theta, double ion_width, double min_ratio,
                     double drop_exponent, double phase, bool reverse,
                     int spoke_count) noexcept
{
    const double period = two_pi / static_cast<double>(spoke_count);
    double phi = periodicPhase(theta, phase, reverse);
    // Fold the phase into a single structure period so the depletion
    // profile repeats spoke_count times around the circumference.
    phi = std::fmod(phi, period);
    if (phi < ion_width) {
        const double s = phi / ion_width;
        return min_ratio + (1.0 - min_ratio) * std::pow(1.0 - s, drop_exponent);
    }

    const double s = (phi - ion_width) / (period - ion_width);
    return min_ratio + (1.0 - min_ratio) * s;
}

// Mirrors the pdf lambda of MakeMultiSpokeSampler() in HallDistribution1D.cpp
[[nodiscard]]
double
multiSpokeDensity (double theta, int spoke_count, double sigma,
                   double phase) noexcept
{
    double density = 0.0;
    const double interval = two_pi / static_cast<double>(spoke_count);
    for (int i = 0; i < spoke_count; ++i) {
        const double center = phase + static_cast<double>(i) * interval;
        density += gaussianDensity(theta, center, sigma);
    }
    return density;
}

} // namespace

int
main (int argc, char* argv[])
{
    const std::string output_file =
        (argc > 1) ? argv[1] : "spoke_azimuthal_density.dat";
    const int num_points = (argc > 2) ? std::atoi(argv[2]) : 3601;
    if (num_points < 2) {
        std::cerr << "num_points must be >= 2\n";
        return 1;
    }

    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Cannot open output file: " << output_file << "\n";
        return 1;
    }

    // Header: column names as a comment (Origin imports them as Long Name)
    out << "# theta[deg]";
    for (int n = 1; n <= max_spoke_count; ++n) {
        out << "\tneutral_n" << n << "\tplasma_n" << n;
    }
    out << "\n";

    // multi_spoke may exceed 1 from peak overlap; normalize each curve to
    // its own maximum (neutral_spoke already peaks at 1 by construction)
    double plasma_max[max_spoke_count + 1];
    for (int n = 1; n <= max_spoke_count; ++n) {
        plasma_max[n] = 0.0;
        for (int j = 0; j < num_points; ++j) {
            const double tj =
                two_pi * static_cast<double>(j) /
                static_cast<double>(num_points - 1);
            plasma_max[n] = std::max(
                plasma_max[n], multiSpokeDensity(tj, n, spoke_sigma,
                                                 plasma_phase));
        }
    }

    out << std::scientific << std::setprecision(8);
    for (int i = 0; i < num_points; ++i) {
        const double theta_deg =
            360.0 * static_cast<double>(i) / static_cast<double>(num_points - 1);
        const double theta = theta_deg * pi / 180.0;
        out << theta_deg;
        for (int n = 1; n <= max_spoke_count; ++n) {
            const double neutral = neutralSpokeDensity(
                theta, spoke_ion_width, spoke_min_ratio, spoke_drop_exponent,
                spoke_phase, spoke_reverse, n);
            const double plasma =
                multiSpokeDensity(theta, n, spoke_sigma, plasma_phase) /
                plasma_max[n];
            out << "\t" << neutral << "\t" << plasma;
        }
        out << "\n";
    }

    std::cout << "Wrote " << num_points << " points to " << output_file << "\n";
    return 0;
}
