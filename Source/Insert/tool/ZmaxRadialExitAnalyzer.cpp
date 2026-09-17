/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

/**
 * Analyze the zmax exit radial statistics written by the WarpX
 * zmax_radial_exit_diag diagnostic (zmax_radial/density.dat and
 * zmax_radial/vdist.dat).
 *
 * Usage:
 *     ZmaxRadialExitAnalyzer <zmax_radial_dir> [mass_kg] [output_file]
 *
 * Arguments:
 *     zmax_radial_dir  Directory containing density.dat and vdist.dat
 *     mass_kg          Particle mass (default: 2.18017e-25, Xe atom)
 *     output_file      Output path (default: <zmax_radial_dir>/radial_analysis.dat)
 *
 * For each radial bin the tool reports:
 *   - total exiting weight and macro-particle count (from density.dat);
 *   - mean velocity per component <vx>, <vy>, <vz> [m/s] and the drift
 *     speed magnitude;
 *   - directional temperatures T_x, T_y, T_z = m sigma_i^2 / k_B [K] and
 *     their average T_avg;
 *   - how well each component histogram matches a Maxwellian (Gaussian)
 *     with the same mean and variance:
 *       TV  = 1/2 integral |p - q| dv   (total variation distance, 0 = perfect)
 *       KL  = integral p ln(p/q) dv     (KL divergence in nats, 0 = perfect)
 *     where p is the normalized measured histogram and q the Maxwellian.
 *
 * In addition, weighted (by num_macro) least-squares fits of the radial
 * profiles flux(r), T_x/T_y/T_z(r) and mean_vz(r) are performed with
 * polynomials in u = r/R (degrees 0-3) and, for the flux, the profile
 * families A(1-u^2)^p and A cos(pi u/2)^p. Fit coefficients and errors are
 * printed to stdout and measured-vs-fit values are written to
 * <zmax_radial_dir>/radial_fit.dat.
 *
 * Post-simulation workflow (updating the fit parameters after a DSMC run):
 *   1. Run the channel DSMC simulation (Script/3d_xe_dsmc) to completion.
 *      At the final step it writes the exit statistics directory
 *      zmax_radial/ (density.dat and vdist.dat) in the run directory.
 *   2. Run this tool on that directory:
 *          cd Source/Insert/tool
 *          ./bin/ZmaxRadialExitAnalyzer <run_dir>/zmax_radial
 *      An absolute path bypasses the bin/paths.conf redirection; a relative
 *      path is resolved against the base path configured there (build/bin).
 *      The default particle mass is already the Xe atom, so mass_kg can be
 *      omitted.
 *   3. From stdout, take the poly3 (cubic in u = r/r_inj) coefficients of
 *      flux, T_x, T_y, T_z and mean_vz, and update the 20 my_constants in
 *      the downstream input files (both must be kept in sync):
 *          Script/3d_xe_exit_inlet   (jf_c0..c3, tx_c0..c3, ty_c0..c3,
 *                                     tz_c0..c3, mvz_c0..c3)
 *          Script/3d_hall_exit_inlet (same parameter set)
 *   4. (Optional) Run Script/3d_xe_exit_inlet (1 step, collisions off) to
 *      validate that injection sampling of the new profiles is correct
 *      before using Script/3d_hall_exit_inlet for production runs.
 *   Note: the fit coefficients are condition-dependent (e.g. wall
 *   temperature). After changing the DSMC setup, the coefficients in both
 *   downstream files must be re-fitted from the new run output.
 *
 * Input file formats (tab-separated, '#' comment lines, one column-header
 * row, vdist.dat has one blank line between radial-bin blocks):
 *   density.dat: r_lo_m  r_hi_m  r_center_m  weight  num_macro  flux_per_m2_s
 *   vdist.dat:   r_center_m  v_lo  v_hi  v_center  w_vx  w_vy  w_vz
 */

namespace {

constexpr double k_B = 1.380649e-23; // Boltzmann constant [J/K]

std::string
trim (const std::string& str) {
    auto start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return "";
    }
    auto end = str.find_last_not_of(" \t\r\n");
    return str.substr(start, end - start + 1);
}

// Same paths.conf mechanism as the other tools in this directory: a relative
// zmax_radial_dir argument is resolved against the base path configured in
// bin/paths.conf (located next to the executable).
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

struct DensityRow {
    double r_hi;
    double r_center;
    double weight;
    double num_macro;
    double flux;
};

struct VDistRow {
    double r_center;
    double v_lo;
    double v_hi;
    double v_center;
    double w[3]; // weights of the vx, vy, vz histograms
};

/** Split a line on whitespace; returns false when the line does not hold
 *  exactly n_expected numeric fields (used to skip column-header rows). */
template <size_t N>
bool
parseDoubles (std::string const& line, double (&out)[N]) {
    if (line.empty() || line[0] == '#') {
        return false;
    }
    std::istringstream iss(line);
    for (size_t i = 0; i < N; ++i) {
        if (!(iss >> out[i])) {
            return false;
        }
    }
    return true;
}

std::vector<DensityRow>
readDensity (std::string const& path) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Error: cannot open " << path << "\n";
        std::exit(1);
    }
    std::vector<DensityRow> rows;
    std::string line;
    while (std::getline(in, line)) {
        double f[6];
        if (!parseDoubles(line, f)) {
            continue;
        }
        rows.push_back(DensityRow{f[1], f[2], f[3], f[4], f[5]});
    }
    return rows;
}

std::vector<VDistRow>
readVDist (std::string const& path) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Error: cannot open " << path << "\n";
        std::exit(1);
    }
    std::vector<VDistRow> rows;
    std::string line;
    while (std::getline(in, line)) {
        double f[7];
        if (!parseDoubles(line, f)) {
            continue;
        }
        rows.push_back(VDistRow{f[0], f[1], f[2], f[3], {f[4], f[5], f[6]}});
    }
    return rows;
}

/** Statistics of one velocity-component histogram of one radial bin. */
struct ComponentStats {
    double weight = 0.0;   // total weight
    double mean = 0.0;     // <v> [m/s]
    double sigma = 0.0;    // std deviation [m/s]
    double temperature = 0.0; // m sigma^2 / k_B [K]
    double tv = 0.0;       // total variation distance to Maxwellian
    double kl = 0.0;       // KL divergence to Maxwellian [nats]
};

ComponentStats
analyzeComponent (std::vector<double> const& v_center,
                  std::vector<double> const& w, double const dv,
                  double const mass) {
    ComponentStats s;
    for (double const wi : w) {
        s.weight += wi;
    }
    if (s.weight <= 0.0) {
        return s;
    }
    double m1 = 0.0;
    for (size_t i = 0; i < w.size(); ++i) {
        m1 += v_center[i] * w[i];
    }
    s.mean = m1 / s.weight;
    double m2 = 0.0;
    for (size_t i = 0; i < w.size(); ++i) {
        double const d = v_center[i] - s.mean;
        m2 += d * d * w[i];
    }
    double const var = m2 / s.weight;
    s.sigma = std::sqrt(var);
    s.temperature = mass * var / k_B;
    if (s.sigma <= 0.0) {
        return s;
    }

    // Compare the normalized histogram p with a Gaussian q of the same
    // mean and variance.
    double const inv_norm = 1.0 / (std::sqrt(2.0 * M_PI) * s.sigma);
    for (size_t i = 0; i < w.size(); ++i) {
        double const p = w[i] / (s.weight * dv);
        double const d = (v_center[i] - s.mean) / s.sigma;
        double const q = inv_norm * std::exp(-0.5 * d * d);
        s.tv += 0.5 * std::abs(p - q) * dv;
        if (p > 0.0 && q > 0.0) {
            s.kl += p * std::log(p / q) * dv;
        }
    }
    return s;
}

/** Result of a weighted least-squares fit of a radial profile y(u),
 *  u = r/R in [0,1]. Models:
 *    "polyN": y = c0 + c1 u + ... + cN u^N
 *    "pow":   y = A (1 - u^2)^p      (coeff = {A, p})
 *    "cos":   y = A cos(pi u / 2)^p  (coeff = {A, p})
 */
struct FitResult {
    std::string model;
    std::vector<double> coeff;
    double wrms_rel = 0.0; // weighted RMS relative error
    double max_rel = 0.0;  // maximum relative error
};

double
evalModel (FitResult const& fit, double const u) {
    if (fit.model == "pow") {
        return fit.coeff[0] * std::pow(1.0 - u * u, fit.coeff[1]);
    }
    if (fit.model == "cos") {
        return fit.coeff[0] * std::pow(std::cos(M_PI_2 * u), fit.coeff[1]);
    }
    // Horner evaluation of the polynomial coefficients
    double f = 0.0;
    for (int j = static_cast<int>(fit.coeff.size()) - 1; j >= 0; --j) {
        f = f * u + fit.coeff[static_cast<size_t>(j)];
    }
    return f;
}

void
setFitErrors (FitResult& fit, std::vector<double> const& u,
              std::vector<double> const& y, std::vector<double> const& w) {
    double sse = 0.0;
    double sw = 0.0;
    fit.max_rel = 0.0;
    for (size_t k = 0; k < u.size(); ++k) {
        double const f = evalModel(fit, u[k]);
        double const rel = std::abs(y[k] - f) / std::abs(y[k]);
        sse += w[k] * rel * rel;
        sw += w[k];
        fit.max_rel = std::max(fit.max_rel, rel);
    }
    fit.wrms_rel = std::sqrt(sse / sw);
}

/** Weighted polynomial fit via the normal equations, solved with
 *  Gauss-Jordan elimination and partial pivoting. Returns false when the
 *  system is numerically singular. */
bool
polyFit (std::vector<double> const& u, std::vector<double> const& y,
         std::vector<double> const& w, int const deg,
         std::vector<double>& coeff) {
    int const n = deg + 1;
    std::vector<std::vector<double>> a(static_cast<size_t>(n),
                                       std::vector<double>(n + 1, 0.0));
    for (size_t k = 0; k < u.size(); ++k) {
        std::vector<double> up(2 * deg + 1, 1.0);
        for (int j = 1; j <= 2 * deg; ++j) {
            up[static_cast<size_t>(j)] = up[static_cast<size_t>(j - 1)] * u[k];
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                a[static_cast<size_t>(i)][static_cast<size_t>(j)] +=
                    w[k] * up[static_cast<size_t>(i + j)];
            }
            a[static_cast<size_t>(i)][static_cast<size_t>(n)] +=
                w[k] * y[k] * up[static_cast<size_t>(i)];
        }
    }
    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int row = col + 1; row < n; ++row) {
            if (std::abs(a[static_cast<size_t>(row)][static_cast<size_t>(col)]) >
                std::abs(a[static_cast<size_t>(piv)][static_cast<size_t>(col)])) {
                piv = row;
            }
        }
        if (std::abs(a[static_cast<size_t>(piv)][static_cast<size_t>(col)]) <
            1e-30) {
            return false;
        }
        std::swap(a[static_cast<size_t>(piv)], a[static_cast<size_t>(col)]);
        double const d = a[static_cast<size_t>(col)][static_cast<size_t>(col)];
        for (int j = col; j <= n; ++j) {
            a[static_cast<size_t>(col)][static_cast<size_t>(j)] /= d;
        }
        for (int row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            double const f =
                a[static_cast<size_t>(row)][static_cast<size_t>(col)];
            for (int j = col; j <= n; ++j) {
                a[static_cast<size_t>(row)][static_cast<size_t>(j)] -=
                    f * a[static_cast<size_t>(col)][static_cast<size_t>(j)];
            }
        }
    }
    coeff.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        coeff[static_cast<size_t>(i)] =
            a[static_cast<size_t>(i)][static_cast<size_t>(n)];
    }
    return true;
}

/** Fit y(u) = A g(u; p) with g = (1-u^2)^p ("pow") or cos(pi u/2)^p
 *  ("cos"). p is found by grid search; for fixed p the amplitude is the
 *  analytic weighted least-squares solution A = sum(w y g) / sum(w g^2). */
FitResult
shapeFit (std::string const& model, std::vector<double> const& u,
          std::vector<double> const& y, std::vector<double> const& w) {
    FitResult best;
    best.model = model;
    best.wrms_rel = 1e300;
    for (double p = 0.1; p <= 4.0001; p += 0.02) {
        double num = 0.0;
        double den = 0.0;
        std::vector<double> g(u.size());
        for (size_t k = 0; k < u.size(); ++k) {
            double const base = (model == "pow")
                                    ? (1.0 - u[k] * u[k])
                                    : std::cos(M_PI_2 * u[k]);
            g[k] = std::pow(std::max(base, 0.0), p);
            num += w[k] * y[k] * g[k];
            den += w[k] * g[k] * g[k];
        }
        if (den <= 0.0) {
            continue;
        }
        FitResult cand;
        cand.model = model;
        cand.coeff = {num / den, p};
        setFitErrors(cand, u, y, w);
        if (cand.wrms_rel < best.wrms_rel) {
            best = cand;
        }
    }
    return best;
}

/** Fit all candidate models to the profile y(u) weighted by w and return
 *  them sorted by weighted RMS relative error (best first). */
std::vector<FitResult>
fitAll (std::vector<double> const& u, std::vector<double> const& y,
        std::vector<double> const& w, bool const with_shapes) {
    std::vector<FitResult> results;
    for (int deg = 0; deg <= 3; ++deg) {
        FitResult r;
        r.model = "poly" + std::to_string(deg);
        if (!polyFit(u, y, w, deg, r.coeff)) {
            continue;
        }
        setFitErrors(r, u, y, w);
        results.push_back(r);
    }
    if (with_shapes) {
        results.push_back(shapeFit("pow", u, y, w));
        results.push_back(shapeFit("cos", u, y, w));
    }
    std::sort(results.begin(), results.end(),
              [](FitResult const& a, FitResult const& b) {
                  return a.wrms_rel < b.wrms_rel;
              });
    return results;
}

void
printFit (std::ostream& os, FitResult const& fit) {
    os << fit.model << "  ";
    if (fit.model == "pow") {
        os << "A=" << fit.coeff[0] << " p=" << fit.coeff[1];
    } else if (fit.model == "cos") {
        os << "A=" << fit.coeff[0] << " p=" << fit.coeff[1];
    } else {
        os << "c=";
        for (size_t j = 0; j < fit.coeff.size(); ++j) {
            if (j > 0) {
                os << ",";
            }
            os << fit.coeff[j];
        }
    }
    os << "  wrms_rel=" << fit.wrms_rel * 100.0
       << "%  max_rel=" << fit.max_rel * 100.0 << "%\n";
}

} // namespace

int
main (int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr
            << "Usage: ZmaxRadialExitAnalyzer <zmax_radial_dir> [mass_kg] "
               "[output_file]\n";
        return 1;
    }
    std::string dir = resolvePath(readPathConfig(), argv[1]);
    {   // paths.conf may be stale; fall back to the literal argument when
        // the resolved directory does not contain the expected files.
        std::ifstream probe(dir + "/density.dat");
        if (!probe && dir != argv[1]) {
            dir = argv[1];
        }
    }
    double const mass = (argc > 2) ? std::atof(argv[2]) : 2.18017e-25;
    std::string const out_path =
        (argc > 3) ? resolvePath(readPathConfig(), argv[3])
                   : dir + "/radial_analysis.dat";

    auto const density = readDensity(dir + "/density.dat");
    auto const vdist = readVDist(dir + "/vdist.dat");
    if (vdist.empty()) {
        std::cerr << "Error: no data rows found in " << dir << "/vdist.dat\n";
        return 1;
    }

    // Group vdist rows into radial bins: consecutive rows share r_center,
    // blocks are separated by blank lines in the file.
    struct Bin {
        double r_center;
        double dv;
        std::vector<double> v_center;
        std::vector<double> w[3];
    };
    std::vector<Bin> bins;
    for (auto const& row : vdist) {
        if (bins.empty() ||
            std::abs(row.r_center - bins.back().r_center) > 1e-15) {
            Bin b;
            b.r_center = row.r_center;
            b.dv = row.v_hi - row.v_lo;
            bins.push_back(std::move(b));
        }
        Bin& b = bins.back();
        b.v_center.push_back(row.v_center);
        for (int c = 0; c < 3; ++c) {
            b.w[c].push_back(row.w[c]);
        }
    }

    if (!density.empty() && density.size() != bins.size()) {
        std::cerr << "Warning: density.dat has " << density.size()
                  << " rows but vdist.dat has " << bins.size()
                  << " radial bins; joining by row index anyway.\n";
    }

    std::ofstream out(out_path);
    if (!out) {
        std::cerr << "Error: cannot open " << out_path << " for writing\n";
        return 1;
    }
    out << "# zmax exit radial analysis\n"
        << "# mass_kg " << mass << "\n"
        << "# T_i = m sigma_i^2 / k_B; T_avg = (T_x+T_y+T_z)/3\n"
        << "# TV = total variation distance to Maxwellian (0 = perfect)\n"
        << "# KL = KL divergence to Maxwellian [nats] (0 = perfect)\n"
        << "r_center_m\tweight\tnum_macro\tmean_vx\tmean_vy\tmean_vz\t"
           "drift_speed\tT_x\tT_y\tT_z\tT_avg\t"
           "TV_x\tTV_y\tTV_z\tKL_x\tKL_y\tKL_z\n";

    // Radial profiles collected for the least-squares fits below
    double R = 0.0;
    std::vector<double> v_r, v_nw, v_flux, v_mv[3], v_temp[3];
    for (auto const& d : density) {
        R = std::max(R, d.r_hi);
    }

    out.precision(6);
    for (size_t ib = 0; ib < bins.size(); ++ib) {
        Bin const& b = bins[ib];
        ComponentStats cs[3];
        for (int c = 0; c < 3; ++c) {
            cs[c] = analyzeComponent(b.v_center, b.w[c], b.dv, mass);
        }
        double const weight = cs[0].weight;
        double const num_macro =
            (ib < density.size()) ? density[ib].num_macro : 0.0;
        // Cross-check against density.dat when available. The two files
        // accumulate the same weights through independent atomic sums, so
        // allow for floating-point accumulation noise (observed ~1e-6).
        if (ib < density.size() && density[ib].weight > 0.0 &&
            std::abs(weight - density[ib].weight) /
                    std::max(density[ib].weight, 1.0) >
                1e-4) {
            std::cerr << "Warning: weight mismatch in bin " << ib
                      << ": density.dat " << density[ib].weight
                      << " vs vdist.dat " << weight << "\n";
        }
        double const drift = std::sqrt(cs[0].mean * cs[0].mean +
                                       cs[1].mean * cs[1].mean +
                                       cs[2].mean * cs[2].mean);
        double const t_avg =
            (cs[0].temperature + cs[1].temperature + cs[2].temperature) / 3.0;

        out << b.r_center << "\t" << weight << "\t" << num_macro;
        for (int c = 0; c < 3; ++c) {
            out << "\t" << cs[c].mean;
        }
        out << "\t" << drift;
        for (int c = 0; c < 3; ++c) {
            out << "\t" << cs[c].temperature;
        }
        out << "\t" << t_avg;
        for (int c = 0; c < 3; ++c) {
            out << "\t" << cs[c].tv;
        }
        for (int c = 0; c < 3; ++c) {
            out << "\t" << cs[c].kl;
        }
        out << "\n";

        if (ib < density.size() && density[ib].num_macro > 0.0) {
            v_r.push_back(b.r_center);
            v_nw.push_back(density[ib].num_macro);
            v_flux.push_back(density[ib].flux);
            for (int c = 0; c < 3; ++c) {
                v_mv[c].push_back(cs[c].mean);
                v_temp[c].push_back(cs[c].temperature);
            }
        }
    }
    out.close();

    std::cout << "Wrote " << out_path << " (" << bins.size()
              << " radial bins)\n";

    // ---- Weighted least-squares fits of the radial profiles ----
    if (R <= 0.0 || v_r.empty()) {
        std::cerr << "Warning: no density.dat data, skipping radial fits\n";
        return 0;
    }
    std::vector<double> u(v_r.size());
    for (size_t k = 0; k < v_r.size(); ++k) {
        u[k] = v_r[k] / R;
    }

    std::cout.precision(4);
    std::cout << std::scientific
              << "\nRadial profile fits (u = r/R, R = " << R
              << " m, weighted by num_macro)\n";

    struct QuantityFit {
        std::string name;
        std::vector<double> const* y;
        bool with_shapes;
        FitResult best;
    };
    QuantityFit fits[] = {
        {"flux", &v_flux, true, {}},     {"T_x", &v_temp[0], false, {}},
        {"T_y", &v_temp[1], false, {}},  {"T_z", &v_temp[2], false, {}},
        {"mean_vz", &v_mv[2], false, {}},
    };
    for (auto& q : fits) {
        std::cout << "\n" << q.name << ":\n";
        auto const results = fitAll(u, *q.y, v_nw, q.with_shapes);
        for (auto const& r : results) {
            std::cout << "  ";
            printFit(std::cout, r);
        }
        q.best = results.front();
        std::cout << "  -> best: " << q.best.model << "\n";
    }
    // <vx>, <vy> are expected to vanish by symmetry; relative errors are
    // meaningless, so report the weighted absolute RMS instead.
    for (int c = 0; c < 2; ++c) {
        double sv = 0.0;
        double sw = 0.0;
        for (size_t k = 0; k < u.size(); ++k) {
            sv += v_nw[k] * v_mv[c][k] * v_mv[c][k];
            sw += v_nw[k];
        }
        std::cout << "\nmean_v" << (c == 0 ? "x" : "y")
                  << ": weighted RMS = " << std::sqrt(sv / sw)
                  << " m/s (consistent with 0 by symmetry)\n";
    }

    // Measured vs best-fit values, for plotting
    std::string const fit_path = dir + "/radial_fit.dat";
    std::ofstream fit_out(fit_path);
    if (fit_out) {
        fit_out << "# zmax exit radial fits (u = r/R, R = " << R << " m)\n"
                << "# models:";
        for (auto const& q : fits) {
            fit_out << " " << q.name << "=" << q.best.model;
        }
        fit_out << "\n"
                << "r_center_m\tflux\tflux_fit\tT_x\tT_x_fit\tT_y\tT_y_fit\t"
                   "T_z\tT_z_fit\tmean_vz\tmean_vz_fit\n";
        fit_out.precision(6);
        for (size_t k = 0; k < u.size(); ++k) {
            fit_out << v_r[k];
            for (auto const& q : fits) {
                fit_out << "\t" << (*q.y)[k] << "\t"
                        << evalModel(q.best, u[k]);
            }
            fit_out << "\n";
        }
        std::cout << "\nWrote " << fit_path << "\n";
    }
    return 0;
}
