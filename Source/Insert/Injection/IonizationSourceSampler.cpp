#include "IonizationSourceSampler.h"

#include "Insert/Utils/InsertUtils.h"
#include "Utils/TextMsg.H"

#include <AMReX_Math.H>
#include <AMReX_Random.H>

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <fstream>
#include <sstream>
#include <utility>

namespace Insert {
namespace {

constexpr amrex::Real Pi =
    amrex::Real(3.141592653589793238462643383279502884);

struct ArrayData
{
    int nrow = 0;
    int ncol = 0;
    std::vector<amrex::Real> values;
};

ArrayData
ReadArray (std::string const& path)
{
    std::ifstream is(path);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        is.good(), "Could not open ionization source array: " + path);

    ArrayData data;
    is >> data.nrow >> data.ncol;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        data.nrow > 0 && data.ncol > 0,
        "Invalid ionization source array shape in: " + path);

    data.values.resize(static_cast<std::size_t>(data.nrow * data.ncol));
    for (auto& value : data.values) {
        is >> value;
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        is.good() || is.eof(),
        "Could not read complete ionization source array: " + path);
    return data;
}

void
RequireEqual (int actual, int expected, std::string const& description)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        actual == expected,
        "Unexpected ionization source " + description + ": got " +
            std::to_string(actual) + ", expected " + std::to_string(expected));
}

amrex::Real
InvertQuadraticCdf (
    amrex::Real target, amrex::Real p0, amrex::Real p1,
    amrex::Real norm)
{
    if (norm <= 0.0) {
        return target;
    }

    amrex::Real lo = 0.0;
    amrex::Real hi = 1.0;
    for (int iter = 0; iter < 12; ++iter) {
        const amrex::Real mid = 0.5 * (lo + hi);
        const amrex::Real cdf = (p0 * mid + 0.5 * p1 * mid * mid) / norm;
        if (cdf < target) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

amrex::Real
InvertCubicCdf (
    amrex::Real target, amrex::Real q0, amrex::Real q1,
    amrex::Real q2, amrex::Real norm)
{
    if (norm <= 0.0) {
        return target;
    }

    amrex::Real lo = 0.0;
    amrex::Real hi = 1.0;
    for (int iter = 0; iter < 12; ++iter) {
        const amrex::Real mid = 0.5 * (lo + hi);
        const amrex::Real cdf =
            (q0 * mid + 0.5 * q1 * mid * mid +
             (amrex::Real(1.0) / amrex::Real(3.0)) * q2 * mid * mid * mid) /
            norm;
        if (cdf < target) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

} // namespace

IonizationSourceSampler::IonizationSourceSampler (std::string source_directory)
    : m_source_directory(std::move(source_directory))
{
    std::ifstream metadata(PathJoin(m_source_directory, "metadata.txt"));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        metadata.good(),
        "Could not open ionization source metadata in " + m_source_directory);

    std::string line;
    while (std::getline(metadata, line)) {
        std::istringstream is(line);
        std::string key;
        std::string equals;
        is >> key >> equals;
        if (key.empty()) {
            continue;
        }
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            equals == "=",
            "Invalid ionization source metadata line: " + line);

        if (key == "center") {
            is >> m_x_center >> m_y_center;
        } else if (key == "r_min") {
            is >> m_r_min;
        } else if (key == "r_max") {
            is >> m_r_max;
        } else if (key == "z_min") {
            is >> m_z_min;
        } else if (key == "z_max") {
            is >> m_z_max;
        } else if (key == "nr") {
            is >> m_nr;
        } else if (key == "nz") {
            is >> m_nz;
        } else if (key == "A_tot") {
            is >> m_total_rate;
        } else if (key == "elec_weight") {
            is >> m_elec_weight;
        }
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_nr > 0 && m_nz > 0,
        "Ionization source metadata requires positive nr and nz.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_r_max > m_r_min && m_z_max > m_z_min,
        "Ionization source metadata requires non-empty r-z bounds.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_total_rate >= 0.0,
        "Ionization source metadata requires non-negative A_tot.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_elec_weight > amrex::ParticleReal(0.0),
        "Ionization source metadata requires positive elec_weight.");

    m_dr = (m_r_max - m_r_min) / static_cast<amrex::Real>(m_nr);
    m_dz = (m_z_max - m_z_min) / static_cast<amrex::Real>(m_nz);

    auto node_rate = ReadArray(PathJoin(m_source_directory, "node_source_rate"));
    RequireEqual(node_rate.nrow, m_nz + 1, "node_source_rate row count");
    RequireEqual(node_rate.ncol, m_nr + 1, "node_source_rate column count");
    m_node_rate = std::move(node_rate.values);

    auto cell_cdf = ReadArray(PathJoin(m_source_directory, "cell_cdf"));
    RequireEqual(cell_cdf.nrow, m_nz, "cell_cdf row count");
    RequireEqual(cell_cdf.ncol, m_nr, "cell_cdf column count");
    m_cell_cdf = std::move(cell_cdf.values);

    if (m_total_rate > 0.0) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            !m_cell_cdf.empty() && m_cell_cdf.back() > 0.0,
            "Ionization source with positive A_tot requires a non-empty CDF.");
    }
}

EmissionSample
IonizationSourceSampler::samplePosition (
    amrex::RandomEngine const& engine) const
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_total_rate > 0.0,
        "Cannot sample from an ionization source with zero total rate.");

    const amrex::Real u0 = amrex::Random(engine);
    const amrex::Real u1 = amrex::Random(engine);
    const amrex::Real u2 = amrex::Random(engine);
    const amrex::Real u3 = amrex::Random(engine);

    const auto cell_iter =
        std::upper_bound(m_cell_cdf.begin(), m_cell_cdf.end(), u0);
    const int flat_index = static_cast<int>(
        std::min<std::ptrdiff_t>(
            std::distance(m_cell_cdf.begin(), cell_iter),
            static_cast<std::ptrdiff_t>(m_cell_cdf.size() - 1)));
    const int iz = flat_index / m_nr;
    const int ir = flat_index - iz * m_nr;

    const amrex::Real s00 = nodeRate(iz, ir);
    const amrex::Real s10 = nodeRate(iz + 1, ir);
    const amrex::Real s01 = nodeRate(iz, ir + 1);
    const amrex::Real s11 = nodeRate(iz + 1, ir + 1);
    const amrex::Real a = s00;
    const amrex::Real b = s10 - s00;
    const amrex::Real c = s01 - s00;
    const amrex::Real d = s11 - s10 - s01 + s00;

    const amrex::Real alpha = a + 0.5 * b;
    const amrex::Real beta = c + 0.5 * d;
    const amrex::Real r0 = m_r_min + static_cast<amrex::Real>(ir) * m_dr;
    const amrex::Real q0 = alpha * r0;
    const amrex::Real q1 = alpha * m_dr + beta * r0;
    const amrex::Real q2 = beta * m_dr;
    const amrex::Real radial_norm =
        q0 + 0.5 * q1 + (amrex::Real(1.0) / amrex::Real(3.0)) * q2;
    const amrex::Real eta = InvertCubicCdf(u1, q0, q1, q2, radial_norm);

    const amrex::Real p0 = a + c * eta;
    const amrex::Real p1 = b + d * eta;
    const amrex::Real axial_norm = p0 + 0.5 * p1;
    const amrex::Real xi = InvertQuadraticCdf(u2, p0, p1, axial_norm);

    const amrex::Real r = r0 + eta * m_dr;
    const amrex::Real theta = 2.0 * Pi * u3;

    EmissionSample sample;
    sample.x = static_cast<amrex::ParticleReal>(m_x_center + r * std::cos(theta));
    sample.y = static_cast<amrex::ParticleReal>(m_y_center + r * std::sin(theta));
    sample.z = static_cast<amrex::ParticleReal>(
        m_z_min + (static_cast<amrex::Real>(iz) + xi) * m_dz);
    return sample;
}

amrex::Real
IonizationSourceSampler::totalRate () const noexcept
{
    return m_total_rate;
}

amrex::ParticleReal
IonizationSourceSampler::electronWeight () const noexcept
{
    return m_elec_weight;
}

amrex::Real
IonizationSourceSampler::nodeRate (int iz, int ir) const
{
    return m_node_rate[static_cast<std::size_t>(nodeIndex(iz, ir))];
}

int
IonizationSourceSampler::nodeIndex (int iz, int ir) const noexcept
{
    return iz * (m_nr + 1) + ir;
}

} // namespace Insert
