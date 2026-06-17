#include "HallDistribution1D.h"

#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"

#include <AMReX.H>
#include <AMReX_Math.H>
#include <AMReX_Parser.H>
#include <AMReX_Random.H>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <utility>

namespace Insert {
namespace {

std::string
ToLower (std::string value)
{
    std::transform(
        value.begin(), value.end(), value.begin(),
        [] (unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

amrex::ParticleReal
TwoPi ()
{
    return amrex::ParticleReal(2.0) * amrex::Math::pi<amrex::ParticleReal>();
}

amrex::ParticleReal
PeriodicDistance (amrex::ParticleReal x, amrex::ParticleReal center) noexcept
{
    const amrex::ParticleReal pi = amrex::Math::pi<amrex::ParticleReal>();
    amrex::ParticleReal distance = std::fmod(x - center + pi, TwoPi());
    if (distance < amrex::ParticleReal(0.0)) {
        distance += TwoPi();
    }
    return distance - pi;
}

amrex::ParticleReal
PeriodicPhase (amrex::ParticleReal theta, amrex::ParticleReal phase,
               bool reverse) noexcept
{
    amrex::ParticleReal local_phase = reverse ? phase - theta : theta - phase;
    local_phase = std::fmod(local_phase, TwoPi());
    if (local_phase < amrex::ParticleReal(0.0)) {
        local_phase += TwoPi();
    }
    return local_phase;
}

amrex::ParticleReal
GaussianDensity (amrex::ParticleReal x, amrex::ParticleReal center,
                 amrex::ParticleReal sigma) noexcept
{
    const amrex::ParticleReal normalized = PeriodicDistance(x, center) / sigma;
    return std::exp(amrex::ParticleReal(-0.5) * normalized * normalized);
}

amrex::ParticleReal
NeutralSpokeDensity (amrex::ParticleReal theta, amrex::ParticleReal ion_width,
                     amrex::ParticleReal min_ratio,
                     amrex::ParticleReal drop_exponent,
                     amrex::ParticleReal phase, bool reverse) noexcept
{
    const amrex::ParticleReal phi = PeriodicPhase(theta, phase, reverse);
    if (phi < ion_width) {
        const amrex::ParticleReal s = phi / ion_width;
        return min_ratio + (amrex::ParticleReal(1.0) - min_ratio) *
                               std::pow(amrex::ParticleReal(1.0) - s,
                                        drop_exponent);
    }

    const amrex::ParticleReal s = (phi - ion_width) / (TwoPi() - ion_width);
    return min_ratio + (amrex::ParticleReal(1.0) - min_ratio) * s;
}

amrex::ParticleReal
InterpolatePdf (std::vector<amrex::ParticleReal> const& values,
                std::vector<amrex::ParticleReal> const& pdf,
                amrex::ParticleReal x) noexcept
{
    if (x <= values.front()) {
        return pdf.front();
    }
    if (x >= values.back()) {
        return pdf.back();
    }

    auto const upper = std::upper_bound(values.begin(), values.end(), x);
    const auto hi = static_cast<std::size_t>(upper - values.begin());
    const auto lo = hi - 1;
    const amrex::ParticleReal fraction =
        (x - values[lo]) / (values[hi] - values[lo]);
    return pdf[lo] + fraction * (pdf[hi] - pdf[lo]);
}

void ValidatePositiveSigma (amrex::ParticleReal sigma, std::string const& name);

void
ValidateNeutralSpokeParameters (amrex::ParticleReal ion_width,
                                amrex::ParticleReal min_ratio,
                                amrex::ParticleReal drop_exponent)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        ion_width > amrex::ParticleReal(0.0) && ion_width < TwoPi(),
        "HallNeutralSpokeDistribution1D requires 0 < ion_width < 2*pi.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        min_ratio >= amrex::ParticleReal(0.0) &&
            min_ratio < amrex::ParticleReal(1.0),
        "HallNeutralSpokeDistribution1D requires 0 <= min_ratio < 1.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        drop_exponent > amrex::ParticleReal(0.0),
        "HallNeutralSpokeDistribution1D requires drop_exponent > 0.");
}

NumericalInverseCDFSampler1D
MakeSingleSpokeSampler (amrex::ParticleReal center, amrex::ParticleReal sigma,
                        int num_bins)
{
    ValidatePositiveSigma(sigma, "HallSingleSpokeDistribution1D");
    return NumericalInverseCDFSampler1D(
        amrex::ParticleReal(0.0), TwoPi(), num_bins,
        [=] (amrex::ParticleReal theta) {
            return GaussianDensity(theta, center, sigma);
        });
}

NumericalInverseCDFSampler1D
MakeMultiSpokeSampler (int spoke_count, amrex::ParticleReal sigma,
                       amrex::ParticleReal phase, int num_bins)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        spoke_count > 0,
        "HallMultiSpokeDistribution1D requires spoke_count > 0.");
    ValidatePositiveSigma(sigma, "HallMultiSpokeDistribution1D");
    return NumericalInverseCDFSampler1D(
        amrex::ParticleReal(0.0), TwoPi(), num_bins,
        [=] (amrex::ParticleReal theta) {
            amrex::ParticleReal density = amrex::ParticleReal(0.0);
            const amrex::ParticleReal interval =
                TwoPi() / static_cast<amrex::ParticleReal>(spoke_count);
            for (int i = 0; i < spoke_count; ++i) {
                const auto center =
                    phase + static_cast<amrex::ParticleReal>(i) * interval;
                density += GaussianDensity(theta, center, sigma);
            }
            return density;
        });
}

NumericalInverseCDFSampler1D
MakeNeutralSpokeSampler (amrex::ParticleReal ion_width,
                         amrex::ParticleReal min_ratio,
                         amrex::ParticleReal drop_exponent,
                         amrex::ParticleReal phase, bool reverse, int num_bins)
{
    ValidateNeutralSpokeParameters(ion_width, min_ratio, drop_exponent);
    return NumericalInverseCDFSampler1D(
        amrex::ParticleReal(0.0), TwoPi(), num_bins,
        [=] (amrex::ParticleReal theta) {
            return NeutralSpokeDensity(theta, ion_width, min_ratio,
                                       drop_exponent, phase, reverse);
        });
}

NumericalInverseCDFSampler1D
MakeTabulatedSampler (std::vector<amrex::ParticleReal> const& values,
                      std::vector<amrex::ParticleReal> const& pdf,
                      int num_bins)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(values.size() == pdf.size() &&
                                         values.size() >= 2,
                                     "HallTabulatedDistribution1D requires "
                                     "values and pdf arrays of size >= 2.");
    for (std::size_t i = 1; i < values.size(); ++i) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            values[i] > values[i - 1],
            "HallTabulatedDistribution1D values must be strictly increasing.");
    }

    return NumericalInverseCDFSampler1D(values.front(), values.back(), num_bins,
                                        [values, pdf] (amrex::ParticleReal x) {
                                            return InterpolatePdf(values, pdf,
                                                                  x);
                                        });
}

NumericalInverseCDFSampler1D
MakeParserSampler (amrex::ParticleReal min, amrex::ParticleReal max,
                   int num_bins, std::string const& variable_name,
                   std::string const& expression)
{
    auto parser = utils::parser::makeParser(expression, {variable_name});
    auto const executor = parser.compile<1>();
    return NumericalInverseCDFSampler1D(
        min, max, num_bins,
        [executor] (amrex::ParticleReal x) { return executor(x); });
}

void
ValidateRange (amrex::ParticleReal min, amrex::ParticleReal max,
               std::string const& name)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(max > min, name + " requires max > min.");
}

void
ValidateSigma (amrex::ParticleReal sigma, std::string const& name)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(sigma >= amrex::ParticleReal(0.0),
                                     name + " requires sigma >= 0.");
}

void
ValidatePositiveSigma (amrex::ParticleReal sigma, std::string const& name)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(sigma > amrex::ParticleReal(0.0),
                                     name + " requires sigma > 0.");
}

amrex::ParticleReal
GetReal (amrex::ParmParse const& pp, std::string const& prefix,
         char const* name)
{
    amrex::ParticleReal value = amrex::ParticleReal(0.0);
    utils::parser::getWithParser(pp, prefix, name, value);
    return value;
}

amrex::ParticleReal
QueryReal (amrex::ParmParse const& pp, std::string const& prefix,
           char const* name, amrex::ParticleReal default_value)
{
    auto value = default_value;
    utils::parser::queryWithParser(pp, prefix, name, value);
    return value;
}

int
QueryInt (amrex::ParmParse const& pp, std::string const& prefix,
          char const* name, int default_value)
{
    auto value = default_value;
    utils::parser::queryWithParser(pp, prefix, name, value);
    return value;
}

bool
QueryBoolWithAliases (amrex::ParmParse const& pp, std::string const& prefix,
                      std::initializer_list<char const*> names,
                      bool default_value)
{
    bool value = default_value;
    for (char const* name : names) {
        if (pp.query(prefix + "." + name, value)) {
            return value;
        }
    }
    return value;
}

amrex::ParticleReal
GetSigmaWithFallback (amrex::ParmParse const& pp, std::string const& prefix)
{
    amrex::ParticleReal sigma = amrex::ParticleReal(0.0);
    if (utils::parser::queryWithParser(pp, prefix, "sigma", sigma)) {
        return sigma;
    }
    utils::parser::getWithParser(pp, prefix, "spoke_sigma", sigma);
    return sigma;
}

amrex::ParticleReal
QueryPhaseWithFallback (amrex::ParmParse const& pp, std::string const& prefix)
{
    amrex::ParticleReal phase = amrex::ParticleReal(0.0);
    if (utils::parser::queryWithParser(pp, prefix, "phase", phase)) {
        return phase;
    }
    if (utils::parser::queryWithParser(pp, prefix, "phase_offset", phase)) {
        return phase;
    }
    utils::parser::queryWithParser(pp, prefix, "spoke_phase", phase);
    return phase;
}

std::string
AxisNameFromPrefix (std::string const& prefix)
{
    const auto pos = prefix.find_last_of('.');
    if (pos == std::string::npos) {
        return prefix;
    }
    return prefix.substr(pos + 1);
}

} // namespace

HallConstantDistribution1D::HallConstantDistribution1D (
    amrex::ParticleReal value) noexcept
    : m_value(value)
{}

amrex::ParticleReal
HallConstantDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    amrex::ignore_unused(engine);
    return m_value;
}

amrex::ParticleReal
HallConstantDistribution1D::min () const noexcept
{
    return m_value;
}

amrex::ParticleReal
HallConstantDistribution1D::max () const noexcept
{
    return m_value;
}

HallUniformDistribution1D::HallUniformDistribution1D (amrex::ParticleReal min,
                                                      amrex::ParticleReal max)
    : m_min(min), m_max(max)
{
    ValidateRange(m_min, m_max, "HallUniformDistribution1D");
}

amrex::ParticleReal
HallUniformDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    return m_min + (m_max - m_min) * amrex::Random(engine);
}

amrex::ParticleReal
HallUniformDistribution1D::min () const noexcept
{
    return m_min;
}

amrex::ParticleReal
HallUniformDistribution1D::max () const noexcept
{
    return m_max;
}

HallAreaUniformDistribution1D::HallAreaUniformDistribution1D (
    amrex::ParticleReal min, amrex::ParticleReal max)
    : m_min(min), m_max(max)
{
    ValidateRange(m_min, m_max, "HallAreaUniformDistribution1D");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_min >= amrex::ParticleReal(0.0),
        "HallAreaUniformDistribution1D requires min >= 0.");
}

amrex::ParticleReal
HallAreaUniformDistribution1D::sample (
    amrex::RandomEngine const& engine) const
{
    const auto r2_min = m_min * m_min;
    const auto r2_max = m_max * m_max;
    return std::sqrt(r2_min + (r2_max - r2_min) * amrex::Random(engine));
}

amrex::ParticleReal
HallAreaUniformDistribution1D::min () const noexcept
{
    return m_min;
}

amrex::ParticleReal
HallAreaUniformDistribution1D::max () const noexcept
{
    return m_max;
}

HallGaussianDistribution1D::HallGaussianDistribution1D (
    amrex::ParticleReal mean, amrex::ParticleReal sigma)
    : m_mean(mean), m_sigma(sigma)
{
    ValidateSigma(m_sigma, "HallGaussianDistribution1D");
}

amrex::ParticleReal
HallGaussianDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    if (m_sigma == amrex::ParticleReal(0.0)) {
        return m_mean;
    }
    return amrex::RandomNormal(m_mean, m_sigma, engine);
}

amrex::ParticleReal
HallGaussianDistribution1D::min () const noexcept
{
    return std::numeric_limits<amrex::ParticleReal>::lowest();
}

amrex::ParticleReal
HallGaussianDistribution1D::max () const noexcept
{
    return std::numeric_limits<amrex::ParticleReal>::max();
}

HallPositiveGaussianDistribution1D::HallPositiveGaussianDistribution1D (
    amrex::ParticleReal mean, amrex::ParticleReal sigma)
    : m_mean(mean), m_sigma(sigma)
{
    ValidateSigma(m_sigma, "HallPositiveGaussianDistribution1D");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE((m_sigma > amrex::ParticleReal(0.0)) ||
                                         (m_mean >= amrex::ParticleReal(0.0)),
                                     "HallPositiveGaussianDistribution1D with "
                                     "sigma == 0 requires mean >= 0.");
}

amrex::ParticleReal
HallPositiveGaussianDistribution1D::sample (
    amrex::RandomEngine const& engine) const
{
    if (m_sigma == amrex::ParticleReal(0.0)) {
        return m_mean;
    }

    amrex::ParticleReal value;
    do {
        value = amrex::RandomNormal(m_mean, m_sigma, engine);
    } while (value < amrex::ParticleReal(0.0));
    return value;
}

amrex::ParticleReal
HallPositiveGaussianDistribution1D::min () const noexcept
{
    return amrex::ParticleReal(0.0);
}

amrex::ParticleReal
HallPositiveGaussianDistribution1D::max () const noexcept
{
    return std::numeric_limits<amrex::ParticleReal>::max();
}

HallSingleSpokeDistribution1D::HallSingleSpokeDistribution1D (
    amrex::ParticleReal center, amrex::ParticleReal sigma, int num_bins)
    : m_sampler(MakeSingleSpokeSampler(center, sigma, num_bins))
{}

amrex::ParticleReal
HallSingleSpokeDistribution1D::sample (
    amrex::RandomEngine const& engine) const
{
    return m_sampler.sample(engine);
}

amrex::ParticleReal
HallSingleSpokeDistribution1D::min () const noexcept
{
    return m_sampler.xmin();
}

amrex::ParticleReal
HallSingleSpokeDistribution1D::max () const noexcept
{
    return m_sampler.xmax();
}

amrex::ParticleReal
HallSingleSpokeDistribution1D::integral () const noexcept
{
    return m_sampler.integral();
}

HallMultiSpokeDistribution1D::HallMultiSpokeDistribution1D (
    int spoke_count, amrex::ParticleReal sigma, amrex::ParticleReal phase,
    int num_bins)
    : m_sampler(MakeMultiSpokeSampler(spoke_count, sigma, phase, num_bins))
{}

amrex::ParticleReal
HallMultiSpokeDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    return m_sampler.sample(engine);
}

amrex::ParticleReal
HallMultiSpokeDistribution1D::min () const noexcept
{
    return m_sampler.xmin();
}

amrex::ParticleReal
HallMultiSpokeDistribution1D::max () const noexcept
{
    return m_sampler.xmax();
}

amrex::ParticleReal
HallMultiSpokeDistribution1D::integral () const noexcept
{
    return m_sampler.integral();
}

HallNeutralSpokeDistribution1D::HallNeutralSpokeDistribution1D (
    amrex::ParticleReal ion_width, amrex::ParticleReal min_ratio,
    amrex::ParticleReal drop_exponent, amrex::ParticleReal phase,
    bool reverse, int num_bins)
    : m_sampler(MakeNeutralSpokeSampler(
          ion_width, min_ratio, drop_exponent, phase, reverse, num_bins))
{}

amrex::ParticleReal
HallNeutralSpokeDistribution1D::sample (
    amrex::RandomEngine const& engine) const
{
    return m_sampler.sample(engine);
}

amrex::ParticleReal
HallNeutralSpokeDistribution1D::min () const noexcept
{
    return m_sampler.xmin();
}

amrex::ParticleReal
HallNeutralSpokeDistribution1D::max () const noexcept
{
    return m_sampler.xmax();
}

amrex::ParticleReal
HallNeutralSpokeDistribution1D::integral () const noexcept
{
    return m_sampler.integral();
}

HallDiscreteDistribution1D::HallDiscreteDistribution1D (
    std::vector<amrex::ParticleReal> values,
    std::vector<amrex::ParticleReal> weights)
    : m_values(std::move(values))
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_values.empty(), "HallDiscreteDistribution1D requires values.");
    if (weights.empty()) {
        weights.assign(m_values.size(), amrex::ParticleReal(1.0));
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(weights.size() == m_values.size(),
                                     "HallDiscreteDistribution1D requires "
                                     "values and weights with equal sizes.");

    m_min = *std::min_element(m_values.begin(), m_values.end());
    m_max = *std::max_element(m_values.begin(), m_values.end());
    m_cdf.assign(weights.size(), amrex::ParticleReal(0.0));

    amrex::ParticleReal cumulative = amrex::ParticleReal(0.0);
    for (std::size_t i = 0; i < weights.size(); ++i) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            weights[i] >= amrex::ParticleReal(0.0),
            "HallDiscreteDistribution1D weights must be non-negative.");
        cumulative += weights[i];
        m_cdf[i] = cumulative;
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        cumulative > amrex::ParticleReal(0.0),
        "HallDiscreteDistribution1D requires a positive total weight.");
    m_integral = cumulative;
    for (auto& value : m_cdf) {
        value /= cumulative;
    }
    m_cdf.back() = amrex::ParticleReal(1.0);
}

amrex::ParticleReal
HallDiscreteDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    const amrex::ParticleReal selector = amrex::Random(engine);
    auto const upper = std::lower_bound(m_cdf.begin(), m_cdf.end(), selector);
    const auto index = static_cast<std::size_t>(upper - m_cdf.begin());
    return m_values[std::min(index, m_values.size() - 1)];
}

amrex::ParticleReal
HallDiscreteDistribution1D::min () const noexcept
{
    return m_min;
}

amrex::ParticleReal
HallDiscreteDistribution1D::max () const noexcept
{
    return m_max;
}

amrex::ParticleReal
HallDiscreteDistribution1D::integral () const noexcept
{
    return m_integral;
}

HallTabulatedDistribution1D::HallTabulatedDistribution1D (
    std::vector<amrex::ParticleReal> values,
    std::vector<amrex::ParticleReal> pdf, int num_bins)
    : m_min(values.empty() ? amrex::ParticleReal(0.0) : values.front()),
      m_max(values.empty() ? amrex::ParticleReal(0.0) : values.back()),
      m_sampler(MakeTabulatedSampler(values, pdf, num_bins))
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_max > m_min,
        "HallTabulatedDistribution1D requires increasing values.");
}

amrex::ParticleReal
HallTabulatedDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    return m_sampler.sample(engine);
}

amrex::ParticleReal
HallTabulatedDistribution1D::min () const noexcept
{
    return m_min;
}

amrex::ParticleReal
HallTabulatedDistribution1D::max () const noexcept
{
    return m_max;
}

amrex::ParticleReal
HallTabulatedDistribution1D::integral () const noexcept
{
    return m_sampler.integral();
}

HallParserDistribution1D::HallParserDistribution1D (
    amrex::ParticleReal min, amrex::ParticleReal max, int num_bins,
    std::string const& variable_name, std::string const& expression)
    : m_min(min), m_max(max),
      m_sampler(
          MakeParserSampler(min, max, num_bins, variable_name, expression))
{
    ValidateRange(m_min, m_max, "HallParserDistribution1D");
}

amrex::ParticleReal
HallParserDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    return m_sampler.sample(engine);
}

amrex::ParticleReal
HallParserDistribution1D::min () const noexcept
{
    return m_min;
}

amrex::ParticleReal
HallParserDistribution1D::max () const noexcept
{
    return m_max;
}

amrex::ParticleReal
HallParserDistribution1D::integral () const noexcept
{
    return m_sampler.integral();
}

std::unique_ptr<HallDistribution1D>
MakeHallDistribution1D (amrex::ParmParse const& pp, std::string const& prefix)
{
    std::string distribution;
    utils::parser::get(pp, prefix, "distribution", distribution);
    distribution = ToLower(distribution);

    if (distribution == "constant") {
        return std::make_unique<HallConstantDistribution1D>(
            GetReal(pp, prefix, "value"));
    }
    if (distribution == "uniform") {
        return std::make_unique<HallUniformDistribution1D>(
            GetReal(pp, prefix, "min"), GetReal(pp, prefix, "max"));
    }
    if (distribution == "area_uniform") {
        return std::make_unique<HallAreaUniformDistribution1D>(
            GetReal(pp, prefix, "min"), GetReal(pp, prefix, "max"));
    }
    if (distribution == "gaussian") {
        return std::make_unique<HallGaussianDistribution1D>(
            QueryReal(pp, prefix, "mean", amrex::ParticleReal(0.0)),
            GetReal(pp, prefix, "sigma"));
    }
    if (distribution == "positive_gaussian") {
        return std::make_unique<HallPositiveGaussianDistribution1D>(
            QueryReal(pp, prefix, "mean", amrex::ParticleReal(0.0)),
            GetReal(pp, prefix, "sigma"));
    }
    if (distribution == "single_spoke") {
        const auto center =
            QueryReal(pp, prefix, "center",
                      QueryReal(pp, prefix, "spoke_center",
                                amrex::Math::pi<amrex::ParticleReal>()));
        return std::make_unique<HallSingleSpokeDistribution1D>(
            center, GetSigmaWithFallback(pp, prefix),
            QueryInt(pp, prefix, "num_bins", 1024));
    }
    if (distribution == "multi_spoke") {
        return std::make_unique<HallMultiSpokeDistribution1D>(
            QueryInt(pp, prefix, "spoke_count", 1),
            GetSigmaWithFallback(pp, prefix),
            QueryReal(pp, prefix, "spoke_phase", amrex::ParticleReal(0.0)),
            QueryInt(pp, prefix, "num_bins", 1024));
    }
    if (distribution == "neutral_spoke" ||
        distribution == "neutral_spoke_depletion") {
        return std::make_unique<HallNeutralSpokeDistribution1D>(
            GetReal(pp, prefix, "ion_width"),
            QueryReal(pp, prefix, "min_ratio", amrex::ParticleReal(0.25)),
            QueryReal(pp, prefix, "drop_exponent", amrex::ParticleReal(4.0)),
            QueryPhaseWithFallback(pp, prefix),
            QueryBoolWithAliases(
                pp, prefix, {"reverse", "reverse_phase", "phase_reverse"}, false),
            QueryInt(pp, prefix, "num_bins", 1024));
    }
    if (distribution == "discrete") {
        std::vector<amrex::ParticleReal> values;
        std::vector<amrex::ParticleReal> weights;
        utils::parser::getArrWithParser(pp, prefix, "values", values);
        utils::parser::queryArrWithParser(pp, prefix, "weights", weights);
        return std::make_unique<HallDiscreteDistribution1D>(std::move(values),
                                                            std::move(weights));
    }
    if (distribution == "tabulated") {
        std::vector<amrex::ParticleReal> values;
        std::vector<amrex::ParticleReal> pdf;
        utils::parser::getArrWithParser(pp, prefix, "values", values);
        utils::parser::getArrWithParser(pp, prefix, "pdf", pdf);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(values.size() == pdf.size() &&
                                             values.size() >= 2,
                                         "tabulated distributions require "
                                         "values and pdf arrays of size >= 2.");
        for (std::size_t i = 1; i < values.size(); ++i) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                values[i] > values[i - 1],
                "tabulated distribution values must be strictly increasing.");
        }
        return std::make_unique<HallTabulatedDistribution1D>(
            std::move(values), std::move(pdf),
            QueryInt(pp, prefix, "num_bins", 1024));
    }
    if (distribution == "parser") {
        const auto axis_name = AxisNameFromPrefix(prefix);
        std::string expression;
        const auto axis_function = "function(" + axis_name + ")";
        if (!utils::parser::Query_parserString(pp, prefix + "." + axis_function,
                                               expression)) {
            utils::parser::Store_parserString(pp, prefix + ".function",
                                              expression);
        }
        return std::make_unique<HallParserDistribution1D>(
            GetReal(pp, prefix, "min"), GetReal(pp, prefix, "max"),
            QueryInt(pp, prefix, "num_bins", 1024), axis_name, expression);
    }

    WARPX_ABORT_WITH_MESSAGE("Unknown HallDistribution1D type: " +
                             distribution);
    return nullptr;
}

} // namespace Insert
