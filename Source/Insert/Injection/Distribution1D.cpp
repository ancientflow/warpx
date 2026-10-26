#include "Distribution1D.h"

#include "Insert/Math/CdfUtils.h"
#include "Insert/Math/Sampling.h"
#include "Insert/Utils/InsertUtils.h"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"

#include <AMReX.H>
#include <AMReX_Math.H>
#include <AMReX_Parser.H>
#include <AMReX_Random.H>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <utility>

namespace Insert {
namespace {

amrex::ParticleReal
NeutralSpokeDensity (amrex::ParticleReal theta, amrex::ParticleReal ion_width,
                     amrex::ParticleReal min_ratio,
                     amrex::ParticleReal drop_exponent,
                     amrex::ParticleReal phase, bool reverse,
                     int spoke_count) noexcept
{
    const amrex::ParticleReal period =
        TwoPi() / static_cast<amrex::ParticleReal>(spoke_count);
    const amrex::ParticleReal local_phase =
        reverse ? phase - theta : theta - phase;
    amrex::ParticleReal phi = Math::WrapToPeriod(local_phase, TwoPi());
    // Fold the phase into a single structure period so the depletion
    // profile repeats spoke_count times around the circumference.
    phi = Math::WrapToPeriod(phi, period);
    if (phi < ion_width) {
        const amrex::ParticleReal s = phi / ion_width;
        return min_ratio + (amrex::ParticleReal(1.0) - min_ratio) *
                               std::pow(amrex::ParticleReal(1.0) - s,
                                        drop_exponent);
    }

    const amrex::ParticleReal s = (phi - ion_width) / (period - ion_width);
    return min_ratio + (amrex::ParticleReal(1.0) - min_ratio) * s;
}

void ValidatePositiveSigma (amrex::ParticleReal sigma, std::string const& name);

void
ValidateNeutralSpokeParameters (amrex::ParticleReal ion_width,
                                amrex::ParticleReal min_ratio,
                                amrex::ParticleReal drop_exponent,
                                int spoke_count)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        spoke_count > 0,
        "NeutralSpokeDistribution1D requires spoke_count > 0.");
    const amrex::ParticleReal period =
        TwoPi() / static_cast<amrex::ParticleReal>(spoke_count);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        ion_width > amrex::ParticleReal(0.0) && ion_width < period,
        "NeutralSpokeDistribution1D requires "
        "0 < ion_width < 2*pi/spoke_count.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        min_ratio >= amrex::ParticleReal(0.0) &&
            min_ratio < amrex::ParticleReal(1.0),
        "NeutralSpokeDistribution1D requires 0 <= min_ratio < 1.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        drop_exponent > amrex::ParticleReal(0.0),
        "NeutralSpokeDistribution1D requires drop_exponent > 0.");
}

NumericalInverseCDFSampler1D
MakeSingleSpokeSampler (amrex::ParticleReal center, amrex::ParticleReal sigma,
                        int num_bins)
{
    ValidatePositiveSigma(sigma, "SingleSpokeDistribution1D");
    return NumericalInverseCDFSampler1D(
        amrex::ParticleReal(0.0), TwoPi(), num_bins,
        [=] (amrex::ParticleReal theta) {
            return Math::WrappedGaussianDensity(theta, center, sigma);
        });
}

NumericalInverseCDFSampler1D
MakeMultiSpokeSampler (int spoke_count, amrex::ParticleReal sigma,
                       amrex::ParticleReal phase, int num_bins)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        spoke_count > 0,
        "MultiSpokeDistribution1D requires spoke_count > 0.");
    ValidatePositiveSigma(sigma, "MultiSpokeDistribution1D");
    return NumericalInverseCDFSampler1D(
        amrex::ParticleReal(0.0), TwoPi(), num_bins,
        [=] (amrex::ParticleReal theta) {
            amrex::ParticleReal density = amrex::ParticleReal(0.0);
            const amrex::ParticleReal interval =
                TwoPi() / static_cast<amrex::ParticleReal>(spoke_count);
            for (int i = 0; i < spoke_count; ++i) {
                const auto center =
                    phase + static_cast<amrex::ParticleReal>(i) * interval;
                density += Math::WrappedGaussianDensity(theta, center, sigma);
            }
            return density;
        });
}

NumericalInverseCDFSampler1D
MakeNeutralSpokeSampler (amrex::ParticleReal ion_width,
                         amrex::ParticleReal min_ratio,
                         amrex::ParticleReal drop_exponent,
                         amrex::ParticleReal phase, bool reverse,
                         int spoke_count, int num_bins)
{
    ValidateNeutralSpokeParameters(ion_width, min_ratio, drop_exponent,
                                   spoke_count);
    return NumericalInverseCDFSampler1D(
        amrex::ParticleReal(0.0), TwoPi(), num_bins,
        [=] (amrex::ParticleReal theta) {
            return NeutralSpokeDensity(theta, ion_width, min_ratio,
                                       drop_exponent, phase, reverse,
                                       spoke_count);
        });
}

NumericalInverseCDFSampler1D
MakeTabulatedSampler (std::vector<amrex::ParticleReal> const& values,
                      std::vector<amrex::ParticleReal> const& pdf,
                      int num_bins)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(values.size() == pdf.size() &&
                                         values.size() >= 2,
                                     "TabulatedDistribution1D requires "
                                     "values and pdf arrays of size >= 2.");
    for (std::size_t i = 1; i < values.size(); ++i) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            values[i] > values[i - 1],
            "TabulatedDistribution1D values must be strictly increasing.");
    }

    return NumericalInverseCDFSampler1D(values.front(), values.back(), num_bins,
                                        [values, pdf] (amrex::ParticleReal x) {
                                            return Math::InterpolateTableClamped(
                                                values, pdf, x);
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

/** True when <prefix>.<name>(r) is given, i.e. the parameter is a parser
 *  function of the local radius evaluated in CoordinateDistribution;
 *  the constant <prefix>.<name> is then unused. */
bool
HasRadialParam (amrex::ParmParse const& pp, std::string const& prefix,
                std::string const& name)
{
    std::string expression;
    return utils::parser::Query_parserString(
        pp, prefix + "." + name + "(r)", expression);
}

} // namespace

ConstantDistribution1D::ConstantDistribution1D (
    amrex::ParticleReal value) noexcept
    : m_value(value)
{}

amrex::ParticleReal
ConstantDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    amrex::ignore_unused(engine);
    return m_value;
}

amrex::ParticleReal
ConstantDistribution1D::min () const noexcept
{
    return m_value;
}

amrex::ParticleReal
ConstantDistribution1D::max () const noexcept
{
    return m_value;
}

UniformDistribution1D::UniformDistribution1D (amrex::ParticleReal min,
                                                      amrex::ParticleReal max)
    : m_min(min), m_max(max)
{
    ValidateRange(m_min, m_max, "UniformDistribution1D");
}

amrex::ParticleReal
UniformDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    return Math::SampleUniform(m_min, m_max, engine);
}

amrex::ParticleReal
UniformDistribution1D::min () const noexcept
{
    return m_min;
}

amrex::ParticleReal
UniformDistribution1D::max () const noexcept
{
    return m_max;
}

AreaUniformDistribution1D::AreaUniformDistribution1D (
    amrex::ParticleReal min, amrex::ParticleReal max)
    : m_min(min), m_max(max)
{
    ValidateRange(m_min, m_max, "AreaUniformDistribution1D");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_min >= amrex::ParticleReal(0.0),
        "AreaUniformDistribution1D requires min >= 0.");
}

amrex::ParticleReal
AreaUniformDistribution1D::sample (
    amrex::RandomEngine const& engine) const
{
    return Math::SampleAreaUniformRadius(m_min, m_max, engine);
}

amrex::ParticleReal
AreaUniformDistribution1D::min () const noexcept
{
    return m_min;
}

amrex::ParticleReal
AreaUniformDistribution1D::max () const noexcept
{
    return m_max;
}

GaussianDistribution1D::GaussianDistribution1D (
    amrex::ParticleReal mean, amrex::ParticleReal sigma)
    : m_mean(mean), m_sigma(sigma)
{
    ValidateSigma(m_sigma, "GaussianDistribution1D");
}

amrex::ParticleReal
GaussianDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    if (m_sigma == amrex::ParticleReal(0.0)) {
        return m_mean;
    }
    return amrex::RandomNormal(m_mean, m_sigma, engine);
}

amrex::ParticleReal
GaussianDistribution1D::min () const noexcept
{
    return std::numeric_limits<amrex::ParticleReal>::lowest();
}

amrex::ParticleReal
GaussianDistribution1D::max () const noexcept
{
    return std::numeric_limits<amrex::ParticleReal>::max();
}

PositiveGaussianDistribution1D::PositiveGaussianDistribution1D (
    amrex::ParticleReal mean, amrex::ParticleReal sigma)
    : m_mean(mean), m_sigma(sigma)
{
    ValidateSigma(m_sigma, "PositiveGaussianDistribution1D");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE((m_sigma > amrex::ParticleReal(0.0)) ||
                                         (m_mean >= amrex::ParticleReal(0.0)),
                                     "PositiveGaussianDistribution1D with "
                                     "sigma == 0 requires mean >= 0.");
}

amrex::ParticleReal
PositiveGaussianDistribution1D::sample (
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
PositiveGaussianDistribution1D::min () const noexcept
{
    return amrex::ParticleReal(0.0);
}

amrex::ParticleReal
PositiveGaussianDistribution1D::max () const noexcept
{
    return std::numeric_limits<amrex::ParticleReal>::max();
}

FluxGaussianDistribution1D::FluxGaussianDistribution1D (
    amrex::ParticleReal sigma)
    : m_sigma(sigma)
{
    ValidateSigma(m_sigma, "FluxGaussianDistribution1D");
}

amrex::ParticleReal
FluxGaussianDistribution1D::sample (
    amrex::RandomEngine const& engine) const
{
    // Flux-weighted half-Maxwellian (Rayleigh) normal component,
    // p(v) ∝ v exp(-v^2 / (2 sigma^2)) for v >= 0, via exact inverse CDF.
    // Subtract before converting so a double random draw cannot round to
    // one in single particle precision and produce log(0).
    amrex::ParticleReal const uniform = static_cast<amrex::ParticleReal>(
        amrex::Real(1.0) - amrex::Random(engine));
    return m_sigma *
        std::sqrt(-amrex::ParticleReal(2.0) * std::log(uniform));
}

amrex::ParticleReal
FluxGaussianDistribution1D::min () const noexcept
{
    return amrex::ParticleReal(0.0);
}

amrex::ParticleReal
FluxGaussianDistribution1D::max () const noexcept
{
    return std::numeric_limits<amrex::ParticleReal>::max();
}

SingleSpokeDistribution1D::SingleSpokeDistribution1D (
    amrex::ParticleReal center, amrex::ParticleReal sigma, int num_bins)
    : m_sampler(MakeSingleSpokeSampler(center, sigma, num_bins))
{}

amrex::ParticleReal
SingleSpokeDistribution1D::sample (
    amrex::RandomEngine const& engine) const
{
    return m_sampler.sample(engine);
}

amrex::ParticleReal
SingleSpokeDistribution1D::min () const noexcept
{
    return m_sampler.xmin();
}

amrex::ParticleReal
SingleSpokeDistribution1D::max () const noexcept
{
    return m_sampler.xmax();
}

amrex::ParticleReal
SingleSpokeDistribution1D::integral () const noexcept
{
    return m_sampler.integral();
}

MultiSpokeDistribution1D::MultiSpokeDistribution1D (
    int spoke_count, amrex::ParticleReal sigma, amrex::ParticleReal phase,
    int num_bins)
    : m_sampler(MakeMultiSpokeSampler(spoke_count, sigma, phase, num_bins))
{}

amrex::ParticleReal
MultiSpokeDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    return m_sampler.sample(engine);
}

amrex::ParticleReal
MultiSpokeDistribution1D::min () const noexcept
{
    return m_sampler.xmin();
}

amrex::ParticleReal
MultiSpokeDistribution1D::max () const noexcept
{
    return m_sampler.xmax();
}

amrex::ParticleReal
MultiSpokeDistribution1D::integral () const noexcept
{
    return m_sampler.integral();
}

NeutralSpokeDistribution1D::NeutralSpokeDistribution1D (
    amrex::ParticleReal ion_width, amrex::ParticleReal min_ratio,
    amrex::ParticleReal drop_exponent, amrex::ParticleReal phase,
    bool reverse, int spoke_count, int num_bins)
    : m_sampler(MakeNeutralSpokeSampler(
          ion_width, min_ratio, drop_exponent, phase, reverse,
          spoke_count, num_bins))
{}

amrex::ParticleReal
NeutralSpokeDistribution1D::sample (
    amrex::RandomEngine const& engine) const
{
    return m_sampler.sample(engine);
}

amrex::ParticleReal
NeutralSpokeDistribution1D::min () const noexcept
{
    return m_sampler.xmin();
}

amrex::ParticleReal
NeutralSpokeDistribution1D::max () const noexcept
{
    return m_sampler.xmax();
}

amrex::ParticleReal
NeutralSpokeDistribution1D::integral () const noexcept
{
    return m_sampler.integral();
}

DiscreteDistribution1D::DiscreteDistribution1D (
    std::vector<amrex::ParticleReal> values,
    std::vector<amrex::ParticleReal> weights)
    : m_values(std::move(values))
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_values.empty(), "DiscreteDistribution1D requires values.");
    if (weights.empty()) {
        weights.assign(m_values.size(), amrex::ParticleReal(1.0));
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(weights.size() == m_values.size(),
                                     "DiscreteDistribution1D requires "
                                     "values and weights with equal sizes.");

    m_min = *std::min_element(m_values.begin(), m_values.end());
    m_max = *std::max_element(m_values.begin(), m_values.end());
    m_integral = Math::BuildNormalizedCumulativeWeights(weights, m_cdf);
}

amrex::ParticleReal
DiscreteDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    const amrex::ParticleReal selector = amrex::Random(engine);
    auto const upper = std::lower_bound(m_cdf.begin(), m_cdf.end(), selector);
    const auto index = static_cast<std::size_t>(upper - m_cdf.begin());
    return m_values[std::min(index, m_values.size() - 1)];
}

amrex::ParticleReal
DiscreteDistribution1D::min () const noexcept
{
    return m_min;
}

amrex::ParticleReal
DiscreteDistribution1D::max () const noexcept
{
    return m_max;
}

amrex::ParticleReal
DiscreteDistribution1D::integral () const noexcept
{
    return m_integral;
}

TabulatedDistribution1D::TabulatedDistribution1D (
    std::vector<amrex::ParticleReal> values,
    std::vector<amrex::ParticleReal> pdf, int num_bins)
    : m_min(values.empty() ? amrex::ParticleReal(0.0) : values.front()),
      m_max(values.empty() ? amrex::ParticleReal(0.0) : values.back()),
      m_sampler(MakeTabulatedSampler(values, pdf, num_bins))
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_max > m_min,
        "TabulatedDistribution1D requires increasing values.");
}

amrex::ParticleReal
TabulatedDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    return m_sampler.sample(engine);
}

amrex::ParticleReal
TabulatedDistribution1D::min () const noexcept
{
    return m_min;
}

amrex::ParticleReal
TabulatedDistribution1D::max () const noexcept
{
    return m_max;
}

amrex::ParticleReal
TabulatedDistribution1D::integral () const noexcept
{
    return m_sampler.integral();
}

ParserDistribution1D::ParserDistribution1D (
    amrex::ParticleReal min, amrex::ParticleReal max, int num_bins,
    std::string const& variable_name, std::string const& expression)
    : m_min(min), m_max(max),
      m_sampler(
          MakeParserSampler(min, max, num_bins, variable_name, expression))
{
    ValidateRange(m_min, m_max, "ParserDistribution1D");
}

amrex::ParticleReal
ParserDistribution1D::sample (amrex::RandomEngine const& engine) const
{
    return m_sampler.sample(engine);
}

amrex::ParticleReal
ParserDistribution1D::min () const noexcept
{
    return m_min;
}

amrex::ParticleReal
ParserDistribution1D::max () const noexcept
{
    return m_max;
}

amrex::ParticleReal
ParserDistribution1D::integral () const noexcept
{
    return m_sampler.integral();
}

std::unique_ptr<Distribution1D>
MakeDistribution1D (amrex::ParmParse const& pp, std::string const& prefix)
{
    std::string distribution;
    utils::parser::get(pp, prefix, "distribution", distribution);
    distribution = ToLower(distribution);

    if (distribution == "constant") {
        return std::make_unique<ConstantDistribution1D>(
            GetWithParser<amrex::ParticleReal>(pp, prefix, "value"));
    }
    if (distribution == "uniform") {
        return std::make_unique<UniformDistribution1D>(
            GetWithParser<amrex::ParticleReal>(pp, prefix, "min"),
            GetWithParser<amrex::ParticleReal>(pp, prefix, "max"));
    }
    if (distribution == "area_uniform") {
        return std::make_unique<AreaUniformDistribution1D>(
            GetWithParser<amrex::ParticleReal>(pp, prefix, "min"),
            GetWithParser<amrex::ParticleReal>(pp, prefix, "max"));
    }
    if (distribution == "gaussian") {
        if (HasRadialParam(pp, prefix, "sigma")) {
            // sigma/mean are functions of r, evaluated per particle in
            // CoordinateDistribution; the constant base is unused.
            return std::make_unique<GaussianDistribution1D>(
                amrex::ParticleReal(0.0), amrex::ParticleReal(0.0));
        }
        return std::make_unique<GaussianDistribution1D>(
            QueryWithParser<amrex::ParticleReal>(
                pp, prefix, "mean", amrex::ParticleReal(0.0)),
            GetWithParser<amrex::ParticleReal>(pp, prefix, "sigma"));
    }
    if (distribution == "positive_gaussian") {
        if (HasRadialParam(pp, prefix, "sigma")) {
            return std::make_unique<PositiveGaussianDistribution1D>(
                amrex::ParticleReal(0.0), amrex::ParticleReal(0.0));
        }
        return std::make_unique<PositiveGaussianDistribution1D>(
            QueryWithParser<amrex::ParticleReal>(
                pp, prefix, "mean", amrex::ParticleReal(0.0)),
            GetWithParser<amrex::ParticleReal>(pp, prefix, "sigma"));
    }
    if (distribution == "flux_gaussian") {
        return std::make_unique<FluxGaussianDistribution1D>(
            GetWithParser<amrex::ParticleReal>(pp, prefix, "sigma"));
    }
    if (distribution == "single_spoke") {
        const auto center =
            QueryWithParser<amrex::ParticleReal>(
                pp, prefix, "center",
                QueryWithParser<amrex::ParticleReal>(
                    pp, prefix, "spoke_center",
                    amrex::Math::pi<amrex::ParticleReal>()));
        return std::make_unique<SingleSpokeDistribution1D>(
            center, GetSigmaWithFallback(pp, prefix),
            QueryWithParser<int>(pp, prefix, "num_bins", 1024));
    }
    if (distribution == "multi_spoke") {
        return std::make_unique<MultiSpokeDistribution1D>(
            QueryWithParser<int>(pp, prefix, "spoke_count", 1),
            GetSigmaWithFallback(pp, prefix),
            QueryWithParser<amrex::ParticleReal>(
                pp, prefix, "spoke_phase", amrex::ParticleReal(0.0)),
            QueryWithParser<int>(pp, prefix, "num_bins", 1024));
    }
    if (distribution == "neutral_spoke" ||
        distribution == "neutral_spoke_depletion") {
        return std::make_unique<NeutralSpokeDistribution1D>(
            GetWithParser<amrex::ParticleReal>(pp, prefix, "ion_width"),
            QueryWithParser<amrex::ParticleReal>(
                pp, prefix, "min_ratio", amrex::ParticleReal(0.25)),
            QueryWithParser<amrex::ParticleReal>(
                pp, prefix, "drop_exponent", amrex::ParticleReal(4.0)),
            QueryPhaseWithFallback(pp, prefix),
            QueryBoolWithAliases(
                pp, prefix, {"reverse", "reverse_phase", "phase_reverse"}, false),
            QueryWithParser<int>(pp, prefix, "spoke_count", 1),
            QueryWithParser<int>(pp, prefix, "num_bins", 1024));
    }
    if (distribution == "discrete") {
        std::vector<amrex::ParticleReal> values;
        std::vector<amrex::ParticleReal> weights;
        utils::parser::getArrWithParser(pp, prefix, "values", values);
        utils::parser::queryArrWithParser(pp, prefix, "weights", weights);
        return std::make_unique<DiscreteDistribution1D>(std::move(values),
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
        return std::make_unique<TabulatedDistribution1D>(
            std::move(values), std::move(pdf),
            QueryWithParser<int>(pp, prefix, "num_bins", 1024));
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
        return std::make_unique<ParserDistribution1D>(
            GetWithParser<amrex::ParticleReal>(pp, prefix, "min"),
            GetWithParser<amrex::ParticleReal>(pp, prefix, "max"),
            QueryWithParser<int>(pp, prefix, "num_bins", 1024),
            axis_name, expression);
    }

    WARPX_ABORT_WITH_MESSAGE("Unknown Distribution1D type: " +
                             distribution);
    return nullptr;
}

} // namespace Insert
