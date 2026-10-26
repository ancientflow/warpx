#pragma once

#include "Insert/Injection/DistributionSampler1D.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_RandomEngine.H>

#include <memory>
#include <string>
#include <vector>

namespace Insert {

class Distribution1D
{
public:
    virtual ~Distribution1D () = default;

    [[nodiscard]] virtual amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const = 0;

    [[nodiscard]] virtual amrex::ParticleReal min () const noexcept = 0;
    [[nodiscard]] virtual amrex::ParticleReal max () const noexcept = 0;

    [[nodiscard]] virtual amrex::ParticleReal
    integral () const noexcept {
        return amrex::ParticleReal(1.0);
    }
};

class ConstantDistribution1D final : public Distribution1D
{
public:
    explicit ConstantDistribution1D (amrex::ParticleReal value) noexcept;

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;

private:
    amrex::ParticleReal m_value;
};

class UniformDistribution1D final : public Distribution1D
{
public:
    UniformDistribution1D (amrex::ParticleReal min,
                               amrex::ParticleReal max);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;

private:
    amrex::ParticleReal m_min;
    amrex::ParticleReal m_max;
};

class AreaUniformDistribution1D final : public Distribution1D
{
public:
    AreaUniformDistribution1D (amrex::ParticleReal min,
                                   amrex::ParticleReal max);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;

private:
    amrex::ParticleReal m_min;
    amrex::ParticleReal m_max;
};

class GaussianDistribution1D final : public Distribution1D
{
public:
    GaussianDistribution1D (amrex::ParticleReal mean,
                                amrex::ParticleReal sigma);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;

private:
    amrex::ParticleReal m_mean;
    amrex::ParticleReal m_sigma;
};

class PositiveGaussianDistribution1D final : public Distribution1D
{
public:
    PositiveGaussianDistribution1D (amrex::ParticleReal mean,
                                        amrex::ParticleReal sigma);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;

private:
    amrex::ParticleReal m_mean;
    amrex::ParticleReal m_sigma;
};

class FluxGaussianDistribution1D final : public Distribution1D
{
public:
    explicit FluxGaussianDistribution1D (amrex::ParticleReal sigma);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;

private:
    amrex::ParticleReal m_sigma;
};

class SingleSpokeDistribution1D final : public Distribution1D
{
public:
    SingleSpokeDistribution1D (amrex::ParticleReal center,
                                   amrex::ParticleReal sigma, int num_bins);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;
    [[nodiscard]] amrex::ParticleReal integral () const noexcept override;

private:
    NumericalInverseCDFSampler1D m_sampler;
};

class MultiSpokeDistribution1D final : public Distribution1D
{
public:
    MultiSpokeDistribution1D (int spoke_count, amrex::ParticleReal sigma,
                                  amrex::ParticleReal phase, int num_bins);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;
    [[nodiscard]] amrex::ParticleReal integral () const noexcept override;

private:
    NumericalInverseCDFSampler1D m_sampler;
};

class NeutralSpokeDistribution1D final : public Distribution1D
{
public:
    NeutralSpokeDistribution1D (
        amrex::ParticleReal ion_width, amrex::ParticleReal min_ratio,
        amrex::ParticleReal drop_exponent, amrex::ParticleReal phase,
        bool reverse, int spoke_count, int num_bins);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;
    [[nodiscard]] amrex::ParticleReal integral () const noexcept override;

private:
    NumericalInverseCDFSampler1D m_sampler;
};

class DiscreteDistribution1D final : public Distribution1D
{
public:
    DiscreteDistribution1D (std::vector<amrex::ParticleReal> values,
                                std::vector<amrex::ParticleReal> weights);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;
    [[nodiscard]] amrex::ParticleReal integral () const noexcept override;

private:
    std::vector<amrex::ParticleReal> m_values;
    std::vector<amrex::ParticleReal> m_cdf;
    amrex::ParticleReal m_min;
    amrex::ParticleReal m_max;
    amrex::ParticleReal m_integral;
};

class TabulatedDistribution1D final : public Distribution1D
{
public:
    TabulatedDistribution1D (std::vector<amrex::ParticleReal> values,
                                 std::vector<amrex::ParticleReal> pdf,
                                 int num_bins);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;
    [[nodiscard]] amrex::ParticleReal integral () const noexcept override;

private:
    amrex::ParticleReal m_min;
    amrex::ParticleReal m_max;
    NumericalInverseCDFSampler1D m_sampler;
};

class ParserDistribution1D final : public Distribution1D
{
public:
    ParserDistribution1D (amrex::ParticleReal min, amrex::ParticleReal max,
                              int num_bins, std::string const& variable_name,
                              std::string const& expression);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;
    [[nodiscard]] amrex::ParticleReal integral () const noexcept override;

private:
    amrex::ParticleReal m_min;
    amrex::ParticleReal m_max;
    NumericalInverseCDFSampler1D m_sampler;
};

std::unique_ptr<Distribution1D>
MakeDistribution1D (amrex::ParmParse const& pp, std::string const& prefix);

} // namespace Insert
