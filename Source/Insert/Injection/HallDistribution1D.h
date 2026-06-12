#pragma once

#include "Insert/Injection/DistributionSampler1D.H"

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_RandomEngine.H>

#include <memory>
#include <string>
#include <vector>

namespace Insert {

class HallDistribution1D
{
public:
    virtual ~HallDistribution1D () = default;

    [[nodiscard]] virtual amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const = 0;

    [[nodiscard]] virtual amrex::ParticleReal min () const noexcept = 0;
    [[nodiscard]] virtual amrex::ParticleReal max () const noexcept = 0;

    [[nodiscard]] virtual amrex::ParticleReal
    integral () const noexcept {
        return amrex::ParticleReal(1.0);
    }
};

class HallConstantDistribution1D final : public HallDistribution1D
{
public:
    explicit HallConstantDistribution1D (amrex::ParticleReal value) noexcept;

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;

private:
    amrex::ParticleReal m_value;
};

class HallUniformDistribution1D final : public HallDistribution1D
{
public:
    HallUniformDistribution1D (amrex::ParticleReal min,
                               amrex::ParticleReal max);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;

private:
    amrex::ParticleReal m_min;
    amrex::ParticleReal m_max;
};

class HallAreaUniformDistribution1D final : public HallDistribution1D
{
public:
    HallAreaUniformDistribution1D (amrex::ParticleReal min,
                                   amrex::ParticleReal max);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;

private:
    amrex::ParticleReal m_min;
    amrex::ParticleReal m_max;
};

class HallGaussianDistribution1D final : public HallDistribution1D
{
public:
    HallGaussianDistribution1D (amrex::ParticleReal mean,
                                amrex::ParticleReal sigma);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;

private:
    amrex::ParticleReal m_mean;
    amrex::ParticleReal m_sigma;
};

class HallPositiveGaussianDistribution1D final : public HallDistribution1D
{
public:
    HallPositiveGaussianDistribution1D (amrex::ParticleReal mean,
                                        amrex::ParticleReal sigma);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;

private:
    amrex::ParticleReal m_mean;
    amrex::ParticleReal m_sigma;
};

class HallSingleSpokeDistribution1D final : public HallDistribution1D
{
public:
    HallSingleSpokeDistribution1D (amrex::ParticleReal center,
                                   amrex::ParticleReal sigma, int num_bins);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;
    [[nodiscard]] amrex::ParticleReal integral () const noexcept override;

private:
    NumericalInverseCDFSampler1D m_sampler;
};

class HallMultiSpokeDistribution1D final : public HallDistribution1D
{
public:
    HallMultiSpokeDistribution1D (int spoke_count, amrex::ParticleReal sigma,
                                  amrex::ParticleReal phase, int num_bins);

    [[nodiscard]] amrex::ParticleReal
    sample (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::ParticleReal min () const noexcept override;
    [[nodiscard]] amrex::ParticleReal max () const noexcept override;
    [[nodiscard]] amrex::ParticleReal integral () const noexcept override;

private:
    NumericalInverseCDFSampler1D m_sampler;
};

class HallDiscreteDistribution1D final : public HallDistribution1D
{
public:
    HallDiscreteDistribution1D (std::vector<amrex::ParticleReal> values,
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

class HallTabulatedDistribution1D final : public HallDistribution1D
{
public:
    HallTabulatedDistribution1D (std::vector<amrex::ParticleReal> values,
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

class HallParserDistribution1D final : public HallDistribution1D
{
public:
    HallParserDistribution1D (amrex::ParticleReal min, amrex::ParticleReal max,
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

std::unique_ptr<HallDistribution1D>
MakeHallDistribution1D (amrex::ParmParse const& pp, std::string const& prefix);

} // namespace Insert
