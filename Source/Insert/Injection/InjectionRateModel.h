#pragma once

#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <memory>
#include <string>

namespace Insert {

class InjectionRateModel
{
public:
    virtual ~InjectionRateModel () = default;

    [[nodiscard]] virtual amrex::Real
    expectedMacroParticles (amrex::Real dt) const = 0;
};

class FixedCountRateModel final : public InjectionRateModel
{
public:
    explicit FixedCountRateModel (amrex::Real count);

    [[nodiscard]] amrex::Real
    expectedMacroParticles (amrex::Real dt) const override;

private:
    amrex::Real m_count;
};

class DensityVolumeRateModel final : public InjectionRateModel
{
public:
    DensityVolumeRateModel (
        amrex::Real density, amrex::Real volume, amrex::Real macro_weight);

    [[nodiscard]] amrex::Real
    expectedMacroParticles (amrex::Real dt) const override;

private:
    amrex::Real m_density;
    amrex::Real m_volume;
    amrex::Real m_macro_weight;
};

class CurrentRateModel final : public InjectionRateModel
{
public:
    CurrentRateModel (
        amrex::Real current, amrex::Real macro_weight,
        amrex::Real l_factor);

    [[nodiscard]] amrex::Real
    expectedMacroParticles (amrex::Real dt) const override;

private:
    amrex::Real m_current;
    amrex::Real m_macro_weight;
    amrex::Real m_l_factor;
};

class MassFlowRateModel final : public InjectionRateModel
{
public:
    MassFlowRateModel (
        amrex::Real mass_flow, amrex::Real particle_mass,
        amrex::Real macro_weight, amrex::Real l_factor);

    [[nodiscard]] amrex::Real
    expectedMacroParticles (amrex::Real dt) const override;

private:
    amrex::Real m_mass_flow;
    amrex::Real m_particle_mass;
    amrex::Real m_macro_weight;
    amrex::Real m_l_factor;
};

class FractionalParticleAccumulator
{
public:
    [[nodiscard]] int
    consume (amrex::Real expected_macro_particles);

    [[nodiscard]] amrex::Real remainder () const noexcept;

private:
    amrex::Real m_remainder = 0.0;
};

std::unique_ptr<InjectionRateModel>
MakeInjectionRateModel (
    amrex::ParmParse const& pp, std::string const& prefix);

} // namespace Insert
