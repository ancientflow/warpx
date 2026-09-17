#pragma once

#include "Insert/Injection/HallCoordinateTransform.h"
#include "Insert/Injection/HallDistribution1D.h"
#include "Insert/Injection/HallPositionSampler.h"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_RandomEngine.H>

#include <array>
#include <memory>
#include <optional>

namespace Insert {

/** Optional radial (r) dependence of a Gaussian velocity axis. When
 *  sigma_of_r is set, the axis is not sampled from its constant base
 *  distribution but as a (drift) Maxwellian whose mean and sigma are
 *  functions of the local radius r about the injection axis, e.g.
 *  reconstructed from exit-plane radial statistics. Only valid for the
 *  gaussian and positive_gaussian distribution types. */
struct HallRadialGaussianOverride
{
    bool positive = false; // truncate to v >= 0 by rejection sampling
    amrex::ParticleReal mean_const = amrex::ParticleReal(0.0);
    // The Parser objects own the compiled-expression data (shared_ptr)
    // backing the executors; they must live as long as the executors.
    amrex::Parser sigma_parser;
    amrex::Parser mean_parser;
    std::optional<amrex::ParserExecutor<1>> sigma_of_r;
    std::optional<amrex::ParserExecutor<1>> mean_of_r;
};

class HallCoordinateDistribution : public HallPositionSampler
{
public:
    HallCoordinateDistribution (
        HallCoordinateSystem coordinate_system,
        std::unique_ptr<HallDistribution1D> q0,
        std::unique_ptr<HallDistribution1D> q1,
        std::unique_ptr<HallDistribution1D> q2,
        HallRotatingAxisFrame rotating_axis_frame = {},
        std::array<std::optional<HallRadialGaussianOverride>, 3>
            radial_velocity = {});

    [[nodiscard]] HallCoordinateSystem coordinateSystem () const noexcept;

    [[nodiscard]] amrex::XDim3
    sampleCoordinates (amrex::RandomEngine const& engine) const;

    [[nodiscard]] EmissionSample
    samplePosition (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::XDim3
    sampleVelocity (
        amrex::RandomEngine const& engine,
        EmissionSample const& position) const;

private:
    [[nodiscard]] amrex::XDim3
    sampleVelocityRadial (
        amrex::RandomEngine const& engine,
        EmissionSample const& position) const;

    HallCoordinateSystem m_coordinate_system;
    std::unique_ptr<HallDistribution1D> m_q0;
    std::unique_ptr<HallDistribution1D> m_q1;
    std::unique_ptr<HallDistribution1D> m_q2;
    HallRotatingAxisFrame m_rotating_axis_frame;
    std::array<std::optional<HallRadialGaussianOverride>, 3>
        m_radial_velocity;
    bool m_has_radial_velocity = false;
};

std::unique_ptr<HallCoordinateDistribution>
MakeHallCoordinateDistribution (
    amrex::ParmParse const& pp, std::string const& prefix,
    HallCoordinateSpace space,
    std::string const& allowed_coupled_distribution = "");

} // namespace Insert
