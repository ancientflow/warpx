#pragma once

#include "Insert/Injection/HallCoordinateTransform.h"
#include "Insert/Injection/HallDistribution1D.h"
#include "Insert/Injection/HallPositionSampler.h"

#include <AMReX_ParmParse.H>
#include <AMReX_RandomEngine.H>

#include <memory>

namespace Insert {

class HallCoordinateDistribution : public HallPositionSampler
{
public:
    HallCoordinateDistribution (
        HallCoordinateSystem coordinate_system,
        std::unique_ptr<HallDistribution1D> q0,
        std::unique_ptr<HallDistribution1D> q1,
        std::unique_ptr<HallDistribution1D> q2);

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
    HallCoordinateSystem m_coordinate_system;
    std::unique_ptr<HallDistribution1D> m_q0;
    std::unique_ptr<HallDistribution1D> m_q1;
    std::unique_ptr<HallDistribution1D> m_q2;
};

std::unique_ptr<HallCoordinateDistribution>
MakeHallCoordinateDistribution (
    amrex::ParmParse const& pp, std::string const& prefix,
    HallCoordinateSpace space);

} // namespace Insert
