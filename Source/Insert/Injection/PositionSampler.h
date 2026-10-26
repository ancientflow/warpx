#pragma once

#include "Insert/Injection/CoordinateTransform.h"

#include <AMReX_RandomEngine.H>

namespace Insert {

class PositionSampler
{
public:
    virtual ~PositionSampler () = default;

    [[nodiscard]] virtual EmissionSample
    samplePosition (amrex::RandomEngine const& engine) const = 0;
};

} // namespace Insert
