#pragma once

#include "Insert/Injection/HallCoordinateTransform.h"

#include <AMReX_RandomEngine.H>

namespace Insert {

class HallPositionSampler
{
public:
    virtual ~HallPositionSampler () = default;

    [[nodiscard]] virtual EmissionSample
    samplePosition (amrex::RandomEngine const& engine) const = 0;
};

} // namespace Insert
