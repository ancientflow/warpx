#include "Insert/Boundary/CeramicInteraction.h"

#include "Insert/Math/ThermalVelocity.h"
#include "Utils/WarpXConst.H"

namespace Insert {

CeramicInteraction::CeramicInteraction (amrex::ParticleReal const secondary_temperature_eV)
    : m_thermal_velocity(static_cast<amrex::ParticleReal>(
          Math::ThermalVelocityFromEV(secondary_temperature_eV, PhysConst::m_e)))
{}

} // namespace Insert
