#ifndef WARPX_INSERT_MATH_THERMALVELOCITY_H_
#define WARPX_INSERT_MATH_THERMALVELOCITY_H_

#include "Utils/WarpXConst.H"

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

/** \brief Temperature-to-thermal-velocity conversions shared by the Insert
 *         boundary-interaction and injection code paths.
 *
 *  The returned value is the one-dimensional standard deviation of a
 *  Maxwellian velocity distribution; it is meant to be combined with a
 *  Gaussian (or Gaussian flux) sampler. The return type is the promoted
 *  type of the constant-times-argument expression, i.e. amrex::Real for
 *  the physical constants used here.
 */
namespace Insert::Math {

/** \brief Thermal velocity for a temperature in Kelvin:
 *         vth = sqrt(kb * T / m). */
template <typename T, typename U>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto
ThermalVelocityFromTemperature (T temperature_K, U mass) noexcept
    -> decltype(PhysConst::kb * temperature_K / mass)
{
    using std::sqrt;
    return sqrt(PhysConst::kb * temperature_K / mass);
}

/** \brief Thermal velocity for a temperature in electron-volts:
 *         vth = sqrt(T_eV * q_e / m). */
template <typename T, typename U>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto
ThermalVelocityFromEV (T temperature_eV, U mass) noexcept
    -> decltype(temperature_eV * PhysConst::q_e / mass)
{
    using std::sqrt;
    return sqrt(temperature_eV * PhysConst::q_e / mass);
}

} // namespace Insert::Math

#endif // WARPX_INSERT_MATH_THERMALVELOCITY_H_
