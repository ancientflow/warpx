#pragma once

#include <AMReX_GpuContainers.H>
#include <AMReX_Vector.H>

namespace Insert {

void ParticleNumber ();
void ShowAndWriteIonzationNum (amrex::Vector<int> num);
void AnodeCurrentCalc ();
extern amrex::Gpu::DeviceVector<amrex::Real> g_accumulated_wall_charge_density;
void InitializeAccumulatedZMinWallChargeDensity (long data_size);
void ZMinWallChargeDeposit ();
void ThrustCalc ();
void BeamDivergenceCalc ();
void IEDFCalc ();
void ClearHallBoundaryParticleCache ();

} // namespace Insert
