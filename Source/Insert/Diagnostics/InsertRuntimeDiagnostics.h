#pragma once

#include "Insert/Boundary/ZMinWallCharge.h"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <memory>

namespace amrex
{
class MultiFab;
}

namespace Insert {

void ParticleNumber ();
void ShowAndWriteIonzationNum (amrex::Vector<int> num);
void AnodeCurrentCalc ();
extern std::unique_ptr<amrex::MultiFab> g_accumulated_wall_charge_density;
void InitializeAccumulatedZMinWallChargeDensity (ZMinWallChargeGrid const& grid);
void ZMinWallChargeDeposit ();
void ThrustCalc ();
void BeamDivergenceCalc ();
void IEDFCalc ();
void ClearHallBoundaryParticleCache ();

} // namespace Insert
