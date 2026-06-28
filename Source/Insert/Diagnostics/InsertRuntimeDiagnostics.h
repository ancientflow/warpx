#pragma once

#include <AMReX_Vector.H>

namespace Insert {

void ParticleNumber ();
void ShowAndWriteIonzationNum (amrex::Vector<int> num);
void AnodeCurrentCalc ();
void ZMinWallChargeDeposit ();
void ThrustCalc ();
void BeamDivergenceCalc ();
void IEDFCalc ();
void ClearHallBoundaryParticleCache ();

} // namespace Insert
