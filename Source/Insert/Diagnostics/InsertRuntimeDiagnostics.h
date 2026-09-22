#pragma once

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

namespace Insert {

void ParticleNumber ();
void ShowAndWriteIonzationNum (amrex::Vector<int> num);
void ThrustCalc ();
void BeamDivergenceCalc ();
void IEDFCalc ();
void ZmaxRadialExitStatsCalc ();
void ClearHallBoundaryParticleCache ();

} // namespace Insert
