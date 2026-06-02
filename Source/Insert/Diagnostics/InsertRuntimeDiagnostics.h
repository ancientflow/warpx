#pragma once

#include <AMReX_Vector.H>

namespace Insert {

void ParticleNumber ();
void ShowAndWriteIonzationNum (amrex::Vector<int> num);
void AnodeCurrentCalc ();

} // namespace Insert
