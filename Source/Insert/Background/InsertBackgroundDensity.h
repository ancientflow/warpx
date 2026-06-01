#pragma once

#include "BackgroundCoupledDensity.h"

#include <AMReX_Vector.H>

#ifdef MCC_DENSITY
extern amrex::Vector<BackgroundCoupledDensity> global_background_density;
#endif

namespace Insert {

void GlobalBackgroundDensityInit ();
void GlobalBackgroundDensityUpdate (int step);
void GlobalBackgroundDensityClean (int step);

} // namespace Insert
