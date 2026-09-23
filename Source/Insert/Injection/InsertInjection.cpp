#include "InsertInjection.h"

#include "WarpX.H"

#include "Insert/Config/WarpXSimulationConfig.h"
#include "Insert/Injection/HallInjector.h"

#include <AMReX_ParmParse.H>

namespace Insert {

void
InitializeHallInjection ()
{
#if defined(HALL3D) || defined(HALL3D_INIT)
    HallInjector::GetInstance().InitializePlasma(WarpX::GetInstance());
#endif
}

void
InjectHallParticles ()
{
#if defined(HALL3D) || defined(HALL3D_INIT)
    amrex::Real dt = 0.0;
    amrex::ParmParse pp_mc("my_constants");
    pp_mc.query("dt", dt);
    HallInjector::GetInstance().InjectParticles(WarpX::GetInstance(), dt, 0);
#endif
}

} // namespace Insert
