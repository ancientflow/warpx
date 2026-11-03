#ifndef WARPX_INSERT_BOUNDARY_PHI_H
#define WARPX_INSERT_BOUNDARY_PHI_H

#include <ablastr/fields/MultiFabRegister.H>

namespace Insert {

class MaterialPoissonSolver;

/** Fill potential ghosts needed by the common valid-region electric-field operator. */
void FillPhiGhosts (
    ablastr::fields::MultiLevelScalarField const& phi,
    MaterialPoissonSolver* material_solver = nullptr);

void AnodeVoltage (int lev = 0);
void DirichletPhiGuardSet (int lev = 0);
void VoltageAdjustment ();

} // namespace Insert

#endif // WARPX_INSERT_BOUNDARY_PHI_H
