#pragma once

#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX_iMultiFab.H>
#include <AMReX_Vector.H>

#include <memory>

namespace Insert {

void AnodeVoltage ();
amrex::Vector<std::unique_ptr<amrex::iMultiFab> > const&
BuildPhiOversetMasks (ablastr::fields::MultiLevelScalarField const& phi);
void DirichletPhiGuardSet ();
void VoltageAdjustment ();

} // namespace Insert
