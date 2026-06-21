#pragma once

#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX_iMultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <memory>
#include <string>

namespace amrex
{
class ParmParse;
}

namespace Insert
{

void Initialize ();
void BeforeStep ();
void ParticleInjection ();
void ReadParticleSubcycling (
    std::string const& species_name, amrex::ParmParse const& pp_species);
int ParticleSubcyclingNdt (std::string const& species_name);
void ApplyParticleSubcycling (
    std::string const& species_name, int step, amrex::Real& dt, bool& do_not_push);
void SetBoundaryPhi ();
amrex::Vector<std::unique_ptr<amrex::iMultiFab> >
BuildPhiOversetMasks (ablastr::fields::MultiLevelScalarField const& phi);
void SetPhiGuards ();
void BeforeCollision (int step);
void AfterCollision (int step);
void AfterDiagnostics ();
void Finalize ();

} // namespace Insert
