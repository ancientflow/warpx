#pragma once

#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX_REAL.H>

namespace amrex
{
class Geometry;
class MultiFab;
}

namespace Insert
{

void ApplyElectrostaticBoundaryCorrection (
    ablastr::fields::MultiLevelScalarField const& phi);

namespace SpectralBoundarySchur
{
    bool Enabled ();

    void ApplyZMinCorrection (
        ablastr::fields::MultiLevelScalarField const& phi);

    void SolveAndReconstruct (
        amrex::Geometry const& geom,
        amrex::MultiFab const& sigma_s_mf,
        amrex::MultiFab const& base_phi);

    amrex::MultiFab const* PhiCorrection ();
}

} // namespace Insert
