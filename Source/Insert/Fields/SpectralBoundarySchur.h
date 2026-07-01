#pragma once

#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX_GpuContainers.H>
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
        amrex::Gpu::DeviceVector<amrex::Real> const& sigma_s_device,
        int nx_face,
        int ny_face);

    amrex::MultiFab const* PhiCorrection ();
}

} // namespace Insert
