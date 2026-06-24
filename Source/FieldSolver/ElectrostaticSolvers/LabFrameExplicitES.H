/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald, Arianna Formenti, Revathi Jambunathan
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_LABFRAMEEXPLICITES_H_
#define WARPX_LABFRAMEEXPLICITES_H_

#include "ElectrostaticSolver.H"

#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <memory>

class LabFrameExplicitES final : public ElectrostaticSolver
{
public:

    explicit LabFrameExplicitES (int nlevs_max) : ElectrostaticSolver (nlevs_max) {
        ReadParameters();
    }

    void InitData () override;

    void ComputeSpaceChargeField (
        ablastr::fields::MultiFabRegister& fields,
        MultiParticleContainer& mpc,
        MultiFluidContainer* mfl,
        int max_level) override;

    void computePhiTriDiagonal (
        const ablastr::fields::MultiLevelScalarField& rho,
        const ablastr::fields::MultiLevelScalarField& phi
    );

    void computePhiTriDiagonal_periodic (
        const ablastr::fields::MultiLevelScalarField& rho,
        const ablastr::fields::MultiLevelScalarField& phi
    );

private:

    void preparePhiExtrapolatedInitialGuess (
        const ablastr::fields::MultiLevelScalarField& phi
    );

    void initializePhiExtrapolationHistory (
        const ablastr::fields::MultiLevelScalarField& phi
    );

    void definePhiExtrapolationCache (
        const ablastr::fields::MultiLevelScalarField& phi
    );

    amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_phi_extrapolation_cache;
    bool m_phi_extrapolation_cache_initialized = false;

};

#endif  // WARPX_LABFRAMEEXPLICITES_H_
