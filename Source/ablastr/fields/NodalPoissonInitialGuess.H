#ifndef ABLASTR_FIELDS_NODAL_POISSON_INITIAL_GUESS_H_
#define ABLASTR_FIELDS_NODAL_POISSON_INITIAL_GUESS_H_

#include "MultiFabRegister.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <array>
#include <memory>

namespace ablastr::fields {

/** Initial-guess orchestration for Cartesian nodal Poisson solves.
 *
 * Owns only history and the Poisson unknown-node mask.
 * Prepares phi in place after native coarse-to-fine interpolation.
 * No particle gather/deposition masks or electric fields are accessed here.
 * The fixed-node rule matches MLNodeLinOp's cell-coverage/Dirichlet mask.
 */
class NodalPoissonInitialGuess
{
public:
    using BC = amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>;

    void beginSolve (MultiLevelScalarField const& phi,
                     amrex::Vector<amrex::Geometry> const& geom,
                     BC const& lobc, BC const& hibc,
                     amrex::Vector<amrex::iMultiFab const*> const& overset,
                     amrex::Real time, int order, amrex::Real alpha, amrex::Real beta);
    void prepare (int lev, amrex::MultiFab& phi, int verbosity);
    void commit (MultiLevelScalarField const& phi);

private:
    struct Prediction {
        int order = -1; // -1: no history, keep the caller's initial guess
        bool changed_interval = false;
        bool finite_fallback = false;
    };
    struct Level {
        // Independent valid-region snapshots, newest first; never aliases the solution.
        std::array<std::unique_ptr<amrex::MultiFab>, 3> history;
        std::array<amrex::Real, 3> time{};
        int depth = 0;
        std::unique_ptr<amrex::iMultiFab> unknown;
        amrex::Geometry geom;
        BC lobc{};
        BC hibc{};
        amrex::iMultiFab const* overset = nullptr;
        bool zero_mean = false;
    };
    void defineHistory (int lev, amrex::MultiFab const& phi);
    void commitHistory (int lev, amrex::MultiFab const& phi, amrex::Real gauge_offset);

    amrex::Vector<Level> m_levels;
    amrex::Real m_time = 0;
    int m_order = 0;
    amrex::Real m_alpha = 0;
    amrex::Real m_beta = 0;
};

} // namespace ablastr::fields
#endif
