#ifndef WARPX_INSERT_HALLFARFIELDPOISSON_H_
#define WARPX_INSERT_HALLFARFIELDPOISSON_H_

#include "Insert/Config/HallElectrostaticConfig.H"

#ifdef WARPX_USE_HALL_FAR_FIELD

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <array>
#include <memory>
#include <string>

class WarpX;

namespace Insert {

/** Extend physical ghosts from synchronized interior samples. */
void ExtendPhysicalGhosts (
    amrex::MultiFab& field, amrex::MultiFab& work,
    amrex::Geometry const& geom, int order);

/** Fixed auxiliary Poisson grid, owned by the single-level electrostatic solver.
 *
 * Transfers charge and solves the coarse material problem with a temporal predictor.
 * Supplies the fine-grid initial guess and Dirichlet data from the coarse solution.
 */
class HallFarFieldPoisson
{
public:
    struct SolverSettings
    {
        amrex::Real m_relative_tolerance = 1.e-11;
        amrex::Real m_absolute_tolerance = 0;
        int m_max_iters = 200;
        int m_verbosity = 2;
        int m_max_coarsening_level = 4;
        int m_num_final_sweeps = 8;
        int m_extrapolation_order = 0;
        amrex::Real m_extrapolation_alpha = 0;
        amrex::Real m_extrapolation_beta = 0;
    };

    void initialize (WarpX const& warpx, amrex::MultiFab const& phi,
                     amrex::MultiFab const& rho, SolverSettings const& defaults);

    /** Restrict synchronized total fine charge, then zero the coarse anode RHS.
     * The input includes wall charge once and has not undergone PEC reflection.
     */
    void restrictCharge (amrex::MultiFab const& rho);

    /** Solve the current coarse source and retain only converged potential history. */
    void solve (amrex::Real time, bool verbose_step);

    /** Interpolate the converged coarse solution onto all valid fine nodes. */
    void interpolateToFine (amrex::MultiFab& phi);

    /** Synchronize E and extend its physical ghosts for particle gathering. */
    void fillElectricFieldGhosts (std::array<amrex::MultiFab*, 3> const& fields);

private:
    void readSolverSettings (SolverSettings const& defaults);
    void defineGeometry (amrex::Geometry const& fine_geom);
    void defineBoxes (amrex::BoxArray const& fine_boxes,
                     amrex::DistributionMapping const& fine_dmap);
    void allocateFields (amrex::MultiFab const& phi, amrex::MultiFab const& rho,
                         int history_capacity);
    void defineChargeTransfer (amrex::MultiFab const& rho);
    void defineSolver ();
    void prepareInitialGuess (amrex::Real time);
    void updateHistory (amrex::Real time);
    void reportLayout () const;

    SolverSettings m_settings;
    amrex::Geometry m_fine_geom;
    amrex::BoxArray m_fine_boxes;
    amrex::DistributionMapping m_dmap;

    amrex::Geometry m_geom;
    amrex::BoxArray m_boxes;
    amrex::IntVect m_ratio;
    // Fine-domain lower corner relative to the expanded coarse-domain lower corner.
    amrex::IntVect m_fine_offset;
    std::string m_layout_file;

    amrex::MultiFab m_rho;
    amrex::MultiFab m_phi;
    amrex::MultiFab m_rhs;
    amrex::MultiFab m_epsilon_r;
    amrex::MultiFab m_anode_phi;
    // Partial restrictions of uniquely owned fine nodes.
    amrex::MultiFab m_rho_patch;
    // Coarsened fine Boxes retain the natural coarse indices, including z offsets.
    amrex::MultiFab m_phi_patch;
    std::array<std::unique_ptr<amrex::MultiFab>, 3> m_fine_e_work;
    // Expanded at physical faces to retain all synchronized deposition ghosts.
    amrex::MultiFab m_rho_source;
    amrex::iMultiFab m_source_anode_mask;
    amrex::iMultiFab m_anode_mask;
    std::unique_ptr<amrex::iMultiFab> m_source_owner;
    std::array<std::unique_ptr<amrex::MultiFab>, 3> m_phi_history;
    std::array<amrex::Real, 3> m_history_times{};
    int m_history_depth = 0;

    // The cached operator borrows the static material fields above. MLMG dies first.
    std::unique_ptr<amrex::MLNodeLaplacian> m_linop;
    std::unique_ptr<amrex::MLMG> m_mlmg;
};

} // namespace Insert

#endif // WARPX_USE_HALL_FAR_FIELD
#endif // WARPX_INSERT_HALLFARFIELDPOISSON_H_
