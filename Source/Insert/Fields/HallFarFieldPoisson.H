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

#include <array>
#include <memory>
#include <string>

class WarpX;

namespace Insert {

/** Fixed auxiliary Poisson grid, owned by the single-level electrostatic solver.
 *
 * Stage one only constructs the layout and storage. Charge transfer, material
 * initialization and the coarse solve will be added in subsequent stages.
 */
class HallFarFieldPoisson
{
public:
    void initialize (WarpX const& warpx, amrex::MultiFab const& phi,
                     amrex::MultiFab const& rho, int history_capacity);

    /** Reject changes to the fine layout instead of retaining stale MPI mappings. */
    void assertLayout (amrex::Geometry const& geom, amrex::MultiFab const& phi,
                       amrex::MultiFab const& rho) const;

private:
    void defineGeometry (amrex::Geometry const& fine_geom);
    void defineBoxes (amrex::BoxArray const& fine_boxes,
                     amrex::DistributionMapping const& fine_dmap);
    void allocateFields (amrex::MultiFab const& phi, amrex::MultiFab const& rho,
                         int history_capacity);
    void reportLayout () const;

    bool m_defined = false;
    amrex::Geometry m_fine_geom;
    amrex::BoxArray m_fine_boxes;
    amrex::DistributionMapping m_dmap;
    amrex::IntVect m_fine_phi_ngrow;
    amrex::IntVect m_fine_rho_ngrow;

    amrex::Geometry m_geom;
    amrex::BoxArray m_boxes;
    amrex::IntVect m_ratio;
    // Fine-domain lower corner relative to the expanded coarse-domain lower corner.
    amrex::IntVect m_fine_offset;
    std::string m_layout_file;

    amrex::MultiFab m_rho;
    amrex::MultiFab m_phi;
    // Coarsened fine Boxes retain their fine-grid ownership and natural indices.
    amrex::MultiFab m_rho_patch;
    amrex::MultiFab m_phi_patch;
    std::array<std::unique_ptr<amrex::MultiFab>, 3> m_phi_history;
    std::array<amrex::Real, 3> m_history_times{};
    int m_history_depth = 0;

    // Populated with the coarse material operator in stage three. MLMG dies first.
    std::unique_ptr<amrex::MLNodeLaplacian> m_linop;
    std::unique_ptr<amrex::MLMG> m_mlmg;
};

} // namespace Insert

#endif // WARPX_USE_HALL_FAR_FIELD
#endif // WARPX_INSERT_HALLFARFIELDPOISSON_H_
