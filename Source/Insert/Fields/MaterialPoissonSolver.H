#ifndef WARPX_INSERT_MATERIAL_POISSON_SOLVER_H_
#define WARPX_INSERT_MATERIAL_POISSON_SOLVER_H_

#include <ablastr/fields/MLMGOptions.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <memory>
#include <vector>

namespace amrex { template <typename MF> class MLMGT; }

namespace Insert {

/** Material Poisson solver with one immutable operator/MLMG cache per AMR level.
 * Each operator solves one level; only level zero has physical monopole faces.
 * No time history, particle data or interlevel interpolation is owned here.
 */
class MaterialPoissonSolver
{
public:
    MaterialPoissonSolver ();
    ~MaterialPoissonSolver ();
    MaterialPoissonSolver (MaterialPoissonSolver const&) = delete;
    MaterialPoissonSolver& operator= (MaterialPoissonSolver const&) = delete;

    void readParameters (amrex::Geometry const& geom);
    [[nodiscard]] bool isMonopole () const noexcept { return m_monopole; }

    /** Common hierarchy callback: RHS is already scaled by -1/epsilon0.
     * Normalize Robin physical-face dual volumes only on level zero, in a private RHS.
     * Preserve interpolated coarse/fine values and ghosts; mask only that private RHS.
     * The returned MLMG remains owned by this solver.
     */
    amrex::MLMGT<amrex::MultiFab>& solveLevel (
        int lev, amrex::MultiFab const& rhs, amrex::MultiFab& phi,
        amrex::MultiFab const& epsilon, amrex::iMultiFab const& mask,
        amrex::MultiFab const& voltage, amrex::Geometry const& geom,
        ablastr::fields::MLMGOptions const& options, amrex::Real absolute_tolerance);

    /** Monopole only: fill one exterior phi layer for centered differences. */
    void fillPhiGhosts (amrex::MultiFab& phi, amrex::Geometry const& geom) const;

    /** Convert an assembled q/dV copy to density per nodal boundary dual volume. */
    static void normalizeBoundarySource (amrex::MultiFab& source, amrex::Geometry const& geom);

private:
    struct Cache;
    std::vector<std::unique_ptr<Cache>> m_caches;
    bool m_monopole = false;
    bool m_diagnostics = false;
    amrex::GpuArray<amrex::Real, 3> m_center{};
    amrex::Real m_phi_infinity = 0;
};

} // namespace Insert
#endif
