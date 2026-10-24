#ifndef WARPX_INSERT_MATERIAL_POISSON_SOLVER_H_
#define WARPX_INSERT_MATERIAL_POISSON_SOLVER_H_

#include <ablastr/fields/MLMGOptions.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <memory>

namespace Insert {

/** Single-level material Poisson solver for Dirichlet and monopole boundaries.
 * Owns one immutable material operator/MLMG cache, allocated on the first solve.
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

    /** Keep rho physical; scale a private RHS and reuse the static operator.
     * In Dirichlet mode the caller must set current face potentials in phi
     * before calling. Voltage may vary at each node and between solves.
     */
    void solve (amrex::MultiFab const& rho, amrex::MultiFab& phi,
                amrex::MultiFab const& epsilon, amrex::iMultiFab const& mask,
                amrex::MultiFab const& voltage, amrex::Geometry const& geom,
                ablastr::fields::MLMGOptions const& options);

    /** Monopole only: fill one exterior phi layer for centered differences. */
    void fillPhiGhosts (amrex::MultiFab& phi, amrex::Geometry const& geom) const;

    /** Convert an assembled q/dV copy to density per nodal boundary dual volume. */
    static void normalizeBoundarySource (amrex::MultiFab& source, amrex::Geometry const& geom);

private:
    struct Cache;
    std::unique_ptr<Cache> m_cache;
    bool m_monopole = false;
    bool m_diagnostics = false;
    amrex::GpuArray<amrex::Real, 3> m_center{};
    amrex::Real m_phi_infinity = 0;
};

} // namespace Insert
#endif
