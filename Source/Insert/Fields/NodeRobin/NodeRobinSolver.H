#ifndef WARPX_INSERT_NODE_ROBIN_SOLVER_H_
#define WARPX_INSERT_NODE_ROBIN_SOLVER_H_

#include <ablastr/fields/MLMGOptions.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <memory>

namespace Insert {

/** Single-level material Poisson solve with a monopole outer boundary. */
class NodeRobinSolver
{
public:
    NodeRobinSolver ();
    ~NodeRobinSolver ();
    NodeRobinSolver (NodeRobinSolver const&) = delete;
    NodeRobinSolver& operator= (NodeRobinSolver const&) = delete;

    void readParameters (amrex::Geometry const& geom);
    [[nodiscard]] bool enabled () const noexcept { return m_enabled; }

    /** Keep rho physical; scale a private RHS and reuse the static operator. */
    void solve (amrex::MultiFab const& rho, amrex::MultiFab& phi,
                amrex::MultiFab const& epsilon, amrex::iMultiFab const& mask,
                amrex::MultiFab const& voltage, amrex::Geometry const& geom,
                ablastr::fields::MLMGOptions const& options);

private:
    struct Cache;
    std::unique_ptr<Cache> m_cache;
    bool m_enabled = false;
    bool m_diagnostics = false;
    amrex::GpuArray<amrex::Real, 3> m_center{};
    amrex::Real m_phi_infinity = 0;
};

} // namespace Insert
#endif
