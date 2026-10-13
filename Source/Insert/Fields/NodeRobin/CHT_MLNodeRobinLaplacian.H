#ifndef CHT_ML_NODE_ROBIN_LAPLACIAN_H_
#define CHT_ML_NODE_ROBIN_LAPLACIAN_H_

#include <AMReX_MLNodeLinOp.H>
#include "CHT_MonopoleBoundary.H"
#include <array>
#include <memory>

namespace cht {

/** Full-domain div(sigma grad(phi)), with a*phi+b*dphi/dn=f on six faces.
 *
 * Direct AMReX nodal-base implementation, usable with amrex::MLMG. Currently
 * supports one AMR level, any number of factor-two MG levels, Cartesian 3D,
 * positive cell-centered scalar sigma, and nonperiodic domain faces.
 * Face order: xlo,xhi,ylo,yhi,zlo,zhi. setRobinBC uses finite nonnegative
 * constants a,b per face, not both zero; f may be a nodal face field.
 * setMonopoleBC supplies analytic spatial coefficients on all six faces.
 * b=0 is Dirichlet; a=0 is Neumann. n always points out of the domain.
 *
 * Pass the PHYSICAL nodal RHS to MLMG::solve (e.g. -rho/epsilon0).
 * MLMG adds the boundary source to its own RHS copy via the virtual hook.
 * Do not add it yourself. sigma=epsilon_r. Boundary terms use lumped nodal
 * surface quadrature and the inherited Neumann half-volume node weights.
 * Internal fixed-potential nodes use AMReX's overset mask. No EB geometry,
 * AMReX patches, or second AMReX library is required.
 */
class MLNodeRobinLaplacian final : public amrex::MLNodeLinOp
{
public:
    struct FaceBC {
        amrex::Real a = 0;
        amrex::Real b = 1;
        amrex::Real f = 0;
    };
    using Faces = std::array<FaceBC, 6>;

    MLNodeRobinLaplacian(const amrex::Vector<amrex::Geometry>& geom,
                        const amrex::Vector<amrex::BoxArray>& grids,
                        const amrex::Vector<amrex::DistributionMapping>& dmap,
                        const amrex::LPInfo& info = amrex::LPInfo());

    std::string name() const override { return "cht::MLNodeRobinLaplacian"; }
    void setSigma(int amrlev, const amrex::MultiFab& sigma);
    void setSigma(amrex::Real sigma);

    /// Resets any face fields supplied by setRobinForcing. Boundary types may
    /// not change after the first prepare/solve; a,b,f may change within a type.
    void setRobinBC(const Faces& faces);
    /// Six open faces with spatial a=(n.r)/|r|^2, b=1, f=a*phi_infinity.
    /// The expansion center must lie strictly inside the physical domain.
    void setMonopoleBC(amrex::GpuArray<amrex::Real,3> const& center,
                       amrex::Real phi_infinity = 0);
    /// AMReX nodal overset convention: 1=unknown, 0=internal fixed potential.
    /// Configure before prepare/solve; the electrode must not touch outer faces.
    /// The caller must cap LPInfo coarsening to resolve the electrode first.
    void setOversetMask(int amrlev, const amrex::iMultiFab& mask, amrex::Real value);
    void setInternalDirichletValue(amrex::Real value);
    amrex::Long internalDirichletNodes(int mglev) const;
    /// Copy f on one face from a one-component nodal MF with the solution
    /// layout. Only surface storage is retained; input may be destroyed.
    /// For b=0 this is a*phi, not phi. Adjacent Dirichlet values and copies of
    /// shared nodal input values must agree. Interior input values are ignored.
    void setRobinForcing(int face, const amrex::MultiFab& f);

    void prepareForSolve() override;
    bool needsUpdate() const override { return m_dirty; }
    void update() override { prepareForSolve(); }
    bool isSingular(int) const override { return m_singular; }
    bool isBottomSingular() const override { return m_singular; }
    amrex::BottomSolver getDefaultBottomSolver() const override { return amrex::BottomSolver::bicgstab; }
    bool supportsAlgMG() const override { return false; }

    void Fapply(int amrlev, int mglev, amrex::MultiFab& out,
                const amrex::MultiFab& in) const override;
    void Fsmooth(int amrlev, int mglev, amrex::MultiFab& sol,
                 const amrex::MultiFab& rhs) const override;
    void normalize(int amrlev, int mglev, amrex::MultiFab& mf) const override;
    void restriction(int amrlev, int cmglev, amrex::MultiFab& coarse,
                     amrex::MultiFab& fine) const override;
    void interpolation(int amrlev, int fmglev, amrex::MultiFab& fine,
                       const amrex::MultiFab& coarse) const override;
    void applyBC(int amrlev, int mglev, amrex::MultiFab& phi, BCMode bc_mode,
                 StateMode state_mode, bool skip_fillboundary = false) const override;
    void applyInhomogNeumannTerm(int amrlev, amrex::MultiFab& rhs) const override;

    /// Diagnostic residual from the physical RHS, including Robin forcing.
    /// The referenced AMReX's 3D MLMG::compResidual skips its inhomogeneous
    /// RHS hook, unlike solve(). Use this helper for physical-source input.
    void compPhysicalResidual(amrex::MultiFab& residual, amrex::MultiFab& phi,
                              const amrex::MultiFab& physical_rhs);

    /// Number of material hierarchy builds; useful for repeated-solve checks.
    int coefficientBuilds() const noexcept { return m_coefficient_builds; }

    // Implementation helpers are public because nvcc requires the enclosing
    // member of an extended device lambda to have public access.
    void fill_sigma_boundary(int mglev);
    void validate_dirichlet_values() const;
    void validate_internal_mask();

protected:
    bool supportRobinBC() const noexcept override { return true; }
    void resizeMultiGrid(int new_size) override;

private:
    // Configuration is owned by setRobinBC; inherited setLevelBC is a final
    // no-op. Do not configure this object through MLLinOp::setDomainBC.
    using amrex::MLLinOp::setDomainBC;

    Faces m_faces{};
    MonopoleBoundary m_monopole;
    std::unique_ptr<amrex::iMultiFab> m_internal_mask;
    amrex::Vector<amrex::Long> m_internal_counts;
    amrex::Real m_internal_value = 0;
    amrex::GpuArray<amrex::Real, 6> m_a_over_b{};
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_sigma;
    // One thin FAB per face/locally-owned box, no full-domain boundary MFs.
    std::array<amrex::Vector<std::unique_ptr<amrex::FArrayBox>>, 6> m_forcing;
    bool m_dirty = true;
    bool m_sigma_dirty = true;
    bool m_prepared = false;
    bool m_singular = true;
    bool m_has_robin_diagonal = false;
    mutable bool m_boundary_values_dirty = true;
    int m_coefficient_builds = 0;
};
} // namespace cht
#endif
