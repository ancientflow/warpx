#include "CHT_MLNodeRobinLaplacian.H"
#include "CHT_NodeRobin_K.H"
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Reduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <cmath>
#include <limits>

namespace cht {
using namespace amrex;
namespace {

Box face_box(Box box, const Box& domain, int face)
{
    int const d = face/2;
    int const p = face%2 ? domain.bigEnd(d) : domain.smallEnd(d);
    Box plane = domain;
    plane.setSmall(d, p);
    plane.setBig(d, p);
    return box & plane;
}

// The volume operator uses reflected coefficient ghosts. The four tangential
// cell contributions give the lumped surface coefficient, including edges.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real face_sigma(int i, int j, int k, int face, Array4<Real const> const& s) noexcept
{
    int const hi = face%2;
    if (face/2 == 0) {
        int const x = i-hi;
        return Real(0.25)*(s(x,j-1,k-1)+s(x,j,k-1)+s(x,j-1,k)+s(x,j,k));
    }
    if (face/2 == 1) {
        int const y = j-hi;
        return Real(0.25)*(s(i-1,y,k-1)+s(i,y,k-1)+s(i-1,y,k)+s(i,y,k));
    }
    int const z = k-hi;
    return Real(0.25)*(s(i-1,j-1,z)+s(i,j-1,z)+s(i-1,j,z)+s(i,j,z));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real robin_diagonal(int i, int j, int k, Array4<Real const> const& sigma,
                    Box const& domain, GpuArray<Real,3> const& dxinv,
                    GpuArray<Real,6> const& ab, MonopoleBoundary const& monopole,
                    GpuArray<Real,3> const& problo) noexcept
{
    int const p[3] = {i,j,k};
    Real d = 0;
    for (int dir=0; dir<3; ++dir) {
        if (p[dir] == domain.smallEnd(dir) && ab[2*dir] != 0) {
            Real const a = monopole.enabled ? monopole.coefficient(i,j,k,2*dir,domain,problo,dxinv) : ab[2*dir];
            d += 2*dxinv[dir]*a*face_sigma(i,j,k,2*dir,sigma);
        }
        if (p[dir] == domain.bigEnd(dir) && ab[2*dir+1] != 0) {
            Real const a = monopole.enabled ? monopole.coefficient(i,j,k,2*dir+1,domain,problo,dxinv) : ab[2*dir+1];
            d += 2*dxinv[dir]*a*face_sigma(i,j,k,2*dir+1,sigma);
        }
    }
    return d;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real diffusion_diagonal(int i, int j, int k, Array4<Real const> const& s,
                        GpuArray<Real,3> const& dxinv) noexcept
{
    Real const fac = -(dxinv[0]*dxinv[0]+dxinv[1]*dxinv[1]+dxinv[2]*dxinv[2])/9;
    return fac*(s(i-1,j-1,k-1)+s(i,j-1,k-1)+s(i-1,j,k-1)+s(i,j,k-1)
               +s(i-1,j-1,k)+s(i,j-1,k)+s(i-1,j,k)+s(i,j,k));
}

void finish_kernel()
{
    if (!Gpu::inNoSyncRegion()) { Gpu::streamSynchronize(); }
}
} // namespace

MLNodeRobinLaplacian::MLNodeRobinLaplacian(
    const Vector<Geometry>& geom, const Vector<BoxArray>& grids,
    const Vector<DistributionMapping>& dmap, const LPInfo& info)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.size()==1 && grids.size()==1 && dmap.size()==1,
        "CHT Robin currently supports one AMR level (multiple MG levels are supported)");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom[0].Coord()==0 && !geom[0].isAnyPeriodic(),
        "CHT Robin requires nonperiodic Cartesian geometry");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(grids[0].ixType().cellCentered(),
        "Pass cell-centered grids to CHT Robin constructor");
    // Explicitly override the EB-enabled AMReX build's default RAP strategy.
    m_coarsening_strategy = CoarseningStrategy::Sigma;
    MLNodeLinOp::define(geom, grids, dmap, info, {}, 0);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_domain_covered[0], "CHT Robin requires full domain coverage");
    for (const auto& ratio : mg_coarsen_ratio_vec) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ratio == IntVect(2),
            "CHT Robin currently requires factor-two coarsening in all directions");
    }
    m_sigma.resize(m_num_mg_levels[0]);
    for (int m=0; m<m_num_mg_levels[0]; ++m) {
        m_sigma[m] = std::make_unique<MultiFab>(m_grids[0][m], m_dmap[0][m], 1, 1);
        m_sigma[m]->setVal(1);
    }
    setRobinBC(Faces{});
}

void MLNodeRobinLaplacian::setSigma(int amrlev, const MultiFab& sigma)
{
    AMREX_ALWAYS_ASSERT(amrlev == 0);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sigma.nComp()==1 &&
        sigma.boxArray()==m_grids[0][0] && sigma.DistributionMap()==m_dmap[0][0],
        "sigma must be a scalar cell-centered MultiFab with the operator layout");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!sigma.contains_nan() && !sigma.contains_inf() && sigma.min(0)>0,
        "sigma must be finite and strictly positive");
    MultiFab::Copy(*m_sigma[0], sigma, 0, 0, 1, 0);
    m_sigma_dirty = m_dirty = true;
}

void MLNodeRobinLaplacian::setSigma(Real sigma)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(sigma) && sigma>0,
        "sigma must be finite and strictly positive");
    m_sigma[0]->setVal(sigma);
    m_sigma_dirty = m_dirty = true;
}

void MLNodeRobinLaplacian::setRobinBC(const Faces& faces)
{
    Array<LinOpBCType,3> lo, hi;
    bool operator_changed = !m_prepared;
    bool const was_singular = m_singular;
    m_singular = !m_internal_mask;
    m_has_robin_diagonal = false;
    operator_changed |= m_monopole.enabled;
    m_monopole.enabled = false;
    for (int f=0; f<6; ++f) {
        const auto& c = faces[f];
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(c.a) && std::isfinite(c.b) &&
            std::isfinite(c.f) && c.a>=0 && c.b>=0 && (c.a>0 || c.b>0),
            "Robin requires finite a,b,f; a,b nonnegative and not both zero");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_prepared || (c.b==0)==(m_faces[f].b==0),
            "Changing Dirichlet/Robin boundary type requires a new operator/MLMG");
        auto const type = c.b==0 ? LinOpBCType::Dirichlet : LinOpBCType::Robin;
        (f%2 ? hi : lo)[f/2] = type;
        operator_changed |= c.a!=m_faces[f].a || c.b!=m_faces[f].b;
        m_a_over_b[f] = c.b>0 ? c.a/c.b : 0;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(m_a_over_b[f]) &&
            std::isfinite(c.f/(c.b>0 ? c.b:c.a)), "Robin coefficient ratio overflow");
        if (c.a>0) { m_singular = false; }
        m_has_robin_diagonal |= m_a_over_b[f]>0;
        m_forcing[f].clear();
    }
    m_faces = faces;
    MLLinOp::setDomainBC(lo, hi);
    // MLNodeLinOp only allocates the finest-level solvability weights when
    // isBottomSingular() is true. A Robin-to-Neumann coefficient update must
    // create those weights before MLMG starts projecting the new RHS.
    if (m_prepared && was_singular != m_singular) {
        m_masks_built = false;
        m_bottom_dot_mask = MultiFab();
        m_coarse_dot_mask = MultiFab();
    }
    m_dirty |= operator_changed;
    m_boundary_values_dirty = true;
}

void MLNodeRobinLaplacian::setMonopoleBC(GpuArray<Real,3> const& center, Real phi_infinity)
{
    auto const lo=m_geom[0][0].ProbLoArray(), hi=m_geom[0][0].ProbHiArray();
    AMREX_ALWAYS_ASSERT(std::isfinite(phi_infinity));
    for (int d=0; d<3; ++d) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(center[d]) && center[d]>lo[d] && center[d]<hi[d],
            "Six-face monopole expansion center must be strictly inside the domain");
    }
    bool const was_enabled=m_monopole.enabled;
    auto const old_center=m_monopole.center;
    bool const old_dirty=m_dirty;
    Faces bc;
    for (auto& f:bc) { f={1,1,phi_infinity}; }
    setRobinBC(bc);
    m_monopole={true,center,phi_infinity};
    // Infinity-potential changes affect only the boundary RHS.
    bool same_center=true;
    for (int d=0; d<3; ++d) { same_center &= old_center[d]==center[d]; }
    m_dirty=old_dirty || !was_enabled || !same_center;
}

void MLNodeRobinLaplacian::setOversetMask(int amrlev, const iMultiFab& mask, Real value)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrlev==0 && !m_prepared && !m_internal_mask,
        "Set the internal electrode mask once, before prepare/solve");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mask.nComp()==1 &&
        mask.boxArray()==amrex::convert(m_grids[0][0],IntVect(1)) &&
        mask.DistributionMap()==m_dmap[0][0], "Internal mask must use the scalar nodal solution layout");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mask.min(0)==0 && mask.max(0)==1,
        "Internal mask requires binary values, at least one fixed node and at least one unknown");
    Box const domain=amrex::surroundingNodes(m_geom[0][0].Domain());
    for (int face=0; face<6; ++face) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mask.min(face_box(domain,domain,face),0)==1,
            "The internal electrode must not touch the external boundary");
    }
    m_internal_mask=std::make_unique<iMultiFab>(mask.boxArray(),mask.DistributionMap(),1,0);
    iMultiFab::Copy(*m_internal_mask,mask,0,0,1,0);
    // Deterministic owner synchronization for Box-shared electrode nodes.
    m_internal_mask->OverrideSync(m_geom[0][0].periodicity());
    MLNodeLinOp::setOversetMask(0,*m_internal_mask);
    setInternalDirichletValue(value);
    m_singular=false;
    m_dirty=true;
}

void MLNodeRobinLaplacian::setInternalDirichletValue(Real value)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_internal_mask && std::isfinite(value),
        "A finite internal potential and an existing electrode mask are required");
    m_internal_value=value;
}

void MLNodeRobinLaplacian::validate_internal_mask()
{
    if (!m_internal_mask || m_internal_counts.size()==static_cast<std::size_t>(m_num_mg_levels[0])) { return; }
    m_internal_counts.clear();
    iMultiFab const* current=m_internal_mask.get();
    std::unique_ptr<iMultiFab> coarse;
    for (int m=0; m<m_num_mg_levels[0]; ++m) {
        if (m>0) {
            auto next=std::make_unique<iMultiFab>(amrex::convert(m_grids[0][m],IntVect(1)),m_dmap[0][m],1,0);
            amrex::average_down_nodal(*current,*next,IntVect(2));
            coarse=std::move(next);
            current=coarse.get();
        }
        auto owner=makeOwnerMask(m_grids[0][m],m_dmap[0][m],m_geom[0][m]);
        iMultiFab fixed(current->boxArray(),current->DistributionMap(),1,0);
        auto const c=current->const_arrays(); auto const o=owner->const_arrays(); auto const f=fixed.arrays();
        ParallelFor(fixed,[=] AMREX_GPU_DEVICE(int box,int i,int j,int k) noexcept {
            f[box](i,j,k)=!c[box](i,j,k) && o[box](i,j,k);
        });
        finish_kernel();
        Long const count=fixed.sum(0);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(count>0,
            "Coarsening erased the internal electrode; reduce LPInfo max coarsening level");
        m_internal_counts.push_back(count);
    }
}

Long MLNodeRobinLaplacian::internalDirichletNodes(int mglev) const
{
    AMREX_ALWAYS_ASSERT(mglev>=0 && mglev<static_cast<int>(m_internal_counts.size()));
    return m_internal_counts[mglev];
}

void MLNodeRobinLaplacian::setRobinForcing(int face, const MultiFab& f)
{
    AMREX_ALWAYS_ASSERT(face>=0 && face<6);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_monopole.enabled, "Use setMonopoleBC to set the far-field potential in monopole mode");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(f.nComp()==1 &&
        f.boxArray()==amrex::convert(m_grids[0][0],IntVect(1)) &&
        f.DistributionMap()==m_dmap[0][0], "Robin forcing must use the scalar nodal solution layout");
    auto& stored = m_forcing[face];
    stored.clear();
    stored.resize(f.local_size());
    Box const domain = amrex::surroundingNodes(m_geom[0][0].Domain());
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<int> reduce_data(reduce_op);
    using Tuple = ReduceData<int>::Type;
    Real const divisor = m_faces[face].b>0 ? m_faces[face].b : m_faces[face].a;
    for (MFIter mfi(f); mfi.isValid(); ++mfi) {
        Box const bx = face_box(mfi.validbox(), domain, face);
        if (!bx.ok()) { continue; }
        auto fab = std::make_unique<FArrayBox>(bx,1,The_Arena());
        auto dest = fab->array();
        auto src = f.const_array(mfi);
        reduce_op.eval(bx,reduce_data,[=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept -> Tuple {
            return {std::isfinite(src(i,j,k)/divisor) ? 0:1};
        });
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept { dest(i,j,k)=src(i,j,k); });
        stored[mfi.LocalIndex()] = std::move(fab);
    }
    Gpu::streamSynchronize(); // owned boundary copy must outlive the caller's input
    int invalid = amrex::get<0>(reduce_data.value());
    ParallelDescriptor::ReduceIntMax(invalid);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(invalid<=0,"Nonfinite Robin forcing on domain face");
    m_boundary_values_dirty = true;
}

void MLNodeRobinLaplacian::validate_dirichlet_values() const
{
    if (!m_boundary_values_dirty) { return; }
    Box const domain = amrex::surroundingNodes(m_geom[0][0].Domain());
    ReduceOps<ReduceOpMax> op;
    ReduceData<Real> data(op);
    using Tuple = ReduceData<Real>::Type;
    for (MFIter mfi(*m_sigma[0]); mfi.isValid(); ++mfi) {
        for (int f=0; f<6; ++f) {
            if (m_faces[f].b!=0) { continue; }
            for (int g=f+1; g<6; ++g) {
                if (m_faces[g].b!=0 || f/2==g/2) { continue; }
                Box const bx = face_box(face_box(amrex::surroundingNodes(mfi.validbox()),domain,f),domain,g);
                if (!bx.ok()) { continue; }
                auto const& fs = m_forcing[f]; auto const& gs = m_forcing[g];
                Array4<Real const> fa = fs.empty() ? Array4<Real const>{} : fs[mfi.LocalIndex()]->const_array();
                Array4<Real const> ga = gs.empty() ? Array4<Real const>{} : gs[mfi.LocalIndex()]->const_array();
                Real const ff=m_faces[f].f, gf=m_faces[g].f, af=m_faces[f].a, ag=m_faces[g].a;
                op.eval(bx,data,[=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept -> Tuple {
                    Real const u=(fa ? fa(i,j,k) : ff)/af;
                    Real const v=(ga ? ga(i,j,k) : gf)/ag;
                    return {std::abs(u-v)/(1+std::abs(u)+std::abs(v))};
                });
            }
        }
    }
    Real error=amrex::get<0>(data.value());
    ParallelDescriptor::ReduceRealMax(error);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(error<=100*std::numeric_limits<Real>::epsilon(),
        "Conflicting Dirichlet values at an edge/corner");
    m_boundary_values_dirty=false;
}

void MLNodeRobinLaplacian::fill_sigma_boundary(int mglev)
{
    MultiFab& s=*m_sigma[mglev];
    s.FillBoundary();
    Box const domain=m_geom[0][mglev].Domain();
    for (MFIter mfi(s); mfi.isValid(); ++mfi) {
        auto a=s.array(mfi);
        Box const bx=mfi.fabbox();
        ParallelFor(bx,[=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept {
            if (!domain.contains(IntVect(i,j,k))) {
                int const ii=amrex::min(amrex::max(i,domain.smallEnd(0)),domain.bigEnd(0));
                int const jj=amrex::min(amrex::max(j,domain.smallEnd(1)),domain.bigEnd(1));
                int const kk=amrex::min(amrex::max(k,domain.smallEnd(2)),domain.bigEnd(2));
                a(i,j,k)=a(ii,jj,kk);
            }
        });
    }
    finish_kernel();
}

void MLNodeRobinLaplacian::prepareForSolve()
{
    BL_PROFILE("cht::MLNodeRobinLaplacian::prepareForSolve");
    MLNodeLinOp::prepareForSolve();
    validate_internal_mask();
    buildMasks();
    if (m_sigma_dirty) {
        for (int m=1; m<m_num_mg_levels[0]; ++m) {
            amrex::average_down(*m_sigma[m-1],*m_sigma[m],m_geom[0][m-1],m_geom[0][m],0,1,IntVect(2));
        }
        for (int m=0; m<m_num_mg_levels[0]; ++m) { fill_sigma_boundary(m); }
        ++m_coefficient_builds;
        m_sigma_dirty=false;
    }
    validate_dirichlet_values();
    m_prepared=true;
    m_dirty=false;
}

void MLNodeRobinLaplacian::applyBC(int amrlev,int mglev,MultiFab& phi,BCMode bc_mode,
                                  StateMode state_mode,bool skip_fillboundary) const
{
    bool const physical = bc_mode==BCMode::Inhomogeneous && state_mode==StateMode::Solution && !m_precond_mode;
    if (physical) { validate_dirichlet_values(); }
    Box const domain=amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());
    // MLMG initializes corrections with zero Dirichlet nodes; the inherited
    // masks, our smoother and interpolation preserve them. As in the native
    // nodal operator, do not re-launch boundary kernels on every MG correction.
    for (int f=0; physical && f<6; ++f) {
        if (m_faces[f].b!=0) { continue; }
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            Box const bx=face_box(mfi.validbox(),domain,f);
            if (!bx.ok()) { continue; }
            auto dest=phi.array(mfi);
            Array4<Real const> forcing;
            if (physical && mglev==0 && !m_forcing[f].empty()) {
                forcing=m_forcing[f][mfi.LocalIndex()]->const_array();
            }
            Real const a=m_faces[f].a, fv=m_faces[f].f;
            ParallelFor(bx,[=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept {
                dest(i,j,k)=physical ? (forcing ? forcing(i,j,k) : fv)/a : Real(0);
            });
        }
    }
    if (physical && mglev==0 && m_internal_mask) {
        auto const mask=m_internal_mask->const_arrays();
        auto const u=phi.arrays();
        Real const value=m_internal_value;
        ParallelFor(phi,[=] AMREX_GPU_DEVICE(int box,int i,int j,int k) noexcept {
            if (!mask[box](i,j,k)) { u[box](i,j,k)=value; }
        });
        finish_kernel();
    }
    MLNodeLinOp::applyBC(amrlev,mglev,phi,bc_mode,state_mode,skip_fillboundary);
}

void MLNodeRobinLaplacian::applyInhomogNeumannTerm(int amrlev,MultiFab& rhs) const
{
    AMREX_ALWAYS_ASSERT(amrlev==0);
    if (m_precond_mode) { return; }
    Box const domain=amrex::surroundingNodes(m_geom[0][0].Domain());
    auto const dxinv=m_geom[0][0].InvCellSizeArray();
    auto const problo=m_geom[0][0].ProbLoArray();
    auto const monopole=m_monopole;
    for (int f=0; f<6; ++f) {
        if (m_faces[f].b==0 || (m_faces[f].f==0 && m_forcing[f].empty())) { continue; }
        for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
            Box const bx=face_box(mfi.validbox(),domain,f);
            if (!bx.ok()) { continue; }
            auto dest=rhs.array(mfi);
            auto mask=m_dirichlet_mask[0][0]->const_array(mfi);
            auto s=m_sigma[0]->const_array(mfi);
            Array4<Real const> forcing = m_forcing[f].empty() ? Array4<Real const>{}
                : m_forcing[f][mfi.LocalIndex()]->const_array();
            Real const b=m_faces[f].b, fv=m_faces[f].f;
            ParallelFor(bx,[=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept {
                if (!mask(i,j,k)) {
                    Real const source=monopole.enabled ? monopole.coefficient(i,j,k,f,domain,problo,dxinv)*monopole.phi_infinity
                        : (forcing ? forcing(i,j,k) : fv)/b;
                    dest(i,j,k)-=2*dxinv[f/2]*face_sigma(i,j,k,f,s)*source;
                }
            });
        }
    }
    finish_kernel();
}

void MLNodeRobinLaplacian::Fapply(int amrlev,int mglev,MultiFab& out,const MultiFab& in) const
{
    BL_PROFILE("cht::MLNodeRobinLaplacian::Fapply");
    auto const x=in.const_arrays(); auto const y=out.arrays();
    auto const s=m_sigma[mglev]->const_arrays();
    auto const mask=m_dirichlet_mask[amrlev][mglev]->const_arrays();
    auto const dxinv=m_geom[amrlev][mglev].InvCellSizeArray();
    Box const domain=amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());
    auto const monopole=m_monopole;
    auto const problo=m_geom[amrlev][mglev].ProbLoArray();
    auto const ab=m_a_over_b;
    bool const robin=m_has_robin_diagonal;
    ParallelFor(out,[=] AMREX_GPU_DEVICE(int box,int i,int j,int k) noexcept {
        Real val=detail::mlndlap_adotx_aa(i,j,k,x[box],s[box],mask[box],dxinv);
        if (robin && !mask[box](i,j,k)) {
            val-=robin_diagonal(i,j,k,s[box],domain,dxinv,ab,monopole,problo)*x[box](i,j,k);
        }
        y[box](i,j,k)=val;
    });
    finish_kernel();
}

void MLNodeRobinLaplacian::compPhysicalResidual(MultiFab& residual, MultiFab& phi,
                                               const MultiFab& physical_rhs)
{
    if (!m_prepared || m_dirty) { prepareForSolve(); }
    MultiFab assembled(physical_rhs.boxArray(),physical_rhs.DistributionMap(),1,0,
                       MFInfo().SetArena(The_Async_Arena()));
    MultiFab::Copy(assembled,physical_rhs,0,0,1,0);
    applyInhomogNeumannTerm(0,assembled);
    solutionResidual(0,residual,phi,assembled);
}

void MLNodeRobinLaplacian::Fsmooth(int amrlev,int mglev,MultiFab& sol,const MultiFab& rhs) const
{
    BL_PROFILE("cht::MLNodeRobinLaplacian::Fsmooth");
    auto const u=sol.arrays(); auto const r=rhs.const_arrays();
    auto const s=m_sigma[mglev]->const_arrays();
    auto const mask=m_dirichlet_mask[amrlev][mglev]->const_arrays();
    auto const dxinv=m_geom[amrlev][mglev].InvCellSizeArray();
    Box const domain=amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());
    auto const monopole=m_monopole;
    auto const problo=m_geom[amrlev][mglev].ProbLoArray();
    auto const ab=m_a_over_b;
    bool const robin=m_has_robin_diagonal;
    // Same eight-color stencil ordering as AMReX's GPU nodal smoother.
    for (int color=0; color<8; ++color) {
        ParallelForStrided(sol,IntVect(2),multicolor_offset(color),
        [=] AMREX_GPU_DEVICE(int box,int i,int j,int k) noexcept {
            if (mask[box](i,j,k)) { u[box](i,j,k)=0; return; }
            Real const d=robin ? robin_diagonal(i,j,k,s[box],domain,dxinv,ab,monopole,problo) : 0;
            Real const ax=detail::mlndlap_adotx_aa(i,j,k,u[box],s[box],mask[box],dxinv)-d*u[box](i,j,k);
            Real const diag=diffusion_diagonal(i,j,k,s[box],dxinv)-d;
            u[box](i,j,k)+=(r[box](i,j,k)-ax)/diag;
        });
    }
    finish_kernel();
    nodalSync(amrlev,mglev,sol);
}

void MLNodeRobinLaplacian::normalize(int amrlev,int mglev,MultiFab& mf) const
{
    auto const x=mf.arrays(); auto const s=m_sigma[mglev]->const_arrays();
    auto const mask=m_dirichlet_mask[amrlev][mglev]->const_arrays();
    auto const dxinv=m_geom[amrlev][mglev].InvCellSizeArray();
    Box const domain=amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());
    auto const monopole=m_monopole;
    auto const problo=m_geom[amrlev][mglev].ProbLoArray();
    auto const ab=m_a_over_b; bool const robin=m_has_robin_diagonal;
    ParallelFor(mf,[=] AMREX_GPU_DEVICE(int box,int i,int j,int k) noexcept {
        if (!mask[box](i,j,k)) {
            Real const d=robin ? robin_diagonal(i,j,k,s[box],domain,dxinv,ab,monopole,problo) : 0;
            x[box](i,j,k)/=diffusion_diagonal(i,j,k,s[box],dxinv)-d;
        }
    });
    finish_kernel();
}

void MLNodeRobinLaplacian::restriction(int amrlev,int cmglev,MultiFab& coarse,MultiFab& fine) const
{
    applyBC(amrlev,cmglev-1,fine,BCMode::Homogeneous,StateMode::Correction);
    bool const copy=!amrex::isMFIterSafe(coarse,fine);
    MultiFab temp;
    if (copy) {
        temp.define(amrex::coarsen(fine.boxArray(),2),fine.DistributionMap(),1,0,
                    MFInfo().SetArena(The_Async_Arena()));
    }
    MultiFab& c=copy ? temp : coarse;
    auto const ca=c.arrays(); auto const fa=fine.const_arrays();
    auto const mask=m_dirichlet_mask[amrlev][cmglev-1]->const_arrays();
    ParallelFor(c,[=] AMREX_GPU_DEVICE(int box,int i,int j,int k) noexcept {
        detail::mlndlap_restriction(i,j,k,ca[box],fa[box],mask[box]);
    });
    finish_kernel();
    if (copy) { coarse.ParallelCopy(temp); }
}

void MLNodeRobinLaplacian::interpolation(int amrlev,int fmglev,MultiFab& fine,const MultiFab& coarse) const
{
    bool const copy=!amrex::isMFIterSafe(coarse,fine);
    MultiFab temp;
    if (copy) {
        temp.define(amrex::coarsen(fine.boxArray(),2),fine.DistributionMap(),1,0,
                    MFInfo().SetArena(The_Async_Arena()));
        temp.ParallelCopy(coarse);
    }
    MultiFab const& c=copy ? temp : coarse;
    auto const ca=c.const_arrays(); auto const fa=fine.arrays();
    auto const mask=m_dirichlet_mask[amrlev][fmglev]->const_arrays();
    auto const s=m_sigma[fmglev]->const_arrays();
    ParallelFor(fine,[=] AMREX_GPU_DEVICE(int box,int i,int j,int k) noexcept {
        detail::mlndlap_interpadd_aa(i,j,k,fa[box],ca[box],s[box],mask[box]);
    });
    finish_kernel();
}

void MLNodeRobinLaplacian::resizeMultiGrid(int new_size)
{
    MLNodeLinOp::resizeMultiGrid(new_size);
    if (new_size>0 && new_size<static_cast<int>(m_sigma.size())) { m_sigma.resize(new_size); }
}
} // namespace cht
