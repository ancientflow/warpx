#include "NodalPoissonInitialGuess.H"

#include <AMReX_GpuLaunch.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>
#include <limits>

namespace ablastr::fields {
namespace {
bool sameInterval (amrex::Real a, amrex::Real b, amrex::Real time)
{
    auto const scale = std::max({std::abs(a), std::abs(b), std::abs(time)});
    return std::abs(a-b) <= 64*std::numeric_limits<amrex::Real>::epsilon()*scale;
}

bool sameGeometry (amrex::Geometry const& a, amrex::Geometry const& b)
{
    if (a.Domain() != b.Domain() || a.Coord() != b.Coord()) { return false; }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (a.ProbLo(d) != b.ProbLo(d) || a.ProbHi(d) != b.ProbHi(d) ||
            a.isPeriodic(d) != b.isPeriodic(d)) { return false; }
    }
    return true;
}

amrex::Real gaugeOffset (amrex::MultiFab const& phi, amrex::Geometry const& geom)
{
    auto owner = amrex::OwnerMask(phi, geom.periodicity());
    auto const domain = amrex::surroundingNodes(geom.Domain());
    auto const periodic = geom.isPeriodicArray();
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> ops;
    amrex::ReduceData<amrex::Real, amrex::Real> data(ops);
    using Tuple = decltype(data)::Type;
    for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto const p = phi.const_array(mfi);
        auto const o = owner->const_array(mfi);
        ops.eval(mfi.validbox(), data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> Tuple {
                amrex::IntVect const iv(AMREX_D_DECL(i,j,k));
                amrex::Real weight = o(i,j,k);
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    if (!periodic[d] &&
                        (iv[d] == domain.smallEnd(d) || iv[d] == domain.bigEnd(d)))
                    { weight *= amrex::Real(0.5); }
                }
                return {weight*p(i,j,k), weight};
            });
    }
    auto const sums = data.value();
    amrex::Real total[2] = {amrex::get<0>(sums), amrex::get<1>(sums)};
    amrex::ParallelDescriptor::ReduceRealSum(total, 2);
    return total[1] > 0 ? total[0]/total[1] : 0;
}
}

void NodalPoissonInitialGuess::beginSolve (
    MultiLevelScalarField const& phi, amrex::Vector<amrex::Geometry> const& geom,
    BC const& lobc, BC const& hibc, amrex::Vector<amrex::iMultiFab const*> const& overset,
    amrex::Real time, int order, amrex::Real alpha, amrex::Real beta)
{
    AMREX_ALWAYS_ASSERT(order >= 0 && order <= 2 && std::isfinite(time) &&
                        std::isfinite(alpha) && std::isfinite(beta));
    AMREX_ALWAYS_ASSERT(geom.size() >= phi.size() &&
                        (overset.empty() || overset.size() == phi.size()));
    m_time = time;
    m_order = order;
    m_alpha = alpha;
    m_beta = beta;
    m_levels.resize(phi.size());
    for (int lev = 0; lev < static_cast<int>(phi.size()); ++lev) {
        AMREX_ALWAYS_ASSERT(phi[lev] && phi[lev]->is_nodal() && phi[lev]->nComp() == 1 &&
                            geom[lev].Coord() == 0);
        auto const& field = *phi[lev];
        auto const* mask = overset.empty() ? nullptr : overset[lev];
        if (mask) {
            AMREX_ALWAYS_ASSERT(mask->boxArray() == field.boxArray() && mask->nComp() == 1 &&
                                mask->DistributionMap() == field.DistributionMap());
        }
        auto& l = m_levels[lev];
        bool const rebuild = !l.unknown || l.unknown->boxArray() != field.boxArray() ||
            l.unknown->DistributionMap() != field.DistributionMap() ||
            !sameGeometry(l.geom, geom[lev]) || l.lobc != lobc || l.hibc != hibc ||
            l.overset != mask;
        if (rebuild) {
            l.geom = geom[lev];
            l.lobc = lobc;
            l.hibc = hibc;
            l.overset = mask;
            l.depth = 0;
            l.unknown = std::make_unique<amrex::iMultiFab>(
                field.boxArray(), field.DistributionMap(), 1, 0);
            amrex::iMultiFab coverage(amrex::convert(field.boxArray(), amrex::IntVect(0)),
                                     field.DistributionMap(), 1, 1);
            // Same values and periodic treatment as MLNodeLinOp::buildMasks:
            // covered=0, uncovered interior=1, outside physical domain=2.
            coverage.BuildMask(geom[lev].Domain(), geom[lev].periodicity(), 0, 1, 2, 0);
            auto const domain = amrex::surroundingNodes(geom[lev].Domain());
            amrex::GpuArray<amrex::LinOpBCType, AMREX_SPACEDIM> lo{}, hi{};
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { lo[d] = lobc[d]; hi[d] = hibc[d]; }
            for (amrex::MFIter mfi(*l.unknown); mfi.isValid(); ++mfi) {
                auto const cells = coverage.const_array(mfi);
                auto const dst = l.unknown->array(mfi);
                auto const fixed = mask ? mask->const_array(mfi) : amrex::Array4<int const>{};
                amrex::ParallelFor(mfi.validbox(),
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        amrex::IntVect const node(AMREX_D_DECL(i,j,k));
                        bool unknown = !fixed || fixed(i,j,k) != 0;
                        for (int c = 0; c < (1 << AMREX_SPACEDIM); ++c) {
                            auto cell = node;
                            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                                cell[d] -= (c >> d) & 1;
                            }
                            unknown &= cells(cell) != 1;
                        }
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            unknown &= !(node[d] == domain.smallEnd(d) &&
                                         lo[d] == amrex::LinOpBCType::Dirichlet);
                            unknown &= !(node[d] == domain.bigEnd(d) &&
                                         hi[d] == amrex::LinOpBCType::Dirichlet);
                        }
                        dst(i,j,k) = int(unknown);
                    });
            }
            l.unknown->OverrideSync(geom[lev].periodicity());
            bool all_neumann = true;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                auto natural = [](amrex::LinOpBCType b) {
                    return b == amrex::LinOpBCType::Periodic ||
                           b == amrex::LinOpBCType::Neumann ||
                           b == amrex::LinOpBCType::inhomogNeumann;
                };
                all_neumann &= natural(lobc[d]) && natural(hibc[d]);
            }
            l.zero_mean = all_neumann && l.unknown->min(0) == 1;
        }
        defineHistory(lev, field);
    }
}

void NodalPoissonInitialGuess::prepare (int lev, amrex::MultiFab& phi, int verbosity)
{
    auto const time = m_time;
    auto const order = m_order;
    auto const alpha = m_alpha;
    auto const beta = m_beta;
    auto const& unknown = *m_levels.at(lev).unknown;
    AMREX_ALWAYS_ASSERT(std::isfinite(alpha) && std::isfinite(beta));
    auto& h = m_levels.at(lev);
    Prediction result;
    // CUDA extended device lambdas must be enclosed by a public member function.
    if (h.depth > 0) {
        bool const same_time = time == h.time[0];
        if (!same_time && h.depth > 1 &&
            !sameInterval(time-h.time[0], h.time[0]-h.time[1], time))
        {
            // Keep the latest solution, then accumulate history at the new interval.
            h.depth = 1;
            result.changed_interval = true;
        }
        result.order = same_time ? 0 : std::min(order, h.depth-1);
        if (result.order == 2 && beta == 0) { result.order = 1; }
        if (result.order == 1 && alpha == 0) { result.order = 0; }
        amrex::GpuArray<amrex::Real, 3> weights{1, 0, 0};
        if (result.order == 1) { weights = {1+alpha, -alpha, 0}; }
        if (result.order == 2) { weights = {1+alpha+beta, -(alpha+2*beta), beta}; }

        // Check the entire candidate before writing anything: a rejected prediction
        // must not leave a mixture of extrapolated and fallback values on a level.
        int bad = 0;
        if (weights[1] != 0 || weights[2] != 0) {
            amrex::ReduceOps<amrex::ReduceOpMax> ops;
            amrex::ReduceData<int> data(ops);
            using Tuple = decltype(data)::Type;
            for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi) {
                auto const mask = unknown.const_array(mfi);
                auto const p0 = h.history[0]->const_array(mfi);
                auto const p1 = result.order >= 1 ? h.history[1]->const_array(mfi) : p0;
                auto const p2 = result.order >= 2 ? h.history[2]->const_array(mfi) : p0;
                int const ncomp = phi.nComp();
                ops.eval(mfi.validbox(), data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> Tuple {
                        int bad = 0;
                        if (mask(i,j,k)) {
                            for (int n = 0; n < ncomp; ++n) {
                                auto const value = weights[0]*p0(i,j,k,n) +
                                    weights[1]*p1(i,j,k,n) + weights[2]*p2(i,j,k,n);
                                bad |= !std::isfinite(value);
                            }
                        }
                        return {bad};
                    });
            }
            bad = amrex::get<0>(data.value());
            amrex::ParallelDescriptor::ReduceIntMax(bad);
        }
        if (bad > 0) {
            weights = {1, 0, 0};
            result.order = 0;
            result.finite_fallback = true;
        }
        for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi) {
            auto const mask = unknown.const_array(mfi);
            auto const dst = phi.array(mfi);
            auto const p0 = h.history[0]->const_array(mfi);
            auto const p1 = result.order >= 1 ? h.history[1]->const_array(mfi) : p0;
            auto const p2 = result.order >= 2 ? h.history[2]->const_array(mfi) : p0;
            amrex::ParallelFor(mfi.validbox(), phi.nComp(),
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                    if (mask(i,j,k)) {
                        dst(i,j,k,n) = weights[0]*p0(i,j,k,n) +
                            weights[1]*p1(i,j,k,n) + weights[2]*p2(i,j,k,n);
                    }
                });
        }
    }
    if (verbosity >= 2) {
        amrex::Print() << "Phi initial guess: level=" << lev << " time=" << m_time
                       << " history_order=" << result.order
                       << " changed_interval=" << result.changed_interval
                       << " finite_fallback=" << result.finite_fallback << '\n';
    }
}

void NodalPoissonInitialGuess::commit (MultiLevelScalarField const& phi)
{
    for (int lev = 0; lev < static_cast<int>(phi.size()); ++lev) {
        auto const& l = m_levels.at(lev);
        auto const offset = l.zero_mean ? gaugeOffset(*phi[lev], l.geom) : 0;
        commitHistory(lev, *phi[lev], offset);
    }
}

void NodalPoissonInitialGuess::defineHistory (int lev, amrex::MultiFab const& phi)
{
    auto const time = m_time;
    auto const order = m_order;
    AMREX_ALWAYS_ASSERT(order >= 0 && order <= 2 && std::isfinite(time));
    auto& h = m_levels.at(lev);
    for (int n = 0; n <= order; ++n) {
        auto& f = h.history[n];
        if (!f || f->boxArray() != phi.boxArray() ||
            f->DistributionMap() != phi.DistributionMap() || f->nComp() != phi.nComp())
        {
            f = std::make_unique<amrex::MultiFab>(
                phi.boxArray(), phi.DistributionMap(), phi.nComp(), 0);
            f->setVal(0);
            h.depth = 0;
        }
    }
    for (int n = order+1; n < 3; ++n) { h.history[n].reset(); }
    h.depth = std::min(h.depth, order+1);
    if (h.depth && time < h.time[0]) { h.depth = 0; }
}

void NodalPoissonInitialGuess::commitHistory (
    int lev, amrex::MultiFab const& phi, amrex::Real gauge_offset)
{
    auto const time = m_time;
    auto const order = m_order;
    auto& h = m_levels.at(lev);
    if (!std::isfinite(gauge_offset) || !phi.is_finite(0, phi.nComp(), 0)) {
        h.depth = 0;
        return;
    }
    // Re-solving the same time replaces that sample, without inventing a time step.
    bool const replace = h.depth > 0 && time == h.time[0];
    if (!replace) {
        for (int n = std::min(h.depth, order); n > 0; --n) {
            std::swap(h.history[n], h.history[n-1]);
            h.time[n] = h.time[n-1];
        }
        h.depth = std::min(h.depth+1, order+1);
    }
    amrex::MultiFab::Copy(*h.history[0], phi, 0, 0, phi.nComp(), 0);
    if (gauge_offset != 0) {
        h.history[0]->plus(-gauge_offset, 0, phi.nComp(), 0);
        if (!h.history[0]->is_finite(0, phi.nComp(), 0)) { h.depth = 0; }
    }
    h.time[0] = time;
}

} // namespace ablastr::fields
