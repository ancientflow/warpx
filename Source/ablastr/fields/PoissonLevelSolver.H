#ifndef ABLASTR_FIELDS_POISSON_LEVEL_SOLVER_H_
#define ABLASTR_FIELDS_POISSON_LEVEL_SOLVER_H_

#include <AMReX_REAL.H>

#include <functional>

namespace amrex {
    struct LPInfo;
    class MultiFab;
    template <typename MF> class MLMGT;
}

namespace ablastr::fields {

struct MLMGOptions;

/** Optional per-level solve used by computePhi; an empty callback uses the native solver.
 *
 * Called synchronously after coarse interpolation and RHS scaling by -1/epsilon0.
 * The callback restores its fixed potentials, configures its MLMG, synchronizes
 * the RHS, and solves this level. It must not rescale the RHS, interpolate to the
 * next AMR level, or commit time history. The caller still owns those operations.
 *
 * The returned MLMG and its operator must outlive the current level's post-solve
 * callback (at least until the next invocation). Captured state is borrowed only
 * during computePhi. No material state is allocated when this callback is empty.
 * Supported here for non-EB multigrid solves; FFT/IGF retains its native path.
 */
using PoissonLevelSolver = std::function<amrex::MLMGT<amrex::MultiFab>& (
    int, amrex::LPInfo const&, amrex::MultiFab&, amrex::MultiFab&,
    MLMGOptions const&, amrex::Real)>;

} // namespace ablastr::fields

#endif
