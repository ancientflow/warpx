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

/** Optional in-place preparation of phi before each level's multigrid solve.
 *
 * Called synchronously after interpolation from the current coarse solution (if any) and
 * before the native or custom level solver. Only valid unknown nodes may change;
 * fixed potentials and ghosts must be preserved. An empty callback leaves phi as is.
 * History and captured state belong to the caller and must outlive computePhi.
 */
using PoissonLevelPrepare = std::function<void (int, amrex::MultiFab&)>;

/** Optional per-level solve used by computePhi; an empty callback uses the native solver.
 *
 * Called synchronously after coarse interpolation and RHS scaling by -1/epsilon0.
 * RHS is a private working copy; masking it cannot change the deposited rho.
 * The callback restores its fixed potentials, configures its MLMG, synchronizes
 * the RHS, and solves this level. Boundary dual-volume weights may be applied
 * to the private RHS. It must not repeat the -1/epsilon0 conversion, interpolate
 * to the next AMR level, or commit time history. The caller owns those operations.
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
