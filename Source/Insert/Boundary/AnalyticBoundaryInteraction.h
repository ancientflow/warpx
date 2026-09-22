#pragma once

namespace Insert {

/** \brief Driver entry for particle interactions with an analytic boundary.
 *
 *  Reads the insert.analytic_walls registry and the per-species
 *  <species>.analytic_wall.<wall>.behaviors policies. It must be called after
 *  the particle push and before Redistribute, so invalidated particles are
 *  removed by redistribution.
 */
void AnalyticBoundaryInteraction ();

/** \brief Write the wall-current diagnostic and reset its accumulator.
 *
 *  The per-step accumulation happens inside AnalyticBoundaryInteraction,
 *  with one electron/ion/net-charge triplet per analytic wall. This function
 *  runs in the diagnostics phase (AfterDiagnostics); on steps selected by
 *  DoBoundaryParticleDiag it MPI-reduces the accumulated wall currents,
 *  appends one row per wall to <my_constants.anode_current_prefix>_<wall>.dat
 *  (default prefix "anode_current") and restarts the accumulation window.
 *  Enabled by my_constants.anode_current_diag.
 */
void AnodeCurrentDiagOutput ();

} // namespace Insert
