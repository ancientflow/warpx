#pragma once

namespace Insert {

/** \brief Radially binned statistics of particles leaving through zmax.
 *
 *  Reads the my_constants.zmax_radial_* parameters and, when enabled,
 *  accumulates per-radius-bin statistics of the configured species scraped
 *  at the zmax domain boundary (requires <species>.save_particles_at_zhi=1):
 *    - density distribution: exiting weight, macro-particle count and number
 *      flux per radial annulus;
 *    - velocity distribution: per radial bin, three one-dimensional
 *      histograms of vx, vy and vz.
 *  Sampling happens every hall_diag_interval steps (DoBoundaryParticleDiag)
 *  starting from zmax_radial_start_step; the zhi boundary buffer is cleared
 *  at every sampled step so that statistics stay continuous and the buffer
 *  does not grow unbounded. Results are written to files on the final step.
 *
 *  Only available for 3D builds; a no-op otherwise. Must be called after the
 *  boundary buffer has been gathered (e.g. from Insert::AfterDiagnostics).
 */
void ZmaxRadialExitStatsCalc ();

} // namespace Insert
