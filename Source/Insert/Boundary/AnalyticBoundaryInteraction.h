#pragma once

namespace Insert {

/** \brief Driver entry for particle interactions with an analytic boundary.
 *
 *  Reads the insert.analytic_boundary.* parameters and, when enabled,
 *  applies the configured wall behavior (absorb / specular / diffuse) in
 *  place to the live particles of the configured species. Must be called
 *  after the particle push and before Redistribute, so that invalidated
 *  (absorbed) particles are removed by the redistribution.
 */
void AnalyticBoundaryInteraction ();

} // namespace Insert
