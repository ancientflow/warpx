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

} // namespace Insert
