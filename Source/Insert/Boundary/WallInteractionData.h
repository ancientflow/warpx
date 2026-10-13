#pragma once

#include "Insert/Boundary/WallBehavior.h"
#include "Insert/Math/ParticleVector.h"

namespace Insert {

/** Compact preprocessing output shared by all interaction handlers. */
struct WallImpact
{
    ParticleVector position;
    ParticleVector velocity;
    ParticleVector normal;
    amrex::Real remaining_time;
};

/** One completed interaction, indexed by the compact candidate index.
 * Handlers finish reflection/emission without invalidating the primary.
 * outgoing_charge includes macro weights and any retained primary's charge.
 * emitted_count counts newly created macro particles, excluding the primary.
 * The common finalizer deposits q_in*w - outgoing_charge, then invalidates
 * primaries with keep_primary == false. No emission buffer size is imposed here.
 */
struct WallInteractionRecord
{
    WallBehavior event = WallBehavior::none;
    bool keep_primary = false;
    int emitted_count = 0;
    amrex::ParticleReal outgoing_charge = 0.0;
};

} // namespace Insert
