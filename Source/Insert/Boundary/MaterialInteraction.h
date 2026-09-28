#pragma once

#include "Insert/Boundary/CeramicInteraction.h"

#include <variant>

namespace Insert {

/** Selects the material-owned host Process entry after common preprocessing.
 * Each alternative manages its own allocation and device emission kernels.
 * Only its trivially copyable sampler is captured inside those kernels.
 */
using MaterialInteraction = std::variant<CeramicInteraction>;

} // namespace Insert
