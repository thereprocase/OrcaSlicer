// snuggle_constants.hpp -- Shared constants for the Snuggle nester
//
// Single source of truth for values used across the nester, CPU evaluator,
// and GPU evaluator. Avoids silent mismatch from duplicated definitions.

#pragma once

namespace snuggle {

constexpr float PI_F      = 3.14159265358979f;
constexpr float TWO_PI_F  = 2.0f * PI_F;

// Number of rotation cache bins (1-degree increments over full circle).
// Used by: rotation cache builder, CPU evaluator angle-to-bin,
// GPU evaluator angle-to-bin, GPU metadata indexing.
constexpr int ROT_CACHE_BINS = 360;

} // namespace snuggle
