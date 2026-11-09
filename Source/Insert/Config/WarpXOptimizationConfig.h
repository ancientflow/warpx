#pragma once

/*
 * Insert compile-time switches for computational efficiency optimizations.
 *
 * Uncomment a disabled #define to enable the corresponding optional path.
 */

/* Cache ion charge density for subcycled MCC ion species. */
// #define MCC_ION_CACHED_RHO

/*
 * Nonrelativistic explicit Boris push for all massive particle species.
 * Treat the stored u as velocity during the push (gamma = 1). The existing
 * split momentum push, including its half-angle magnetic rotation, is retained.
 * This switch only changes the pusher; coupled MCC is always nonrelativistic.
 * Injection and diagnostics retain their conventions for low-energy simulations.
 * Uncomment and rebuild to enable; leave undefined for the relativistic push.
 */
#define WARPX_NONRELATIVISTIC_PUSH
