#pragma once

/*
 * Insert compile-time feature switches.
 *
 * Uncomment a disabled #define to enable the corresponding optional path.
 */

/* Diagnostics. */
// #define FIELDSIZE
#define NUMP
// #define COLLISION_RECORD
#define BENCHMARK_2D

/* Background MCC density coupling. AVERAGE_USE loads an offline-averaged FAB. */
#define MCC_DENSITY
// #define MCC_DENSITY_AVERAGE_USE

/* Coupled MCC collision behavior. */
#define MCC_DELETE
// #define MCC_EXCITATION

/*
 * Averaged ionization-source workflow.
 *
 * IONIZATION_SOURCE_RECORD writes source tables from a 3D MCC run.
 * IONIZATION_SOURCE_INJECT reads source tables and injects particles from them.
 */
// #define IONIZATION_SOURCE_RECORD
// #define IONIZATION_SOURCE_INJECT

#include "Insert/Config/WarpXOptimizationConfig.h"

/* Configuration checks. */
#ifdef MCC_DENSITY_AVERAGE_CALC
#error "MCC_DENSITY_AVERAGE_CALC has been removed. Average rho_<species> diagnostics offline instead."
#endif

#if defined(MCC_ION_CACHED_RHO) && defined(MCC_DENSITY_AVERAGE_USE)
#error "MCC_ION_CACHED_RHO requires direct coupled MCC density updates; disable MCC_DENSITY_AVERAGE_USE"
#endif

#if defined(IONIZATION_SOURCE_INJECT) && !defined(MCC_DENSITY_AVERAGE_USE)
#error "IONIZATION_SOURCE_INJECT requires MCC_DENSITY_AVERAGE_USE"
#endif

#if defined(IONIZATION_SOURCE_RECORD) && defined(IONIZATION_SOURCE_INJECT)
#error "IONIZATION_SOURCE_RECORD and IONIZATION_SOURCE_INJECT are mutually exclusive"
#endif

#if (defined(IONIZATION_SOURCE_RECORD) || defined(IONIZATION_SOURCE_INJECT)) && \
    !defined(WARPX_DIM_3D)
#error "Average ionization source currently supports WARPX_DIM_3D only"
#endif

#if defined(IONIZATION_SOURCE_RECORD) && !defined(MCC_DENSITY)
#error "IONIZATION_SOURCE_RECORD requires MCC_DENSITY"
#endif
