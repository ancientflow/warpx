#pragma once

/*
 * Insert compile-time feature switches.
 *
 * Uncomment a disabled #define to enable the corresponding optional path.
 */

/* Diagnostics and benchmark hooks. */
// #define FIELDSIZE
#define NUMP
// #define COLLISION_RECORD
// #define BENCHMARK_2D

/* Background MCC density coupling. */
#define MCC_DENSITY
// #define MCC_DENSITY_AVERAGE_CALC
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

/* Particle push / injection helpers. */
#define PUSH_GAP

/* Configuration checks. */
#if defined(MCC_DENSITY_AVERAGE_CALC) && defined(MCC_DENSITY_AVERAGE_USE)
#error "MCC_DENSITY_AVERAGE_CALC and MCC_DENSITY_AVERAGE_USE are mutually exclusive"
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
