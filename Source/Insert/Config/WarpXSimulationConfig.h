#pragma once

/*
 * Insert simulation-mode switches derived from WarpX dimensionality.
 *
 * Uncomment a disabled #define to enable the corresponding optional path.
 */

/* 3D Hall simulation hooks. */
#ifdef WARPX_DIM_3D
#define HALL3D
// #define HALL3D_INIT
#endif

/* type of plane hall thruster */
#ifdef HALL3D
#define PHT        // normal plane hall thruster
#define SLOPE_ANODE // plane hall thruster with slope anode
#endif


/* 1D wave simulation hooks. */
#ifdef WARPX_DIM_1D_Z
#define WAVE1D
#endif
