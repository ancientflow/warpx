#pragma once

/*
 * Insert simulation-mode switches derived from WarpX dimensionality.
 *
 * Uncomment a disabled #define to enable the corresponding optional path.
 */

/* 3D Insert simulation hooks. */
#ifdef WARPX_DIM_3D
#define WARPX_INSERT_3D
// #define WARPX_INSERT_INITIALIZATION_3D
#endif
