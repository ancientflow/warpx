#pragma once
// #define FIELDSIZE
#define MCC_DENSITY
//#define MCC_DENSITY_AVERAGE_CALC
//#define MCC_DENSITY_AVERAGE_USE
//#define IONIZATION_SOURCE_RECORD
//#define IONIZATION_SOURCE_INJECT
#if defined(MCC_DENSITY_AVERAGE_CALC) && defined(MCC_DENSITY_AVERAGE_USE)
#error "MCC_DENSITY_AVERAGE_CALC and MCC_DENSITY_AVERAGE_USE are mutually exclusive"
#endif
#if defined(IONIZATION_SOURCE_RECORD) && defined(IONIZATION_SOURCE_INJECT)
#error "IONIZATION_SOURCE_RECORD and IONIZATION_SOURCE_INJECT are mutually exclusive"
#endif
#if defined(IONIZATION_SOURCE_RECORD) && \
    (defined(MCC_DENSITY_AVERAGE_CALC) || defined(MCC_DENSITY_AVERAGE_USE))
#error "IONIZATION_SOURCE_RECORD replaces the MCC_DENSITY_AVERAGE_* path"
#endif
#if defined(IONIZATION_SOURCE_INJECT) && \
    (defined(MCC_DENSITY_AVERAGE_CALC) || defined(MCC_DENSITY_AVERAGE_USE))
#error "IONIZATION_SOURCE_INJECT replaces the MCC_DENSITY_AVERAGE_* path"
#endif
#if (defined(IONIZATION_SOURCE_RECORD) || defined(IONIZATION_SOURCE_INJECT)) && \
    !defined(WARPX_DIM_3D)
#error "Average ionization source currently supports WARPX_DIM_3D only"
#endif
#if defined(IONIZATION_SOURCE_RECORD) && !defined(MCC_DENSITY)
#error "IONIZATION_SOURCE_RECORD requires MCC_DENSITY"
#endif
#define MCC_DELETE
//#define MCC_EXCITATION
#define PUSH_GAP
#define NUMP
//#define COLLISION_RECORD
// #define BENCHMARK_2D
