#pragma once
// #define FIELDSIZE
#define MCC_DENSITY
//#define MCC_DENSITY_AVERAGE_CALC
//#define MCC_DENSITY_AVERAGE_USE
#if defined(MCC_DENSITY_AVERAGE_CALC) && defined(MCC_DENSITY_AVERAGE_USE)
#error "MCC_DENSITY_AVERAGE_CALC and MCC_DENSITY_AVERAGE_USE are mutually exclusive"
#endif
#define MCC_DELETE
//#define MCC_EXCITATION
#define PUSH_GAP
#define NUMP
//#define COLLISION_RECORD
// #define BENCHMARK_2D
