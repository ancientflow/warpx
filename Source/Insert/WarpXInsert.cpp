#include "WarpXInsert.h"

#include "InsertBackgroundDensity.h"
#include "InsertBoundaryParticles.h"
#include "InsertBoundaryPhi.h"
#include "InsertInjection.h"
#include "InsertRuntimeDiagnostics.h"
#include "WarpXFunctionConfig.h"
#include "WarpXSimulationConfig.h"

/**
 * 粒子注入入口
 */
namespace Insert {

void
ParticleInjection () {
#ifdef HALL3D
    CathodeInjection3D();
    XeInjection();
#endif
#ifdef HALL3D_INIT
    XeFastInjection();
#endif
}

/**
 * 电势修正入口
 */
void
PhiAdjustmentEntrance () {
#ifdef BENCHMARK_2D
    VoltageAdjustment();
#endif
#ifdef HALL3D
    // GetPhiFromFile();
#endif
}

/**
 * 边界电势设置入口
 */
void
SetBoundaryPhi () {
#ifdef HALL3D
    AnodeVoltage();
#endif
}

/**
 * 自定义诊断入口
 */
void
BeforeStep () {
#ifdef NUMP
    ParticleNumber();
#endif
}

/**
 * 自定义初始化入口
 */
void
Initialize () {
#ifdef PUSH_GAP
    // PushGapInit();
#endif
#ifdef MCC_DENSITY
    GlobalBackgroundDensityInit();
#endif
#ifdef HALL3D
    PlasmaInit();
#endif
}

/**
 * 碰撞记录入口
 */
#ifdef COLLISION_RECORD
void
CollisionRecord (amrex::Vector<int> vec) {
    ShowAndWriteIonzationNum(vec);
}
#endif

void
AfterDiagnostics () {
#ifdef HALL3D
    SecondaryEmission();
    AnodeCurrentCalc();
#endif
}

/**
 * 共置网格下，对于第一类边界条件的guard cell设置
 */
void
SetPhiGuards () {
#ifndef WAVE1D
    DirichletPhiGuardSet();
#endif
}

/**
 * 碰撞前
 */
void
BeforeCollision (int step) {
#ifdef MCC_DENSITY
    GlobalBackgroundDensityUpdate(step);
#endif
}

/**
 * 碰撞前
 */
void
AfterCollision (int step) {
#ifdef MCC_DENSITY
    GlobalBackgroundDensityClean(step);
#endif
}

} // namespace Insert
