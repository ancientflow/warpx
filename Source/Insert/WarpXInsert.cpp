#include "WarpXInsert.h"
#include "WarpXInsertFunction.h"
#include "AMReX_Vector.H"
#include "WarpXSimulationFunction.h"

/**
 * 粒子注入入口
 */
void ParticleInjectionEntrance () {
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
void PhiAdjustmentEntrance () {
#ifdef BENCHMAKR2D
    VoltageAdjustment(warpx_instance);
#endif
#ifdef HALL3D
    //GetPhiFromFile();
#endif
}

/**
 * 边界电势设置入口
 */
void BoundaryPhiSetEntrance () {
#ifdef HALL3D
    AnodeVoltage();
#endif
}


/**
 * 自定义诊断入口
 */
void MyDiag () 
{
#ifdef NUMP
    ParticleNumber();
#endif
}

/**
 * 自定义输出入口
 */
void MyOutput()
{
#ifdef HALL3D
    AnodeCurrentCalc();
#endif
}

/**
 * 碰撞记录入口
 */
#ifdef COLLISION_RECORD
void CollisionRecord(amrex::Vector<int> vec)
{
    ShowAndWriteIonzationNum(vec);
}
#endif

/**
 * 数据检验入口
 */
void DataExamine()
{
#ifdef PHIEXAMINE
    PhiBCExamine();
#endif
}

/**
 * 自定义初始化入口
 */
void MyInit()
{
#ifdef PUSH_GAP
    //PushGapInit();
#endif
#ifdef MCC_DENSITY
    GlobalBackgroundDensityInit();
#endif
}

/**
 * 共置网格下，对于第一类边界条件的guard cell设置
 */
void PhiGuardSetEntrance()
{
#ifndef WAVE1D
    DirichletPhiGuardSet();
#endif
}

/**
 * 碰撞前
 */
void
BeforeCollision (int step, bool if_split) {
#ifdef MCC_DENSITY
    GlobalBackgroundDensityUpdate(step, if_split);
#endif
}

/**
 * 碰撞前
 */
void
AfterCollision (int step, bool if_split) {
#ifdef MCC_DENSITY
    GlobalBackgroundDensityClean(step, if_split);
#endif
}