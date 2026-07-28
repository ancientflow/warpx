#include "WarpXInsert.h"

#include "Insert/Background/InsertBackgroundDensity.h"
#include "Insert/Boundary/InsertBoundaryParticles.h"
#include "Insert/Boundary/InsertBoundaryPhi.h"
#include "Insert/Collisions/IonizationSourceTable.h"
#include "Insert/Config/WarpXFunctionConfig.h"
#include "Insert/Config/WarpXSimulationConfig.h"
#include "Insert/Diagnostics/InsertRuntimeDiagnostics.h"
#include "Insert/Injection/InsertInjection.h"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <map>
#include <string>

namespace {
std::map<std::string, int> particle_subcycling_ndt;
}

/**
 * 粒子注入入口
 */
namespace Insert {

void
ParticleInjection () {
#if defined(HALL3D) || defined(HALL3D_INIT)
    InjectHallParticles();
#endif
}

void
ReadParticleSubcycling (std::string const& species_name,
                        amrex::ParmParse const& pp_species) {
    int ndt = 1;
    pp_species.query("ndt", ndt);
    particle_subcycling_ndt[species_name] = ndt;
}

int
ParticleSubcyclingNdt (std::string const& species_name) {
    auto const iter = particle_subcycling_ndt.find(species_name);
    if (iter != particle_subcycling_ndt.end()) {
        return iter->second;
    }
    return 1;
}

void
ApplyParticleSubcycling (std::string const& species_name, int step,
                         amrex::Real& dt, bool& do_not_push) {
    int const ndt = ParticleSubcyclingNdt(species_name);

    if (step % ndt == 0) {
        amrex::Print() << "push species: " << species_name << "\n";
        dt *= ndt;
        do_not_push = false;
    } else {
        do_not_push = true;
    }
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
    InitializeHallInjection();
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
    NeutralAtomEBInteraction();
    // SecondaryEmission();
    // AnodeIonNeutralization();
    AnodeCurrentCalc();
    ZMinWallChargeDeposit();
    ThrustCalc();
    BeamDivergenceCalc();
    IEDFCalc();
    ClearHallBoundaryParticleCache();
#endif
}

/**
 * 共置网格下，对于平板霍尔推力器zmin电势边界的guard cell设置
 */
void
SetPhiGuards () {
#ifdef HALL3D
    HallThrusterPhiGuardSet();
#elif !defined(WAVE1D)
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

void
Finalize () {
#ifdef IONIZATION_SOURCE_RECORD
    IonizationSourceFinalize();
#endif
#ifdef MCC_DENSITY
    GlobalBackgroundDensityFinalize();
#endif
}

} // namespace Insert
