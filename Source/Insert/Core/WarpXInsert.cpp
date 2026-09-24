#include "WarpXInsert.h"

#include "Insert/Background/InsertBackgroundDensity.h"
#include "Insert/Boundary/AnalyticBoundaryInteraction.h"
#include "Insert/Boundary/InsertBoundaryPhi.h"
#include "Insert/Collisions/IonizationSourceTable.h"
#include "Insert/Config/WarpXFunctionConfig.h"
#include "Insert/Config/WarpXSimulationConfig.h"
#include "Insert/Diagnostics/InsertRuntimeDiagnostics.h"
#include "Insert/Diagnostics/ZmaxRadialExitStats.h"
#include "Insert/Injection/InsertInjection.h"
#include "Utils/TextMsg.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <initializer_list>
#include <map>
#include <string>

namespace {
std::map<std::string, int> particle_subcycling_ndt;
}

namespace Insert {

void
BackwardCompatibility ()
{
    amrex::ParmParse const pp_mc("my_constants");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!pp_mc.contains("anode_current_path"),
        "my_constants.anode_current_path has been removed. Use "
        "my_constants.anode_current_prefix with insert.analytic_walls; "
        "the diagnostic writes one file per wall.");

    for (char const* name : {"zmin_wall_charge_diag", "zmin_wall_charge_dir",
                            "zmin_wall_charge_interval", "zmin_wall_charge_write_interval"}) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!pp_mc.contains(name),
            std::string("my_constants.") + name + " has been removed. Use "
            "insert.use_electrostatic_materials and insert.analytic_walls "
            "for persistent wall charge.");
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::ParmParse::getEntries("insert.schur_boundary").empty(),
        "insert.schur_boundary.* has been removed. Use insert.use_electrostatic_materials "
        "with insert.anode_implicit_function and insert.anode_potential_function.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::ParmParse::getEntries("insert.neutral_atom_eb").empty(),
        "insert.neutral_atom_eb.* has been removed. Configure insert.analytic_walls "
        "and <species>.analytic_wall.<wall>.behaviors instead.");
}

/**
 * 粒子注入入口
 */
void
ParticleInjection () {
#if defined(WARPX_DIM_XZ) && defined(BENCHMARK_2D)
    Benchmark2DParticleInjection();
#elif defined(HALL3D) || defined(HALL3D_INIT)
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
#if defined(WARPX_DIM_XZ) && defined(BENCHMARK_2D)
    VoltageAdjustment();
#endif
}

/**
 * 边界电势设置入口
 */
void
SetBoundaryPhi () {
#ifdef PHT
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
    // Must run before any diagnostic that clears the boundary buffer
    // (e.g. ClearHallBoundaryParticleCache), otherwise the zmax exit
    // statistics see an already-cleared buffer and record nothing.
    ZmaxRadialExitStatsCalc();
    // Written outside the HALL3D block: the diagnostic accumulates inside the
    // analytic wall interaction, which is guarded by WARPX_DIM_3D only.
    AnodeCurrentDiagOutput();
#ifdef HALL3D
    ThrustCalc();
    BeamDivergenceCalc();
    IEDFCalc();
    ClearHallBoundaryParticleCache();
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
