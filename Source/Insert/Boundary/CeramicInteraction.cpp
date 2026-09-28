#include "Insert/Boundary/CeramicInteraction.h"

#include "Insert/Math/ThermalVelocity.h"
#include "Particles/ParticleCreation/DefaultInitialization.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartUtils.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Scan.H>

#include <limits>

namespace Insert {

CeramicInteraction::CeramicInteraction (amrex::ParticleReal const secondary_temperature_eV)
    : m_thermal_velocity(static_cast<amrex::ParticleReal>(
          Math::ThermalVelocityFromEV(secondary_temperature_eV, PhysConst::m_e)))
{}

#if defined(WARPX_DIM_3D)
void CeramicInteraction::Process (
    WarpXParticleContainer& source, WarpXParticleContainer& destination,
    WarpXParIter const& pti, amrex::Gpu::DeviceVector<int> const& indices,
    amrex::Gpu::DeviceVector<WallImpact> const& impacts,
    amrex::Gpu::DeviceVector<WallInteractionRecord>& records) const
{
    int const nhits = static_cast<int>(indices.size());
    if (nhits == 0) { return; }
    int const lev = pti.GetLevel();
    auto const model = *this;
    auto const* const source_index = indices.dataPtr();
    auto const* const hit = impacts.dataPtr();
    auto* const record = records.dataPtr();
    auto const input = source.ParticlesAt(lev, pti).getParticleTileData();
    amrex::ParticleReal const source_charge = source.getCharge();
    amrex::ParticleReal const target_charge = destination.getCharge();
    amrex::Gpu::DeviceVector<amrex::Long> counts(nhits), offsets(nhits);
    auto* const count = counts.dataPtr();
    auto* const offset = offsets.dataPtr();

    amrex::ParallelForRNG(nhits, [=] AMREX_GPU_DEVICE (
        int const k, amrex::RandomEngine const& engine) noexcept {
        int const i = source_index[k];
        auto const impact = hit[k];
        WallInteractionRecord result;
        result.event = model.SelectEvent(impact.velocity, engine);
        result.keep_primary = result.event == WallBehavior::specular;
        count[k] = WallEmissionCount(result.event);
        if (result.keep_primary) {
            auto const u = model.SampleReflection(impact.velocity, impact.normal);
            input.m_rdata[PIdx::ux][i] = u.x;
            input.m_rdata[PIdx::uy][i] = u.y;
            input.m_rdata[PIdx::uz][i] = u.z;
            input.m_rdata[PIdx::x][i] = impact.position.x + u.x * impact.remaining_time;
            input.m_rdata[PIdx::y][i] = impact.position.y + u.y * impact.remaining_time;
            input.m_rdata[PIdx::z][i] = impact.position.z + u.z * impact.remaining_time;
            result.outgoing_charge = source_charge * input.m_rdata[PIdx::w][i];
        }
        record[k] = result;
    });
    amrex::Long const added = amrex::Scan::ExclusiveSum(nhits, count, offsets.data());
    amrex::Gpu::streamSynchronize(); // finish source access before resizing
    if (added == 0) { return; }

    SmartCopyFactory const factory(source, destination);
    auto const Copy = factory.getSmartCopy();
    auto& target_tile = destination.DefineAndReturnParticleTile(
        lev, pti.index(), pti.LocalTileIndex());
    amrex::Long const old_size = target_tile.numParticles();
    if (old_size + added > std::numeric_limits<int>::max()) {
        amrex::Abort("Ceramic emission exceeds the particle tile index range.");
    }
    target_tile.resize(old_size + added);
    auto const src = source.ParticlesAt(lev, pti).getParticleTileData();
    auto const dst = target_tile.getParticleTileData();
    amrex::ParallelForRNG(nhits, [=] AMREX_GPU_DEVICE (
        int const k, amrex::RandomEngine const& engine) noexcept {
        if (count[k] == 0) { return; }
        int const i = source_index[k];
        auto const impact = hit[k];
        // One joint sampling call per event; SEE2 retains independent draws.
        auto const velocities = model.SampleEmission(
            record[k].event, impact.velocity, impact.normal, engine);
        for (int j = 0; j < count[k]; ++j) {
            int const target = static_cast<int>(old_size + offset[k] + j);
            auto const& u = velocities[j];
            Copy(dst, src, i, target, engine);
            dst.m_rdata[PIdx::ux][target] = u.x;
            dst.m_rdata[PIdx::uy][target] = u.y;
            dst.m_rdata[PIdx::uz][target] = u.z;
            dst.m_rdata[PIdx::x][target] = impact.position.x + u.x * impact.remaining_time;
            dst.m_rdata[PIdx::y][target] = impact.position.y + u.y * impact.remaining_time;
            dst.m_rdata[PIdx::z][target] = impact.position.z + u.z * impact.remaining_time;
            record[k].outgoing_charge += target_charge * dst.m_rdata[PIdx::w][target];
            ++record[k].emitted_count;
        }
    });
    ParticleCreation::DefaultInitializeRuntimeAttributes(
        target_tile, destination, old_size, old_size + added);
    setNewParticleIDs(target_tile, old_size, added);
    amrex::Gpu::streamSynchronize();
}
#endif

} // namespace Insert
