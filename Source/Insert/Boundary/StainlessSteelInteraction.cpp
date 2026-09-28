#include "Insert/Boundary/StainlessSteelInteraction.h"

#include "Particles/ParticleCreation/DefaultInitialization.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartUtils.H"
#include "Particles/WarpXParticleContainer.H"

#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Reduce.H>
#include <AMReX_Scan.H>

#include <limits>
#include <sstream>

namespace Insert {

StainlessSteelInteraction::StainlessSteelInteraction (int const binomial_trials)
    : m_trials(binomial_trials)
{
    if (m_trials < 1 || m_trials > spectrum_count) {
        amrex::Abort("Stainless steel binomial_trials must be between 1 and 10.");
    }
}

#if defined(WARPX_DIM_3D)
void StainlessSteelInteraction::Process (
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
    amrex::ParticleReal const target_charge = destination.getCharge();
    amrex::Gpu::DeviceVector<amrex::Long> counts(nhits), offsets(nhits);
    auto* const count = counts.dataPtr();
    auto* const offset = offsets.dataPtr();

    amrex::Gpu::DeviceVector<Parameters> parameters(nhits);
    amrex::Gpu::DeviceVector<Event> events(nhits);
    auto* const parameter = parameters.dataPtr();
    auto* const event = events.dataPtr();
    amrex::ParallelFor(nhits, [=] AMREX_GPU_DEVICE (int const k) noexcept {
        parameter[k] = model.Evaluate(hit[k].velocity, hit[k].normal);
    });
    // 有效项返回 -1，无效项返回紧凑索引；最大值定位一个无效撞击。
    // 归约结果返回主机后再决定是否中止，检查发生在当前 tile 的粒子修改之前。
    int const bad = amrex::Reduce::Max<int>(nhits,
        [=] AMREX_GPU_DEVICE (int const k) noexcept {
            return parameter[k].status == Status::valid ? -1 : k;
        });
    if (bad >= 0) {
        Parameters p;
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, parameters.begin()+bad,
            parameters.begin()+bad+1, &p);
        std::ostringstream message;
        message << "Invalid stainless-steel wall impact: ";
        switch (p.status) {
        case Status::probabilities:
            message << "elastic + rediffused probabilities exceed one";
            break;
        case Status::binomial_mean:
            message << "conditional mean must be less than binomial_trials";
            break;
        default: break;
        }
        message << "; E=" << p.energy << " eV, cos(theta)=" << p.cosine
                << ", delta_e=" << p.elastic << ", delta_r=" << p.rediffused
                << ", delta_ts=" << p.true_yield << ", mu=" << p.mean
                << ", M=" << m_trials << ", level=" << lev << ", compact index=" << bad;
        amrex::Abort(message.str());
    }
    amrex::ParallelForRNG(nhits, [=] AMREX_GPU_DEVICE (
        int const k, amrex::RandomEngine const& engine) noexcept {
        event[k] = model.SelectEvent(parameter[k], engine);
        count[k] = event[k].count;
        WallInteractionRecord result;
        result.event = event[k].kind;
        // Common cleanup deletes the incident after all products exist.
        result.keep_primary = false;
        record[k] = result;
    });
    // 前缀和为每个事件分配互不重叠的写入区间，added 是实际新增粒子总数。
    amrex::Long const added = amrex::Scan::ExclusiveSum(nhits, count, offsets.data());
    amrex::Gpu::streamSynchronize(); // finish source access before resizing
    if (added == 0) { return; }

    SmartCopyFactory const factory(source, destination);
    auto const Copy = factory.getSmartCopy();
    auto& target_tile = destination.DefineAndReturnParticleTile(
        lev, pti.index(), pti.LocalTileIndex());
    amrex::Long const old_size = target_tile.numParticles();
    if (old_size + added > std::numeric_limits<int>::max()) {
        amrex::Abort("Stainless-steel emission exceeds the particle tile index range.");
    }
    target_tile.resize(old_size + added);
    // 源和目标可能相同；resize 后重新获取数据指针，避免使用失效地址。
    auto const src = source.ParticlesAt(lev, pti).getParticleTileData();
    auto const dst = target_tile.getParticleTileData();
    amrex::ParallelForRNG(nhits, [=] AMREX_GPU_DEVICE (
        int const k, amrex::RandomEngine const& engine) noexcept {
        if (count[k] == 0) { return; }
        int const i = source_index[k];
        auto const impact = hit[k];
        // Sample the joint spectrum without changing the selected multiplicity.
        auto const velocities = model.SampleEmission(
            parameter[k], event[k], engine);
        for (int j = 0; j < count[k]; ++j) {
            int const target = static_cast<int>(old_size + offset[k] + j);
            auto const& u = velocities[j];
            // 继承宏粒子权重及共享属性，再覆盖出射速度和剩余时间推进后的坐标。
            Copy(dst, src, i, target, engine);
            dst.m_rdata[PIdx::ux][target] = u.x;
            dst.m_rdata[PIdx::uy][target] = u.y;
            dst.m_rdata[PIdx::uz][target] = u.z;
            dst.m_rdata[PIdx::x][target] = impact.position.x + u.x * impact.remaining_time;
            dst.m_rdata[PIdx::y][target] = impact.position.y + u.y * impact.remaining_time;
            dst.m_rdata[PIdx::z][target] = impact.position.z + u.z * impact.remaining_time;
            // 记录带宏粒子权重的实际出射电荷，供公共流程计算净沉积电荷。
            record[k].outgoing_charge += target_charge * dst.m_rdata[PIdx::w][target];
            ++record[k].emitted_count;
        }
    });
    // 补齐运行时属性并赋予新 ID；同步保证临时数组释放前设备读写已完成。
    ParticleCreation::DefaultInitializeRuntimeAttributes(
        target_tile, destination, old_size, old_size + added);
    setNewParticleIDs(target_tile, old_size, added);
    amrex::Gpu::streamSynchronize();
}
#endif

} // namespace Insert
