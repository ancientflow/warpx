#pragma once

#include "Insert/Boundary/WallBehavior.h"
#include "Insert/Boundary/WallInteractionData.h"
#include "Insert/Boundary/WallInteractionOperators.h"
#include "Insert/Math/ParticleVector.h"
#include "Utils/WarpXConst.H"

#include "Particles/WarpXParticleContainer_fwd.H"

#include <AMReX_Array.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Random.H>

#include <cmath>
#include <type_traits>

namespace Insert {

/** Electron/ceramic interaction from Script/3d_hall_dielectric_slot.
 * Owns event selection, reflection, allocation and secondary emission.
 * Deposition, primary deletion and diagnostics remain in the common finalizer.
 * Input velocities are physical m/s.
 */
class CeramicInteraction
{
public:
    /** 构造陶瓷预设，缓存 sqrt(kT/m_e) 对应的热速度（m/s）。
     * secondary_temperature_eV 是 kT（eV），不是出射电子的平均动能。
     */
    explicit CeramicInteraction (amrex::ParticleReal secondary_temperature_eV = 3.0);

    static constexpr int max_emitted = 2;

#if defined(WARPX_DIM_3D)
    /** 主机入口：处理当前粒子 tile 中已经完成碰撞点反推的粒子。
     * indices[k] 是入射粒子在 source tile 中的索引；impacts[k] 包含撞击点、
     * 撞击速度、指向计算域的法向和碰撞后剩余时间。records 须按 indices 等长分配。
     * 镜面反射直接更新原粒子；SEE 在 destination 中追加粒子，允许源和目标相同。
     * 返回前完成记录与粒子生成；电荷沉积、入射粒子失效标记和诊断由公共流程完成。
     */
    void Process (
        WarpXParticleContainer& source, WarpXParticleContainer& destination,
        WarpXParIter const& pti, amrex::Gpu::DeviceVector<int> const& indices,
        amrex::Gpu::DeviceVector<WallImpact> const& impacts,
        amrex::Gpu::DeviceVector<WallInteractionRecord>& records) const;
#endif

    /** 根据撞击速度（m/s）计算非相对论动能（eV），一次抽样确定事件。
     * 按吸收、镜面反射、双电子发射的累积概率判断，剩余概率为单电子发射。
     * 本预设概率只依赖能量；返回值的发射数由 WallEmissionCount 给出。
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    WallBehavior SelectEvent (
        ParticleVector const& hit_velocity, amrex::RandomEngine const& engine) const noexcept;

    /** Sample the whole emission event, retaining model-specific spectra and
     * correlations. Only WallEmissionCount(event) entries are used. This ceramic
     * preset draws each secondary independently from a Maxwellian flux at kT,
     * with mean kinetic energy 2 kT and no event-wise energy truncation.
     * normal 指向计算域，内部归一化；返回实验室坐标系速度（m/s）。
     * hit_velocity 传给漫反射算子，本模型的热出射速度不继承其动能。
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuArray<ParticleVector, max_emitted> SampleEmission (
        WallBehavior event, ParticleVector const& hit_velocity,
        ParticleVector const& normal, amrex::RandomEngine const& engine) const noexcept;

    /** 对撞击点速度作镜面反射：v_out = v_in - 2(v_in dot n)n。
     * normal 指向计算域，内部归一化；保持速率与切向速度，不移动粒子。
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleVector SampleReflection (
        ParticleVector const& hit_velocity, ParticleVector const& normal) const noexcept;

private:
    amrex::ParticleReal m_thermal_velocity; ///< 热速度 sqrt(kT/m_e)，单位 m/s。
};

static_assert(std::is_trivially_copyable_v<CeramicInteraction>);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
WallBehavior CeramicInteraction::SelectEvent (
    ParticleVector const& hit_velocity, amrex::RandomEngine const& engine) const noexcept
{
    amrex::ParticleReal const speed_squared = Math::NormSquared(hit_velocity);

    amrex::ParticleReal const energy_eV = PhysConst::m_e * speed_squared / (2.0 * PhysConst::q_e);
    amrex::ParticleReal const a = energy_eV / 43.4592;
    amrex::ParticleReal const r = energy_eV / 30.0;
    amrex::ParticleReal const s = energy_eV / 127.8958;
    amrex::ParticleReal const p_absorb = 0.5 * std::exp(-a*a);
    amrex::ParticleReal const p_specular = 0.5 * std::exp(-r*r);
    // expm1 避免低能量时直接计算 1-exp(-s*s) 的相消误差。
    amrex::ParticleReal const p_see2 = -std::expm1(-s*s);
    amrex::ParticleReal const draw = amrex::Random(engine);
    if (draw < p_absorb) {
        return WallBehavior::absorb;
    }
    if (draw < p_absorb + p_specular) {
        return WallBehavior::specular;
    }
    return draw < p_absorb + p_specular + p_see2 ?
        WallBehavior::secondary_electron_2 : WallBehavior::secondary_electron_1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<ParticleVector, CeramicInteraction::max_emitted>
CeramicInteraction::SampleEmission (
    WallBehavior const event, ParticleVector const& hit_velocity,
    ParticleVector const& normal_to_domain, amrex::RandomEngine const& engine) const noexcept
{
    ParticleVector normal = normal_to_domain;
    Math::Normalize(normal_to_domain, normal);
    amrex::GpuArray<ParticleVector, max_emitted> velocities{};
    int const count = WallEmissionCount(event);
    for (int j = 0; j < count; ++j) {
        DiffuseReemissionOperator{m_thermal_velocity}(
            normal, hit_velocity, velocities[j], engine);
    }
    return velocities;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ParticleVector CeramicInteraction::SampleReflection (
    ParticleVector const& hit_velocity, ParticleVector const& normal_to_domain) const noexcept
{
    ParticleVector normal = normal_to_domain;
    Math::Normalize(normal_to_domain, normal);
    return Math::Reflect(hit_velocity, normal);
}

} // namespace Insert
