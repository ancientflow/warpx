#pragma once

#include "Insert/Boundary/WallInteractionData.h"
#include "Insert/Math/ParticleVector.h"
#include "Utils/WarpXConst.H"

#include "Particles/WarpXParticleContainer_fwd.H"

#include <AMReX_Array.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Random.H>

#include <cmath>
#include <limits>
#include <type_traits>

namespace Insert {

/** Furman--Pivi unconditioned SLAC 304 preset (2002), Tables I/II.
 * Uses the penetrated-electron binomial model, Eqs. (39)--(46).
 * M defaults to 10; the tabulated spectra only support multiplicities 1--10.
 * Impacts with invalid probabilities or conditional mean are fully absorbed.
 * No probability clipping, angular saturation or Poisson truncation.
 */
class StainlessSteelInteraction
{
public:
    using Scalar = amrex::ParticleReal;
    static constexpr int spectrum_count = 10;

    enum class Status : int {
        valid, probabilities, binomial_mean
    };

    /** 单次撞击的模型参数；产额为每个入射电子对应的平均电子数。
     * elastic/rediffused 也是各自单电子背散射分支的概率。
     * mean 是进入剩余分支后的条件均值，不是总产额；energy 单位为 eV。
     */
    struct Parameters
    {
        ParticleVector normal{};
        Scalar energy = 0;
        Scalar cosine = 1;
        Scalar elastic = 0;
        Scalar rediffused = 0;
        Scalar true_yield = 0;
        Scalar mean = 0;
        Status status = Status::valid;
    };

    /** 已选定的机制及出射数量；count=0 表示吸收，能谱抽样不再改变 count。 */
    struct Event
    {
        WallBehavior kind = WallBehavior::absorb;
        int count = 0;
    };

    /** 设置二项分布试验次数 M；1<=M<=10，受现有能谱表覆盖范围限制。
     * 条件均值不满足 mu<M 的撞击回退为吸收；这里不实现泊松分布或其截断。
     */
    explicit StainlessSteelInteraction (int binomial_trials = spectrum_count);

#if defined(WARPX_DIM_3D)
    /** 主机入口：参数计算与检查、事件抽样、精确分配及出射粒子生成。
     * indices[k] 映射到 source 的当前 tile；impacts[k] 是对应预处理结果，
     * records 须预先分配同样长度。允许 source 与 destination 为同一容器。
     * 所有出射分支均创建新粒子；返回前同步，原粒子仍有效，供公共流程
     * 按入射电荷减出射电荷沉积，再标记删除并汇总诊断。
     */
    void Process (
        WarpXParticleContainer& source, WarpXParticleContainer& destination,
        WarpXParIter const& pti, amrex::Gpu::DeviceVector<int> const& indices,
        amrex::Gpu::DeviceVector<WallImpact> const& impacts,
        amrex::Gpu::DeviceVector<WallInteractionRecord>& records) const;
#endif

    /** 从物理速度（m/s）及指向计算域的法向计算 E0、入射角和三个产额。
     * 内部归一化法向，cos(theta0)=-v dot n/|v|；假定输入状态有效。
     * 不限制能量或角度的拟合范围，仅通过 status 报告概率及二项均值问题。
     * 此函数不抽随机数、不修改粒子，也不在设备端中止程序。
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Parameters Evaluate (
        ParticleVector const& velocity, ParticleVector const& normal) const noexcept;

    /** p.status 无效时直接吸收，不抽随机数；有效时先选择两个背散射分支。
     * 剩余概率 A 内采样 Binomial(M,mu/M)，其中 n=0 为吸收。
     * A 不是吸收概率：实际吸收概率为 A*(1-mu/M)^M。
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Event SelectEvent (Parameters const& p, amrex::RandomEngine const& engine) const noexcept;

    /** 对 SelectEvent 返回的固定事件采样，返回实验室坐标系速度（m/s）。
     * 仅前 event.count 项有效；弹性为有限宽度能谱，再扩散为幂律能谱，
     * 真二次为给定总能量上限 E0 的联合 Gamma 能谱。所有机制使用余弦角分布。
     * 要求有效的 Parameters 和对应 Event；本函数不分配粒子、不改变发射数。
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuArray<ParticleVector, spectrum_count> SampleEmission (
        Parameters const& p, Event const& event, amrex::RandomEngine const& engine) const noexcept;

private:
    int m_trials; ///< 条件二项分布的试验次数 M，也是单次发射数量上限。

    /** 生成严格位于 (0,1) 的粒子精度随机数，供对数和逆分布变换使用。
     * 转换为 ParticleReal 后再次排除端点，避免低精度舍入得到 1。
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Scalar UniformOpen (amrex::RandomEngine const& engine) noexcept
    {
        Scalar u;
        do {
            u = static_cast<Scalar>(amrex::Random(engine));
        } while (u <= Scalar(0) || u >= Scalar(1));
        return u;
    }

    /** Marsaglia--Tsang 拒绝采样：Gamma(shape,1)，要求 shape>1。
     * 用于总能量及 Dirichlet 权重；现有能谱表的所有 shape 均满足条件。
     * 拒绝只重抽连续随机量，不重新选择碰撞事件或发射数量。
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Scalar Gamma (Scalar shape, amrex::RandomEngine const& engine) noexcept;

    /** 抽样 S~Gamma(shape,scale) 在 S<=limit 下的条件分布。
     * shape>1、scale>0、limit>=0；scale 和 limit 的单位均为 eV。
     * 高上限直接拒绝超限 Gamma；低上限用幂律提议分布避免极低接受率。
     * limit=0 返回零；此处截断的是事件总能量，不是电子数分布。
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Scalar TruncatedGamma (
        Scalar shape, Scalar scale, Scalar limit, amrex::RandomEngine const& engine) noexcept;

    /** 抽样弹性出射能量（eV）：以入射 energy 为中心、sigma=1.9 eV 的
     * 高斯分布限制在 [0,energy] 上；等价于截断半正态能量损失。
     * energy=0 取零能量极限，不用直接返回入射能量的单能近似。
     */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Scalar ElasticEnergy (Scalar energy, amrex::RandomEngine const& engine) noexcept;
};

static_assert(std::is_trivially_copyable_v<StainlessSteelInteraction>);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
StainlessSteelInteraction::Parameters StainlessSteelInteraction::Evaluate (
    ParticleVector const& velocity, ParticleVector const& normal) const noexcept
{
    Parameters p;
    Scalar const v2 = Math::NormSquared(velocity);
    // Preprocessing supplies a valid incident velocity and inward surface normal.
    Math::Normalize(normal, p.normal);
    p.energy = (PhysConst::m_e_v<Scalar> / (Scalar(2)*PhysConst::q_e_v<Scalar>)) * v2;
    // A stationary limiting incident has E=0; take normal incidence for its yield.
    p.cosine = v2 > Scalar(0) ? -Math::Dot(velocity, p.normal)/std::sqrt(v2) : Scalar(1);
    // Only roundoff in a normalized dot product is clamped, never a physical angle.
    p.cosine = p.cosine < Scalar(0) ? Scalar(0) :
               p.cosine > Scalar(1) ? Scalar(1) : p.cosine;
    // 表 I：两个背散射产额共享角度修正，真二次峰值及峰值能量各自修正。
    Scalar const factor = Scalar(1) + Scalar(0.26)*(Scalar(1)-p.cosine*p.cosine);
    p.elastic = (Scalar(0.07) + Scalar(0.43)*
        std::exp(-std::pow(p.energy/Scalar(100), Scalar(0.9))/Scalar(0.9))) * factor;
    p.rediffused = Scalar(0.74)*(-std::expm1(-p.energy/Scalar(40))) * factor;
    Scalar const peak_yield = Scalar(1.22)*
        (Scalar(1)+Scalar(0.66)*(Scalar(1)-std::pow(p.cosine, Scalar(0.8))));
    Scalar const peak_energy = Scalar(310)*(Scalar(1)+Scalar(0.7)*(Scalar(1)-p.cosine));
    Scalar const x = p.energy/peak_energy;
    p.true_yield = peak_yield*Scalar(1.813)*x/(Scalar(0.813)+std::pow(x, Scalar(1.813)));
    // 剩余分支同时包含吸收和真二次；A=0 仅在真二次产额也为零时可用。
    Scalar const a = Scalar(1)-p.elastic-p.rediffused;
    if (!std::isfinite(a) || !std::isfinite(p.true_yield) ||
        a < Scalar(0) || (a == Scalar(0) && p.true_yield > Scalar(0))) {
        p.status = Status::probabilities;
        return p;
    }
    p.mean = a > Scalar(0) ? p.true_yield/a : Scalar(0);
    if (!std::isfinite(p.mean) || p.mean < Scalar(0) || p.mean >= Scalar(m_trials)) {
        p.status = Status::binomial_mean;
    }
    return p;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
StainlessSteelInteraction::Event StainlessSteelInteraction::SelectEvent (
    Parameters const& p, amrex::RandomEngine const& engine) const noexcept
{
    // Outside the sampler's admissible domain, absorb the whole incident.
    // The common finalizer handles its charge deposition and deletion.
    if (p.status != Status::valid) { return {WallBehavior::absorb, 0}; }
    Scalar const u = UniformOpen(engine);
    if (u < p.elastic) { return {WallBehavior::elastic_backscatter, 1}; }
    if (u < p.elastic+p.rediffused) { return {WallBehavior::rediffused_backscatter, 1}; }
    int n = 0;
    Scalar const b = p.mean/Scalar(m_trials);
    for (int j = 0; j < m_trials; ++j) {
        n += UniformOpen(engine) < b;
    }
    return {n == 0 ? WallBehavior::absorb : WallBehavior::true_secondary, n};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
StainlessSteelInteraction::Scalar StainlessSteelInteraction::Gamma (
    Scalar const shape, amrex::RandomEngine const& engine) noexcept
{
    Scalar const d = shape-Scalar(1)/Scalar(3);
    Scalar const c = Scalar(1)/std::sqrt(Scalar(9)*d);
    for (;;) {
        Scalar const z = static_cast<Scalar>(
            amrex::RandomNormal(amrex::Real(0), amrex::Real(1), engine));
        Scalar const t = Scalar(1)+c*z;
        if (t <= Scalar(0)) { continue; }
        Scalar const v = t*t*t;
        if (std::log(UniformOpen(engine)) <
            Scalar(0.5)*z*z + d*(Scalar(1)-v+std::log(v))) {
            return d*v;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
StainlessSteelInteraction::Scalar StainlessSteelInteraction::TruncatedGamma (
    Scalar const shape, Scalar const scale, Scalar const limit,
    amrex::RandomEngine const& engine) noexcept
{
    if (limit == Scalar(0)) { return Scalar(0); }
    Scalar const x = limit/scale;
    // 上限不低于未截断分布的均值时，直接 Gamma 拒绝抽样效率足够。
    if (x >= shape) {
        Scalar s;
        do { s = Gamma(shape, engine); } while (s > x);
        Scalar const result = scale*s;
        return result > limit ? limit : result; // multiplication roundoff only
    }
    // t=S/limit has density t^(shape-1)*exp(-x*t), 0<t<1.
    // Propose Beta(k,1), k=max(shape-x,1); divide by the exact rejection envelope.
    Scalar const k = shape-x > Scalar(1) ? shape-x : Scalar(1);
    Scalar const mode = x > Scalar(0) && (shape-k) < x ? (shape-k)/x : Scalar(1);
    Scalar const log_envelope = (shape-k)*std::log(mode)-x*mode;
    for (;;) {
        Scalar const log_t = std::log(UniformOpen(engine))/k;
        Scalar const t = std::exp(log_t);
        Scalar const log_accept = mode == Scalar(1) && k > Scalar(1) ?
            x*(log_t-std::expm1(log_t)) : (shape-k)*log_t-x*t-log_envelope;
        if (std::log(UniformOpen(engine)) <= log_accept) { return limit*t; }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
StainlessSteelInteraction::Scalar StainlessSteelInteraction::ElasticEnergy (
    Scalar const energy, amrex::RandomEngine const& engine) noexcept
{
    if (energy == Scalar(0)) { return Scalar(0); }
    constexpr Scalar sigma = Scalar(1.9);
    // At low energy, uniform rejection avoids the vanishing acceptance rate of
    // rejecting an untruncated half-normal. Both branches sample Eq. (26).
    if (energy < sigma) {
        for (;;) {
            Scalar const loss = energy*UniformOpen(engine);
            if (std::log(UniformOpen(engine)) <= -Scalar(0.5)*(loss/sigma)*(loss/sigma)) {
                return energy-loss;
            }
        }
    }
    Scalar loss;
    do {
        loss = sigma*std::abs(static_cast<Scalar>(
            amrex::RandomNormal(amrex::Real(0), amrex::Real(1), engine)));
    } while (loss > energy);
    return energy-loss;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<ParticleVector, StainlessSteelInteraction::spectrum_count>
StainlessSteelInteraction::SampleEmission (
    Parameters const& p, Event const& event, amrex::RandomEngine const& engine) const noexcept
{
    amrex::GpuArray<Scalar, spectrum_count> energies{};
    if (event.kind == WallBehavior::elastic_backscatter) {
        energies[0] = ElasticEnergy(p.energy, engine);
    } else if (event.kind == WallBehavior::rediffused_backscatter) {
        energies[0] = p.energy*std::pow(UniformOpen(engine), Scalar(1)/Scalar(1.4));
    } else if (event.kind == WallBehavior::true_secondary) {
        // 表 II 按发射数量 n 选第 n-1 项；同一次事件中所有电子共用 p_n、epsilon_n。
        constexpr amrex::GpuArray<Scalar, spectrum_count> shapes{
            Scalar(1.6), Scalar(2), Scalar(1.8), Scalar(4.7), Scalar(1.8),
            Scalar(2.4), Scalar(1.8), Scalar(1.8), Scalar(2.3), Scalar(1.8)};
        constexpr amrex::GpuArray<Scalar, spectrum_count> scales{
            Scalar(3.9), Scalar(6.2), Scalar(13), Scalar(8.8), Scalar(6.25),
            Scalar(2.25), Scalar(9.2), Scalar(5.3), Scalar(17.8), Scalar(10)};
        Scalar const shape = shapes[event.count-1];
        Scalar remaining = TruncatedGamma(
            Scalar(event.count)*shape, scales[event.count-1], p.energy, engine);
        // Independent Gamma weights yield Dirichlet fractions, independent of S.
        // Truncating S alone gives exactly the document's conditional joint law.
        amrex::GpuArray<Scalar, spectrum_count> weights{};
        Scalar weight_sum = 0;
        for (int j = 0; j < event.count; ++j) {
            weights[j] = Gamma(shape, engine);
            weight_sum += weights[j];
        }
        // 逐项分配剩余总能量，最后一项取余量，减小独立乘法造成的总和误差。
        for (int j = 0; j+1 < event.count; ++j) {
            Scalar const fraction = weights[j]/weight_sum;
            energies[j] = remaining*(fraction > Scalar(1) ? Scalar(1) : fraction);
            remaining -= energies[j];
            weight_sum -= weights[j];
        }
        energies[event.count-1] = remaining;
    }

    ParticleVector tangent1, tangent2;
    Math::BuildOrthonormalBasis(p.normal, tangent1, tangent2);
    amrex::GpuArray<ParticleVector, spectrum_count> velocities{};
    constexpr Scalar energy_to_v2 = Scalar(2)*PhysConst::q_e_v<Scalar>/PhysConst::m_e_v<Scalar>;
    Scalar energy_sum = 0;
    // dP/dOmega=cos(theta)/pi：sin(theta)=sqrt(U)，方位角均匀。
    // 包括弹性背散射在内都使用该角分布，不采用镜面方向。
    for (int j = 0; j < event.count; ++j) {
        Scalar const u = UniformOpen(engine);
        Scalar const phi = Scalar(6.2831853071795864769)*UniformOpen(engine);
        Scalar const speed = std::sqrt(energy_to_v2*energies[j]);
        Scalar const tangent = speed*std::sqrt(u);
        velocities[j] = Math::ExpandInBasis(
            ParticleVector{tangent*std::cos(phi), tangent*std::sin(phi),
                           speed*std::sqrt(Scalar(1)-u)}, tangent1, tangent2, p.normal);
        energy_sum += Math::NormSquared(velocities[j])/energy_to_v2;
    }
    if (energy_sum > p.energy) {
        // Correct only floating-point overshoot after basis rotation/storage.
        Scalar const scale = std::sqrt(p.energy/energy_sum) *
            (Scalar(1)-Scalar(16)*std::numeric_limits<Scalar>::epsilon());
        for (int j = 0; j < event.count; ++j) {
            velocities[j].x *= scale;
            velocities[j].y *= scale;
            velocities[j].z *= scale;
        }
    }
    return velocities;
}

} // namespace Insert
