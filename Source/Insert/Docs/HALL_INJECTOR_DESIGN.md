# Hall Injector Design

## 背景和目标

`Source/Insert` 中的粒子注入逻辑只服务当前霍尔推力器模型，不再以合并回
WarpX 官方 `PlasmaInjector` 为设计约束。设计目标是把现有电子、离子、原子
注入逻辑整理成可组合、可扩展、便于维护的私有注入框架，并为后续真实空心
阴极点发射、孔阵列发射和 spoke 初始化分布预留接口。

当前实现主要分散在 `InsertInjection.cpp`：

- `PlasmaInit()`：初始化环形等离子体，同时生成电子和 Xe 离子。
- `CathodeInjection3D()`：按阴极电流持续注入电子。
- `XeInjection()` / `XeFastInjection()`：按质量流量从平面或孔阵列注入 Xe 原子。
- spoke 相关角向分布函数已经存在，但当前 `MakePlasmaThetaSampler()` 实际使用
  `uniform_distribution`。

本设计不要求兼容 WarpX 官方输入参数格式，也不要求复用 `PlasmaInjector` 的
`injection_style`。WarpX 官方 `PlasmaInjector` 可以作为参考，但不是继承对象。

## 设计原则

1. `Insert` 注入只表达霍尔推力器物理，不为通用 PIC 注入抽象付出额外复杂度。
2. 位置空间、速度空间和注入速率必须可组合。
3. 位置空间和速度空间分别选择坐标系；霍尔推力器默认位置采用柱坐标，速度采用
   笛卡尔坐标。
4. 每个坐标方向独立配置分布函数并采样，后续再扩展多坐标耦合分布。
5. spoke 是柱坐标角向分布的一种，不应硬编码在 `PlasmaInit()` 中。
6. 空心阴极点发射需要位置、法向和权重因子，不能只返回坐标。
7. 多物种同步注入应由源配置描述，避免多个函数各自硬编码 species 名称。
8. 先保留 CPU/host 侧 `AddNParticles()` 路径，后续确有性能需求再迁移到 tile
   级 GPU 注入。

## 总体结构

建议在 `Source/Insert/Injection` 中逐步整理为：

```text
Source/Insert/Injection/
  HallInjector.h
  HallInjector.cpp
  HallInjectionConfig.h
  HallInjectionConfig.cpp
  HallInjectionSource.h
  HallInjectionSource.cpp
  HallCoordinateDistribution.h
  HallCoordinateDistribution.cpp
  HallCoordinateTransform.h
  HallCoordinateTransform.cpp
  HallDistribution1D.h
  HallDistribution1D.cpp
  HallRateModel.h
  HallRateModel.cpp
  DistributionSampler1D.H
  DistributionSampler1D.cpp
```

其中 `NumericalInverseCDFSampler1D` 继续复用，用于 spoke、任意表格化方向分布和
其他无法方便写出解析反 CDF 的一维概率密度采样。`DistributionSampler1D` 仅作为
兼容旧代码的类型别名保留。

## 核心数据模型

### EmissionSample

位置分布不只返回坐标，而返回完整发射样本：

```cpp
struct EmissionSample
{
    amrex::ParticleReal x;
    amrex::ParticleReal y;
    amrex::ParticleReal z;
    amrex::ParticleReal nx;
    amrex::ParticleReal ny;
    amrex::ParticleReal nz;
    amrex::ParticleReal weight_factor = 1.0_prt;
};
```

字段含义：

- `x/y/z`：粒子初始位置。
- `nx/ny/nz`：局部发射法向。体积初始化可使用零向量或默认轴向。
- `weight_factor`：几何或发射点权重，例如不同孔面积、不同点源强度。

### HallInjectionSource

每个源负责一个 species 的一次注入配置：

```cpp
class HallInjectionSource
{
public:
    void Initialize();
    void Inject(WarpX& warpx, amrex::Real dt, int step);

private:
    std::string m_species_name;
    amrex::ParticleReal m_macro_weight;
    HallRateModel m_rate;
    HallCoordinateDistribution m_position_space;
    HallCoordinateDistribution m_velocity_space;
    amrex::Real m_fractional_particles = 0.0;
};
```

它负责：

1. 根据 `RateModel` 计算本步应注入的宏粒子数。
2. 调用位置空间分布采样坐标、法向和几何权重。
3. 调用速度空间分布采样速度坐标。
4. 调用对应 species 的 `AddNParticles()`。

### HallInjector

全局入口只做调度：

```cpp
class HallInjector
{
public:
    void ReadParameters();
    void InitializePlasma(WarpX& warpx);
    void InjectParticles(WarpX& warpx, amrex::Real dt, int step);

private:
    amrex::Vector<HallInjectionSource> m_initial_sources;
    amrex::Vector<HallInjectionSource> m_continuous_sources;
};
```

`Insert::Initialize()` 调用 `InitializePlasma()`；`Insert::ParticleInjection()` 调用
`InjectParticles()`。

## 注入速率模型

`HallRateModel` 负责把物理量转成宏粒子数，并维护小数残差。

建议支持：

```text
fixed_count
density_volume
current
mass_flow
```

### fixed_count

直接每次注入固定宏粒子数。适用于调试、基准测试和快速初始化。

### density_volume

用于初始环形等离子体：

```text
macro_particles = density * volume / macro_weight
```

spoke 分布默认只改变角向分布，不改变总粒子数。若后续需要把角向分布积分用于
缩放总数，可以显式增加 `scale_count_by_distribution_integral`。

### current

用于阴极电子：

```text
real_particles = current * dt / e
macro_particles = real_particles / macro_weight
```

当前代码中还有 `l_factor` 缩放，应作为私有参数保留：

```text
real_particles = current * dt / l_factor^2 / e
```

### mass_flow

用于 Xe 原子注入：

```text
real_particles = mass_flow * dt / atom_mass
macro_particles = real_particles / macro_weight
```

同样保留当前 `l_factor^2` 缩放。

## 坐标空间分布

`HallCoordinateDistribution` 是位置和速度共用的分布描述。每个空间先选择坐标系，
再为每个坐标方向选择一维分布函数。第一阶段只支持方向独立采样，即联合分布为：

```text
P(q0, q1, q2) = P0(q0) * P1(q1) * P2(q2)
```

后续需要孔阵列局部圆斑、倾斜喷口、速度角散度等多坐标相关分布时，再加入
`coupled_distribution`。

支持的坐标系：

```text
cartesian:   x y z      或 vx vy vz
cylindrical: r theta z  或 vr vtheta vz
```

霍尔推力器默认约定：

```text
position.coordinate_system = cylindrical
velocity.coordinate_system = cartesian
```

也就是说，粒子位置通常在柱坐标中描述环形通道、spoke 和轴向范围；速度直接在
WarpX `AddNParticles()` 需要的笛卡尔分量中采样。

### 采样流程

每个源的采样流程为：

1. 在位置坐标系中采样 `q0/q1/q2`。
2. 将位置坐标转换为笛卡尔 `x/y/z`。
3. 根据位置样本给出默认法向和几何权重。
4. 在速度坐标系中采样 `u0/u1/u2`。
5. 若速度坐标系不是笛卡尔，则根据位置角度和局部基矢转换为 `vx/vy/vz`。
6. 调用 `AddNParticles()`。

### 方向分布函数

每个方向支持以下一维分布：

```text
constant
uniform
area_uniform
gaussian
positive_gaussian
single_spoke
multi_spoke
discrete
tabulated
parser
```

`constant` 用于固定平面、固定点源或固定速度分量。

`uniform` 在 `[min, max]` 上均匀采样。

`area_uniform` 用于柱坐标半径方向，使横截面积均匀：

```cpp
r = sqrt((r_max*r_max - r_min*r_min) * random + r_min*r_min);
```

当前阴极注入中的等价形式可以保留：

```cpp
r = sqrt((r1 + r2) * (random * (r2-r1) + r1) - r1*r2);
```

`gaussian` 用于普通正态速度或坐标扰动。

`positive_gaussian` 用于法向正向发射。第一阶段用于笛卡尔 `vz`：

```cpp
do {
    vz = RandomNormal(mean, sigma);
} while (vz < 0);
```

如果未来速度空间选择局部坐标或柱坐标，可把正向判定推广到指定方向。

`single_spoke` 和 `multi_spoke` 是 `theta` 方向分布。当前
`single_spoke_distribution()` 中固定为：

```text
center = pi
sigma = pi / 4
```

迁移后应改为输入参数：

```text
spoke_center
spoke_sigma
spoke_count
spoke_phase
```

`multi_spoke` 保留当前思想：把每个 spoke 区间映射回单 spoke 分布。

`discrete` 用于从一组离散坐标值中按权重选择，例如多个点发射器或离散孔中心。

`tabulated` 使用表格概率密度构造 CDF；`NumericalInverseCDFSampler1D` 已经提供
基本采样能力。

`parser` 为后续解析一维函数预留，例如 `f(theta)` 或 `f(z)`。

### 一维分布类层次

一维方向分布统一用类实现，并继承共同基类。第一阶段注入走 host 侧
`AddNParticles()` 批量路径，使用虚函数多态即可；后续如果需要 tile 级 GPU 注入，
再把这一层改成 `std::variant` 或 WarpX 风格的 `enum + union`。

基类建议为：

```cpp
class HallDistribution1D
{
public:
    virtual ~HallDistribution1D() = default;

    [[nodiscard]] virtual amrex::ParticleReal
    sample(amrex::RandomEngine const& engine) const = 0;

    [[nodiscard]] virtual amrex::ParticleReal
    min() const noexcept = 0;

    [[nodiscard]] virtual amrex::ParticleReal
    max() const noexcept = 0;

    [[nodiscard]] virtual amrex::ParticleReal
    integral() const noexcept
    {
        return 1.0_prt;
    }
};
```

各分布对应派生类：

```text
HallConstantDistribution1D
HallUniformDistribution1D
HallAreaUniformDistribution1D
HallGaussianDistribution1D
HallPositiveGaussianDistribution1D
HallSingleSpokeDistribution1D
HallMultiSpokeDistribution1D
HallDiscreteDistribution1D
HallTabulatedDistribution1D
HallParserDistribution1D
```

`HallCoordinateDistribution` 持有三个方向分布：

```cpp
class HallCoordinateDistribution
{
private:
    CoordinateSystem m_coordinate_system;
    std::unique_ptr<HallDistribution1D> m_q0;
    std::unique_ptr<HallDistribution1D> m_q1;
    std::unique_ptr<HallDistribution1D> m_q2;
};
```

分布解析通过工厂函数集中处理：

```cpp
std::unique_ptr<HallDistribution1D>
MakeHallDistribution1D(amrex::ParmParse const& pp,
                       std::string const& prefix);
```

例如 `prefix = "channel_plasma.position.theta"` 时，工厂读取：

```text
channel_plasma.position.theta.distribution = multi_spoke
channel_plasma.position.theta.spoke_count = 4
channel_plasma.position.theta.spoke_sigma = pi / 4
```

`single_spoke`、`multi_spoke`、`tabulated` 和 `parser` 内部可以继续复用
`NumericalInverseCDFSampler1D` 进行数值 CDF 反演采样。`constant`、`uniform`、
`area_uniform`、`gaussian` 和 `positive_gaussian` 可直接实现采样，不需要预构造
CDF。

### 位置空间示例

环形体积初始化：

```text
position.coordinate_system = cylindrical
position.r.distribution = area_uniform
position.r.min = 0.0105 / l_factor
position.r.max = 0.0155 / l_factor
position.theta.distribution = multi_spoke
position.theta.spoke_count = 4
position.theta.spoke_sigma = pi / 4
position.z.distribution = uniform
position.z.min = 0.001 / l_factor
position.z.max = 0.004 / l_factor
```

环形平面注入：

```text
position.coordinate_system = cylindrical
position.r.distribution = area_uniform
position.r.min = r1
position.r.max = r2
position.theta.distribution = uniform
position.z.distribution = constant
position.z.value = 0
```

点发射：

```text
position.coordinate_system = cartesian
position.x.distribution = constant
position.x.value = cathode_x
position.y.distribution = constant
position.y.value = cathode_y
position.z.distribution = constant
position.z.value = cathode_z
```

多个点发射器可先用 `discrete` 选择发射点编号，再查表得到 `x/y/z/normal/weight`。
这属于轻量耦合；若需要在 `x/y/z` 上同时采样局部圆斑，应进入
`coupled_distribution`。

### 速度空间示例

电子和初始离子的各向同性高斯速度：

```text
velocity.coordinate_system = cartesian
velocity.vx.distribution = gaussian
velocity.vx.mean = 0
velocity.vx.sigma = sigma
velocity.vy.distribution = gaussian
velocity.vy.mean = 0
velocity.vy.sigma = sigma
velocity.vz.distribution = gaussian
velocity.vz.mean = 0
velocity.vz.sigma = sigma
```

Xe 原子平面正向注入：

```text
velocity.coordinate_system = cartesian
velocity.vx.distribution = gaussian
velocity.vx.mean = 0
velocity.vx.sigma = sqrt(kb * Tx / mass)
velocity.vy.distribution = gaussian
velocity.vy.mean = 0
velocity.vy.sigma = sqrt(kb * Ty / mass)
velocity.vz.distribution = positive_gaussian
velocity.vz.mean = vz0
velocity.vz.sigma = sqrt(kb * Tz / mass)
```

真实空心阴极若需要沿局部法向发射，建议新增速度局部坐标系：

```text
velocity.coordinate_system = local_normal
velocity.vnormal.distribution = positive_gaussian
velocity.vt1.distribution = gaussian
velocity.vt2.distribution = gaussian
```

`local_normal` 不作为第一阶段要求，但 `EmissionSample::normal` 已经为它预留了
数据。

### 多坐标耦合分布

第一阶段的独立方向分布可以覆盖：

- 环形体积初始化。
- 环形平面注入。
- spoke 角向初始化。
- 固定点源。
- 笛卡尔速度高斯、各向异性高斯和正向截断高斯。

以下功能需要耦合分布：

- 孔阵列平面注入：先选孔中心，再在孔局部圆斑内采样。
- 多点源带有限发射斑：先选点，再在点局部坐标中采样。
- 速度锥角发射：速度方向和法向、角度、方位角耦合。

建议预留接口：

```text
position.coupled_distribution = hole_array_plane
velocity.coupled_distribution = beam_cone
```

迁移期可以保留当前 `HoleArrayPlane` 的 named sampler，等独立方向分布稳定后再
统一到 `coupled_distribution`。

## 多物种初始化

当前 `PlasmaInit()` 对同一组位置先写电子，再用相同位置写 Xe 离子，只改变速度
sigma。这个功能应显式建模为共享位置样本。

建议提供：

```cpp
class HallCoupledInitialSource
{
public:
    void InitializePair(WarpX& warpx);

private:
    HallCoordinateDistribution m_position_space;
    amrex::Vector<SpeciesVelocityConfig> m_species;
};
```

配置上允许：

```text
shared_position = true
species = electrons xe_ions
```

这样可以保证初始准中性等离子体中电子和离子位置一一对应。

## 参数组织建议

建议先使用 `my_constants` 或 `insert.injection.*` 命名空间，不追求 WarpX 官方
参数兼容。

示例：

```text
insert.initial_sources = channel_plasma

channel_plasma.type = initial
channel_plasma.position.coordinate_system = cylindrical
channel_plasma.position.r.distribution = area_uniform
channel_plasma.position.r.min = 0.0105 / l_factor
channel_plasma.position.r.max = 0.0155 / l_factor
channel_plasma.position.theta.distribution = multi_spoke
channel_plasma.position.theta.spoke_count = 4
channel_plasma.position.theta.spoke_sigma = pi / 4
channel_plasma.position.z.distribution = uniform
channel_plasma.position.z.min = 0.001 / l_factor
channel_plasma.position.z.max = 0.004 / l_factor
channel_plasma.density = 1e18
channel_plasma.shared_species = electrons xe_ions
channel_plasma.electrons.weight = elec_weight
channel_plasma.electrons.velocity.coordinate_system = cartesian
channel_plasma.electrons.velocity.vx.distribution = gaussian
channel_plasma.electrons.velocity.vx.sigma = 592982
channel_plasma.electrons.velocity.vy.distribution = gaussian
channel_plasma.electrons.velocity.vy.sigma = 592982
channel_plasma.electrons.velocity.vz.distribution = gaussian
channel_plasma.electrons.velocity.vz.sigma = 592982
channel_plasma.xe_ions.weight = elec_weight
channel_plasma.xe_ions.velocity.coordinate_system = cartesian
channel_plasma.xe_ions.velocity.vx.distribution = gaussian
channel_plasma.xe_ions.velocity.vx.sigma = 1212.41
channel_plasma.xe_ions.velocity.vy.distribution = gaussian
channel_plasma.xe_ions.velocity.vy.sigma = 1212.41
channel_plasma.xe_ions.velocity.vz.distribution = gaussian
channel_plasma.xe_ions.velocity.vz.sigma = 1212.41
```

Xe 原子孔阵列注入示例：

```text
insert.continuous_sources = xe_neutral_inlet

xe_neutral_inlet.species = xe_netural
xe_neutral_inlet.rate = mass_flow
xe_neutral_inlet.mass_flow = m_dot
xe_neutral_inlet.macro_weight = xe_weight
xe_neutral_inlet.mass = 2.179e-25
xe_neutral_inlet.position.coordinate_system = cylindrical
xe_neutral_inlet.position.coupled_distribution = hole_array_plane
xe_neutral_inlet.hole_count = 48
xe_neutral_inlet.hole_radius = 0.001 / l_factor
xe_neutral_inlet.hole_ring_radius = (0.021 + 0.031) / 4 / l_factor
xe_neutral_inlet.z = 0
xe_neutral_inlet.velocity.coordinate_system = cartesian
xe_neutral_inlet.velocity.vx.distribution = gaussian
xe_neutral_inlet.velocity.vx.sigma = sqrt(kb * Tx / mass)
xe_neutral_inlet.velocity.vy.distribution = gaussian
xe_neutral_inlet.velocity.vy.sigma = sqrt(kb * Ty / mass)
xe_neutral_inlet.velocity.vz.distribution = positive_gaussian
xe_neutral_inlet.velocity.vz.mean = vz0
xe_neutral_inlet.velocity.vz.sigma = sqrt(kb * Tz / mass)
```

真实空心阴极点发射示例：

```text
cathode_electron.species = electrons
cathode_electron.rate = current
cathode_electron.current = Ic
cathode_electron.macro_weight = elec_weight
cathode_electron.position.coordinate_system = cartesian
cathode_electron.position.x.distribution = constant
cathode_electron.position.x.value = cathode_x
cathode_electron.position.y.distribution = constant
cathode_electron.position.y.value = cathode_y
cathode_electron.position.z.distribution = constant
cathode_electron.position.z.value = cathode_z
cathode_electron.normal = nx ny nz
cathode_electron.spot_radius = 0.0001
cathode_electron.velocity.coordinate_system = local_normal
cathode_electron.velocity.coupled_distribution = beam_cone
cathode_electron.mean_speed = cathode_v0
cathode_electron.angle_sigma = cathode_angle_sigma
```

## 迁移计划

### 阶段 1：无行为改变的工具提取

1. 提取 `HallCoordinateDistribution` 和一维方向分布函数，保持当前
   `PlasmaInit()` 和 `CathodeInjection3D()` 的采样结果一致。
2. 提取柱坐标到笛卡尔坐标的转换，位置空间默认使用 `cylindrical`，速度空间默认
   使用 `cartesian`。
3. 用方向分布表达各向同性高斯、各向异性高斯和 `vz` 正向截断高斯。
4. `HoleArrayPlane` 暂时作为 named coupled sampler 保留，保持当前
   `XeInjection()` 行为一致。
5. 保留旧函数入口，只把内部实现替换为新类。

验收标准：

- 当前 `HALL3D` 和 `HALL3D_INIT` 宏路径仍可编译。
- 相同随机种子下统计分布一致，粒子数量一致。

### 阶段 2：统一源配置

1. 新增 `HallInjectionSource` 和 `HallInjector`。
2. `Insert::Initialize()` 通过 `HallInjector` 执行初始源。
3. `Insert::ParticleInjection()` 通过 `HallInjector` 执行连续源。
4. 旧函数变为兼容包装或删除。

验收标准：

- `PlasmaInit()`、`CathodeInjection3D()`、`XeInjection()` 的物理配置能用统一源
  表达。
- 源数量、species 名称和权重不再硬编码在多个函数中。

### 阶段 3：启用 spoke 参数化

1. 把 `single_spoke_distribution` 和 `multi_spoke_distribution` 移入
   `theta` 方向分布。
2. `MakePlasmaThetaSampler()` 改为根据 `position.theta.distribution` 选择分布。
3. 保留 `uniform` 作为默认值，避免默认行为改变。

验收标准：

- 默认配置仍为均匀角向分布。
- 设置 `position.theta.distribution = multi_spoke` 后，初始化粒子角向分布呈指定
  spoke 数。

### 阶段 4：真实空心阴极扩展

1. 新增 `cartesian + constant` 点源配置和 `discrete` 多点源配置。
2. 新增 `local_normal` 速度坐标系和 `beam_cone` 耦合速度采样。
3. 支持多个发射点及发射权重。

验收标准：

- 单点源和多点源能复用同一套 rate/velocity 组合。
- 点源法向变化后速度方向随之变化。

## 与 PlasmaInjector 的边界

WarpX 官方 `PlasmaInjector` 仍负责官方通用注入路径。`Insert` 的 Hall 注入器只在
`Insert::Initialize()` 和 `Insert::ParticleInjection()` 中调用，不新增
`PlasmaInjector` 派生类，也不修改官方 `injection_style`。

可借鉴但不直接依赖的设计：

- `InjectorPosition` 的“运行时选择具体位置分布”思想。
- `InjectorMomentum` 的“速度分布独立于位置分布”思想。
- `NumericalInverseCDFSampler1D` 形式的数值 CDF 反演采样。
- `PlasmaInjector` 的参数解析分层思想。

不采用的部分：

- 不使用 `PlasmaInjector` 的 `enum + union` GPU 可拷贝限制作为第一阶段约束；
  一维分布优先使用共同基类和虚函数实现。
- 不追求输入文件参数与 WarpX 官方注入参数一致。
- 不把 Hall 专用概念塞入 WarpX 通用 `Source/Initialization`。

## 风险和注意事项

- 当前 `MakePlasmaThetaSampler()` 虽然有 spoke 函数，但默认使用均匀分布；迁移时
  必须保持默认行为不变。
- 当前 Xe species 名称为 `xe_netural`，疑似拼写历史遗留。迁移时不要擅自改名，
  除非同步修改输入文件和下游逻辑。
- 当前使用 `std::rand()` 选择孔阵列起点，应迁移为 AMReX 随机数，避免并行可重复性
  问题。
- 如果未来需要在 GPU 上大规模注入，应把 `AddNParticles()` 批量 host 路径改为
  tile 级直接写入；这不作为第一阶段目标。
- 小数粒子残差必须作为源对象成员保存，不能继续依赖函数内 `static` 状态。
