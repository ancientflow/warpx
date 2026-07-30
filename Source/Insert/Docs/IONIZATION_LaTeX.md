# 自洽电离源统计与双线性采样设计

## 1. 目标和边界

当前用平均中性密度替代瞬时中性粒子反馈的方案仍然保留了电子碰撞电离的瞬时反馈链：
电子密度、电子能量和局部电场的低频振荡仍会调制电离率，因此不能稳定地抑制呼吸振荡。

本方案改为两阶段闭环：

```text
完整自洽放电: CoupledBackgroundMCCCollision 产生真实电离事件
        -> 统计加权电离源 S_bar(r,z)
        -> 构造双线性采样表和采样器
        -> 合并注入重构分支后，由其接口按采样结果注入电子-离子对
```

设计目标：

1. 在正常放电阶段统计自洽 MCC 电离源，而不是统计平均中性密度。
2. 在抑制呼吸振荡的控制算例中，电离对数量和空间分布只由统计表决定。
3. 当前分支只实现到电离源统计、采样表构造和双线性取点采样；实际粒子注入等待从 `my-dev` 出发的注入重构分支合并。
4. 继续使用三维粒子、三维网格、三维边界和现有静电/电磁求解路径。
5. 初期仅支持当前 `Source/Insert` 私有三维 Hall 算例路径，即 `WARPX_DIM_3D`、`HALL3D`、`MCC_DENSITY`。
6. 第一版直接使用网格内双线性源项采样，不实现分片常数源。
7. 第一版不考虑 MPI 和自适应网格，只支持单 MPI rank、单 level。
8. 所有功能控制都由编译宏完成；相关宏关闭时，代码行为必须与修改前一致。

`MCC_DENSITY_AVERAGE_CALC` 和 `MCC_DENSITY_AVERAGE_USE` 对本目标已经没有保留意义：它们应被平均电离源统计、双线性采样以及未来注入路径完全替代。新的控制算例不得再依赖平均中性密度文件，也不得同时启用平均中性密度和平均电离源路径。

## 2. 当前 WarpX 代码事实

相关代码路径：

```text
Source/Insert/Core/WarpXInsert.cpp
Source/Insert/Background/BackgroundCoupledDensity.cpp
Source/Insert/Collisions/CoupledBackgroundMCCCollision.cpp
Source/Insert/Collisions/FilterCopyTransformCoupled.H
Source/Insert/Injection/InsertInjection.cpp
Source/Insert/Diagnostics/InsertRuntimeDiagnostics.cpp
Source/Insert/Config/WarpXFunctionConfig.h
```

当前私有 MCC 路径的关键行为：

1. `Insert::BeforeCollision(step)` 更新 `global_background_density`。
2. `CoupledBackgroundMCCCollision::doCollisions()` 从 `BackgroundCoupledDensity::m_background_density_fabs` 读取局部中性密度。
3. `doBackgroundIonizationCouple<depos_order>()` 使用 `CoupledImpactIonizationFilterFunc` 判定电离事件。
4. `filterCopyTransformParticles<1, depos_order>()` 对满足条件的入射电子生成一个新电子和一个离子。
5. `ImpactIonizationTransformFunc` 修改入射电子能量，并设置新电子和新离子的速度。
6. 完整自洽统计算例中，还会按 `elec_weight` 消耗中性粒子权重并更新背景密度。
7. 旧的 `MCC_DENSITY_AVERAGE_CALC/USE` 只平均背景中性密度，不能切断瞬时电离反馈；新实现应停止使用这两个宏作为控制算例路径。

电离源统计和采样属于碰撞模块，应放在 `Source/Insert/Collisions`。当前分支只复用：

1. `doBackgroundIonizationCouple<depos_order>()` 中已经生成的电离事件和新增电子权重。
2. `my_constants.elec_weight` 作为统计权重一致性的检查依据。
3. `my_constants.l_factor` 和当前 Hall3D 坐标约定。
4. AMReX `MultiFab` 和 plotfile 输出接口。

当前 Hall3D 自定义注入代码采用的坐标约定是：

```text
z       : 推进轴向
x, y    : 横截面方向
center  : (L/2/l_factor, L/2/l_factor)
r       : sqrt((x-center_x)^2 + (y-center_y)^2)
theta   : atan2(y-center_y, x-center_x)
```

因此文档中的轴向变量记为 `z_axis`，三维 Cartesian 变量仍为 `(x,y,z)`。

## 3. 物理量定义

一次电离事件记为：

$$
(t_e, r_e, z_e, W_e, \mathbf{u}_{e,b}, \mathbf{u}_{i,b}) .
$$

其中：

1. `r_e` 是相对通道中心的径向位置。
2. `z_e` 是轴向位置，即代码中的 `PIdx::z`。
3. `W_e` 是该事件代表的真实电离对数。当前控制算例要求完整统计算例和未来平均源项注入算例严格使用相同的电子、离子宏粒子权重，因此第一版固定取 `W_e = my_constants.elec_weight`。
4. `u_{e,b}` 是新生电子速度分布，可选统计。
5. `u_{i,b}` 是新生离子速度分布，可选统计。

平均体源项定义为：

$$
\bar S(r,z) \quad [\mathrm{m^{-3}s^{-1}}] .
$$

在轴对称统计表中，单元权重不是单纯的 `S_bar`，而是体积分后的真实电离率：

$$
A_{i,k}
= \int_{z_i}^{z_{i+1}} \int_{r_k}^{r_{k+1}}
  \bar S(r,z) 2\pi r\,dr\,dz
\quad [\mathrm{s^{-1}}] .
$$

第一版不使用分片常数近似，而是在 `r-z` 源项网格节点上保存 `\bar S^{node}`，并在每个源项单元内做双线性插值：

$$
S(\xi,\eta)=a+b\xi+c\eta+d\xi\eta,
\qquad
\xi=\frac{z-z_i}{\Delta z},
\quad
\eta=\frac{r-r_k}{\Delta r}.
$$

单元权重必须由双线性源项积分得到：

$$
A_{i,k}=2\pi\Delta z_i\Delta r_k
\int_0^1\int_0^1 S(\xi,\eta)(r_k+\eta\Delta r_k)\,d\eta\,d\xi .
$$

总电离率为：

$$
A_{\mathrm{tot}}=\sum_{i,k}A_{i,k}.
$$

## 4. 正常放电统计阶段

### 4.1 统计窗口

当前算例已经通过 `PlasmaInit()` 初始化等离子体，本身避免了点火过程，因此第一版不再设置统计起止步。只要定义 `IONIZATION_SOURCE_RECORD`，就从第一个碰撞步开始全量统计，直到 `Insert::Finalize()` 输出结果。

这意味着统计时长为完整运行时间内实际发生 MCC 电离统计的时间：

$$
T_{\mathrm{avg}} = N_{\mathrm{collision\ sample}}\Delta t .
$$

如果后续需要跳过非稳态阶段，应重新论证物理必要性；当前设计不提供 start/stop 运行时参数。

### 4.2 事件位置

当前 `doBackgroundIonizationCouple()` 的电离事件由入射电子触发，`SmartCopy` 会把产物位置复制到入射电子位置。因此第一版统计事件位置应取入射电子或新增产物电子的位置，两者一致。

对当前 Hall3D 约定：

$$
r_e = \sqrt{(x_e-x_c)^2+(y_e-y_c)^2},
\qquad
z_e = z_e .
$$

中心位置固定为计算域横截面中心，不做运行时配置：

$$
x_c = \frac{x_{\min}+x_{\max}}{2},
\qquad
y_c = \frac{y_{\min}+y_{\max}}{2}.
$$

径向范围和轴向范围也由当前 `Geometry` 推导。对当前 Hall3D 轴向为 `z` 的约定：

$$
r_{\min}=0,
\qquad
r_{\max}=\min(x_c-x_{\min}, x_{\max}-x_c, y_c-y_{\min}, y_{\max}-y_c),
$$

$$
z_{\min}=\mathrm{ProbLo}(2),
\qquad
z_{\max}=\mathrm{ProbHi}(2).
$$

其他中心设置原则上不合理，第一版不提供。`nr` 和 `nz` 也由编译期常量或源项表类内部常量控制。

### 4.3 事件沉积

第一版使用节点型 `r-z` 源项网格。事件落入单元 `(i,k)` 后，按双线性权重沉积到四个节点：

$$
C^{node}_{i,k} \mathrel{+}= W_e(1-\xi)(1-\eta),
\qquad
C^{node}_{i+1,k} \mathrel{+}= W_e\xi(1-\eta),
$$

$$
C^{node}_{i,k+1} \mathrel{+}= W_e(1-\xi)\eta,
\qquad
C^{node}_{i+1,k+1} \mathrel{+}= W_e\xi\eta .
$$

统计时长为：

$$
T_{\mathrm{avg}}=(N_{\mathrm{sample}}\Delta t).
$$

节点源项由节点累计量除以统计时间和节点控制体积得到：

$$
\bar S^{node}_{i,k}=\frac{C^{node}_{i,k}}{T_{\mathrm{avg}}V^{node}_{i,k}}.
$$

随后每个单元用四个节点源项构造双线性函数 `S(\xi,\eta)`，并计算单元权重 `A_{i,k}`。由于两次计算的电子、离子宏粒子权重严格一致，事件权重固定为 `W_e = elec_weight`。

### 4.4 推荐实现接入点

推荐在 `doBackgroundIonizationCouple()` 中，在 `filterCopyTransformParticles<1, depos_order>()` 返回 `num_added` 后，对新增粒子区间统计：

```cpp
const auto np_elec = elec_tile.numParticles();
const auto np_ion = ion_tile.numParticles();
const auto num_added = filterCopyTransformParticles(...);

// 新增电子范围: [np_elec, np_elec + num_added)
// 新增离子范围: [np_ion, np_ion + num_added)
// 对新增电子位置和权重做 r-z 沉积。
```

该方式的优点是：

1. 不需要改 `FilterCopyTransformCoupled.H` 的 mask 内核接口。
2. 统计对象就是实际产生的粒子，和 `num_added`、`setNewParticleIDs()` 一致。
3. 可以同时统计新生电子速度分布，用于后续重放。

注意：新增粒子统计必须发生在 `setNewParticleIDs()` 前后均可，但不要依赖粒子 ID；应依赖新增区间 `[old_np, old_np + num_added)`。

### 4.5 单进程和单层级限制

第一版不考虑 MPI 并行统计，也不考虑自适应网格。定义 `IONIZATION_SOURCE_RECORD` 或 `IONIZATION_SOURCE_INJECT` 时，应在初始化阶段检查：

```cpp
WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
    amrex::ParallelDescriptor::NProcs() == 1,
    "Average ionization source only supports one MPI rank in the first version.");
```

同时只支持单 level：

```cpp
WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
    warpx_instance.maxLevel() == 0,
    "Average ionization source only supports one mesh level in the first version.");
```

这与当前 `CoupledBackgroundMCCCollision` 已经只支持 `flvl == 0` 保持一致。

## 5. 源项文件格式

必须新增独立目录，不复用 `background_density_fab`，也不再输出或读取平均中性密度作为控制算例输入：

```text
ionization_source_fab/
  metadata.txt
  node_source_rate
  node_source_count
  cell_rate
  cell_cdf
ionization_source_plt/
  Header
  Level_0/
```

`node_source_rate` 保存节点源项 `\bar S^{node}_{i,k}`，单位为 `m^-3 s^-1`；`node_source_count` 保存统计窗口内沉积到节点的累计真实电离对数，用于检查统计噪声。`cell_rate` 保存由双线性节点源项积分得到的单元总率 `A_{i,k}`，单位为 `s^-1`；`cell_cdf` 保存采样器使用的一维累计分布。

运行完成后必须额外输出 ParaView 可直接打开的 AMReX plotfile 目录 `ionization_source_plt/`。该 plotfile 用于检查和展示统计得到的二维轴对称电离源网格，不作为控制算例的高性能读取格式。

plotfile 输出要求：

1. 使用 AMReX 提供的 plotfile 写出接口，不手写 VTK/XML 文件。ParaView/VisIt 读取 AMReX plotfile 时会把该数据视为三维数据，但 `y` 方向只有一个 cell。
2. 推荐包含 `<AMReX_PlotFileUtil.H>` 并调用 `amrex::WriteMultiLevelPlotfile`；第一版 `nlev = 1`。
3. 在 `WARPX_DIM_3D` 构建中，将二维 `r-z` 源项网格构造成带单层 dummy 维度的 cell-centered `amrex::MultiFab`，`BoxArray` 尺寸为 `nr x 1 x nz`。
4. `Geometry` 使用源项网格物理范围：`r_min/r_max`、dummy `y` 方向 `[0, 1]`、`z_min/z_max`。
5. plotfile 变量至少包含 `source_rate`、`source_density`、`source_count`、`cell_volume`。
6. plotfile 中的 `source_rate` 对应 `cell_rate`，`source_count` 对应单元累计计数，均为 cell-centered 可视化变量。
7. `source_density = source_rate / cell_volume`，单位为 `m^-3 s^-1`。
8. 输出目录名固定为 `ionization_source_plt`，不提供运行时改名参数。

`metadata.txt` 至少包含：

```text
version = 1
axis = z
center = <x_c> <y_c>
r_min = <...>
r_max = <...>
z_min = <...>
z_max = <...>
nr = <...>
nz = <...>
T_avg = <...>
A_tot = <...>
elec_weight = <my_constants.elec_weight used during recording and required during injection>
coordinate_units = simulation_units
node_source_units = real_pairs_per_m3_per_second
cell_rate_units = real_pairs_per_second_per_rz_cell
plotfile = ionization_source_plt
```

如果后续需要统计速度分布，可额外保存：

```text
electron_birth_energy_hist
ion_birth_velocity_moments
```

第一版不建议把速度分布作为必须项，否则会拖慢验证闭环。

## 6. 采样器与未来注入阶段

### 6.1 当前分支实现范围

当前分支只实现电离源统计、源项表输出和双线性位置采样器，不创建粒子，不接入 `AddNParticles()`。

当前分支使用的宏：

```cpp
#define IONIZATION_SOURCE_RECORD
```

定义 `IONIZATION_SOURCE_RECORD` 时，完整自洽放电会统计电离源，输出 `ionization_source_fab/` 和 `ionization_source_plt/`，并可构造 `IonizationSourceSampler` 做采样验证。

两个宏都未定义时，不执行任何平均源项统计、读取或注入代码，功能必须与修改前完全一致。

### 6.2 合并注入重构分支后的接口约束

从 `my-dev` 出发的另一个分支已经对粒子注入进行了重构。实际电子-离子对创建应等待该分支合并后，使用它提供的注入接口。

合并后可恢复如下宏语义：

```cpp
#define IONIZATION_SOURCE_INJECT
```

定义 `IONIZATION_SOURCE_INJECT` 时，代码读取固定目录 `ionization_source_fab/` 中的源项表，使用 `IonizationSourceSampler` 给出新生电子-离子对位置，并通过注入重构分支提供的接口创建粒子。同时禁用 MCC 电离产物生成。

`IONIZATION_SOURCE_RECORD` 和未来的 `IONIZATION_SOURCE_INJECT` 必须互斥：

```cpp
#if defined(IONIZATION_SOURCE_RECORD) && defined(IONIZATION_SOURCE_INJECT)
#   error "IONIZATION_SOURCE_RECORD and IONIZATION_SOURCE_INJECT are mutually exclusive"
#endif
```

本分支不得新增 `IonizationSourceInjection` 文件，也不得直接调用旧的 `AddNParticles()` 注入路径。

### 6.3 未来注入数量逻辑

读入源项表后：

$$
A_{\mathrm{tot}}=\sum_{i,k}A_{i,k}.
$$

一个时间步的期望真实电离对数：

$$
\Lambda=A_{\mathrm{tot}}\Delta t.
$$

固定宏粒子权重 `w_p = my_constants.elec_weight` 时，期望宏粒子对数：

$$
\lambda=\frac{A_{\mathrm{tot}}\Delta t}{w_p}.
$$

为了减少源项噪声，使用 deterministic reservoir：

$$
R\leftarrow R+\lambda,
\qquad
N_{\mathrm{pair}}=\lfloor R\rfloor,
\qquad
R\leftarrow R-N_{\mathrm{pair}}.
$$

长期平均满足：

$$
\left\langle N_{\mathrm{pair}} w_p\right\rangle
=A_{\mathrm{tot}}\Delta t.
$$

第一版固定使用 `reservoir`，不提供 Poisson 计数模型开关。控制源项的目的就是降低电离源自身噪声，因此不应重新引入泊松涨落。

## 7. 空间采样

先用双线性源项积分得到每个单元的真实电离率 `A_{i,k}`，再构造归一化概率、CDF 和采样用 alias table：

$$
P_{i,k}=\frac{A_{i,k}}{A_{\mathrm{tot}}},
\qquad
F_m=\sum_{q=0}^{m}P_q .
$$

`cell_cdf` 保留为可读输出和调试路径。性能路径应在读入源项表后用 `P_{i,k}` 构造 Vose alias table，使单元选择为 `O(1)`。为保持四随机数接口不变，alias 采样可使用 `j=floor(u0*N_active)` 选列，并用 `v=fract(u0*N_active)` 做 alias 判别。

二维单元展平顺序固定为：

$$
m = i_z N_r + i_r .
$$

对每个采样点使用四个随机数：

$$
u_0,u_1,u_2,u_3\in[0,1).
$$

含义：

1. `u0` 通过预构造的 alias table 确定 `r-z` 源项单元；`cell_cdf` 只作为可读输出和调试路径保留。
2. `u1` 通过径向边缘分布反解 `\eta`。
3. `u2` 通过给定 `\eta` 后的轴向条件分布反解 `\xi`。
4. `u3` 确定方位角。

选中单元后，使用四个节点源项构造：

$$
S(\xi,\eta)=a+b\xi+c\eta+d\xi\eta,
\qquad
r(\eta)=r_k+\eta\Delta r.
$$

单元内概率密度为：

$$
p(\xi,\eta)\propto S(\xi,\eta)r(\eta).
$$

性能优先的采样顺序是先径向、后轴向。径向边缘分布为：

$$
p(\eta)\propto \left(a+\frac{b}{2}+\left(c+\frac{d}{2}\right)\eta\right)
\left(r_k+\eta\Delta r\right).
$$

该边缘分布是二次多项式，CDF 是单调三次多项式。`u1` 反解得到 `\eta` 后，轴向条件分布为：

$$
p(\xi|\eta)\propto a+c\eta+(b+d\eta)\xi .
$$

`u2` 也使用固定 8 步二分反解。这样轴向和径向使用同一类固定迭代热路径，不需要解析退化分支。

最终：

$$
z=z_i+\xi\Delta z_i,
\qquad
r=r_k+\eta\Delta r_k,
\qquad
\theta=2\pi u_3.
$$

三维坐标：

$$
x=x_c+r\cos\theta,
\qquad
y=y_c+r\sin\theta,
\qquad
z=z .
$$

第一版不实现分片常数源项，也不提供体积均匀采样 fallback。

## 8. 未来注入阶段的粒子属性约束

当前分支不创建粒子。本节只记录合并注入重构分支后应保持的平均源项注入逻辑，并作为采样器返回权重信息的设计约束。

### 8.1 权重

未来注入阶段固定：

$$
w_e=w_i=w_p=\texttt{my_constants.elec_weight}.
$$

这和当前 `doBackgroundIonizationCouple()` 消耗中性粒子权重时使用 `elec_weight` 的路径一致。

因此源项表只需要保存统计时使用的 `elec_weight`，注入阶段读取后校验其与当前输入文件中的 `my_constants.elec_weight` 一致；若不一致应直接报错，而不是自动重标定源项。

### 8.2 位置

未来注入接口应保证电子和离子同位生成：

$$
\mathbf{x}_e=\mathbf{x}_i=\mathbf{x}_{\mathrm{birth}}.
$$

这样保持局部准中性，避免人为注入强空间电荷噪声。

### 8.3 新生离子速度

当前 MCC 电离中，`ImpactIonizationTransformFunc` 先用背景中性温度计算 `ion_vel_std`，再给新生离子三个速度分量分别抽取 Maxwellian 随机数。未来注入阶段应复用当前 Hall3D 中 Xe 离子热速度尺度，例如 `PlasmaInit()` 中的 `sigma_xe` 思路：

$$
\mathbf{u}_i \sim \mathcal{N}(0,\sigma_{i,x})
\times \mathcal{N}(0,\sigma_{i,y})
\times \mathcal{N}(0,\sigma_{i,z}).
$$

未来注入阶段不提供运行时速度参数，直接复用当前 `PlasmaInit()` 中 Xe 离子热速度常量。若后续需要改变该模型，应通过代码常量或新的编译宏显式修改。

### 8.4 新生电子速度

当前 `ImpactIonizationTransformFunc` 的能量分配是：先计算入射电子碰撞前动能 `E_coll`，扣除电离能 `m_energy_cost`，再把剩余能量二等分给散射后的入射电子和新生电子：

$$
E_{e,1}^{out}=E_{e,2}^{out}=\frac{E_{coll}-E_{ion}}{2}.
$$

两个电子随后使用相同动量幅值做各向同性散射。未来平均源项注入若完全固定源项，则入射电子能量损失消失，这是一个重要物理差异。

未来注入阶段固定使用低温 Maxwellian 新生电子模型，例如 `T_e,birth = 1 eV`，该值作为代码常量实现，不提供运行时参数。

并在验证中重点比较平均电子温度。如果控制算例电子温度偏高，再通过新的编译宏升级为统计正常放电中新生电子能量直方图，而不是在同一套运行时参数中混合多种模型。

## 9. 几何和边界合法性

当前分支构造采样表前应通过源项表预先排除非法区域，而不是在采样时大量拒绝。未来注入阶段应继承同一合法区域。

构造 `A_{i,k}` 时应置零的区域：

1. 计算域外。
2. 固体壁面或嵌入边界内部。
3. 非放电区、阳极/阴极实体区域。
4. `r < r_min` 或 `r > r_max`。
5. `z < z_min` 或 `z > z_max`。

采样器仍应保留安全检查。如果采样点非法，允许有限次数重采样；超过次数后返回失败状态并在诊断中计数。

## 10. 代码结构建议

新增文件：

```text
Source/Insert/Collisions/IonizationSourceTable.h
Source/Insert/Collisions/IonizationSourceTable.cpp
Source/Insert/Collisions/IonizationSourceRecorder.h
Source/Insert/Collisions/IonizationSourceRecorder.cpp
Source/Insert/Collisions/IonizationSourceSampler.h
Source/Insert/Collisions/IonizationSourceSampler.cpp
Source/Insert/Collisions/IonizationSourceOutput.h
Source/Insert/Collisions/IonizationSourceOutput.cpp
```

修改文件：

```text
Source/Insert/CMakeLists.txt
Source/Insert/Config/WarpXFunctionConfig.h
Source/Insert/Core/WarpXInsert.h
Source/Insert/Core/WarpXInsert.cpp
Source/Insert/Collisions/CoupledBackgroundMCCCollision.cpp
```

所有功能控制都使用编译宏完成：

```cpp
// #define IONIZATION_SOURCE_RECORD
// #define IONIZATION_SOURCE_INJECT  // reserved for the injection-refactor branch
```

宏语义：

1. 两个宏都未定义：代码必须与修改前功能一致，不统计、不读取、不注入平均电离源。
2. `IONIZATION_SOURCE_RECORD`：完整自洽放电，统计电离源，构建采样表，并在结束时输出源项表和 plotfile。
3. `IONIZATION_SOURCE_INJECT`：当前分支不实现；等待注入重构分支合并后，用其接口读取固定源项表并注入平均电离源。
4. 两个宏互斥，不能同时定义。

实现任务应删除或永久停用 `MCC_DENSITY_AVERAGE_CALC` 和 `MCC_DENSITY_AVERAGE_USE` 控制路径。

## 11. 推荐最小实现顺序

### 阶段 A：统计源项并构造采样器

1. 新增 `IonizationSourceTable`，维护 `r-z` 节点源项网格、累计数组、单元权重、CDF 和 alias table 输入概率。
2. 新增 `IonizationSourceSampler`，实现双线性源项的 `samplePosition(u0, u1, u2, u3)`。
3. 新增 `IonizationSourceOutput`，负责写 `metadata.txt`、内部读取用 `node_source_rate/node_source_count/cell_rate/cell_cdf` 和 ParaView 用 AMReX plotfile `ionization_source_plt/`。
4. 在 `doBackgroundIonizationCouple()` 的新增电子区间沉积 `W_e` 到节点源项。
5. 在 `Insert::Finalize()` 输出 `node_source_rate`、`node_source_count`、`cell_rate`、`cell_cdf`、`metadata.txt` 和 AMReX plotfile `ionization_source_plt/`。
6. 增加诊断输出 `A_tot(t)` 或每步 `num_added * elec_weight / dt`。

验收标准：统计得到的：

$$
A_{\mathrm{tot}}^{\mathrm{table}}
=\sum_{i,k}A_{i,k}
$$

与 `COLLISION_RECORD` 或新增诊断中的平均电离率一致，并且运行结束后 AMReX plotfile `ionization_source_plt/` 能被 ParaView 直接打开。

### 阶段 B：只实现按表采样器

1. 新增 `IonizationSourceSampler`，读入 `node_source_rate` 或已生成的 `cell_rate/cell_cdf`，构建双线性单元 alias table 和必要的反解系数。
2. 实现 `samplePosition(u0, u1, u2, u3)`，返回双线性源项采样得到的 `(x,y,z)`。
3. 采样器只返回位置和源项权重信息，不创建粒子，不调用 `AddNParticles()`。
4. 通过离线或测试调用采样器，将大量采样点重新统计回同一 `r-z` 网格。

验收标准：采样点恢复双线性源项表的单元权重：

$$
\frac{C_{i,k}^{\mathrm{sample}}}{N_{\mathrm{sample}}}
\approx
\frac{A_{i,k}}{\sum_{j,l}A_{j,l}}.
$$

实际电子-离子对注入等待从 `my-dev` 出发的注入重构分支合并后，通过该分支提供的注入接口接入；本计划保持源项数量、位置采样和权重逻辑不变。

### 阶段 C：合并注入重构分支后的控制算例接入

该阶段不属于当前分支实现范围，只保留接口计划：

1. 删除或停用 `MCC_DENSITY_AVERAGE_CALC` 和 `MCC_DENSITY_AVERAGE_USE` 相关控制路径。
2. 输入文件移除 MCC 电离过程或禁用电离产物生成，避免和平均源项注入双计数。
3. 保留其他非电离碰撞、边界、阴极电子注入和场求解。
4. 通过注入重构分支提供的接口按 `IonizationSourceSampler` 返回的位置、权重和未来速度模型创建电子-离子对。
5. 运行同样时间长度，对比低频谱和平均放电量。

未来验收标准：控制算例不再读取 `background_density_fab` 或任何平均中性密度文件；低频呼吸振荡明显降低，而 ECDI 高频结构仍存在。

## 12. 双线性源项实现细节

第一版直接实现双线性源项采样，不保留分片常数源项路径。

单元内双线性源项：

$$
S(\xi,\eta)=a+b\xi+c\eta+d\xi\eta,
\qquad
\xi=\frac{z-z_i}{\Delta z},
\quad
\eta=\frac{r-r_k}{\Delta r}.
$$

单元权重：

$$
A_{i,k}=2\pi\Delta z\Delta r
\int_0^1\int_0^1 S(\xi,\eta)(r_k+\eta\Delta r)
\,d\eta\,d\xi .
$$

为避免采样时重复计算，每个非零单元预计算并缓存：

$$
\alpha=a+\frac{b}{2},
\qquad
\beta=c+\frac{d}{2},
$$

$$
q_0=\alpha r_k,
\qquad
q_1=\alpha\Delta r+\beta r_k,
\qquad
q_2=\beta\Delta r.
$$

径向边缘 CDF 写成：

$$
F_\eta(\eta)=
\frac{q_0\eta+\frac{q_1}{2}\eta^2+\frac{q_2}{3}\eta^3}
{q_0+\frac{q_1}{2}+\frac{q_2}{3}}.
$$

反解 `F_\eta(\eta)=u_1` 时，性能路径使用固定 8 步二分，不做 Newton、不做收敛判断、不写 Cardano 三次闭式根。定义

$$
D=q_0+\frac{q_1}{2}+\frac{q_2}{3},
\qquad
T_\eta=u_1D,
$$

$$
H_\eta(x)=q_0x+\frac{q_1}{2}x^2+\frac{q_2}{3}x^3.
$$

在 `[0,1]` 上做 8 次固定二分：每次取 `mid=(lo+hi)/2`，比较 `H_\eta(mid)` 与 `T_\eta`，用条件选择更新 `lo/hi`。采样热路径不提前退出；退化单元必须在构造 active-cell 表时剔除。

给定 `\eta` 后定义：

$$
A=a+c\eta,
\qquad
B=b+d\eta.
$$

轴向条件 CDF 为：

$$
F_\xi(\xi|\eta)=
\frac{A\xi+\frac{B}{2}\xi^2}{A+\frac{B}{2}}.
$$

同样使用固定 8 步二分反解 `F_\xi(\xi|\eta)=u_2`。定义

$$
T_\xi=u_2\left(A+\frac{B}{2}\right),
\qquad
H_\xi(x)=Ax+\frac{B}{2}x^2.
$$

在 `[0,1]` 上做 8 次固定二分。由于 active-cell 构造阶段已经保证节点源项非负且单元权重大于阈值，热路径只需要固定循环和条件选择，不需要 `|B|` 退化分支或判别式修正。

8 步二分后的归一化区间宽度为 `1/256`，若取区间中点，归一化位置误差不超过 `1/512`。对于 `\Delta r` 或 `\Delta z` 约为 `4e-5` 的网格，物理位置误差约为 `7.8e-8`，低于 `1e-6` 的物理误差目标。

实现要求：

1. 节点源项必须非负；如果数值误差产生极小负值，应在构造采样表时报错或截断并记录诊断。
2. `A_{i,k}` 必须由双线性函数积分得到，不允许用单元平均值替代。
3. 单元内取点必须按 `S(\xi,\eta)(r_k+\eta\Delta r)` 分布采样。
4. 径向和轴向反解都使用固定 8 步二分；不写 Cardano 三次闭式根，也不使用轴向二次解析分支。
5. 若单元四个节点源项全为零，则该单元权重为零，不进入 active-cell alias table；`cell_cdf` 输出中可保留零增量。
6. `IonizationSourceSampler` 使用 active-cell SoA 缓存 `cell_index`、`alias_prob`、`alias_index`、`q0/q1/q2/D`、`a/b/c/d` 和几何系数，采样热路径不得做动态分配。
7. 采样函数应可内联，避免虚调用；全局 `cell_cdf` 二分查找只作为调试或 alias table 构造失败时的 fallback。单元内固定 8 步二分是正式热路径。
8. 不为常数、线性、纯二次等退化单元在采样热路径增加专门分支；这类情况由同一固定 8 步二分处理，非法或零权重单元只在 active-cell 构造阶段处理。

## 13. 验证指标

### 13.1 总电离率

统计表总电离率：

$$
A_{\mathrm{tot}}^{\mathrm{table}}=\sum_{i,k}A_{i,k}.
$$

要求它与碰撞模块记录的累计新增电子权重一致：

$$
A_{\mathrm{tot}}^{\mathrm{table}}
\approx \frac{1}{T}\sum_{n}N_n^{\mathrm{ion}} w_p .
$$

### 13.2 空间源项分布与采样闭环

对 `IonizationSourceSampler` 做大量离线采样，并把采样点重新统计到同一 `r-z` 网格：

$$
\frac{C_{i,k}^{\mathrm{sample}}}{N_{\mathrm{sample}} V_{i,k}}
\propto \bar S_{i,k}.
$$

归一化后，采样分布应恢复双线性源项表给出的单元权重 `A_{i,k}`。该验证不创建粒子，也不依赖注入接口。

### 13.3 呼吸振荡抑制

检查：

1. 每步电离率 `N_pair(t) w_p / dt`。
2. 总电子数和总离子数。
3. 放电电流 `I_d(t)`。
4. 低频谱峰值。

### 13.4 平均放电状态

与完整自洽放电比较：

1. 平均电子密度。
2. 平均离子密度。
3. 平均电势分布。
4. 平均放电电流。
5. 离子束流。
6. 电子温度。

如果平均电子温度显著偏高，说明关闭 MCC 电离能量损失后需要补偿能量通道。

### 13.5 ECDI 保留情况

目标不是消除所有非轴对称结构，而是移除低频电离反馈。需要确认：

1. 低频呼吸振荡削弱。
2. 高频 ECDI 谱仍存在。
3. 高频谱不被平均源项采样噪声淹没。

## 14. 关键结论

当前 WarpX 私有代码已经具备完整自洽电离路径，但平均中性密度方案仍让电离率由瞬时电子动力学决定，不能稳定切断呼吸振荡反馈。

更合适的控制变量是电离源本身：

$$
\text{自洽 MCC 电离事件}
\rightarrow \bar S^{node}_{i,k}
\rightarrow \text{双线性源项单元权重}
\rightarrow \text{CDF}
\rightarrow \text{双线性取点采样器}.
$$

第一版应采用：

1. `r-z` 节点源项表和单元内双线性插值。
2. 固定宏粒子权重 `my_constants.elec_weight`，并利用当前算例两次计算权重一致的前提简化设计。
3. 全量统计，从第一个 MCC 电离步累计到 `Insert::Finalize()`。
4. 四随机数采样：`u0` 选源项单元，`u1,u2,u3` 在单元内按双线性源项和周向角采样三维坐标。
5. AMReX plotfile `ionization_source_plt/` 输出，使结果可直接由 ParaView 打开。
6. 本分支不创建粒子，不接入旧 `AddNParticles()` 注入路径；实际电子-离子对注入等待 `my-dev` 注入重构分支合并后接入其接口。

这样能最大程度复用当前碰撞模块结构，同时把呼吸振荡控制变量从“不稳定的平均中性密度”改为“直接固定的平均电离源统计与采样”。
