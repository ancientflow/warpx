# Insert 诊断与注入参数

`Source/Insert` 按功能组织材料静电、壁面交互、粒子注入、碰撞和诊断，
不按 CHT/PHT 或推力器类型划分。具体材料、Hall 物理项和装置算例仍使用其物理名称。
[命名与模块职责说明](ModuleNaming.md) 列出本次迁移及保留的专用模型。

本文档整理 `Source/Insert` 当前运行时诊断和粒子注入输入参数。参数名按
AMReX `ParmParse` 写法列出，例如 `my_constants.foo` 和 `insert.foo`。

## 已移除的旧边界参数

初始化时的 `BackwardCompatibility()` 会检查以下废弃参数，无论诊断是否开启：

- `my_constants.anode_current_path`：改用 `my_constants.anode_current_prefix`。
- `my_constants.zmin_wall_charge_diag`、`zmin_wall_charge_dir`、
  `zmin_wall_charge_interval` 和 `zmin_wall_charge_write_interval`：旧 zmin 面壁面
  电荷实现已移除，持久壁面电荷由静电介质与解析壁面处理。
- `insert.schur_boundary.*`：旧 zmin 混合边界 Schur 修正已移除，改用体阳极和
  静电介质求解。
- `insert.neutral_atom_eb.*`：改用解析壁面几何及物种对应的壁面行为。

这些参数即使设置为 `0` 也会报出迁移提示。旧输入中的 `my_constants.voltage`
仅是表达式常量；zmin 环形阳极可用 `insert.zmin_ring_anode`（见下节）显式打开，
或在体阳极电势表达式中引用。

## zmin 面环形阳极电势

无陶瓷的简化 Hall 模型（阳极环位于 zmin 域边界面上）可用以下开关：

```text
insert.zmin_ring_anode = 1
```

每次 Poisson 求解前，在 `boundary.potential_lo_z` 写入整面常量之后，把 zmin
面内阳极环覆盖到的节点改写为 `my_constants.voltage`；环外节点保持
`boundary.potential_lo_z` 的值（通常设 `0`）。环几何由 `my_constants` 的
`l_factor`、`L` 决定：环心为 `(prob_lo_x, prob_lo_y) + L/l_factor/2`，内、外半径
分别为 `0.021/(2*l_factor)` 和 `0.031/(2*l_factor)`。仅限 3D，且当
`insert.use_electrostatic_materials = 1`（体阳极 overset 路径）启用时自动跳过。
该设置只改电势，不含壁面电荷沉积、Schur 修正或粒子边界相互作用。

## 静电介质与体阳极

当前统一材料后端用于 3D Cartesian、lab-frame、MLMG、非周期静态网格，不支持 EB、IGF、
relativistic electrostatic solver 或 Python `poissonsolver` callback。

官方 WarpX/ablastr 求解路径由唯一的编译宏
`WARPX_USE_ELECTROSTATIC_MATERIALS` 控制。宏在
`Source/Insert/Config/ElectrostaticMaterialConfig.H` 中只为 3D target 定义。注释对应的
`#define` 行并重新编译，即可使官方目录编译原有求解路径；
Insert 下的材料实现仍可保留在源文件列表中。

启用后，`LabFrameExplicitES` 持有统一材料 Poisson 求解器，其缓存按需初始化。
静态网格、介电常数和阳极 mask 在首次调用时初始化，后续求解复用多重网格层级。
每步求解都恢复阳极指定电势。Dirichlet 已支持局部多层网格，每层持有独立算子和 MLMG；
monopole 仍限完整单层域。细层必须严格位于物理域内部，各方向细化比均为 2。
不检查阳极与粗细接口之间的模板缓冲，也不以此限制 MG 粗化。
局部细层可以不含阳极，此时只使用物理/覆盖约束，不人为添加固定节点。
同层 Box 接缝仍由邻接网格提供数据，不作为粗细固定边界。
详细实现、运行记录和精度限制见 [阶段 4 记录](AMRStage4.md)。
该路径不支持运行期间 regrid、moving window 或重新分配网格，
需设置 `warpx.regrid_int = -1`、`algo.load_balance_intervals = 0`。

同一路径还会分配持久的、与 `rho` 同布局的 `wall_charge` 网格源，初始为零。
已实现的多层壁电荷基础设施在粒子沉积及 `afterdeposition` 回调之后，将原始副本加入本步 `rho`，
再共同执行一次 `SyncRho`、PEC 电荷边界处理及材料约束。
单层保留原有顺序：粒子 `rho` 完成同步、滤波及 PEC 处理之后，才加入独立同步并
清除 ghost 的壁电荷副本；单层壁电荷不参与粒子的滤波及 PEC 镜像处理。
持久源保留各 FAB 及 ghost 上尚未求和的原始贡献，不接收同步或粗层限制副本。
多层要求 `warpx.use_filter = 0` 且 ECDI 关闭，避免隐式改变持久源的物理含义。

撞壁净电荷及交点由 GPU kernel 直接写入临时粒子容器，不经过主机撞击队列；
随后通过 `Redistribute` 路由到交点的
最细覆盖层和所属 MPI box；不使用入射粒子的旧 tile 作为沉积目标。
壁面产物仍由步末的粒子 `Redistribute` 按最终位置分配层级。

可选 `insert.wall_charge_diagnostic = 1` 写入 `wall_charge.dat`，列为
`step time deposited_C stored_C level0_C ... composite_C composite_minus_deposited_C`。
`deposited_C` 为累计撞击净电荷，`stored_C` 与各层值对未组装的原始沉积项求和，
因而不重复计算共享节点或限制副本。`composite_C` 在独立临时副本上执行 `SyncRho`，
以粗细覆盖掩码排除被细网格覆盖的粗层单元，并用每个未覆盖单元角点的 `dV/8`
权重积分。最后一列显式报告复合节点积分与撞击账目的差异，需在粗细界面及
物理边界处检查；该积分位于 PEC 镜像电荷处理和阳极清零之前。诊断默认关闭。

宏启用后，使用以下输入打开材料功能：

```text
insert.use_electrostatic_materials = 1

my_constants.ceramic_top = 2.0e-3
my_constants.ceramic_epsilon_r = 4.0
my_constants.anode_xmin = -5.0e-3
my_constants.anode_xmax =  5.0e-3
my_constants.anode_ymin = -5.0e-3
my_constants.anode_ymax =  5.0e-3
my_constants.anode_zmin =  1.8e-3
my_constants.anode_zmax =  2.0e-3
my_constants.anode_voltage = 300.0

# Cell-centered、无量纲的相对介电常数。这里 z < ceramic_top 为陶瓷。
insert.relative_permittivity_function(x,y,z) = "if(z < ceramic_top, ceramic_epsilon_r, 1.0)"

# Nodal 隐式函数；小于或等于零的节点属于体阳极。
insert.anode_implicit_function(x,y,z) = "if((x >= anode_xmin) * (x <= anode_xmax) * (y >= anode_ymin) * (y <= anode_ymax) * (z >= anode_zmin) * (z <= anode_zmax), -1.0, 1.0)"

# Nodal 阳极固定电势，单位 V；第一版为静态空间表达式。
insert.anode_potential_function(x,y,z) = "anode_voltage"
```

`relative_permittivity_function` 必须在整个 cell-centered 网格上返回有限且严格大于
零的值。等离子体和真空通常取 `1`，陶瓷可取 `4`。输入值是
`epsilon_r`，不能乘入 `epsilon_0`；Poisson 方程右端仍由 WarpX 除以
`epsilon_0`。

阳极表达式在节点上只解析一次，`anode_implicit_function <= 0` 的节点写入 overset
mask `0`，其他节点写入 `1`。体阳极必须非空并至少覆盖两个 z 向 nodal 层。
每次 Poisson 求解时，mask 为零的节点恢复到阳极电势；仅清除私有求解 RHS
上的对应源项，保留沉积场 `rho`。两种材料边界均支持空间变化阳极电压。
陶瓷和阳极均处于统一场计算域中；粒子吸收、壁面自由电荷和阳极电流仍由独立的
Insert 模块处理。

计算域外边界的电势和 ghost nodes 沿用求解器自身的边界处理。原生共置网格的 MLMG
路径启用 `setFinalFillBC(true)`。注意 MLMG 不会填 Dirichlet 面的 ghost 节点
（保持为 0），而 collocated 网格的 `computeE` 在边界节点的中心差分会读到它们；
因此在未启用材料路径时，求解后由 `Insert::SetPhiGuards()`（`DirichletPhiGuardSet`）
按 `phi_ghost = 2*phi_bnd - phi_int - rho_bnd*dx^2/eps0` 重填六个 PEC 面的 ghost
节点，使中心差分退化为 Poisson 一致的单边通量。材料路径使用逐层的独立处理：
粗层电势插值提供细区 ghost，同层通信覆盖内部 box 边界；物理 Dirichlet 面使用
二阶单边导数对应的二次外推；默认 Dirichlet 材料路径仅使用第一类物理边界。电场 ghost 按分量实际
staggering 从粗层插值并执行同层通信，供粒子 gather 使用。
内部体阳极的固定电势仍由上述 overset mask 和 `anode_potential_function` 提供。

### 单层单极 Robin 求解

`insert.electrostatic_boundary` 默认为 `dirichlet`。材料启用时，两种模式统一使用
`Insert::Poisson::MLNodeRobinLaplacian` 和同一缓存生命周期；非材料路径保留原生算子。
Dirichlet 已接入公共逐层求解和粗到细插值；monopole 仍限完整单层域，多层支持留待阶段 5。
Dirichlet 保留原生 `boundary.potential_*` 时间函数及非零值，每步将已设置的节点面值
复制为 `a=1,b=0` 的 forcing。选择 `monopole` 后，同一算子施加六面单极 Robin 条件，
内部阳极仍固定电势：

```text
max_step = 5
amr.max_level = 0
warpx.grid_type = staggered
algo.particle_shape = 1
algo.field_gathering = energy-conserving
interpolation.galerkin_scheme = 1
warpx.use_filter = 0
insert.ecdi_control.enabled = 0
boundary.field_lo = none none none
boundary.field_hi = none none none
boundary.particle_lo = absorbing absorbing absorbing
boundary.particle_hi = absorbing absorbing absorbing
insert.use_electrostatic_materials = 1
insert.electrostatic_boundary = monopole
# 可省略，默认取物理域中心；坐标必须严格位于域内
insert.robin_center = 0.005 0.005 0.00375
insert.robin_phi_infinity = 0
insert.robin_diagnostics = 1
```

本阶段支持 3D Cartesian、完整非周期静态域、单层、staggered 或 collocated、
一阶粒子形状、
关闭 rho/ECDI 滤波，以及域内固定电势阳极。支持显式 labframe 静电多步 PIC 演化；
静磁耦合与隐式演化尚未接入此边界模式。沿用材料路径对 EB、Python Poisson 回调、重网格、负载均衡和 moving window
的限制。`robin_center` 和 `robin_phi_infinity` 支持参数表达式；后者默认 0 V。
阳极不得接触外域边界；`anode_potential_function(x,y,z)` 可以随空间变化。

入口复用材料 rho/epsilon/mask，但不对外面执行旧 PEC rho 镜像和 Dirichlet phi 设置。
输出 rho 保持原沉积 `q/ΔV` 约定。在独立 RHS 副本上，对外域面/棱/角节点分别
乘 2/4/8，匹配节点算子的半/四分之一/八分之一控制体积，再乘 `-1/epsilon0`；
算子自行加入 Robin 源。此过程不改写 rho 或持久壁电荷。壁电荷复合诊断的临时副本
在 Robin 模式下使用同一控制体积修正；默认 Dirichlet 诊断不变。
算子与 MLMG 持久缓存，网格、材料/mask 对象或层次选项变化时拒绝继续使用旧缓存。
粗化上限逐级检查阳极 z 厚度、细节点类别、间隙和介质跳变，
并遵守 `self_fields_max_coarsening_level`；详见 [阶段 3 记录](AMRStage3.md)。
为确保各 Box 可安全二倍粗化，所选上限可能随 Box 划分而更保守。
关闭半粗化，仅支持原生 smoother/Krylov 底层求解器。

`robin_diagnostics=1` 对两种材料模式都生效，输出实际 MG 层数和固定节点数、
粗化上限、有效 pre/post、迭代次数、材料层次构建次数、`compPhysicalResidual` 最大范数，
以及准备/求解/额外残差计算耗时；提高 MLMG/bottom 日志等级以记录 bottom 失败或回退。
Dirichlet 多层诊断另检查固定节点误差（必须为零）并统计 Box 边界上的未知节点条目；
条目统计包含同层共享节点副本。诊断默认关闭，不额外分配残差场。

Robin 面将边界节点保留为未知数。两种模式均通过统一选项入口应用
`warpx.self_fields_num_pre_smooth` 和 `warpx.self_fields_num_post_smooth`，默认 2/2；
旧材料 Dirichlet 的硬编码 1/2 已移除。增加扫描次数的收益取决于具体离散问题，
应以同一固定源、相同容差的迭代数和耗时验证，不能用单个独立算例的收敛因子概括。

六个场边界必须显式设为 `none`；这使通用场处理不再执行 PEC/PMC 镜像或切向场清零，
真正的静电边界由 Robin 算子负责。粒子吸收/反射仍由 `boundary.particle_*` 独立控制。
外加 E/B 仍按现有顺序叠加，外加 E 不经过 PEC 截断。

Robin 模式只允许 `algo.field_gathering=energy-conserving`（默认值），拒绝
`momentum-conserving`。交错网格必须保持 `interpolation.galerkin_scheme=1`；
共置网格按 WarpX 原有规则使用同阶节点形函数。此设置约束取场算法，不意味着
含吸收边界、碰撞或外加场的整个系统严格守恒能量。

`Insert::FillPhiGhosts` 统一分派旧边界、材料 Dirichlet 和 Robin 的电势 ghost 处理。
Robin 直接使用原 `phi_fp`，执行同层通信，不再分配宽电势副本或计算域外 E。
公共 `ElectrostaticSolver::computeE` 仅计算有效区，E 的同层通信仍保留。

一阶交错 Galerkin 取场在 E 法向使用零阶形函数，不经过交错到节点的场居中；
域内粒子不需要物理域外 E。有效区前向差分只使用物理域内及边界电势节点，
因此交错模式不再填充 Robin 域外电势。
共置网格的一阶取场也不需要物理域外 E，但边界节点的中心差分需要第一层 phi ghost：
`phi_ghost = phi_interior - 2*dx*q*(phi_boundary-phi_infinity)`，其中 q 使用外法向。
只填充各面紧邻的这一层、其余坐标处于物理域内的节点，覆盖有效区面/棱/角求导。
这使中心差分自然满足 `E_n=q*(phi-phi_infinity)`；其他域外 phi ghost 不作物理延拓。

当前仍拒绝高阶粒子形状和 rho 滤波：它们需要独立的外域沉积处理。
单极延拓是远场近似，不保证任意近场分布或延伸至外边界的陶瓷都适合这个模型。
32³/64³ 真空球体验证显示误差随网格细化下降；分段介质算例已验证运行、缓存和 MPI
一致性，但实际羽流仍须分别检查域大小、展开中心及介质外延的敏感性。

## 解析壁面粒子相互作用

解析壁面由全局列表声明，列表顺序为交界处的优先级。无效区域必须由输入保证不相交。
每个物种可独立选择与每个壁面是否交互。配置了
`<species>.analytic_wall.<wall>.material` 时优先使用对应材料模型，忽略同一物种、
同一壁面上的 `behaviors` 和 `p_*` 参数；未配置材料时保持通用解析边界处理。
两者都未配置时，该物种不与此解析壁面交互。未知材料名称会报错，不会退回通用路径。

目前支持 `material = ceramic` 和 `material = stainless_steel`，
要求入射和二次电子物种均使用物理电子质量和电荷。
`secondary_electron_species` 默认与入射物种相同，
陶瓷的 `secondary_electron_temperature_eV` 默认 `3.0`，`deposit_wall_charge` 默认 `1`。
材料路径仍使用同一解析几何、壁面电荷场和电流诊断，不重复执行通用碰撞处理。例如：

```text
electrons.analytic_wall.ceramic.material = ceramic
electrons.analytic_wall.ceramic.secondary_electron_species = electrons
electrons.analytic_wall.ceramic.secondary_electron_temperature_eV = 3.0
electrons.analytic_wall.ceramic.deposit_wall_charge = 1
```

不锈钢采用 Furman–Pivi 的 SLAC 304 未处理表面参数集，独立实现弹性背散射、
再扩散和多电子真二次发射能谱，三类出射均使用余弦角分布：

```text
electrons.analytic_wall.anode.material = stainless_steel
electrons.analytic_wall.anode.secondary_electron_species = electrons
electrons.analytic_wall.anode.binomial_trials = 10
electrons.analytic_wall.anode.deposit_wall_charge = 1
```

`binomial_trials` 为文档二项方案的 M，默认 10，允许 1–10；拟合系数与能谱表编码在模型中，
不使用陶瓷的温度参数。真二次能量按受总能量约束的联合分布采样。
实现假定预处理提供有效入射状态，不检查输入状态异常，也不限制拟合能量或角度范围；
超出拟合范围时继续使用相同公式。仍需满足
`A=1-delta_e-delta_r>=0`（A=0 时需 delta_ts=0）及 `mu=delta_ts/A<M`。
例如 100 eV、30° 已不满足 M=10，不能仅按能量范围判断可用性。
概率或条件均值不满足上述要求时，整个撞击直接回退为完全吸收，不再中止计算：
不产生任何出射电子，删除入射粒子，并在启用壁面电荷沉积时沉积完整入射电荷。
这是一项额外的物理近似，会降低这些撞击的有效二次电子产额；有效撞击仍使用原模型。
不裁剪概率，也不截断泊松尾部。
详见 [材料实现说明](../Boundary/MaterialInteraction.md)。

通用路径中，`behaviors` 的顺序定义累积概率抽样；除
最后一项外，每项以 `p_<behavior>(E_eV,u_n,u_t,x,y,z,t)` 给出概率，最后一项采用余量。
支持 `absorb`、`specular`、`diffuse`、`convert`、`secondary_electron_1` 和
`secondary_electron_2`，最多六项。

```text
insert.analytic_walls = ceramic anode

# F > 0 是有效的等离子体区域，F < 0 是该壁面的无效区域。
analytic_wall.ceramic.signed_value(x,y,z) = "z - ceramic_top"
analytic_wall.ceramic.normal_x(x,y,z) = "0.0"
analytic_wall.ceramic.normal_y(x,y,z) = "0.0"
analytic_wall.ceramic.normal_z(x,y,z) = "1.0"

electrons.analytic_wall.ceramic.behaviors = absorb secondary_electron_1 diffuse
electrons.analytic_wall.ceramic.p_absorb(E_eV,u_n,u_t,x,y,z,t) = "0.8"
electrons.analytic_wall.ceramic.p_secondary_electron_1(E_eV,u_n,u_t,x,y,z,t) = "0.1"
electrons.analytic_wall.ceramic.secondary_electron_species = electrons
electrons.analytic_wall.ceramic.secondary_electron_temperature_eV = 3.0
electrons.analytic_wall.ceramic.wall_temperature = 400.0

ions.analytic_wall.anode.behaviors = convert specular
ions.analytic_wall.anode.p_convert(E_eV,u_n,u_t,x,y,z,t) = "0.5"
ions.analytic_wall.anode.product_species = xe_neutral
ions.analytic_wall.anode.wall_temperature = 400.0
ions.analytic_wall.anode.product_charge = 0.0
ions.analytic_wall.anode.deposit_wall_charge = 1
```

两种路径共享预处理：先收集有效且位于对应边界外的粒子索引，再只对候选粒子计算
撞击点、撞击速度、法向和剩余时间。后续交互、产物创建、沉积及诊断按候选数量遍历；
不再为域内粒子分配碰撞信息数组。索引压缩可能改变随机数与粒子的对应顺序，
相同种子下的逐粒子结果不保证与旧遍历方式完全一致。

预处理后，通用路径和材料路径各自完成事件选择、目标空间分配、反射与发射。
事件只抽样一次，材料路径自行按其能谱生成并写入粒子，不经过通用发射函数。
两者返回是否保留原粒子、实际新增数量和含宏粒子权重的出射总电荷；公共收尾据此
沉积电荷、标记删除并记录诊断，保证在发射及初始化完成之后才删除入射粒子。
陶瓷 SEE2 保留两个独立的发射采样；通用 SEE2 保留原有的共享速度采样。

相互作用仅在物种实际推进的步执行。`specular` 与 `diffuse` 保留入射粒子并从交点推进
剩余子步；`convert` 和二次电子发射使入射粒子失效，并以相同宏粒子权重创建目标
物种。转换产物的温度为 `wall_temperature`（K），二次电子产物的温度为
`secondary_electron_temperature_eV`。壁面获得的电荷为
`(q_in - sum(q_out)) * weight`，以交点形函数沉积到持久 `wall_charge`。交点坐标及
形函数权重使用 `ParticleReal`；只有写入场 FAB 时转换为 `Real`。该功能要求前述材料
Poisson 路径已启用，以分配和求解持久壁面电荷。

两个可选参数控制壁面电荷沉积。`deposit_wall_charge`（默认 `1`）为 `0` 时，该物种在此
壁面的吸收类事件不向 `wall_charge` 沉积电荷（适用于固定电势的导体壁面）。
`product_charge` 覆盖转换产物的电荷（默认为产物物种的 `charge`）：当产物物种
为了经电荷沉积获得数密度场而携带记账用非零电荷（如 `charge = 1`）时，应显式设为
 `0.0`，使转换事件的壁面沉积保持物理正确的 `(q_in - 0) * weight`。

## 运行时诊断

### 总体节奏

壁面电流和边界粒子缓存相关诊断共享同一个间隔：

```text
my_constants.boundary_particle_diag_interval = 10
```

旧参数 `my_constants.hall_diag_interval` 和
`my_constants.clear_hall_boundary_particle_cache_diag` 继续作为兼容别名，分别映射到
`my_constants.boundary_particle_diag_interval` 和
`my_constants.clear_boundary_particle_cache_diag`。新旧名称同时出现时必须取相同值，
否则初始化报错。

默认值为 `10`，内部会限制为至少 `1`。`AnodeCurrentDiagOutput` 使用该间隔
作为阳极电流的采样窗口（写出并清零累加器）。以下路径均使用该间隔读取
`ParticleBoundaryBuffer`：

- `ThrustCalc`
- `BeamDivergenceCalc`
- `IEDFCalc`
- `ClearBoundaryParticleCache`

这意味着边界粒子缓存的读取和清理使用同一时步判断，避免诊断间隔不一致导致
粒子在被统计前被清掉。

### 开关与输出

| 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `my_constants.particle_number_diag` | `0` | 打印每个 species 当前粒子数。需要编译宏 `NUMP`。 |
| `my_constants.wall_interaction_diag` | `0` | 每步按入射 species 打印解析壁面交互的宏粒子计数，同一物种的所有解析壁面合并统计。需要启用 `insert.analytic_walls`。 |
| `my_constants.collision_record_diag` | `0` | 写出碰撞产生的电子和离子宏粒子数到 `collision_record.dat`。需要编译宏 `COLLISION_RECORD`。 |
| `my_constants.anode_current_diag` | `0` | 统计每面解析壁面吸收的电流，每壁面一个输出文件。需要启用 `insert.analytic_walls`。 |
| `my_constants.anode_current_prefix` | `anode_current` | 壁面电流输出文件前缀，文件名为 `<prefix>_<wall>.dat`。 |
| `my_constants.thrust_diag` | `0` | 统计出口离子轴向动量并计算推力。需要 `WARPX_INSERT_3D`。 |
| `my_constants.beam_divergence_diag` | `0` | 统计出口离子束流发散角。需要 `WARPX_INSERT_3D`。 |
| `my_constants.iedf_diag` | `0` | 统计出口离子能量分布函数。需要 `WARPX_INSERT_3D`。 |
| `my_constants.clear_boundary_particle_cache_diag` | `0` | 按 `boundary_particle_diag_interval` 清理边界粒子缓存。需要 `WARPX_INSERT_3D`。 |

`clear_boundary_particle_cache_diag` 应在所有读取边界缓存的诊断之后执行。
当前 `Insert::AfterDiagnostics()` 中的调用顺序已满足这一点。

### 解析壁面交互计数

```text
my_constants.wall_interaction_diag = 1
```

启用后，每步在标准输出中按入射 species 分别列出，每个物种两行：第一行为
物种名，第二行依次为：

- `removed`：被移除的入射宏粒子总数（只给总和，即 `absorb`、`convert`
  与二次电子发射事件之和）。
- `specular` / `diffuse`：镜面反射 / 漫反射事件数。
- `emitted products`：转换产生的产物宏粒子数，等于 `convert`。
- `emitted secondaries`：产生的二次电子宏粒子数，等于 `SEE1 + 2 * SEE2`；
  括号内 `SEE1`、`SEE2` 分别为 `secondary_electron_1` 和
  `secondary_electron_2` 的事件数。

计数跨 MPI ranks 和所有解析壁面求和，只由 IO rank 输出；仅统计配置了解析
壁面策略的物种。没有事件（包括该物种本步未推进）时输出零。每步重新计数，
不受 `boundary_particle_diag_interval` 控制，不包括计算域外边界和 EB 的交互。
所有数值均为宏粒子事件计数，不是权重之和，也不是 species 净增减量。
产物数归属于入射物种：例如 `xe_ions` 的 `emitted products` 表示由离子撞壁
产生的原子，不会记入 `xe_neutral` 的入射事件。

默认关闭，替代此前无条件输出的整体壁面计数。关闭时跳过事件计数的归约，
不影响壁面交互、壁面电荷沉积或独立的 `anode_current_diag` 诊断。

### 阳极电流

```text
my_constants.anode_current_diag = 1
my_constants.anode_current_prefix = "anode_current"
```

壁面电流统计基于解析壁面（`insert.analytic_walls`）：每步在
`AnalyticBoundaryInteraction` 的壁面吸收 kernel 中累加，按
`boundary_particle_diag_interval` 在诊断阶段写出并清零，即每行是一个采样窗口内的
累计量。**每面解析壁面一个输出文件**，文件名为
`<anode_current_prefix>_<wall>.dat`（如 `anode_current_anode_ring.dat`、
`anode_current_ceramic_wall.dat`）；阳极电流即阳极壁面对应的文件。
`anode_current_prefix` 默认值为 `anode_current`。每个文件的输出列为：

```text
time    electron    ion    net_charge
```

统计对象（只计入 absorb / convert / 二次电子等移除入射粒子的行为，
specular / diffuse 反射不产生净电流）：

- `electron` / `ion`：窗口内被该壁面移除的电子/离子宏粒子权重之和
  （按物种电荷正负分类）。
- `net_charge`：窗口内该壁面接收的净电荷，单位 C，口径与壁面电荷沉积
  一致，即 `(charge - q_out) * w`；除以窗口时长即为平均电流。

对任何物种都没有配置 behaviors 的壁面也会生成文件（内容恒为零），
便于后处理脚本统一处理。旧参数 `my_constants.anode_current_path`（单文件
zlo 缓存实现）已移除，设置它会触发报错并提示改用
`anode_current_prefix`。

注意：旧实现从 zlo 域边界粒子缓存统计，新解析壁面吸收的粒子不再进入该
缓存，故旧实现已删除。一步内穿越壁面到 domain zmin 之间区域的粒子
（`v * dt` 大于壁面到 zmin 的距离）会被域边界吸收而不经壁面模型，两种
实现都统计不到。

### 推力

```text
my_constants.thrust_diag = 1
my_constants.thrust_diag_path = "thrust.dat"
```

`thrust_diag_path` 默认值为 `thrust.dat`。统计出口边界：

```text
xlo xhi ylo yhi zhi
```

输出列为：

```text
time    dt    ion_weight    thrust_N    axial_momentum_kg_m_per_s
```

`dt` 为本次诊断与上次诊断之间的物理时间；首次统计时如果无法由时间差得到，
使用 `boundary_particle_diag_interval * WarpX::getdt(0)` 估算。

### 束流发散

```text
my_constants.beam_divergence_diag = 1
my_constants.beam_divergence_path = "beam_divergence.dat"
```

`beam_divergence_path` 默认值为 `beam_divergence.dat`。输出列为：

```text
time    dt    ion_weight    divergence_angle_rad    divergence_angle_deg
axial_momentum_kg_m_per_s    transverse_momentum_kg_m_per_s
```

发散角按 `atan2(transverse_momentum, axial_momentum)` 计算。

### IEDF

```text
my_constants.iedf_diag = 1
my_constants.iedf_path = "iedf.dat"
my_constants.iedf_bins = 200
my_constants.iedf_min_eV = 0
my_constants.iedf_max_eV = 500
```

参数：

| 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `my_constants.iedf_path` | `"iedf.dat"` | 输出文件。 |
| `my_constants.iedf_bins` | `200` | 能量 bin 数，内部限制为至少 `1`。 |
| `my_constants.iedf_min_eV` | `0` | 能量下限。 |
| `my_constants.iedf_max_eV` | `500` | 能量上限；若不大于下限，会设为 `min + 1`。 |

输出列为：

```text
time    bin_center_eV    bin_lo_eV    bin_hi_eV    ion_weight
current_A    pdf_per_eV
```

每次诊断块写完后会额外写一个空行。

## 粒子注入

### 入口参数

当前 `WARPX_INSERT_3D` / `WARPX_INSERT_INITIALIZATION_3D` 注入入口统一走 `ParticleInjector`。全局入口参数为：

```text
insert.initial_sources = source1 source2
insert.continuous_sources = source3 source4
my_constants.dt = 3e-12
```

- `insert.initial_sources` 在初始化阶段执行一次。
- `insert.continuous_sources` 在每个粒子注入步执行。
- 两个数组为空时，不会创建任何默认 source。
- 连续注入使用 `my_constants.dt` 作为传入 source 的时间步。若未设置，该值为 `0`，
  依赖 `dt` 的速率模型不会产生粒子。`dt` 支持数学表达式（含 `my_constants`
  中定义的其他常量），例如 `my_constants.dt = 5e-8 / l_factor`。

### Source 基本结构

普通 source 使用 source 名作为参数前缀：

```text
<source>.rate = current
<source>.species = electrons
<source>.macro_weight = elec_weight

<source>.position.coordinate_system = cylindrical
<source>.position.r.distribution = area_uniform
...

<source>.velocity.coordinate_system = cartesian
<source>.velocity.vx.distribution = gaussian
...
```

多物种 source 使用同一组位置样本：

```text
<source>.shared_species = electrons xe_ions
```

每个物种的粒子权重读取优先级为：

```text
<source>.<species>.weight
<source>.macro_weight
1.0
```

每个物种的速度分布读取优先级为：

```text
<source>.<species>.velocity.*
<source>.velocity.*
```

可选 source 级参数：

| 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `<source>.x_offset` | `0.0` | 对采样位置 x 加偏移。 |
| `<source>.y_offset` | `0.0` | 对采样位置 y 加偏移。 |
| `<source>.batch_multiplier` | 未启用 | 批量注入阈值倍数；设置时必须 `>= 1.0`。 |

### 注入数量模型

`<source>.rate` 支持以下取值：

| rate | 必填参数 | 可选参数 | 宏粒子期望数 |
| --- | --- | --- | --- |
| `fixed_count` | `count` | 无 | `count` |
| `density_volume` | `density`, `volume`, `macro_weight` | 无 | `density * volume / macro_weight` |
| `current` | `current`, `macro_weight` | `l_factor = 1.0` | `current * dt / (l_factor^2 * q_e * macro_weight)` |
| `mass_flow` | `mass_flow`, `mass`, `macro_weight` | `l_factor = 1.0` | `mass_flow * dt / (l_factor^2 * mass * macro_weight)` |

`fixed_count` 和 `density_volume` 不使用 `dt`，适合一次性初始化。
`current` 和 `mass_flow` 用于连续注入。

### 坐标系统

位置分布默认使用 `cylindrical`，速度分布默认使用 `cartesian`。

| coordinate_system | position 轴名 | velocity 轴名 |
| --- | --- | --- |
| `cartesian` | `x`, `y`, `z` | `vx`, `vy`, `vz` |
| `cylindrical` | `r`, `theta`, `z` | `vr`, `vtheta`, `vz` |
| `local_normal` | 不支持 | `vnormal`, `vt1`, `vt2` |

`local_normal` 速度坐标会根据采样位置中的法向量转换到笛卡尔速度。
普通坐标分布产生的位置法向量默认为 `(0, 0, 1)`。

### 一维分布

每个坐标轴通过 `<prefix>.distribution` 指定一维分布：

| distribution | 必填参数 | 可选参数 | 说明 |
| --- | --- | --- | --- |
| `constant` | `value` | 无 | 固定值。 |
| `uniform` | `min`, `max` | 无 | 区间均匀分布，要求 `max > min`。 |
| `area_uniform` | `min`, `max` | 无 | 面积均匀半径分布，要求 `min >= 0` 且 `max > min`。 |
| `gaussian` | `sigma` | `mean = 0` | 高斯分布，要求 `sigma >= 0`。 |
| `positive_gaussian` | `sigma` | `mean = 0` | 重采样直到值非负。 |
| `flux_gaussian` | `sigma` | 无 | 通量加权半麦克斯韦（Rayleigh）法向分量，`p(v) ∝ v exp(-v²/2σ²)`，`v >= 0`，精确逆 CDF 采样；用于平衡态储罐/缓冲腔的孔口入口。 |
| `single_spoke` | `sigma` 或 `spoke_sigma` | `center` / `spoke_center = pi`, `num_bins = 1024` | 单 spoke 角向分布。 |
| `multi_spoke` | `sigma` 或 `spoke_sigma` | `spoke_count = 1`, `spoke_phase = 0`, `num_bins = 1024` | 多 spoke 角向分布。 |
| `neutral_spoke` | `ion_width` | `min_ratio = 0.25`, `drop_exponent = 4.0`, `phase` / `phase_offset` / `spoke_phase = 0`, `reverse = 0`, `spoke_count = 1`, `num_bins = 1024` | 中性气体 spoke 耗尽分布；`spoke_count > 1` 时结构在圆周上重复，要求 `ion_width < 2*pi/spoke_count`。 |
| `neutral_spoke_depletion` | 同 `neutral_spoke` | 同 `neutral_spoke` | `neutral_spoke` 的别名。 |
| `discrete` | `values` | `weights` | 离散采样；未给 `weights` 时等权。 |
| `tabulated` | `values`, `pdf` | `num_bins = 1024` | 表格 PDF，`values` 必须严格递增。 |
| `parser` | `min`, `max`, `function(<axis>)` 或 `function` | `num_bins = 1024` | 解析表达式 PDF。 |

### hole_array_plane 位置耦合分布

孔阵列是更高一层的坐标分派层，本身不实现任何分布，必须与
`<source>.position.*` 的常规分布搭配；底层分布定义单孔内的局部坐标：

```text
<source>.position.coupled_distribution = hole_array_plane
<source>.hole_count = 48
<source>.hole_ring_radius = <ring_radius>

<source>.position.coordinate_system = cylindrical
<source>.position.r.distribution = parser
<source>.position.r.min = 0.0
<source>.position.r.max = <hole_radius>
<source>.position.r.function(r) = r * flux(r)
<source>.position.theta.distribution = uniform
<source>.position.theta.min = 0.0
<source>.position.theta.max = 2*pi
<source>.position.z.distribution = constant
<source>.position.z.value = 0.0
```

参数：

| 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `<source>.hole_count` | `48` | 孔数量，必须 `> 0`。 |
| `<source>.hole_ring_radius` | 若缺省则读取 `<source>.ring_radius` | 孔中心所在环半径，必须 `>= 0`。 |

粒子按序号轮流分派到各孔（每次注入随机起始孔），采样位置连同局部
坐标原点一起偏移到孔中心；`sigma(r)/mean(r)` 的 `r` 与局部方位角均
以孔中心为参考。孔径与孔内剖面由底层 `position.*` 分布定义（例如
`r` 轴用 `function(r) = r * flux(r)` 实现通量加权）。速度仍按普通
velocity 分布配置。

### 平均电离源

编译宏 `IONIZATION_SOURCE_INJECT` 打开时，`insert.initial_sources` 或
`insert.continuous_sources` 可列出特殊 source：

```text
insert.continuous_sources = average_ionization_source
```

该 source 不读取普通 `<source>.*` 参数，而是读取 `ionization_source_fab` 目录：

```text
ionization_source_fab/metadata.txt
ionization_source_fab/node_source_rate
ionization_source_fab/cell_cdf
```

`metadata.txt` 需要包含：

```text
center = <x_center> <y_center>
r_min = <value>
r_max = <value>
z_min = <value>
z_max = <value>
nr = <value>
nz = <value>
A_tot = <total_rate>
elec_weight = <weight>
```

运行时会要求 `metadata.txt` 中的 `elec_weight` 与当前
`my_constants.elec_weight` 一致。注入物种固定为 `electrons` 和 `xe_ions`，
权重均为 `my_constants.elec_weight`。

### 常用配置片段

连续阴极电子注入：

```text
insert.continuous_sources = cathode_electron

cathode_electron.species = electrons
cathode_electron.rate = current
cathode_electron.current = Ic
cathode_electron.macro_weight = elec_weight
cathode_electron.l_factor = l_factor
cathode_electron.x_offset = 0.025 / l_factor
cathode_electron.y_offset = 0.025 / l_factor

cathode_electron.position.coordinate_system = cylindrical
cathode_electron.position.r.distribution = area_uniform
cathode_electron.position.r.min = 0.016 / l_factor
cathode_electron.position.r.max = 0.020 / l_factor
cathode_electron.position.theta.distribution = uniform
cathode_electron.position.theta.min = 0.0
cathode_electron.position.theta.max = 2*pi
cathode_electron.position.z.distribution = uniform
cathode_electron.position.z.min = 0.042 / l_factor
cathode_electron.position.z.max = 0.046 / l_factor

cathode_electron.velocity.coordinate_system = cartesian
cathode_electron.velocity.vx.distribution = gaussian
cathode_electron.velocity.vx.sigma = 592982
cathode_electron.velocity.vy.distribution = gaussian
cathode_electron.velocity.vy.sigma = 592982
cathode_electron.velocity.vz.distribution = gaussian
cathode_electron.velocity.vz.sigma = 592982
```

中性气体质量流量注入：

```text
insert.continuous_sources = xe_neutral_inlet

xe_neutral_inlet.species = xe_neutral
xe_neutral_inlet.rate = mass_flow
xe_neutral_inlet.mass_flow = m_dot
xe_neutral_inlet.mass = m_xe
xe_neutral_inlet.macro_weight = xe_weight
xe_neutral_inlet.l_factor = l_factor
xe_neutral_inlet.x_offset = 0.025 / l_factor
xe_neutral_inlet.y_offset = 0.025 / l_factor

xe_neutral_inlet.position.coupled_distribution = hole_array_plane
xe_neutral_inlet.hole_count = 48
xe_neutral_inlet.hole_ring_radius = (0.021 + 0.031) / 4 / l_factor

xe_neutral_inlet.position.coordinate_system = cylindrical
xe_neutral_inlet.position.r.distribution = area_uniform
xe_neutral_inlet.position.r.min = 0.0
xe_neutral_inlet.position.r.max = 0.001 / l_factor
xe_neutral_inlet.position.theta.distribution = uniform
xe_neutral_inlet.position.theta.min = 0.0
xe_neutral_inlet.position.theta.max = 2*pi
xe_neutral_inlet.position.z.distribution = constant
xe_neutral_inlet.position.z.value = 0.0

xe_neutral_inlet.velocity.coordinate_system = cartesian
xe_neutral_inlet.velocity.vx.distribution = gaussian
xe_neutral_inlet.velocity.vx.sigma = sqrt(kb * Tx / m_xe)
xe_neutral_inlet.velocity.vy.distribution = gaussian
xe_neutral_inlet.velocity.vy.sigma = sqrt(kb * Ty / m_xe)
xe_neutral_inlet.velocity.vz.distribution = positive_gaussian
xe_neutral_inlet.velocity.vz.mean = vz0
xe_neutral_inlet.velocity.vz.sigma = sqrt(kb * Tz / m_xe)
```

圆锥面可将上例的位置配置替换为：

```text
xe_neutral_inlet.position.coupled_distribution = truncated_cone_surface
xe_neutral_inlet.position.slope = eb_k
xe_neutral_inlet.position.r_min = inlet_r_min
xe_neutral_inlet.position.r_max = inlet_r_max
xe_neutral_inlet.position.r_reference = eb_a1
xe_neutral_inlet.position.z_reference = 0.0
xe_neutral_inlet.position.theta_min = 0.0
xe_neutral_inlet.position.theta_max = 2*pi
```

采样器返回
`z = z_reference + slope * (r - r_reference)` 上的坐标；径向采用
`r^2` 均匀采样，方位角采用均匀采样，已包含圆锥面积元中的极坐标权重。
`z` 由采样得到的 `r` 唯一确定，不应再配置独立的
`position.z.distribution`。

若需要绕装置轴线旋转的非各向同性进气速度，可配置：

```text
xe_neutral_inlet.velocity.coordinate_system = rotating_axis
xe_neutral_inlet.velocity.axis_at_theta0 = inlet_axis_r inlet_axis_theta inlet_axis_z
xe_neutral_inlet.velocity.vx.distribution = gaussian
xe_neutral_inlet.velocity.vx.mean = 0.0
xe_neutral_inlet.velocity.vx.sigma = sqrt(kb * Tx / m_xe)
xe_neutral_inlet.velocity.vy.distribution = gaussian
xe_neutral_inlet.velocity.vy.mean = 0.0
xe_neutral_inlet.velocity.vy.sigma = sqrt(kb * Ty / m_xe)
xe_neutral_inlet.velocity.vz.distribution = positive_gaussian
xe_neutral_inlet.velocity.vz.mean = vz0
xe_neutral_inlet.velocity.vz.sigma = sqrt(kb * Tz / m_xe)
```

`axis_at_theta0` 在 `theta=0` 处的局部 `(r, theta, z)` 基中定义并自动
归一化。采样器先生成局部速度，再将局部 `+z` 轴对齐到该向量，最后绕
全局 `z` 轴旋转到粒子相对 source 轴心的方位角。`vz` 分布必须保证采样
值非负；注入流量仍由 source 的 rate model 和宏粒子权重决定。

## 相关源文件

- `Source/Insert/Diagnostics/InsertRuntimeDiagnostics.cpp`
- `Source/Insert/Utils/InsertUtils.cpp`
- `Source/Insert/Injection/ParticleInjector.cpp`
- `Source/Insert/Injection/InjectionRateModel.cpp`
- `Source/Insert/Injection/CoordinateDistribution.cpp`
- `Source/Insert/Injection/Distribution1D.cpp`
- `Source/Insert/Injection/ParticleInjectionSource.cpp`
- `Source/Insert/Geometry/TruncatedConeSurfaceSampler.cpp`
