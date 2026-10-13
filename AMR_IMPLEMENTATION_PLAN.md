# WarpX 原生 AMR 改造实施计划

> 目标：从当前单层 Hall 静电实现出发，引入两层静态网格细化（MR），
> 让 level 0 覆盖扩展域并携带羽流粒子，level 1 保留通道附近的细网格。
>
> 实施基线：分支 `AMR-support`，HEAD `b50110b2b`。
> 阶段 0 实施时实际分支为 `Amr-support`，HEAD 为 `01f1b2bfca2566e56164b8c008636a724885f946`。
> 下列复选框记录实施状态；阶段 0 的外部数据版本仍待补齐。

> 2026-10-13 调整：阶段 0 的配置与阶段 1 的功能代码作为独立提交保存，
> 构建及数值验收状态不因此改为完成。暂缓阶段 2 及后续 AMR 开发，
> 当前优先保留单层求解行为，并规划单层单极 Robin 远场边界；见第 6 节。
> 第 3 节保留改造前的问题清单，实际实现状态以第 4 节记录为准。

## 1. 当前实现基础

- 当前主算例是 `Script/PHT/3d_hall_dielectric_slot`：`labframe`、MLMG、staggered、
  `amr.max_level = 0`，包含介电材料、体阳极、解析墙和耦合 MCC。
  原计划的 collocated 描述与当前输入不一致；阶段 0 沿用 staggered，
  后续场计算、ghost 和 gather 适配需覆盖实际使用的布局。
- 当前域为 `[0, 10] × [0, 10] × [-2.5, 10] mm`，网格 `256×256×320`，
  各方向间距均为 `39.0625 μm`；六面电势为零，粒子边界为 absorbing。
- `Source/FieldSolver/ElectrostaticSolvers/LabFrameExplicitES.cpp` 已使用多层
  `rho_fp`、`phi_fp`，包含粒子沉积、`SyncRho`、逐层求解和 `computeE`。
- `Source/ablastr/fields/PoissonSolver.H::computePhi` 是逐层求解：先解粗层，
  再由 `interpolatePhiBetweenLevels` 为细层提供边界数据及初值。
- 材料存储、壁电荷存储和解析墙遍历有逐层骨架，但仍有明确的单层断言或索引假设。
- `LabFrameExplicitES` 的电势初值外推已逐层存储，改造时保留。

## 2. 几何与物理约定

level 0 为完整扩展域，level 1 为局部细化区；初版使用 `amr.max_level = 1`
和 `amr.ref_ratio_vect = 2 2 2`，不启用动态 regrid 或 AMR subcycling。
层数、Geometry、BoxArray、DistributionMapping、介电系数、阳极 mask、边界类型
及算子构造参数在初始化后固定，材料数据和求解器缓存按此约定持续复用。
静电不支持 `warpx.do_subcycling = 1`；算例的 species `ndt` 按独立的物种推进路径处理。

阶段 0 选定的扩展域保持旧域中心，将每个方向长度加倍：
`[-5, 15] × [-5, 15] × [-8.75, 16.25] mm`。
level 0 取 `256×256×320`，粗层间距为 `78.125 μm`，level 1 恢复原分辨率。
细化区选为原计算域，z 向覆盖体阳极及全部解析注入源；外场覆盖与羽流截断影响
仍须在阶段 1–2 验证，不能由几何配置推断其物理收敛性。

- 保持 `lx/ly/lz`、阳极中心、环缝、注入面等物理参数含义，单独定义扩展域端点，
  不通过放大这些共享常数隐式改变装置尺寸、粒子权重或注入位置。
- 初始 tagging 使用 `warpx.ref_patch_function(x,y,z)`，
  或 `warpx.fine_tag_lo/hi`。
  结合 blocking factor 和网格聚类确定最终 BoxArray 的 patch 覆盖范围。
  各层 Geometry 共享同一物理域，局部细化范围由该层 BoxArray 决定。
- 细区覆盖阳极、主要介质界面、注入区、主要电离区及必要的下游缓冲，
  将主要电离区与细化边界隔开，留出背景密度沉积和采样模板所需的距离。
  粗层主要承载低密度羽流；本工况采用细化界面附近背景密度缺项对电离影响可忽略的近似。
- 当前 `relative_permittivity_function` 对所有 `z < 0` 都给陶瓷介电常数；
  `ceramic_wall.signed_value` 对阳极环之外的整个平面生效。
  扩展域后，同一物理壁面同时穿过粗层独占区和细化区。两层使用相同的物理坐标、
  墙面表达式和材料参数，分别处理各自粒子的撞壁、吸收、反射、转化和二次发射。
  细化边界是数值分区，壁面及其物理响应在此保持连续。
  初版沿用当前平面模型在扩展域内的定义，覆盖新增粗层壁面区域。
- 外边界初版沿用有限距离接地电势和粒子吸收条件。

## 3. 当前阻塞点与必须处理的耦合

### 3.1 明确的硬限制

| 位置（相对仓库根目录） | 当前行为 | 初版处理 |
|---|---|---|
| `Source/Insert/Fields/HallWallCharge.cpp::prepare` | `max_level == 0` | 实现跨层壁电荷语义后解除 |
| `Source/Insert/Fields/HallElectrostaticMaterial.cpp::prepare` | `max_level == 0` | 材料、阳极掩码和求解链路一起适配 |
| `Source/ablastr/fields/PoissonSolver.H::computePhi` | 材料 cache 要求 `finest_level == 0` | 直接实现逐层 linop/MLMG 缓存，按当前求解层复用对应对象 |
| `Source/Insert/Collisions/CoupledBackgroundMCCCollision.cpp::doCollisions` | `flvl == 0` | 属于主算例必改项，不能仅删除断言 |
| `Source/WarpX.cpp` 外部 grid 场读取检查 | 文件 E/B grid 场要求 `max_level == 0` | 主算例确用 `warpx.B_ext_grid_init_style = read_from_file`，必须适配 |
| `Source/Insert/Fields/ECDIChargeFilter.cpp` | 单层、单 rank、单 box 限制 | 初版保持禁用及明确报错 |
| `Source/Insert/Collisions/IonizationSourceTable.cpp` | 单层、单 rank 限制 | 相关可选模式初版保持禁用及明确报错 |

材料求解直接采用 per-level cache，连同 `LabFrameExplicitES`、`ElectrostaticSolver`
和 ablastr 的接口一起改。每层保存一组配对的 `MLNodeLaplacian` 与 `MLMG`，
首次求解时创建，后续时间步复用。

- [ ] 将单组 `m_material_linop/m_material_mlmg` 改为按层管理的缓存，
  在 `computePhi` 的层循环中访问当前 `lev` 的对象。缓存接口约束改为层数匹配、
  每层 linop/MLMG 配对及布局一致。
- [ ] 首次创建时使用该层 Geometry、BoxArray、DistributionMapping、介电系数和阳极 mask，
  配置算子、边界类型及多重网格参数，再构造引用该算子的 MLMG。
- [ ] 每次求解更新 RHS、phi 初值、指定电势以及可变的求解控制参数。
  粗到细插值后恢复该层阳极电势；固定电势数值的更新与算子缓存分开处理。
- [ ] 随时间变化的粒子电荷与壁电荷仅更新 RHS，各层算子及 MLMG 持续复用。
- [ ] 程序结束释放缓存时，每层先销毁 MLMG，再销毁其 linop。
  缓存结构应保证引用的算子对象在 MLMG 使用期间保持有效。

### 3.2 材料与体阳极：插值后必须恢复约束

`HallElectrostaticMaterial::initialize` 写入固定阳极电势，`prepare` 中重新写
`phi` 的语句当前被注释，只清零掩码内 `rho`。粗层求解后的
`interpolatePhiBetweenLevels` 会重写细层初值，可能覆盖细层掩码内的固定电势。

- [ ] 在每层实际求解前、粗到细插值之后，按该层 mask 恢复阳极指定电势。
  仅在整个 `computePhi` 调用之前恢复一次不足以处理后续插值。
- [ ] 适配两层介电系数、阳极 mask、薄结构离散与局部细区边界条件。
- [ ] 完善 collocated 布局下两层 phi ghost 填充、E 计算与粒子 gather，
  处理物理边界和粗细界面的数据衔接。

### 3.3 壁电荷：统一 RHS 同步与两层壁面归属

多层壁电荷采用统一电荷同步路径：各层持久壁电荷以原始沉积状态保存，
每次求解将其副本加入本步粒子 `rho_fp`，共同经过一次 `SyncRho`。
复用其中的 fine-to-coarse 限制、共享节点处理和 MPI 求和；
电势继续由 `interpolatePhiBetweenLevels` 从粗层插值到细层。

单层兼容路径保留原有装配顺序：粒子 rho 完成同步、滤波及 PEC 处理后，
才加入独立同步的壁电荷副本，不将壁电荷送入粒子 PEC 镜像处理。

```text
各层粒子本步沉积 ──────────┐
                          ├→ 本步 rho_fp → SyncRho → 电荷边界处理
各层持久壁电荷的原始贡献副本 ┘                         → 材料约束 → 泊松求解
```

- [ ] 在首次撞壁处理前初始化各层持久壁电荷。将 `ApplyWallPolicy` 中的
  `get(0)`、`MakeWallChargeGrid(Geom(0))` 改为匹配实际沉积层的访问，
  同时匹配粒子 tile、box、DistributionMapping 和该层的体积归一化。
- [ ] 持久存储保留各 FAB 的原始沉积贡献，包括共享节点和 ghost 上尚未求和的贡献。
  调整 `HallWallCharge::addTo/synchronize` 接口，使加入本步 RHS 前的数据状态
  与粒子原始沉积一致；现有 `SumBoundary` 和清 ghost 操作应从持久存储路径移出。
  `SyncRho` 仅作用于本步工作场，其同步、滤波和粗层限制结果均不回写持久存储。
- [ ] 每次求解重新装配粒子电荷与累计壁电荷，核对回调的执行顺序和数据语义。
  粗层 RHS 同时包含粗层真实撞壁贡献和细层壁电荷的限制贡献，
  每项在本步装配中只加入一次。
- [ ] 粗层独占区和细化区均执行完整壁面策略，每次撞壁只生成一份净留存电荷
  （入射电荷减去出射产物电荷），并只计一次壁电流/事件统计。
  持久电荷按撞击位置的最细覆盖层归属；对粒子在推进中跨越细化边界后撞壁的情形，
  核查墙处理与 `Redistribute` 的时序，必要时将撞击沉积路由到对应层和 box。
  产物粒子按产生位置归入对应层级，保证粗细两侧采用相同的材料模型。
- [ ] 核查撞击沉积模板跨越 box、MPI 或细化边界时的归属与传递，保留完整净电荷。
  粗层中由细层限制得到的电荷是求解副本，独立于粗层真实撞壁的持久贡献。
  总壁电荷诊断使用同步后的临时副本、节点积分权重和粗细覆盖掩码，
  排除共享节点与限制副本的重复计数，并与累计撞击净电荷核对。
- [ ] 明确统一装配后普通 rho 滤波、PEC 电荷边界处理及 ECDI 对壁电荷的作用。
  初版保持 `warpx.use_filter = 0` 并禁用 ECDI，处理 PEC 附近的壁电荷边界语义；
  后续启用滤波时明确其适用源项。

### 3.4 外加磁场、MCC、注入与诊断

- **文件磁场**：`Source/Initialization/WarpXInitData.cpp::ReadExternalFieldFromFile`
  内部写死 `Geom(0)`。需传入层号或 Geometry，检查所有调用点和所需 fine/coarse patch
  场的初始化。读取器使用目标 nodal/cell 布局及 ghost 范围，必须核查
  `../../data/B3D_dielectric.bp5` 的实际覆盖域及域外处理，
  确保文件数据覆盖扩展域和所需 ghost 范围后适配读取检查。
- **MCC**：`doCollisions` 的 `base_lev = 0`、`ncell/inv_gap` 和背景密度最大值
  `m_background_density_fabs[0].max(0)` 均需处理。
  碰撞率上界要覆盖被处理层的密度（逐层上界或全层上界），否则可能漏采碰撞。
  同时修正 ionization 路径 `LowerCorner(box, 0, ...)` 和
  `Source/Insert/Core/WarpXFunc.h` 的 `ParticlesAt(0, pti)` / `LowerCorner(...,0,...)`。
- **背景密度**：各层只沉积本层中性粒子，供同层 MCC 使用，
  采用层内密度和 bins 完成电离消耗及产物创建。粗层在细区覆盖下的背景密度允许为零，
  本轮省略背景密度的跨层限制和细化界面补齐，保留同层跨 box/MPI 同步及空 tile 处理。
  粗层仍按本层密度执行 MCC，电荷与壁电荷继续通过统一 `SyncRho` 做跨层同步。
  `MCC_DENSITY_AVERAGE_USE`
  写死层 0 且读单 FAB，初版遇 MR 应明确报错。CALC 的文件输出布局需按层适配。
- **注入**：`WarpXParticleContainer::AddNParticles` 忽略 lev 参数，先写层 0
  的临时 tile，末尾调用 `Redistribute()`。沿用该路径，
  按粒子位置完成细区、粗区及跨 rank 的归属分配。
- **诊断**：`ZmaxRadialExitStats` 的 `Geom(0).ProbHi(2)` 表示整个计算域出口。
  明确需要的是扩展域出口还是原通道截面；后者需要独立的截面穿越诊断，
  径向 bin 间距应显式约定。
  `InsertInjection.cpp` 的 2D 配置及 `InsertBoundaryPhi.cpp::VoltageAdjustment`
  是另行适配项，不算入首版 3D 完成范围。

### 3.5 初始化与时间推进

- 初始建层路径为 `WarpX::InitFromScratch` → `AmrCore::InitFromScratch`
  → `MakeNewLevelFromScratch`。
- 初版显式关闭动态负载均衡（`algo.load_balance_intervals = 0`）。
  按初始化确定的网格布局创建材料、持久壁电荷、求解器和背景 bins。
- 协调当前 species `ndt` 下跨层粒子移动、MCC 和场更新时序。

## 4. 分阶段实施

### 阶段 0：配置与几何准备

- [ ] 记录 HEAD、编译选项、Insert 宏、输入和磁场文件版本。
- [x] 确定扩展域端点、粗层网格及静态细化区，保持装置尺寸和注入参数的物理含义。
- [x] 配置两层网格，关闭 AMR subcycling 和动态负载均衡。

#### 阶段 0 实施记录（2026-10-10）

新增配置入口 [`Script/PHT/3d_hall_dielectric_slot_amr`](Script/PHT/3d_hall_dielectric_slot_amr)，
通过 `FILE` 继承单层主算例及其 MCC 输入，随后覆盖几何和静态 MR 参数。
相对路径以仓库内 `<build>/bin` 为工作目录解析，与原主算例一致。
阶段 2 完成后可从该目录传入 `../../Script/PHT/3d_hall_dielectric_slot_amr`。
当前阶段仅准备配置；材料、壁电荷、外场及 MCC 的单层限制仍然存在。

**代码与构建基线**

- 分支 `Amr-support`，HEAD `01f1b2bfca2566e56164b8c008636a724885f946`。
  主输入最后修改提交同为该 HEAD。实施开始时已有计划文件的未提交修改和
  暂存的 `不锈钢SEE模型_Furman-Pivi.md`，本次保留这些改动。
- CLion Release 配置来自 `.idea/workspace.xml`：
  `-DBUILD_TESTING:BOOL=OFF -DWarpX_COMPUTE:STRING=CUDA`
  `-DCMAKE_CUDA_COMPILER:FILEPATH=/usr/local/cuda-13.2/bin/nvcc -DWarpX_FFT:BOOL=ON`，
  构建参数 `-j 8`，匹配目标为 `app_3d`。
- 下表读取现有 CMake cache，表示配置记录，未声称本轮重建了可执行文件。

| 选项 | `cmake-build-release` | `build-codex` |
|---|---|---|
| `CMAKE_BUILD_TYPE` / `WarpX_DIMS` | Release / 3 | Release / 3 |
| `WarpX_COMPUTE` | CUDA | OMP |
| `WarpX_MPI` / `WarpX_OPENPMD` | ON / ON | ON / ON |
| `WarpX_EB` / `WarpX_QED` | ON / ON | ON / ON |
| `WarpX_FFT` / `WarpX_PYTHON` | ON / OFF | OFF / ON |
| `WarpX_PRECISION` / `WarpX_PARTICLE_PRECISION` | DOUBLE / DOUBLE | DOUBLE / DOUBLE |
| `CMAKE_CXX_COMPILER` | `/usr/bin/c++` | `/usr/bin/c++` |

Release 的 CUDA 架构为 `native`，C++/CUDA Release flags 均为 `-O3 -DNDEBUG`；
`WarpX_FASTMATH=OFF`、`AMReX_CUDA_FASTMATH=ON`、`WarpX_IPO=OFF`。
依赖配置为 AMReX `66028f892b7d0e0902ffe821e10fcc0915ab1cda`、
openPMD `0.17.1`、PICSAR `26.05`。

Insert 启用宏为 `NUMP`、`BENCHMARK_2D`、`MCC_DENSITY`、`MCC_DELETE`，
3D 条件下另启用 `HALL3D` 和 `WARPX_USE_HALL_ELECTROSTATIC_MATERIALS`。
`FIELDSIZE`、`COLLISION_RECORD`、`MCC_DENSITY_AVERAGE_CALC/USE`、
`MCC_EXCITATION`、`IONIZATION_SOURCE_RECORD/INJECT`、`MCC_ION_CACHED_RHO`、
`HALL3D_INIT` 均未启用。

以下为修改前文件内容的 SHA-256，便于后续区分输入或构建配置漂移：

| 文件 | SHA-256 |
|---|---|
| `Script/PHT/3d_hall_dielectric_slot` | `a28a0501973873f5eb7c993631ad524f29a9e71e24106985a92d20300e381853` |
| `Script/3d_hall_mcccoll` | `4884b2f43cd5b26a5b684db0eb6ded7604436ea4a0f7543f34b1213b5a93fc60` |
| `Source/Insert/Config/HallElectrostaticConfig.H` | `93be8137618824a1bc6b9c101d2f321c8632b025699fad208ce41db186c2bcd8` |
| `Source/Insert/Config/WarpXFunctionConfig.h` | `54c66b1e62ca357bb95f1192c9e82d9c95b5b6d54079c665b32b032da6109e8c` |
| `Source/Insert/Config/WarpXOptimizationConfig.h` | `d325dd69b814539312009c7d7e969d784dab3ae851bf5c493e44416e5c127383` |
| `Source/Insert/Config/WarpXSimulationConfig.h` | `f42008ccfaf68423765b4be7fa337e704d6b732bffa2d6ad2e21cad8ff2c01bb` |
| `cmake-build-release/CMakeCache.txt` | `3e62ad4c1616d583c8d884ec22df38ec231f1c97dbde28448813359bc732f125` |
| `build-codex/CMakeCache.txt` | `fc3d405da67462a50c8e175734eba6354e63dfd8bddaed94be43185920b315f2` |

**外部数据待补项**

磁场引用为 `../../data/B3D_dielectric.bp5`，从约定工作目录解析到
`<repo>/data/B3D_dielectric.bp5`；初始中性粒子引用为
`../../data/hall_start_256_600K_slot.bp5`。本工作区没有这两个数据集，
因此不能记录其实际内容版本、openPMD 元数据或确认磁场对扩展域及 ghost 的覆盖。
第一个复选框保留未完成：数据可用后记录 BP5 目录内各文件的内容哈希、
iteration、轴顺序、网格 spacing/offset、单位及各分量 extent/position；
不得仅凭文件名或沿用原路径认定覆盖已经满足。

**几何决策与检查**

| 项目 | 阶段 0 配置 |
|---|---|
| level 0 物理域（mm） | `[-5,15] × [-5,15] × [-8.75,16.25]` |
| level 0 网格 / 间距 | `256×256×320` / `78.125 μm` |
| level 1 tagging 区（mm） | `[0,10] × [0,10] × [-2.5,10]` |
| 细化率 / level 1 间距 | `2×2×2` / `39.0625 μm` |
| 粗层 tagging 半开索引 | `[64,192) × [64,192) × [80,240)` |
| 对应细层半开索引 | `[128,384) × [128,384) × [160,480)` |
| blocking factor / max grid size / tag buffer | `32` / `128` / `0` |
| 静态控制 | `warpx.regrid_int=-1`、`warpx.do_subcycling=0`、`algo.load_balance_intervals=0` |

各层 Geometry 仍共享扩展物理域；表中 level 1 范围是局部细化覆盖目标。
按 `WarpX::ErrorEst` 的 cell-center 判据核对了 tagging 索引，细层各边界与
blocking factor 对齐，目标细区包含 `256×256×320` 个细胞。
实际 BoxArray 的分块及 DistributionMapping 需在后续可初始化的构建中确认。
设置 `n_error_buf=0` 避免默认 tag 膨胀；源区缓冲由显式选定的细区提供。

装置中心仍为 `(5,5) mm`；阳极半径 `2.1–3.1 mm`、z 范围 `[-1.25,0] mm`；
初始等离子体 z 为 `0.2–0.8 mm`，环缝在 z=0，阴极源半径 `3.2–4 mm`、
z 为 `8.4–9.2 mm`，均位于细区内部。阴极源距细区 z 上边界至少 `0.8 mm`，
初始通道等离子体距该边界至少 `9.2 mm`；演化后主要电离区的范围留待阶段 2 核查。
上外边界位于 `16.25 mm`，下外边界位于 `-8.75 mm`。
`lx/ly/lz`、粒子权重、SEE、流量、注入、碰撞以及 species `ndt=1/10/100` 均继承主输入。
陶瓷平面模型延续至新增粗层壁面区域，六面保持接地电势及粒子吸收。
使用 `Geom(0).ProbHi(2)` 的出口统计将对应 `z=16.25 mm`，不是旧的 `z=10 mm` 截面。

验证：使用 Python 递归解析两个 `FILE` include，核对覆盖项只涉及域及网格控制，
检查了粗细间距、细化边界对齐、解析源区包含关系及关闭 subcycling/负载均衡的参数，全部通过。
未运行完整 Hall 算例或声称验证实际 BoxArray；外部数据缺失且后续阶段的单层断言仍在。
本轮没有 C++ 修改，也没有编译；当前会话未提供 CLion MCP 构建工具。

### 阶段 1：材料、阳极、壁电荷和外场

- [ ] 完成逐层材料 linop/MLMG 缓存及全链路接口调整，实现首次创建与后续复用。
- [ ] 完成两层材料准备和插值后阳极电势恢复，将每层材料数据接入对应缓存算子。
- [ ] 完成持久原始壁电荷存储及 RHS 装配接口，将两层壁电荷接入统一 `SyncRho`。
- [ ] 完成粗层、细层和跨细化边界撞壁的电荷归属、产物迁移及统计去重。
- [ ] 将文件磁场读取改为按目标层 Geometry 采样，处理扩展域及 ghost 数据覆盖。

#### 阶段 1 实现记录（2026-10-10，待构建及运行验收）

阶段 1 的功能代码已写入工作区；按用户要求不提供测试相关实现，由用户另行测试。
上面的验收项暂不勾选，
因为本会话没有 CLion MCP 构建工具，未能按仓库要求重建 `app_3d`，
也未运行修改后的 WarpX。旧可执行文件不能验证本次 C++ 修改。
生产磁场文件仍然缺失，因此其实际扩展域/ghost 覆盖尚未验证。

- **材料求解缓存**：新增 `Source/Insert/Fields/MaterialPoissonCache.H`，
  每层拥有独立的 `MLNodeLaplacian`/`MLMG`，结构成员顺序保证先释放 MLMG。
  `LabFrameExplicitES` → `ElectrostaticSolver` → ablastr 接口已贯通，
  复用前核对层数、布局、Geometry、材料/mask 对象、网格布局及边界类型；
  RHS、初值、求解容差和迭代控制逐次更新。
- **材料与阳极**：移除单层限制，保留逐层介电常数正定性、阳极薄结构检查。
  每层求解之前恢复固定电势，位置位于粗到细插值之后，避免细层 mask 中的
  电势被粗层初值覆盖。缓存只允许静态网格，显式拒绝负载均衡、regrid 和 moving window。
- **phi/E ghost**：材料路径逐层填充 collocated 电势 ghost；物理 Dirichlet 面
  使用二次外推，同层通信优先于粗层插值。按用户范围仅处理第一类物理边界，
  已移除材料 ghost 路径的 Neumann 和周期专门处理。
  E ghost 插值按各分量实际 nodal/cell staggering 处理，支持当前主算例的 staggered
  及原计划中的 collocated。静电 aux 路径直接复制各层 fine-patch 场，
  不使用电磁求解器的 coarse-patch 差值修正。
- **壁电荷装配**：持久场保留原始 FAB/ghost 贡献，删除持久场上的 `SumBoundary`
  和 ghost 清零。多层在粒子沉积及 `afterdeposition` 回调后加入本步 rho，
  然后统一 `SyncRho` → PEC 电荷处理 → 材料掩码 → 求解，拒绝普通 rho 滤波和 ECDI。
  单层在粒子 rho 处理后，将独立 `SumBoundary` 并清 ghost 的壁电荷副本加入 RHS，
  保留原有源项处理顺序及滤波适用范围。两条路径均不回写持久原始数据。
- **撞击路由**：每次事件仅记录一份净留存电荷，由 GPU kernel 直接追加到临时
  AMReX 粒子容器，已移除主机撞击队列和对应 OpenMP 并行/临界区处理。容器将交点路由到
  最细覆盖层及所属 MPI box，再使用该层体积归一化沉积。
  净电荷记录与事件/电流统计不重复执行材料策略。
  产物沿用物种创建路径并由 `HandleParticlesAtBoundaries` 中随后的全层
  `Redistribute` 按最终位置分配；壁电荷在首次撞壁之前可以按现有 rho 布局初始化。
- **诊断**：`insert.wall_charge_diagnostic=1` 输出 `wall_charge.dat`。
  包含累计撞击净电荷、持久原始存储总量及各层分量，以及在独立临时副本上
  同步后的复合节点积分和差值。复合积分排除粗层被细网格覆盖的单元，
  用单元角点积分权重避免共享节点重复计数；不对持久数据做同步。
  粗细界面附近复合积分与原始沉积账目的差异会显式输出，仍需运行时定量检查，
  不把原始账目守恒当作复合求解链路已经验收。
- **文件外场**：`ReadExternalFieldFromFile` 改为显式接收 `lev`，
  所有 grid/particle E/B 和 hybrid 外部矢势调用点已更新。
  使用目标层 Geometry、各分量布局及 ghost 范围采样；仅 MR 完整覆盖路径处理
  openPMD `gridUnitSI` 和 component `position`。
  3D MR 要求文件覆盖全部请求的插值模板；不足时明确报错，不在域外默默置零。
  单层路径保留原有坐标解释及域外零值策略。按用户要求不增加场幅值 `unitSI` 换算，
  避免改变既有磁场文件的解释。静电外场仅持有 `*_fp_external`，
  `AddExternalFields(lev)` 将其逐层加回 fine-patch 场，aux 从该场更新。

按用户要求，已移除本次新增的回归输入、分析脚本、合成外场生成脚本、
CTest 注册及临时验证环境。阶段 1 功能实现和壁电荷运行诊断保留，
后续测试由用户自行安排。本次 C++ 修改尚未经编译或运行验收；
生产输入仍需阶段 2 的 MCC 多层适配及阶段 0 所缺的数据集。

### 阶段 2：MCC、粒子跨层与主算例接入（暂缓）

- [ ] 完成按层独立背景密度的 MCC、逐层碰撞率上界、同层消耗/产物和辅助复制路径。
- [ ] 衔接注入、双向跨层迁移、gather/deposit、外边界吸收及墙面产物粒子的处理时序。
- [ ] 在 Hall 主算例中引入扩展域和静态 fine patch，沿用扩展域内的墙/材料定义，
  使粗层新增壁面区域与细区壁面同时参与粒子交互和电荷反馈。
- [ ] 对不支持的可选模式在初始化时明确报错；保持 SEE、注入和碰撞物理参数一致。

待单层 Robin 接入和单层求解验收后，再决定是否恢复本阶段。

### 后续可选范围

ECDI 多层滤波、IonizationSourceTable、AVERAGE_USE/CALC 文件约定、
2D 及其他维度作为后续独立任务。
首版实现范围为上述 3D 直接背景密度模式。

## 5. 实施约定

- 按求解精度和资源预算调整扩展域与细区布局，必要时评估复合多层泊松或粒子重采样方案。
- 移除或重命名用户参数时，补充相应的兼容性检查。
- 遵循仓库对自动生成文件的维护规则。

## 6. 单层兼容性与 Robin 接入计划（2026-10-13）

本节记录接入方案，不表示 Robin 已加入生产求解路径。现有单层主输入保持
`amr.max_level=0`；两层输入仅保留为阶段 0–1 的开发入口，MCC 多层限制继续保留。

### 6.1 本次单层兼容处理

- 材料缓存只有 level 0 时仍使用 `amrex::MLNodeLaplacian` 与持久 `MLMG`，
  方程、介电系数、网格布局和时间外推公式不变。
- 每层求解前恢复指定阳极电势，防止初值修改破坏固定电势。
- 单层壁电荷保留在粒子同步、滤波、PEC 处理之后加入；原始存储在临时副本上
  同步，避免重复累加。同一接口的多层分支才参与统一 `SyncRho`。
- 外场坐标元数据换算及严格覆盖检查只作用于 MR；单层使用既有文件约定。
- 保留已确认的 Dirichlet phi ghost 外推及 E ghost 处理；它们是此前明确的边界修正，
  不宣称相对改造前逐位相同。阶段 1 的壁电荷交点路由也需核对 MPI/多 Box 行为。
- 本次完成源码路径核对和差异检查；未通过 CLion 构建修改后的 `app_3d`，
  未复跑单层数值对照。因此不能把上述隔离措施当作运行验收通过。

### 6.2 依赖与导入范围

- AMReX 已更新至 `47ceffc47cea4b157fa20aed944689557ca72134`（版本标签 26.10）。
  已核对 `supportsAlgMG`、`ParallelForStrided`、`multicolor_offset` 及节点算子虚接口；
  相关 MLMG/节点基类与参考包验证基线 `334d08c6...` 无差异。
- 将 `cht_node_robin/src` 四个核心文件及来源/版权记录导入
  `Source/Insert/Fields/NodeRobin/`，仅编入 `lib_3d`，复用宿主 AMReX/CUDA 配置。
  不复制第二份 AMReX，不直接套用参考包的独立 CMake 最低版本要求。
- 首版仅支持单层、完整域、3D Cartesian、非周期和静态网格，保留多个 Box。
  单层 AMR 与多层 MG 是不同概念；关闭该分支的半粗化，保留各方向二倍 MG 粗化。

### 6.3 求解入口与缓存

- 在 Insert 建立专用单层 Robin 求解入口及缓存，由 `LabFrameExplicitES` 显式选择。
  默认继续使用当前 Dirichlet 材料求解，避免改变既有输入含义。
- Robin 缓存持有 `cht::MLNodeRobinLaplacian` 和 `MLMG`，后者先析构；
  首次设置 `epsilon_r`、单极中心、远场电势及三参数 `setOversetMask(0, mask, voltage)`。
  不经基类两参数 mask 接口绕过新算子的电极管理，也不由旧 `setDomainBC` 覆盖 Robin 配置。
- 当前主算例阳极为域内统一电压，可复用已有几何生成的 mask；首版不承诺支持
  任意空间变化的阳极电势。粗化上限按实际 mask 核对，不能仅沿用固定值 4。
- 每步复用粒子/壁电荷装配、时间外推及阳极约束；物理 RHS 仅缩放一次为
  `-rho/epsilon0`，Robin 源项由算子钩子加入。诊断使用 `compPhysicalResidual`，
  避免直接用三维 `MLMG::compResidual` 漏掉非齐次边界源项。

### 6.4 外边界与电场链路

- 单极条件为 `dphi/dn + q*(phi-phi_infinity)=0`，
  `q=n·(x-center)/|x-center|²`，法向朝外；中心必须位于物理域内部。
- 新模式必须贯通边界识别、`setPhiBC`、电荷边界处理和 E 边界处理，
  避免 Poisson 为 Robin 而外侧仍执行 PEC 电荷镜像或切向 E 清零。
  内部阳极保留 Dirichlet；外侧粒子吸收仍由粒子边界独立控制。
- 首先适配当前 staggered 主输入，复用内部 `computeE` 差分并处理 E ghost/gather。
  若支持 collocated，按 Robin 法向导数补物理 ghost，不能直接把算子内部
  齐次镜像或现有 Dirichlet 二次外推当作已验证的 Robin 取场模板。
- 当前陶瓷表达式延伸到部分外边界；必须评估单极近似与材料外延、域大小及
  展开中心的关系。球体真空验证不代表真实羽流工况已获得同样精度。

### 6.5 验收顺序

1. 通过 CLion 构建匹配的 `app_3d`；先验收当前单层 Dirichlet 路径，检查
   电势、电场、阳极电压、壁电荷账目、迭代次数及多步时间外推。
2. 导入后复核参考球体和阳极问题，包括物理残差、边界 E 和重复求解缓存。
3. 在 WarpX 单层装配中比较单/多 Box，并独立验收 MPI。
4. 最后检查真实工况的边界模型与计算域独立性，再决定是否恢复 AMR 阶段 2。

本次不新增测试工程或注册；运行验收需使用已构建的对应代码及可用的外部数据。
