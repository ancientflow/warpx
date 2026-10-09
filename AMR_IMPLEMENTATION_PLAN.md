# WarpX 原生 AMR 改造实施计划

> 目标：用 WarpX 原生 AMR（2 层静态网格）替代自研双层远场泊松（`HallFarFieldPoisson`），
> 粗层（level 0）覆盖扩展域并携带粒子，解决远场电势自洽问题。
>
> 本文基于 2026-10 对 `Source/Insert/` 及 fork 主代码改动的盘点结论。

## 1. 现状评估摘要

### 1.1 有利条件

- **fork 对主代码污染很轻**：相对上游基线（`36c582069`，BLAST-WarpX development）
  仅改 21 个文件、+735/−22 行，且集中在静电求解器一条线
  （`LabFrameExplicitES.cpp`、`PoissonSolver.H`、`ElectrostaticSolver.*`）。
- **上游原生支持 labframe 静电 + MR**，且有 CI 回归测试：
  `Examples/Tests/electrostatic_sphere/inputs_test_3d_electrostatic_sphere_lab_frame_mr_emass_10`
  （`amr.max_level = 1` + `warpx.do_electrostatic = labframe`）。
- **多层骨架现成可用（零改动）**：
  多层电荷沉积 + average-down（`WarpXParticleContainer.cpp:1825-1898`）、
  `SyncRho` 粗细同步（`Source/Parallelization/WarpXComm.cpp:1354-1420`）、
  逐层 push/gather（`WarpXEvolve.cpp:1408-1428`）、逐层边界/scrape、
  逐层泊松求解 + `interpolatePhiBetweenLevels`
  （粗层先解，插值为细层 Dirichlet 边界 + 初值，`PoissonSolver.H:532-540`）
  ——这正是自研远场方案的官方等价物。
- **Insert 各模块的层级意识比预期好**：壁面交互的粒子遍历已逐层
  （`AnalyticBoundaryInteraction.cpp:668`），背景密度 CALC 主路径已逐层，
  材料介电初始化/prepare 已逐层（`HallElectrostaticMaterial.cpp:204-259`），
  壁电荷存储/synchronize 已 per-level（`HallWallCharge.cpp:59-109`）。
- **网格静态**：WarpX 的 Evolve 从不调 `regrid`，层级由初始 tagging 决定、
  运行期不变。Insert 各模块的一次性缓存（`m_defined`）因此不会因 regrid 失效。

### 1.2 关键设计决策：几何语义反转

AMR 下 level 0 必须是**覆盖扩展域的粗层**：

- `geometry.prob_lo/hi` 从细通道域（10×10×12.5 mm）改为扩展域
  （20×20×25 mm，与现远场粗层范围一致，`ff_xc±lx` 等）；
- level 0 网格 256×256×320、间距为现细层 2 倍；
- 细通道区通过初始 tagging（`warpx.refine_patch` parser 或 `fine_tag_lo/hi`）
  生成 level 1，分辨率与现细层一致；
- **fork 中所有 `Geom(0)`/`rho_fp[0]`/`phi_fp[0]` 写死处的语义随之从
  "细通道网格"变为"扩展粗网格"**——这是改造中最容易出错的全局性变化。

### 1.3 已知限制（与本改造相关的上游约束）

- subcycling 不支持静电（`WarpX.cpp:770-779`）：粗细层同 dt，无影响。
- `MakeNewLevelFromCoarse` 是 abort 桩（`WarpX.cpp:2326-2330`）：
  动态升层不可用，但初始 tagging 路径可用（上游 MR 测试为证）；实施前需在
  阶段 0 验证确认。
- 外部场文件读取要求 `max_level==0`（`WarpX.cpp:745`）——**本算例外加磁场
  从文件读取，这是一个必须解决的硬冲突**（见第 5 节风险 R2）。

## 2. 开启 `max_level=1` 后的故障清单（改造的完整工作集）

### 2.1 立即崩溃（硬断言）

| # | 位置 | 内容 |
|---|---|---|
| A1 | `HallWallCharge.cpp:46-48` | `max_level==0` 断言（材料开关打开即触发） |
| A2 | `ECDIChargeFilter.cpp:360-365` | 单层 + 单 MPI rank + 单 box 断言（仅 `insert.ecdi_control.enabled` 时） |
| A3 | `IonizationSourceTable.cpp:45-51` | 单层 + 单 rank 断言（仅平均电离源表模式） |

### 2.2 不报错但结果错误（隐性单层假设）

| # | 位置 | 内容 |
|---|---|---|
| B1 | `LabFrameExplicitES.cpp:256-263, 298, 311` | 远场耦合点全部只吃 `rho_fp[0]`/`phi_fp[0]`（AMR 下语义倒置） |
| B2 | `PoissonSolver.H:347-358` | 材料 linop cache 用 level-0 算子解所有层（`m_material_linop/m_material_mlmg`，`LabFrameExplicitES.H:78-79`） |
| B3 | `AnalyticBoundaryInteraction.cpp:815-817` | 壁电荷只沉积到 `get(0)`、`MakeWallChargeGrid(Geom(0))`；lev>0 tile 会错位/越界（:710） |
| B4 | `CoupledBackgroundMCCCollision.cpp:256-258, 346-348, 392-395, 620` | `Geom(0)`/`ncell` 写死、`flvl==0` 断言、`LowerCorner(...,0,...)` |
| B5 | `WarpXFunc.h:83, 91` | `ParticlesAt(0, pti)`、`LowerCorner(box, 0, ...)` |
| B6 | `InsertInjection.cpp` / `HallInjectionSource.cpp:204` | 全部 `AddNParticles(0, ...)`（原生实现忽略 lev，先放 level 0 tile(0,0,0) 再 Redistribute） |
| B7 | `BackgroundCoupledDensity.cpp:167, 244-250` | AVERAGE_USE 模式写死 `flvl=0` |
| B8 | `ZmaxRadialExitStats.cpp:60-73`、`InsertInjection.cpp:60` | 诊断/配置用 `Geom(0)` 的 dx/zmax/prob 做分箱 |
| B9 | `LabFrameExplicitES.cpp` 的 `ExtendPhysicalGhosts`/`fillElectricFieldGhosts` | 只处理 `phi_fp[0]`/`Efield_fp[0]` |
| B10 | `WarpXFieldBoundaries.cpp:56-64` | 远场模式下跳过 PEC E 边界的逻辑，需随远场模块退役重新审视 |

## 3. 总体方案

```
level 0（粗层）：扩展域 20×20×25 mm，256×256×320，间距 2×现细层
                 携带粒子（羽流），六面接地 Dirichlet
level 1（细层）：tagging 出的通道区，分辨率与现细层一致
                 解析墙/注入/MCC/材料/壁电荷均在此层为主
求解：原生逐层泊松（粗层先解 → interpolatePhiBetweenLevels 给细层
      Dirichlet 边界 + 初值），完全替代 HallFarFieldPoisson 的职责
```

退役对象：`HallFarFieldPoisson.*` 整体、`DefineExternalDirichletBCs`
（`PoissonBoundaryHandler.cpp:155-161`）、`WARPX_USE_HALL_FAR_FIELD` 相关分支。
粗层电荷现在由粗层粒子自洽沉积产生，`restrictCharge`/frozen-boundary
插值/时序外推历史（`HallFarFieldPoisson.cpp:432-650`）全部不再需要。

## 4. 分阶段计划

### 阶段 0：基线验证（1–2 天）

- [ ] 在本 fork 上运行上游回归测试
  `electrostatic_sphere_lab_frame_mr_emass_10`（CPU 构建即可），
  确认多层沉积/SyncRho/逐层求解/插值链路未被 fork 破坏。
  （忽略 checksum 差异，只看物理合理性与不崩溃。）
- [ ] 确认初始 tagging 生成 level 1 的路径在 `MakeNewLevelFromCoarse`
  abort 桩存在的情况下可用（上游测试通过即说明可用；记录实际调用路径）。
- [ ] 极简对比算例：均匀等离子体块 + 接地边界，分别用
  ① 现有双层远场、② 2 层原生 AMR 求解，对比 phi 场，
  确认两种"粗层边界"语义物理等价后再动主算例。

**验收**：① ② 的 phi 全场相对偏差在求解容差量级。

### 阶段 1：最小可用版（约 1 周）

目标：`3d_hall_dielectric_slot_farfield` 在 2 层 AMR 下跑通，
10 万步内 phi/rho 演化与现方案一致。

- [ ] **几何切换**：`geometry.prob_*` 改扩展域；新增 tagging 配置圈出
  细通道区（x/y 全域中心 ±5 mm、z 从 z_lo 到通道出口稍下游）；
  核对 level 1 box 布局与 MPI 划分（沿用现 2×2 半域约束的对应要求）。
- [ ] **外加磁场冲突（R2）**：`WarpX.cpp:745` 的 `max_level==0` 限制——
  将 B 场文件读取/插值扩展到逐层（先按最简方案：两层都从同一文件
  按各自 Geom 采样；粗层分辨率减半对通道区无影响，通道区在 level 1）。
- [ ] **退役远场**：删除/旁路 B1 全部耦合点与 `DefineExternalDirichletBCs`；
  细层边界改由原生 `setPhiBC`（`ElectrostaticSolver.cpp:98`）
  + `interpolatePhiBetweenLevels` 提供。
  保留 `HallFarFieldPoisson` 源码与编译开关以便回退对比。
- [ ] **B2 材料 linop cache 逐层化**：`m_material_linop/m_material_mlmg`
  改为 per-level vector（或先走非 cache 分支 `PoissonSolver.H:360-368`，
  性能换正确性，后续再优化）。
- [ ] **A1 + B3 壁电荷逐层化**：放开断言；`get(lev)`、
  `MakeWallChargeGrid(Geom(lev))`；`ApplyWallPolicy` 的 lev 贯穿到底；
  定义粗细重叠区壁电荷归属（壁面几何都在 level 1 区域内，
  预期粗层无壁电荷，验证即可）。
- [ ] **B4/B5 MCC 修正**：`Geom(lev)`、`ncell` 按层、放开 `flvl==0` 断言、
  修 `LowerCorner(...,0,...)` 与 `ParticlesAt(0,pti)`；
  背景中性气 bins 在粗层的分布语义按"同一物理数密度、按层网格分箱"实现。
- [ ] **B9**：`ExtendPhysicalGhosts`/`fillElectricFieldGhosts` 改到最细层
  （或逐层），与原生 aux 填充流程对齐。
- [ ] **B8 等诊断**：`Geom(0)` → `Geom(finest_level)`，逐处几行。
- [ ] **ECDI / 电离源表 / AVERAGE_USE（A2/A3/B7）**：本阶段一律在初始化时
  检测 `max_level>0` 给出明确 `amrex::Abort` 提示（不允许静默错误），
  功能本身推迟到阶段 3。
- [ ] **B6 注入验证**：原生 `AddNParticles` 放 level 0 后经 Redistribute
  按位置分配层级——用一个小算例验证生于细区的粒子确实落到 level 1；
  若 WarpX 的 per-level Redistribute 不跨层迁移，则注入时按位置显式
  路由到最细覆盖层。

**验收**：2 层 AMR 下 10 万步，phi max 稳定 ~300–400 V 量级、
通道区密度/电势剖面与现双层方案一致（定性 + 关键截面定量）。

### 阶段 2：物理语义与性能调通（1–2 周）

- [ ] **粗层粒子**：确认粗层 push/gather 使用粗层 E 场（原生
  `Efield_aux` 逐层已有）；检查粗细界面处粒子跨层迁移时 gather/deposit
  无双计、无漏计（SyncRho 的 OwnerMask 路径）。
- [ ] **壁面交互在粗层**：通道壁面几何全在 level 1 内，预期粗层无撞壁；
  验证扩展域外边界的粒子吸收（`boundary.particle_* = absorbing`）
  在粗层正常工作。
- [ ] **阳极电流/壁电荷诊断跨层汇总**：累加器本身层级无关，
  验证数值与单层基线连续。
- [ ] **性能摸底**：粗层新增 ~10⁸ 宏粒子（前期估算：0.3–0.4 μs 渡越
  × 3.4×10¹⁴ 宏粒子/s 注入率，电子+离子合计 ~2.4×10⁸，~15 GB 内存），
  测量每步耗时增量（预期粒子部分 5–6×），评估 arena/显存预算
  （注意本机另一进程长期占 ~78 GB 显存）。
- [ ] **启动暂态拉长**：粗层灌满需 ~10 万步（0.3 μs），确认物理判稳
  时间点相应后移，必要时延长初始化算例或采用密度渐变启动。

**验收**：稳态段（10–15 万步）阳极电流、壁电荷、推力等物理量与
双层方案基线在合理偏差内；无 22.5 万步 6 kV 类失控（配合 SEE 配置）。

### 阶段 3：全功能对齐（可选，2 周以上）

- [ ] **ECDI 滤波器多层重设计**（A2，算法级）：最细层做 (r,z) 投影 →
  限制到粗层，或逐层独立滤波；同时解除单 rank/单 box 限制。
- [ ] **IonizationSourceTable 多层化**（A3）。
- [ ] **AVERAGE_USE 背景模式**（B7）：per-level FAB 输入约定。
- [ ] **regrid 生命周期挂钩**（防御性）：若未来启用动态 regrid，
  材料/壁电荷/背景密度/诊断缓存需挂到 `RemakeLevel` 重建。
  当前静态网格下仅留 TODO 注释。
- [ ] **3D 以外的维度**：`HallElectrostaticConfig.H:5-10` 的编译开关
  与 2D benchmark 路径（`InsertInjection.cpp:60`）适配。

**验收**：全部 Insert 功能在 `max_level=1` 下可用，无 Abort 占位。

## 5. 风险与对策

| # | 风险 | 影响 | 对策 |
|---|---|---|---|
| R1 | `Geom(0)` 语义反转导致隐性错误（注入偏移、诊断分箱、阳极几何判断用了错误层的 dx） | 结果错误且不报错 | 全库 grep `Geom(0)`/`geom(0)`/`[0]` 硬编码清单逐项过；阶段 0 对比算例兜底 |
| R2 | 外加磁场文件读取 `max_level==0` 限制 | 初始化即失败 | 阶段 1 第一项解决；粗层从同文件按粗 Geom 采样 |
| R3 | 注入粒子层级归属不符预期 | 粒子消失或层级错位 | 阶段 1 专项小算例验证；必要时显式路由 |
| R4 | 粗层粒子开销超预期（~10⁸ 宏粒子、显存、每步 5–6×） | 性能不可接受 | 阶段 2 摸底后决策：提高粗层粒子权重（仅粗层）、缩小扩展域、或回退中间路线 |
| R5 | 粗细界面电荷/粒子归属双计或漏计 | 界面处 rho/phi 假结构 | 依赖原生 SyncRho/OwnerMask；界面附近画 rho 剖面专项核查 |
| R6 | 初始 tagging 路径与 abort 桩冲突 | level 1 无法生成 | 阶段 0 用上游 MR 测试先行确认实际调用路径 |
| R7 | 与"去 SEE"等物理实验配置叠加后难以归因 | 物理结论混淆 | AMR 切换与物理配置变更分开验证，每次只动一个变量 |

## 6. 备选中间路线（不改造也可缓解远场问题）

保留现有双层电势结构，**把穿过细层边界的离子通量折算为电荷注入粗层 rho**
（粒子本身不进粗层）。工作量约 3–5 天，可把羽流电荷的主导项补进粗层
泊松解。若阶段 2 性能摸底不通过，回退到此方案。

## 7. 工作量汇总

| 阶段 | 内容 | 估时 |
|---|---|---|
| 0 | 基线验证 + 物理等价性对比 | 1–2 天 |
| 1 | 最小可用版（主算例 2 层 AMR 跑通） | ~1 周 |
| 2 | 物理语义与性能调通 | 1–2 周 |
| 3 | 全功能对齐（ECDI/电离源表/USE 模式，可选） | 2 周+ |
| 合计（最小可用） | 阶段 0–2 | **2–3 周** |
| 合计（全功能） | 阶段 0–3 | **4–6 周** |

主要风险不在代码量（主代码骨架现成、Insert 已半多层级化），
而在物理验证：粗细界面的电荷/粒子归属、粗层分辨率减半对羽流
自洽电势的影响、以及 10⁸ 量级粗层粒子带来的性能与显存压力。
