# AMR 阶段 0：接口与基线记录

源码基线为 `Amr-support` / `e1d320e5c`。本次已完成接口迁移和静态核对；
CLion 编译及原生数值基线仍待运行，不能将阶段 0 整体标为数值验收完成。
构建配置、实际依赖版本、输入 SHA-256 和外部数据存在性见
[AMRStage0Baseline.json](AMRStage0Baseline.json)。

## 已落实的接口边界

`ablastr::fields::computePhi` 和 `ElectrostaticSolver::computePhi` 不再包含
Insert 材料头文件、缓存类型、材料开关或阳极恢复逻辑。
公共的 `PoissonLevelSolver` 是可选同步回调，默认空值仍调用原生算子。
材料启用后，`LabFrameExplicitES` 才构造 `Insert::MaterialPoissonSolver` 适配器。

| 数据/操作 | 所有者及调用约定 |
|---|---|
| rho/phi、Geometry、BoxArray、DistributionMapping | WarpX / fields；回调只借用当前求解期间的数据 |
| 材料 epsilon、mask、指定电势 | `ElectrostaticMaterial`；适配器借用，公共求解器不解释其物理含义 |
| 材料 linop/MLMG | `LabFrameExplicitES` 的逐层材料缓存持有，先析构 MLMG 再析构算子 |
| 回调对象 | 当前 `ComputeSpaceChargeField` 调用的局部变量；不存入 MLMG，不跨时间步保存借用引用 |
| RHS 缩放与还原 | 公共层保留当前 `-rho/epsilon0` 与求解后还原；材料回调不得再次缩放 |
| 逐层求解 | 回调接收层号、LPInfo、phi、已缩放 RHS、MLMGOptions 和有效绝对容差；完成该层求解后返回仍有效的 MLMG 引用 |
| 材料固定电势、RHS 共享节点同步 | 材料适配器在当前层求解前完成；全场粗插值之后恢复阳极电压 |
| 粗到细插值和 post-solve 操作 | 公共层在回调返回后执行；回调不得重复插值、重复提交历史 |
| phi 历史 | 仍由 `LabFrameExplicitES` 持有，本阶段只规定契约，独立组件及逐层深度在阶段 1 实现 |
| 耦合背景密度与碰撞 bins | `global_background_density` 中各背景物种对象；不借用材料算子或 phi 历史 |

回调的 MLMG 与算子至少活到当前层 post-solve 操作结束；本次材料适配器返回的是
持久缓存中的 MLMG。空回调不会构造材料适配器，也不会填充材料缓存。
原生缓存对象/求解器由其自身路径管理。公共回调目前限定非 EB 的 MLMG 求解，
FFT/IGF 和 EB 原生调用仍走既有分支。

本次保持材料原生 `MLNodeLaplacian`、MG 粗化上限 4、pre/post 为 1/2 的既有行为，
只迁移其归属。这些设置及自定义算子统一接入在阶段 3 修改，避免把接口重构与数值
变化混在同一个基线中。Robin 单层入口不变。私有物理 RHS、细层初值保护和逐层历史
仍属后续工作；此回调尚不是完整的 AMR 历史/接口协调器。

## 独立启用及能力范围

| 模块 | 现有启用方式 | 未启用时 | 当前限制 |
|---|---|---|---|
| 材料 | `insert.use_electrostatic_materials=1` | 默认 0；不构造材料求解适配器，材料缓存为空 | 3D lab-frame、MLMG、非 EB、静态网格；保持现有限制 |
| Robin | `insert.electrostatic_boundary=monopole`，同时要求材料启用 | 默认 dirichlet；Robin 内部 cache 为 null | 单个完整 AMR 域、统一阳极电压等原有限制 |
| 时间初值 | `warpx.self_fields_phi_extrapolation_order/alpha/beta` | 无新增总开关；order 0 是沿用上一解，不代表关闭历史 | 当前深度仍全层共享，细层初值仍会被粗插值覆盖 |
| 耦合 MCC | `<collision>.type=coupled_background_mcc`；背景由 `collisions.background_species` 指定 | 未选择该类型时不构造其碰撞对象；背景列表为空时不创建密度 MultiFab/bins，更新/清理直接返回 | 当前仍有 `flvl==0` 限制，阶段 2 适配 |

`WARPX_USE_ELECTROSTATIC_MATERIALS` 当前由 3D 配置头自动定义，不是独立 CMake
功能开关；不能以运行时材料关闭冒充“材料代码未编译”。本次公共求解头已不依赖该宏。
`MCC_DENSITY`、`MCC_DELETE`、`NUMP`、`BENCHMARK_2D` 及 3D `WARPX_INSERT_3D` 的原有宏
未修改。显式指定背景物种仍会创建背景数据，即使未选择耦合碰撞；这不是材料模块隐式启用。

## 时间、提交与失效契约

- 场电势对应 `warpx.gett_new(0)`，当前显式静电在粒子边界处理/重分布之后求解。
  `LabFrameExplicitES` 的顺序为沉积与 SyncRho → 材料源/约束 → 求解 → phi ghost →
  历史提交 → computeE。此轮不重排调用，不让逐层求解回调提交历史。
- 阶段 1 将以实际场求解时间管理历史：同刻重复求解不推进历史，所有层成功后统一提交，
  步长改变按计划降阶。当前代码尚没有时间戳，也没有逐层独立有效深度，不能声称已实现。
- 碰撞接受 `CollisionHandler` 传下来的 `cur_time` 和 `dt`。subcycle 使用
  `dt/ndt` 并逐次更新子步时间；supercycle 自 `start_step` 起每隔 ndt 步使用 `dt*ndt`。
  该间隔与场历史时间、背景物种推进间隔分别管理，不能用 AMR 层号替代。
- 显式静电碰撞前后分别调用 `Insert::BeforeCollision/AfterCollision`；背景更新由
  原有物种推进/排序条件触发。本阶段只记录，不改变背景采样与碰撞时序。
- 新进程 restart 时，求解器缓存和 phi 历史不会作为 MLMG 对象恢复，应重新建立；
  当前历史从深度 0 开始。阶段 1 负责时间戳/候选的冷启动契约。
  `Insert::Initialize` 在进入 Evolve 时初始化背景；checkpoint 对物理背景/持久壁电荷
  的完整恢复需要另行验收，不能把求解缓存可重建等同于物理状态已正确恢复。
- 当前材料布局变化仍明确拒绝动态 regrid、负载均衡及 moving window；未来公共层
  分别通知历史、材料和碰撞缓存失效，各模块只能清理自己拥有的状态。
  背景密度释放继续在 AMReX finalize 前执行，材料适配器本身不延长数据生命周期。

## 原生基线及当前证据

已复用仓库现有输入，未新增测试工程、修改分析容差或重置 checksum。

| 基线 | 现有测试/输入 | 覆盖及状态 |
|---|---|---|
| N0 单层 | `test_3d_electrostatic_sphere_lab_frame` | 无材料/耦合背景配置；输入和 base 文件指纹已记录，数值未运行 |
| N0 多层 | `test_3d_electrostatic_sphere_lab_frame_mr_emass_10` | 两层，fine tagging 覆盖整个物理域；不能代替局部粗细接口验收，数值未运行 |
| N0 局部细区 | 以上 MR 输入的独立运行副本，fine_tag_lo/hi 改为 `-0.25` / `0.25` 三方向 | 待 CLion 构建后运行；原分析阈值能否适用须单独评估，不能伪装成既有测试通过 |

源检查确认这些基线默认材料关闭、背景列表为空；原生调用传入空逐层回调。
静态对照修改前非材料分支，算子构造/边界配置、求解操作、求解后插值/post-solve 回调
和 rho 恢复的源码操作保持一致；不构成编译、GPU/MPI 或数值等价证明。

`ctest --test-dir build-codex -N` 的精确过滤已列出上述两个现有测试各自的
`.run/.analysis/.checksum`，共 6 项。这仅证明旧配置中的注册存在。
`build-codex` 的缓存仍指向旧 AMReX，且没有可用于本次改动验收的已核实二进制；
不能直接将该目录的旧配置作为当前源码验证。

后续执行：在 CLion 中刷新选定的仓库内 3D 配置并启用 `BUILD_TESTING=ON`，记录实际
依赖版本，仅构建 `app_3d`。然后在该构建目录运行上述测试的 run/analysis；checksum
按仓库规则单独忽略，不修改物理分析断言。若选择 `cmake-build-release`，示例为：

```bash
ctest --test-dir cmake-build-release --output-on-failure \
  -R '^test_3d_electrostatic_sphere_lab_frame(_mr_emass_10)?\.(run|analysis)$'
```

局部细区变体需隔离工作目录，并使用相同初始源/容差记录粗细层残差、phi/E 与接口误差；
原生物理精度未通过时，先定位基线问题，再进入后续模块改造。

本会话缺少 CLion 工具，未执行编译和数值运行。两个 Hall 外部 BP5 数据集仍不存在，
已在 JSON 中记录；它们不阻碍原生球体基线，但阻碍完整 Hall 验收。
