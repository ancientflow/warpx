# 阶段 7：Hall 主算例联合接入

日期：2026-10-27。实施起点：`47efa679d4eec4bff70143ead03712169795c25b`。
状态：联合输入、源码修复和回归用例已写入；静态检查通过，编译及运行验收待完成。
阶段 2 的独立数值验收仍未完成；阶段 5–6 已有记录不替代完整 Hall 联合验收。

## 联合输入

`Script/PHT/3d_hall_dielectric_slot_amr` 继续包含原单层主输入，独立覆盖域和网格：

- 扩展域 `[-5,15] × [-5,15] × [-8.75,16.25] mm`；粗层 `256×256×320`。
- 细区 `[0,10] × [0,10] × [-2.5,10] mm`，细化比 2；两层间距分别为
  `78.125 μm` 和 `39.0625 μm`。
- `lx/ly/lz`、阳极/环缝/阴极位置、粒子权重、流量、SEE 及 species ndt 保持主输入值。
- 两层按各自物理坐标执行相同墙面和介电函数。粗层新增的 `z=0` 陶瓷区域
  不被裁剪掉；阳极环仍保持其物理尺寸。
- 显式选择 Dirichlet，保留二阶电势初猜、直接背景 MCC、DSMC 和文件外加磁场。
  静态布局、关闭 AMR subcycling/负载均衡、关闭普通及 ECDI 滤波的限制保持生效。
- 开启 `insert.wall_charge_diagnostic=1`；Full 输出改为当前运行目录下
  `diags/hall_amr`，不再继承主输入中的个人磁盘绝对路径。
- 出口统计若开启，对应扩展域 `z=16.25 mm`，不代表旧域 `z=10 mm` 的截面统计。

原单层输入未修改。把本扩展域输入的 `amr.max_level` 改成 0 只会删除细层，
不会恢复旧域和旧分辨率。继续沿用 staggered 主设置，但其生产接口自力仍待验证；
没有将阶段 6 的 collocated 小问题或 24 格缓冲结论推广到本算例。

## 修复的集成时序问题

此前 `WarpX::InitData` 已求解初始场并提交 phi 历史、写出初始诊断，
随后 `Evolve` 中的 `Insert::Initialize` 才执行 `insert.initial_sources`。
初次半步动量推进因此可能读取尚未包含这些源的场，初始诊断也遗漏这些粒子。
此外，每次分段调用 `Evolve` 都会再次执行初始注入。

现在 `WarpX::InitFromScratch` 在 `mypc->InitData()` 之后调用
`Insert::InitializeParticleInjection()`。该入口内部按维度/编译开关处理，
未配置初始源时不注入。`AddNParticles` 已包含全层 `Redistribute`，
首次沉积前粒子已有正确层及 MPI 归属。移除 `Insert::Initialize` 中的重复入口，
不增加首次推进前的额外求场或历史提交。checkpoint 恢复不会再次注入初始源；
这不代表壁电荷等物理状态的 checkpoint 恢复已通过验收。

连续注入此前单独读取 `my_constants.dt`。现在从 `warpx.getdt(0)` 读取本次推进步长，
并传入当前步号；框架内 source 的显式步长覆盖仍生效。主 Hall 输入中的两者相等，
其流量不因此改变。原来两者不等的配置会按实际推进时间修正注入量。
这两项修复同时适用于单层和多层。

## 数据所有权及消费时刻的源码核对

| 数据 | 生产与消费顺序 | 避免过期/重复数据的方式 |
|---|---|---|
| 初始源与 phi 历史 | 初始注入/重分布 → 初始沉积/求场 → 提交历史 → 首次半步推进 | 本次修复；不在 Evolve 中重新注入 |
| 连续源 | 取本次场 dt → 注入/重分布 → 碰撞 → 推进 | 不再单独读取可能不一致的 dt 常量 |
| MCC 密度、bins、可消耗数 | 每次实际耦合碰撞调用前重建；每个 subcycle 也重建 | 使用当前中性粒子布局、权重和索引，不以 species ndt 判断有效性 |
| 碰撞 dt | CollisionHandler 按 subcycle/supercycle 传入 | 不乘 AMR 细化率，不改物种推进调度 |
| 墙面交点及产物 | 推进 → 墙面处理 → 交点路由/电荷沉积 → 粒子全层 Redistribute | 交点按最细覆盖层归属；产物按最终位置归属；单事件只记录一次净留存电荷 |
| 多层壁电荷 | 持久原始 FAB/ghost 贡献 → 加入重新沉积的工作 rho → 一次 SyncRho | 限制/同层求和不回写持久场；诊断在独立副本上的 SyncRho 不计为再次装配物理源 |
| 单层壁电荷 | 粒子 rho 同步/PEC/滤波 → 独立同步壁电荷副本 → 加入 RHS | 保留原单层源处理顺序 |
| 细层 phi 初猜 | 当前粗层求解 → 粗插值 → 本层未知节点历史外推 → 约束恢复/求解 | 历史独立存储，接口使用当前粗解；所有层成功后才提交历史 |
| 材料缓存 | 初始化固定 epsilon/mask/布局；每次更新 RHS 和指定电势/边界 forcing | 不把算子复用等同于物理源复用；布局变化仍明确拒绝 |
| 外场及 aux/cax | computeE → 各层叠加外场 → 推进/同步前更新 ghost、aux/cax | 沿用阶段 5–6 的公共路径，不在材料模块另加外场或缓存 |

对应实现主要位于 `WarpXEvolve.cpp`、`LabFrameExplicitES.cpp`、`WarpXComm.cpp`、
`CoupledBackgroundMCCCollision.cpp`、`BackgroundCoupledDensity.cpp`、
`AnalyticBoundaryInteraction.cpp` 和 `WallCharge.cpp`。这些是源码证据，
本会话没有断点、日志或数值结果可证明上述联合运行行为。

## 新增回归及其边界

`Examples/Tests/hall_amr` 注册 `test_3d_hall_amr_single`、`multi`、`mpi` 三组
run/analysis，分别使用单层、两层、两层双 MPI rank。小问题采用 `32^3` 粗网格、
4 步、介电跳变、二阶外推、每步两次 MCC 子步、每步排序和恒定外加磁场。
中性物种 ndt=100；第一步在粗区和细区各中性化一个离子。

分析预先检查：

- t=0 已有两电子、两离子、两中性粒子，分离的带电初始源已产生非零场。
- 实际 dt 为 `my_constants.dt/2`，连续源每实际步注入期望 2.3 个宏粒子，
  累计计数为 `2+floor(2.3*step)`；这能识别旧步长和旧初始注入时序。
- 两离子在第一步变为两中性产物；此后粒子数量与壁电荷保持正确。
- 累计壁电荷及独立复合积分均为 `200e`，多层的原始壁电荷各为 `100e`；
  反复装配不允许重复累加。交点刻意远离粗细接口，不能据此证明接口积分精度。
- 各层内部 Bz 始终为 `0.001 T`，不允许累加；phi/E 有限。

截面为独立的非零弹性小表，不使用生产 Xe2X 数据。此用例覆盖 MCC 调用、
排序和墙面产物存在时的组合路径，但不要求发生随机散射事件，
也不检验电离消耗、SEE 产额、文件磁场覆盖或生产自力误差。
不能用它替代阶段 2 独立碰撞验收、阶段 5–6 精度验收和生产统计。
分段 Evolve 不重复注入目前只有调用位置的源码证据，新用例没有直接覆盖 Python 分段调用。

## 本次验证与剩余阻塞

- `git diff --check`、分析脚本 Python 语法检查通过。
- 使用仓库真实 `add_warpx_test` 函数的独立 CMake 注册检查通过：六个 run/analysis、
  单层/多层参数和截面复制正确。注册检查使用不可执行的占位目标，未编译、未运行 WarpX。
- 静态核对联合输入的 include、域尺寸、层间分辨率及物理参数继承。
- 本会话无 CLion 工具。仓库要求只能通过 CLion 编译，所以未重建 `app_3d`；
  现有 Release 可执行文件日期早于当前源码，未将其用作新改动的验收证据。
- 仓库 `data/` 不存在，缺少 `B3D_dielectric.bp5`、`hall_start_256_600K_slot.bp5`；
  约定运行目录中的生产 `Xe2X` 截面目录也不存在。不能确认实际数据版本和外场模板覆盖。
- 未执行新回归、单层原生回归、独立 MCC 或完整 Hall 运行，尚无耗时及守恒实测结果。

具备构建工具后，先通过 CLion 刷新启用测试的 3D 配置，只构建 `app_3d`，再运行：

```bash
ctest --test-dir <build> --output-on-failure \
  -R '^test_3d_(hall_amr_(single|multi|mpi)|particle_rotating_axis|electrostatic_sphere_lab_frame)\.(run|analysis)$'
```

随后补齐阶段 2 独立验收及实际数据，按阶段 7 的联合输入运行短步检查，
记录实际源、壁电荷和碰撞账目。阶段 7 在这些检查完成前保持未整体验收状态。
