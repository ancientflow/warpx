# 阶段 3：单层材料后端统一

日期：2026-10-24。范围：3D、完整单层、静态材料与几何。
实施起点：`Amr-support` / `4d20d4ca6d4d943e9d1327769d9a64f5b78dbe5b`。

## 实现

`Insert::MaterialPoissonSolver` 统一持有 `cht::MLNodeRobinLaplacian` 和 MLMG。
旧的原生材料 `MaterialPoissonCache` 及独立 `NodeRobinSolver` 已移除。
非材料调用仍进入原生 `computePhi`，不分配材料缓存。材料 AMR 在初始化时明确拒绝，
局部覆盖待阶段 4；不回退到原生材料算子。

- Dirichlet 初始化通过 `setRobinBC` 设置六面 `a=1,b=0`，每步在原生 `setPhiBC`
  和合法初值准备后从 `phi` 复制六面 forcing。系数只设一次，不清掉逐步更新的 forcing。
  保留原生面相交处的覆盖顺序、非零值及时间函数。原生 `boundary.potential_*`
  解析器仍只接受 `t`；适配器本身逐节点复制，可接收已填入 `phi` 的空间变化面值。
- monopole 仍通过 `setMonopoleBC` 配置，保留中心、无穷远电势和场边界 `none` 的约定。
- 内部 `setInternalDirichletValues` 持有节点电势场副本，允许非均匀阳极电压及逐步更新。
  仅物理解读取指定电势；修正方程使用原有齐次掩码、平滑和插值约束。
  标量 `setInternalDirichletValue` 接口仍可用，调用后撤销场覆盖。
- 保留现有粒子沉积、SyncRho、壁电荷、PEC 处理次序。固定节点源的清零移入私有 RHS，
  `rho_fp` 不再因阳极掩码被原地清零。Robin 的面/棱/角双体积归一化仅处理私有源，
  再乘 `-1/epsilon_0`；Dirichlet 不施加该归一化。
- 物理残差使用 `compPhysicalResidual` 加入非齐次 Robin 源，固定节点由节点算子排除。
  `computeE`、现有物理面 ghost、aux、gather、沉积及缓冲配置未增加新的实现路径。
- 显式 pre/post 仅经 `applyMLMGOptions` 设置，不再硬编码 1/2。
  未指定时统一采用 `MLMGOptions` 默认 2/2。与旧 Dirichlet 的 1/2 耗时不可直接等同比较。

## 缓存与 MG

算子先于 MLMG 声明，保证逆序析构时 MLMG 先释放。复用时检查布局、几何、材料对象和
层次构造选项，并将实际 epsilon/mask 与算子拥有的副本比较，拒绝原地修改。
只更新 RHS、phi、固定电势和求解控制选项；时间历史仍由上层公共组件持有。

粗化候选按实际节点注入逐级检查：

- 固定节点在 z 向仍跨至少两层，不能仅保留一个阳极节点；
- 粗单元角点不能丢失其内部细节点的固定/未知类别；
- 两条细边的端点类别相同时，中点不能相反，避免粗化抹掉薄结构或连接隔开的区域；
- 同一候选粗单元的八个细胞 epsilon 必须一致，避免直接平均掉介质跳变。

这是保守的离散几何保护；曲面、斜面或连续变化 epsilon 可能使粗化停在很细的层，
甚至只剩一个 MG 层。它不是连续几何拓扑或生产性能已通过的证明。实际层数还受 AMReX
布局与层次选项限制。若性能不足，应另行评估几何粗化方案，不能直接恢复固定上限或
“至少一个阳极节点”的旧判据。

`insert.robin_diagnostics=1` 保留原参数并扩展到两种材料边界：输出有效 pre/post、
实际 MG 层数、粗化上限、各 MG 层唯一固定节点数、初末残差、外层迭代次数、
bottom Krylov 迭代序列和材料系数重建次数。准备、求解、额外物理残差分别计时。
诊断模式提高 MLMG/bottom 日志等级，保留 AMReX 的 bottom 失败与回退信息；
当前 AMReX 没有已存储 bottom 返回码的公开查询接口，不将迭代数冒充成功状态。

核对的 AMReX 提交为 `b552c85c432defab450688581a1c87aead7b461d`。
`MLLinOp::setDomainBC` 将 Robin 存为内部 Neumann，父类已完成反射 ghost 和半体积权重。
删除了检查 `LoBC()/HiBC()==Robin` 的无效重复分支，没有重复增加边界权重。

## 验证记录

CLion Release 已成功构建仓库内 `build` 的 `app_3d`；CUDA 13.2、MPI、DOUBLE、EB、
FFT、openPMD，AMReX 版本如上。最终可执行文件更新时间为 2026-10-24 10:33（本地时间）。
完整重编译约 38 分钟，主要尾部等待来自现有 `ImplicitPushPX.cpp` 的 CUDA 汇编。
CLion 初次工具调用超过 300 秒后，IDE 内构建继续执行，最终编译和链接成功。

运行使用 NVIDIA GeForce RTX 5070 Laptop GPU、1 个 MPI rank。复用现有
`Examples/Tests/electrostatic_sphere/inputs_base_3d`，以命令行覆盖参数；未新增测试输入、
Python 分析文件或 CTest 注册。当前 `BUILD_TESTING=OFF`，`ctest -N` 为 0 项，
以下结果为直接运行新二进制的记录。仓库 `main` 会在输入名前加 `../../Script/`，
运行目录和相对输入名按此约定设置。复杂表达式的引号经 CLion 运行工具传递时丢失，
实际验收使用命令行启动已由 CLion 构建的可执行文件。

固定源基线：`[-0.5,0.5]^3`，`32^3` 单元、8 个同层 Box；原球形电子电荷源，
关闭粒子推进；epsilon 为 `if(z<0,4,1)`，阳极为 `[-0.375,-0.25]^3`，
阳极电势为 `5+x+2*y+3*z`。相对/绝对容差均为 `1e-10`，运行 3 步，共 4 次求解。
两种模式均得到 3 个 MG 层，唯一固定节点数为 `125 / 27 / 8`，
每次运行的系数构建次数始终为 1，无 bottom 失败信息。

| 边界 / pre/post | 首次迭代数 | 首次求解耗时（s） | 第 4 次物理残差 | 第 4 次求解耗时（s） |
|---|---:|---:|---:|---:|
| Dirichlet / 2/2 | 7 | 0.120996 | 7.4351e-11 | 0.002379 |
| Dirichlet / 4/4 | 5 | 0.101479 | 6.3300e-11 | 0.002329 |
| monopole / 2/2 | 10 | 0.280736 | 2.1199e-11 | 0.002142 |
| monopole / 4/4 | 8 | 0.219414 | 2.2721e-11 | 0.001545 |

耗时为本机单次、rank 0 记录，诊断额外开销另计，不作为生产性能结论。
MLMG 使用 `greater` 范数，首次求解的相对阈值由较大的初始残差决定，
后续固定源求解进一步达到绝对阈值；因此表中明确区分首次和第 4 次求解。

另外 5 项运行退出码均为 0：

- `dirichlet-time`：六面 `boundary.potential_*=1+1.e5*t`；每步新面值触发求解，
  最终物理残差 `6.0872e-9`，满足当前相对准则，系数仍只构建一次。
- `dirichlet-staggered`：最终物理残差 `7.3911e-11`。
- `monopole-staggered`：最终物理残差 `1.9976e-11`，启用 Galerkin gathering 配置。
- `dirichlet-volume-source`：电子密度覆盖整个域，包括固定阳极节点；最终物理残差
  `7.6707e-11`，验证私有 RHS 固定节点处理可运行。
- `native`：关闭材料，原生求解正常结束，残差 `7.2579e-12`，无材料缓存诊断输出。

额外压力检查 `monopole-strict` 将相对及绝对容差都设为 `1e-12`：首次求解成功，
下一次求解在 200 次迭代后残差为 `8.6816e-12`，未达到绝对阈值并中止。
这是本次观测到的更严格精度限制，保留失败日志，不计为通过；未修改断言或容差掩盖它。
当前验收范围是上述 `1e-10` 设置的单层基线；CPU、多 rank、生产 Hall 输入及空间误差
收敛性仍需后续验证。阶段 1 各阶历史的独立验收仍按原计划进行。

日志保存在仓库内 `build/amr-stage3-<名称>/run.log`。可从仓库根目录复现 Dirichlet 基线：

```bash
mkdir -p build/amr-stage3-dirichlet
cd build/amr-stage3-dirichlet
../bin/warpx.3d \
    '../Examples/Tests/electrostatic_sphere/inputs_base_3d' \
    'max_step=3' \
    'amr.n_cell=32 32 32' \
    'amr.max_grid_size=16' \
    'amr.blocking_factor=8' \
    'diagnostics.enable=0' \
    'electron.do_not_push=1' \
    'electron.do_not_gather=1' \
    'algo.field_gathering=energy-conserving' \
    'warpx.grid_type=collocated' \
    'warpx.self_fields_required_precision=1.e-10' \
    'warpx.self_fields_absolute_tolerance=1.e-10' \
    'warpx.self_fields_phi_extrapolation_order=2' \
    'warpx.limit_verbose_step=0' \
    'amrex.the_arena_init_size=16777216' \
    'insert.use_electrostatic_materials=1' \
    'insert.relative_permittivity_function(x,y,z)="if(z<0,4,1)"' \
    'insert.anode_implicit_function(x,y,z)="if((x>=-0.375)*(x<=-0.25)*(y>=-0.375)*(y<=-0.25)*(z>=-0.375)*(z<=-0.25),-1,1)"' \
    'insert.anode_potential_function(x,y,z)="5+x+2*y+3*z"' \
    'insert.robin_diagnostics=1' \
    'warpx.do_electrostatic=labframe' > run.log 2>&1
```

monopole 基线在同组参数末尾追加 `insert.electrostatic_boundary=monopole`、
`'boundary.field_lo=none none none'`、`'boundary.field_hi=none none none'`、
`warpx.use_filter=0` 和 `insert.ecdi_control.enabled=0`。4/4 扫描再追加
`warpx.self_fields_num_pre_smooth=4 warpx.self_fields_num_post_smooth=4`。
`git diff --check` 通过。
