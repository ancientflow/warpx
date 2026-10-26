# 阶段 5：多层 Robin、历史与取场

日期：2026-10-26。实施起点：`027ab3ad1`。
范围：3D Cartesian、静态非周期网格、逐层单向 Poisson 求解；内部细区不接触物理边界。

## 实现

- `LabFrameExplicitES` 的两种材料边界都通过公共 `computePhi` 和同一个
  `PoissonLevelSolver` 回调逐层求解。删除独立单层 `MaterialPoissonSolver::solve`
  入口，继续复用公共粗插值、`NodalPoissonInitialGuess` 及求解后统一历史提交。
- monopole 初始化显式设置 boundary handler 的 Robin 类型、非周期标志及非
  Dirichlet 标志，不再留下未初始化的边界数组。六个场边界仍要求 `none`。
- `solveLevel` 只在 `lev==0` 设置物理 monopole 条件。内部细层使用阶段 4 的局部
  Dirichlet 算子，接口电势由公共粗插值提供；各层独立复用算子、MLMG 和材料缓存。
- 公共调度复制物理 rho 并乘一次 `-1/epsilon0`。材料在自己的 RHS 上，仅对
  Robin level 0 的物理面/棱/角再乘 2/4/8；这一无量纲权重不作用于同层 Box 接缝
  或粗细接口，也不改写 rho 或持久壁电荷。物理 Robin forcing 仍由算子添加。
- Robin 的 phi 和 E 同层 ghost 刷新扩展到所有层。只有 level 0 追加物理 Robin
  phi ghost；细层保留公共粗插值提供的取场 ghost。Dirichlet 物理面原有二阶
  单边导数等价 ghost 保留，未增加材料专用 E 插值或 gather/deposition 路径。

## 公共静电 coarse-aux 修复

取场验证发现了计划中记录的原生疑点：公共静电 aux 分支刷新了 fine aux，却没有
刷新 `Efield_cax/Bfield_cax`。缓冲区粒子仍按原生分区读取 cax，因而可能读取初始零场。
关闭材料后同样复现，问题属于公共静电取场路径。

对照使用三个电子探针，分别位于细区内部、细区接口缓冲带、粗区。外加恒定网格电场
`(0.25,-0.5,0.75)` V/m，初速度为零，宏粒子权重 `1e-30`，三步 `dt=1e-10` s。
常电势材料解使自洽 E 接近零；原生对照使用六面同一电势。预期三者的速度均为
`(-13.19115006,26.38230013,-39.57345019)` m/s。

修复前接口粒子的 x 速度约为零（材料 `-2.35e-17`、原生 `-1.90e-11` m/s），
另外两个粒子符合解析值。使用 `clion-debugger` 技能尝试在粗层 gather 处设置
logpoint，Release 二进制报告“没有可执行代码与此行相关联”，未取得断点变量证据。
随后以原生/材料对照、源码数据流及修复前后的解析运动验证定位；已清理 agent 断点。

`WarpXComm.cpp` 的两种静电 aux helper 现在都通过公共辅助函数，将上一层当前 aux
复制到本层的 coarsened cax 布局（含所需 ghost），同时处理 E、B。沿用现有
`ParallelCopy` 和通信精度；电磁分支及粒子分区、gather kernel、缓冲宽度均不改变。
网格外场在 aux 更新之前加入 fp，随后通过同一复制进入 cax，不重复叠加。

修复后，原生、Dirichlet、monopole、collocated/staggered 的探针相对速度误差均
小于 `1.4e-11`，低于预设 `1e-7` 阈值。另以原生 momentum-conserving 配置覆盖
staggered-to-nodal aux 分支，误差 `1.37e-11`。非零网格磁场 `Bz=0.001` T、
初始 `ux/c=1e-4`、三步 `dt=1e-12` s 的 Robin 探针也符合回旋解析值，
相对动量范数误差 `1.36e-12`，三个位置一致。双 MPI rank 的 Robin 电场探针
相对速度误差 `3.37e-12`，验证跨 rank 的粗层场复制。

这证明了本次小问题中的缓冲区粗场刷新及网格外场传递；系统自力、跨接口轨迹、
buffer 宽度扫描和文件外场数据仍由后续阶段验收。

## 历史与约束检查

CLion `app_3d` 构建通过；CUDA 13.2 / DOUBLE / MPI，AMReX 仍为
`b552c85c432defab450688581a1c87aead7b461d`。复用已有输入并用命令行覆盖，
未创建测试输入、分析脚本或 CTest 注册；当前构建 `BUILD_TESTING=OFF`。

基础配置：`[-0.5,0.5]^3`，粗层 `32^3`，三级，max grid size 16、blocking factor 8，
tagging `[-0.25,0.25]^3`，原有球形电子源冻结位置；epsilon 为 `if(z<0,4,1)`，
阳极 `[-0.1875,-0.125]^3`，电势 `0.05+0.01*x+0.02*y+0.03*z` V。
相对/绝对容差均为 `1e-10`，pre/post 2/2，运行 4 步共 5 次逐层求解。
Robin 远场电势为 0.01 V，关闭 rho/ECDI 滤波，energy-conserving gather，
staggered 配置启用 Galerkin。外推 alpha、beta 显式设为 1；默认值 0 不产生高阶预测。

Dirichlet/monopole × collocated/staggered × order 0/1/2，共 12 组运行全部退出 0。
每层历史日志分别为：

| 指定阶数 | 五次求解实际阶数 |
|---|---|
| 0 | `-1,0,0,0,0` |
| 1 | `-1,0,1,1,1` |
| 2 | `-1,0,1,2,2` |

`-1` 表示无历史，保留当前合法初值。所有层固定节点误差为 0，系数构建次数为 1。
以 order 0 为对照，最终 plotfile 所有层/Box 的最大电势差 `8.88e-14` V，
最大电场分量差 `1.94e-12` V/m，低于事先采用的 `1e-8` V / `1e-6` V/m 比较阈值。
MLMG 仍采用 `greater` 相对准则，不能把每次相对收敛都解释为已达到绝对残差阈值。

补充检查：

- Robin 使用有限参数 alpha=0、beta=`1e308`，在历史足够时候选权重溢出；
  三层均记录 `finite_fallback=1` 并退回 order 0，求解正常结束。
- Dirichlet 六面电势 `0.01+1000*t`，三层每步固定节点误差仍为 0。
- Robin 两个细层均无阳极的三级问题正常求解，各细层所有 MG 电极计数为 0。
- Robin 三级 2 MPI ranks 运行通过；细区接触物理域边界仍明确拒绝。
- 三级 Robin 两次求解的 CUDA memcheck 退出 0，`ERROR SUMMARY: 0 errors`。

## 单向场求解空间检查

使用相同外域和各自相同物理边界的均匀 `128^3` 参考解。对照的 AMR 粗层为
`32^3` 或 `64^3`，细化比 2，细区固定为 `[-0.25,0.25]^3`，`n_error_buf=0`。
epsilon=1，阳极 `[-0.375,-0.3125]^3`、0.05 V，细区无阳极，仍使用原球形电子源。
所有比较使用 collocated 场；不同边界的解不互相作为参考。

将参考 plotfile 的 cell-centered 输出按体积平均到被比较的网格，使用 child mask
排除粗层覆盖区。定义 `r=max(|x|,|y|,|z|)`：细区内部 `r<0.1875`，接口带
`0.1875<=r<0.3125`，粗区 `0.3125<=r<=0.4375`，物理边界带 `r>0.4375`。
排除电极邻域 `(-0.4375,-0.25)^3`，避免将电极边缘奇异性混入接口指标。
电势 RMS 除以参考全域电势范围，向量 E RMS 除以参考全域最大 `|E|`；
验收预算预设为 5% / 10%，不是每个接近零场点的局部相对误差。

| 边界 / 区域 | phi 相对 RMS：32→64 | E 相对 RMS：32→64 |
|---|---|---|
| Dirichlet / 粗区 | `1.666e-3 → 4.824e-4` | `1.374e-3 → 3.476e-4` |
| Dirichlet / 细区 | `7.295e-4 → 2.046e-4` | `1.864e-4 → 5.196e-5` |
| Dirichlet / 接口带 | `1.486e-3 → 4.155e-4` | `6.699e-4 → 1.786e-4` |
| Dirichlet / 物理边界 | `1.236e-3 → 3.529e-4` | `1.905e-3 → 5.237e-4` |
| Robin / 粗区 | `4.271e-3 → 1.226e-3` | `1.453e-3 → 3.798e-4` |
| Robin / 细区 | `3.377e-3 → 9.607e-4` | `5.267e-4 → 1.503e-4` |
| Robin / 接口带 | `3.908e-3 → 1.113e-3` | `8.809e-4 → 2.468e-4` |
| Robin / 物理边界 | `5.379e-3 → 1.535e-3` | `1.331e-3 → 3.498e-4` |

各区域均随粗层加密改善。另将 32 粗层的 patch 缩至 `[-0.125,0.125]^3`，在公共
中心区域 `[-0.0625,0.0625]^3` 比较：Dirichlet 电势相对 RMS 为 `6.117e-4 → 6.211e-4`，
Robin 为 `3.264e-3 → 3.275e-3`；对应 E 误差也在预算内。未引入细层到粗层的场反馈。

## 日志、复现与限制

日志位于 `build/amr-stage5-<名称>/run.log`，plotfile 为同目录的 `diag000000` /
`diag000004`。历史矩阵名称为 `<dirichlet|monopole>-<collocated|staggered>-h<0|1|2>`，
空间矩阵为 `space-<dirichlet|monopole>-<32|64|128|small>`。其他名称与上文对应：
`robin-nonfinite`、`dirichlet-time`、`robin-no-anode-valid`、`robin-mpi2`、
`bad-physical`、`memcheck`、`probe-*-fixed`、`probe-*-staggered`、
`probe-native-momentum`、`probe-magnetic`、`probe-mpi2`；Dirichlet collocated
探针名称为 `probe-dirichlet`。修复前探针日志/输出保留在
`probe-robin`、`probe-native`；所有相对路径均以 `build/amr-stage5-` 为前缀。

复现沿用阶段 4 的球形输入命令，按上述基础配置替换阳极电势、max_step 和历史参数。
Robin 加 `insert.electrostatic_boundary=monopole`、六面 `boundary.field_lo/hi=none none none`
和 `insert.robin_phi_infinity=0.01`。输出使用已有 Full/plotfile 诊断，字段
`phi Ex Ey Ez rho`；分析以终端内联 Python/yt 完成，没有写入新的测试脚本。

一次三级无阳极试运行组合 `n_error_buf=0`、同一 tagging 范围和 max grid size 32，
在 AMReX 初始化调整网格时触发 `Bfield_fp[level=1] already exists`，尚未进入 Poisson。
该日志保留在 `robin-no-anode`，不计为通过；随后采用已有可运行的三级布局
（n_error_buf=1、max grid size 16、tagging `[-0.125,0.125]^3`）完成无阳极验收。
此网格重建问题未通过放开材料静态网格限制来绕过。

以上属于阶段 5 的小问题验收。CPU、变时间步/同刻重算的独立历史验收、生产文件场、
壁电荷跨接口守恒、完整缓冲宽度与自力扫描仍按后续/对应模块计划执行。
阶段 3/4 记录的更严格绝对残差停滞没有因本次接入而消失，也未修改容差掩盖它。
