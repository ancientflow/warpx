# 阶段 4：自定义材料算子的局部覆盖

日期：2026-10-24。实施起点：`Amr-support` / `d84181b2c`。
范围：静态 3D Cartesian、非周期 Dirichlet、逐层单向求解，细化比为 2。
monopole 仍限完整单层域，多层 Robin 和取场精度验收留在阶段 5。

## 实现与约束

`MaterialPoissonSolver` 为每个 AMR level 分别缓存一个自定义算子和 MLMG。
每个实例仍只接收一个 AMR level，Geometry 始终为该分辨率的完整物理域；
局部覆盖由 BoxArray 表示。初始化拒绝细区接触物理边界，不通过缩小 Geometry
把粗细接口伪装成物理面。

- 复用 `MLNodeLinOp::buildMasks`：在每个 MG level，根据整层 BoxArray 建立 cell
  coverage，节点的八个相邻单元中只要有域内未覆盖单元，就将该节点设为固定。
  同层 Box 的共享接缝有完整覆盖，保持未知。阳极 overset 与物理 Dirichlet
  约束合并到同一固定掩码，求解、平滑、插值、残差和私有 RHS 使用一致语义。
- 物理解保留公共粗插值给出的接口值；MG 修正采用齐次固定值。restriction 在完成
  可能的布局转移后，再用目标 MG 层掩码清除固定节点；prolongation 跳过固定节点。
- 每个 MG 层的 epsilon ghost 先清零，再用同层 `FillBoundary` 填入有效邻居，
  仅物理域外 ghost 使用原有镜像。局部域外没有虚构介质值；未知节点的体算子
  只读取被覆盖的相邻单元，系数加权 MG 插值可读取边界处的零延拓。
- 按后续要求移除了阳极与粗细接口的模板缓冲检查，初始化和 MG 粗化均不再
  以该缓冲为限制。保留阶段 3 的介质跳变、厚度和离散几何保护。
  阳极与粗细固定节点重合时，沿用现有阳极指定电势覆盖逻辑。
- 细层允许没有阳极：保留全 1 的静态材料 mask，但不安装虚构的 overset 约束，
  MG 电极计数为 0，粗化不再因“无阳极厚度”而提前终止。局部覆盖本身有粗细
  Dirichlet 约束，因此不作奇异系统投影；完整域纯 Neumann 的奇异性仍按实际约束判断。

Dirichlet 通过 `PoissonLevelSolver` 接回公共 `computePhi`。公共调度负责当步粗到细
插值、各层初值准备和求解顺序；上层继续用同一个 `NodalPoissonInitialGuess` 提交历史。
没有复制材料专用的历史或层级插值算法。回调收到已经乘 `-1/epsilon0` 的私有工作 RHS，
材料缓存再持有其受约束的求解副本。源清零覆盖物理面、阳极和粗细固定节点，
不修改沉积 `rho`，也不对 `rho` 作乘除缩放往返。非材料分支保留原生行为。

## ghost 审查

求解器 ghost 和取场 ghost 具有不同用途：

| 环节 | 数据来源与写入范围 |
|---|---|
| 求解内部 phi | 父类刷新同层邻居；粗细固定节点保留粗插值值；修正方程保持齐次 |
| 细层取场 phi | 公共 `interpolatePhiBetweenLevels` 在求解前提供有效节点及 collocated 所需的一层 ghost；MLMG 禁用 final-fill 回拷，求解后只刷新同层邻居 |
| Dirichlet 物理面取场 phi | 保留 `phi_g = 3 phi_b - 3 phi_1 + phi_2`，朴素中心差分等价于二阶单边导数 |
| 单层 monopole 取场 phi | 保留现有 Robin 一层 ghost 公式，未扩展到多层 |
| E | 删除材料专用的粗到细 E 插值；保留同层刷新和既有物理域外延拓，后续仍执行 WarpX 公共场边界策略 |

删除了 `ElectrostaticMaterial::fillPhiGhosts` 中额外分配整层临时 phi、再次粗插值
再复制 ghost 的逻辑。保留由公共调度写入的 coarse/fine ghost。
`computeE`、aux/cax、粒子分区、gather、沉积和缓冲区配置未增加新实现路径。
本阶段运行关闭粒子 gather/push；全链路取场和空间误差收敛不由这些代数检查代替。

## 阶段 4 原始验证记录（移除缓冲检查前）

通过 CLion 构建仓库内 `build` 的 `app_3d`，最终增量构建和链接成功。
CUDA 13.2 / DOUBLE / MPI / EB / FFT / openPMD；AMReX 提交仍为
`b552c85c432defab450688581a1c87aead7b461d`，GPU 为 RTX 5070 Laptop。
运行复用 `Examples/Tests/electrostatic_sphere/inputs_base_3d` 并使用命令行覆盖，
未创建测试输入、Python 分析文件或 CTest 注册。当前 `ctest -N` 为 0 项，
因此以下均为直接运行 CLion 构建二进制的记录。

共同配置：物理域 `[-0.5,0.5]^3`，粗网格 `32^3`，`max_grid_size=16`、
`blocking_factor=8`，epsilon 为 `if(z<0,4,1)`；保留球形电子源，关闭 push/gather、
rho/ECDI 滤波；pre/post 为 2/2，初值历史阶数为 0，相对/绝对容差均为 `1e-10`。
运行 3 步，共 4 次逐层求解。阳极默认为 `[-0.1875,-0.125]^3`，与电子球分离。
tagging 范围 `[-0.25,0.25]^3` 经 AMReX 缓冲和聚类后，各层分别为 8 / 27 / 125 个 Box，
细层覆盖率分别为物理域的 42.1875% / 24.4140625%。

| 运行目录后缀 | 配置 | MG 层数（各 AMR 层） | 第 4 次物理残差（各 AMR 层） |
|---|---|---|---|
| `two-source` | 两层，阳极 `5+x+2*y+3*z` V | 2 / 3 | `9.3551e-12` / `3.9379e-11` |
| `three-low-voltage` | 三级，阳极 `0.5+0.1*x+0.2*y+0.3*z` V | 2 / 3 / 4 | `1.4924e-11` / `1.1359e-11` / `2.6250e-11` |
| `no-anode-interior` | 三级，两个细层都无阳极 | 2 / 4 / 4 | `2.9141e-12` / `4.8224e-11` / `1.4773e-12` |
| `mpi2` | 两层，低电压，2 MPI ranks | 2 / 3 | `1.4773e-11` / `1.0135e-11` |
| `staggered` | 三级，低电压，staggered | 2 / 3 / 4 | `1.4924e-11` / `1.1359e-11` / `2.6250e-11` |
| `single-dirichlet` | 单层 Dirichlet 回归 | 2 | `9.3533e-12` |
| `single-monopole` | 单层 monopole 回归 | 2 | `1.3043e-11` |
| `native-amr` | 两层，关闭材料 | 原生后端 | 正常退出，无材料缓存输出 |

上述运行退出码均为 0。材料每层系数构建次数始终为 1；诊断在求解前复用残差工作场
记录固定电势，求解后检查所有固定节点，误差均为 0。三级含阳极场景每层的未知
Box 边界条目为 5768 / 26372 / 149424（包含共享副本），确认接缝没有整体固定。
无阳极细层的所有 MG 电极计数均为 0。首次求解采用 MLMG `greater` 相对准则，
表格列出重复求解后的绝对残差，不能把首次相对收敛误称为已达到绝对阈值。

最终二进制另以 `compute-sanitizer --tool memcheck --error-exitcode 99` 运行三级无阳极
细区（`max_step=0`，其他设置同 `no-anode-interior`），退出码 0，报告
`ERROR SUMMARY: 0 errors`，日志为 `build/amr-stage4-memcheck/run.log`。
`git diff --check` 通过。当前环境未提供 pre-commit 可执行程序，未运行完整 hook 集合。

非法输入检查均按预期中止：

- `bad-interface`：当时阳极 `[-0.375,-0.25]^3` 接触实际细区边界，触发模板缓冲断言；
  该断言现已移除，此项仅保留为历史记录。
- `bad-physical`：细化整个物理域，触发细区必须严格在域内的断言。
- `bad-monopole`：两层 monopole，触发仍限完整单层的断言。

保留的精度压力失败：`three` 使用与两层相同的约 5 V 阳极，首次三级求解成功，
下一次 level 2 在 200 次迭代后残差为 `1.864958260e-10`，没有达到 `1e-10`
绝对容差。将电压降低十倍、保持容差不变后上述三级用例通过；这一现象与细网格
大电势的残差舍入误差相符，但本阶段没有证明其唯一原因，也没有修改收敛判据。
该失败不计入通过项，不能据此宣称任意电压/网格都能达到同一绝对残差。
另一次无阳极配置实际仅在最细区边界截到一层阳极节点，被厚度断言拒绝；
最终无阳极用例缩小 tagging 范围，使阳极完全位于细区外。

运行日志保存在 `build/amr-stage4-<后缀>/run.log`。从仓库根目录复现三级基线：

```bash
mkdir -p build/amr-stage4-three-low-voltage
cd build/amr-stage4-three-low-voltage
../bin/warpx.3d \
    '../Examples/Tests/electrostatic_sphere/inputs_base_3d' \
    'max_step=3' 'amr.n_cell=32 32 32' 'amr.max_grid_size=16' \
    'amr.blocking_factor=8' 'amr.max_level=2' \
    'warpx.fine_tag_lo=-0.25 -0.25 -0.25' 'warpx.fine_tag_hi=0.25 0.25 0.25' \
    'diagnostics.enable=0' 'electron.do_not_push=1' 'electron.do_not_gather=1' \
    'algo.field_gathering=energy-conserving' 'warpx.grid_type=collocated' \
    'warpx.self_fields_required_precision=1.e-10' \
    'warpx.self_fields_absolute_tolerance=1.e-10' \
    'warpx.self_fields_phi_extrapolation_order=0' 'warpx.limit_verbose_step=0' \
    'amrex.the_arena_init_size=16777216' 'insert.use_electrostatic_materials=1' \
    'insert.relative_permittivity_function(x,y,z)="if(z<0,4,1)"' \
    'insert.anode_implicit_function(x,y,z)="if((x>=-0.1875)*(x<=-0.125)*(y>=-0.1875)*(y<=-0.125)*(z>=-0.1875)*(z<=-0.125),-1,1)"' \
    'insert.anode_potential_function(x,y,z)="0.5+0.1*x+0.2*y+0.3*z"' \
    'insert.robin_diagnostics=1' 'warpx.do_electrostatic=labframe' \
    'warpx.use_filter=0' 'insert.ecdi_control.enabled=0' > run.log 2>&1
```

无阳极细区用例将 tagging 改为 `[-0.125,0.125]^3`，阳极改为
`[-0.375,-0.3125]^3`；双 rank 用例将 `amr.max_level=1` 并用 `mpiexec -n 2` 启动。
单层用例设置 `amr.max_level=0`；monopole 另设置 `insert.electrostatic_boundary=monopole`
和六面 `boundary.field_lo/hi=none none none`。本机以 root 运行 MPI 时增加
`--allow-run-as-root`。实际运行仍遵循仓库 main 的 `../../Script/` 输入路径前缀约定。

## 移除缓冲检查后的验证

CLion `app_3d` 增量构建通过。复用原 `bad-interface` 的两层配置（阳极
`[-0.375,-0.25]^3`，阳极电势 `5+x+2*y+3*z`，`max_step=0`，相对/绝对容差
均为 `1e-10`），现在能够完成求解并正常退出，两层固定节点误差均为 0。
日志位于 `build/amr-stage4-interface-no-buffer/run.log`。
此次检查验证原先被拒绝的配置可以求解，不替代接口空间精度验收。
