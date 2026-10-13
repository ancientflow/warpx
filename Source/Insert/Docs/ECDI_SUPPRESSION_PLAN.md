# ECDI Suppression Control Case Plan

## 目标

本计划用于实现一个与完整三维 ECDI 算例对照的控制算例。控制算例仍使用三维粒子、
三维网格、三维边界条件和三维 Poisson 求解，但在每次静电场求解前，将三维电荷密度
源项投影到轴对称的周向平均分量，再回填到三维网格：

```text
完整算例: rho_3d -> 3D Poisson -> E_3d
控制算例: rho_3d -> m=0 charge filter -> rho_filtered -> 3D Poisson -> E_filtered
```

该功能的物理含义是人为切断周向非轴对称电荷扰动对静电场的反馈，用来量化 ECDI
及相关周向模对完整三维算例的贡献。它不是一个真实的无 ECDI 物理模型。

## 实现边界

1. 只处理三维静电 Poisson 路径，初期限定在 `WARPX_DIM_3D` 和
   `LabFrameExplicitES::ComputeSpaceChargeField`。
2. 不对求解后的电势 `phi_fp` 做二次平均或滤波。
3. 不修改粒子推进、碰撞、电离、注入和边界粒子处理逻辑。
4. 默认只过滤 `rho_fp`，让 Poisson 方程求解
   `L_h phi = -rho_filtered / epsilon_0`。
5. 若后续扩展到电磁求解或依赖电流连续性的算法，需要同时评估 `J` 的一致滤波。

## 建议代码结构

新增文件：

```text
Source/Insert/Fields/ECDIChargeFilter.h
Source/Insert/Fields/ECDIChargeFilter.cpp
```

修改文件：

```text
Source/Insert/CMakeLists.txt
Source/Insert/Config/WarpXFunctionConfig.h
Source/Insert/Core/WarpXInsert.h
Source/Insert/Core/WarpXInsert.cpp
Source/FieldSolver/ElectrostaticSolvers/LabFrameExplicitES.cpp
```

建议新增 hook：

```cpp
namespace Insert
{
void FilterRhoForECDIControl (
    ablastr::fields::MultiLevelScalarField const& rho_fp,
    int max_level);
}
```

建议接入位置在 `LabFrameExplicitES::ComputeSpaceChargeField` 中：

```text
mpc.DepositCharge(...)
mfl->DepositCharge(...)
warpx.SyncRho(...)
warpx.ApplyRhofieldBoundary(...)
Insert::FilterRhoForECDIControl(rho_fp, max_level)
warpx.ApplyRhofieldBoundary(...)
computePhi(rho_fp, phi_fp, ...)
computeE(Efield_fp, phi_fp, ...)
```

第二次 `ApplyRhofieldBoundary` 用于在过滤后重新同步物理边界和 guard cell。过滤算子
本身只应改写有效区域内的 `rho`。

## 运行时参数

建议使用运行时参数而不是只依赖编译期宏：

```text
insert.ecdi_control.enabled = 0
insert.ecdi_control.axis = z
insert.ecdi_control.center = 0.0 0.0
insert.ecdi_control.nr = -1
insert.ecdi_control.filter_interval = 1
insert.ecdi_control.diagnostics = 1
```

参数含义：

- `enabled`：是否启用控制算例过滤。
- `axis`：周向平均轴，初期只支持 `z`。
- `center`：旋转轴在横截面内的位置，`z` 轴时为 `(x_c, y_c)`。
- `nr`：径向 bin 数；`-1` 表示根据三维网格和问题区域自动生成。
- `filter_interval`：每隔多少次 Poisson 求解过滤一次；对照算例建议为 `1`。
- `diagnostics`：是否输出过滤强度和守恒误差。

编译期可加私有保护宏：

```cpp
// #define ECDI_CHARGE_FILTER
```

但实际开启仍由 `insert.ecdi_control.enabled` 控制。

## 连续数学定义

以 `z` 轴为推进轴和周向平均轴，三维 Cartesian 坐标到柱坐标的映射为

$$
r = \sqrt{(x-x_c)^2 + (y-y_c)^2},
\qquad
\theta = \operatorname{atan2}(y-y_c, x-x_c),
\qquad
z = z .
$$

完整三维电荷密度可写为 Fourier 展开：

$$
\rho(r,\theta,z)
= \sum_{m=-\infty}^{\infty} \rho_m(r,z) e^{i m \theta}.
$$

周向平均过滤只保留 `m=0` 分量：

$$
\bar{\rho}(r,z)
= \rho_0(r,z)
= \frac{1}{2\pi}\int_0^{2\pi}\rho(r,\theta,z)\,d\theta .
$$

回填到三维后：

$$
\rho_f(x,y,z) = \bar{\rho}\left(
\sqrt{(x-x_c)^2 + (y-y_c)^2}, z
\right).
$$

三维 Poisson 方程仍在原三维网格上求解：

$$
L_h \phi = -\frac{\rho_f}{\epsilon_0},
\qquad
\mathbf{E} = -G_h \phi .
$$

因此控制算例满足的是

$$
\nabla_h \cdot \mathbf{E} = \frac{\rho_f}{\epsilon_0},
$$

而不是

$$
\nabla_h \cdot \mathbf{E} = \frac{\rho}{\epsilon_0}.
$$

被人为切断反馈的源项为

$$
\rho_\perp = \rho - \rho_f .
$$

## 守恒离散投影算法

### 网格和权重

WarpX 的静电电荷密度 `rho_fp` 定义在节点上，因此过滤对象不是 cell-centered
密度，而是 nodal MultiFab 上的节点值。令三维有效节点的索引为 `p`，节点坐标为
`(x_p,y_p,z_p)`，节点控制体积为 `\Omega_p`，原始节点电荷密度为 `rho_p`。令二维
轴对称投影网格的 bin 索引为 `a = (i_r, i_z)`。三维节点对应的柱坐标为

$$
r_p = \sqrt{(x_p-x_c)^2 + (y_p-y_c)^2},
\qquad
z_p = z_p .
$$

节点控制体积 `\Omega_p` 应与 Poisson 源项所在的 nodal 离散空间一致。均匀 Cartesian
网格内部节点可取

$$
\Omega_p = \Delta x \Delta y \Delta z .
$$

物理边界节点按控制体积截断，例如单个方向边界节点为半控制体积，两个方向边界节点为
四分之一控制体积，角点为八分之一控制体积。实现时应优先复用或仿照 WarpX/AMReX
对 nodal 源项做积分诊断时的体积分权，避免把边界节点当作完整 cell 体积。

定义从三维节点到二维 bin 的插值权重

$$
W_{a p} = W_a(r_p,z_p),
\qquad
\sum_a W_{a p} = 1 .
$$

初期建议使用分片常数或双线性权重。分片常数更简单、更稳；双线性权重更平滑，但需要
更仔细处理边界 bin。

### 三维到二维：电荷守恒求和

先按电荷量而不是密度求和：

$$
Q_a = \sum_p \rho_p \Omega_p W_{a p}.
$$

同时计算同一组权重对应的离散节点控制体积：

$$
D_a = \sum_p \Omega_p W_{a p}.
$$

轴对称平均密度定义为

$$
\bar{\rho}_a =
\begin{cases}
Q_a / D_a, & D_a > 0, \\
0, & D_a = 0 .
\end{cases}
$$

这里 `D_a` 必须使用离散节点控制体积和，而不是直接使用解析环形体积。这样可以自然处理
nodal 边界权重、Cartesian 计算域不是完整圆柱、边界裁剪、并行分块和 AMR 层局部覆盖
带来的体积差异。如果计算域是完整圆柱，`D_a` 的径向求和会逼近解析环形体积

$$
V_a^{\mathrm{annulus}}
= \pi\left(r_{i+1/2}^2-r_{i-1/2}^2\right)\Delta z .
$$

### 二维到三维：守恒回填

使用相同权重回填：

$$
\rho^f_p = \sum_a \bar{\rho}_a W_{a p}.
$$

若使用分片常数权重，上式直接保证全局电荷守恒：

$$
\sum_p \rho^f_p \Omega_p
= \sum_p \Omega_p \sum_a \frac{Q_a}{D_a} W_{a p}
= \sum_a \frac{Q_a}{D_a} \sum_p \Omega_p W_{a p}
= \sum_a Q_a
= \sum_p \rho_p \Omega_p .
$$

若使用双线性权重，只要投影和回填使用同一组 `W_{a p}`，并且 `D_a` 按同一组权重计算，
上述守恒关系仍成立。实现中需要避免回填时改用另一套几何插值权重。

### 并行归约

每个 MPI rank 先在本地累积

$$
Q_a^{(rank)} = \sum_{p \in rank} \rho_p \Omega_p W_{a p},
\qquad
D_a^{(rank)} = \sum_{p \in rank} \Omega_p W_{a p}.
$$

然后做全局求和：

$$
Q_a = \sum_{rank} Q_a^{(rank)},
\qquad
D_a = \sum_{rank} D_a^{(rank)}.
$$

AMReX 实现上可先使用 host 端 `Vector<Real>` 保存 `Q` 和 `D`，通过
`ParallelDescriptor::ReduceRealSum` 做全局归约。若后续性能不足，再改成 GPU
端局部归约加 host 端全局归约。

实现遍历时必须使用 `rho_fp` 的 nodal valid box，而不是 cell-centered valid box。
节点坐标应由 nodal index 和 `Geometry` 的 `ProbLo`、`CellSize` 计算；边界节点的
`\Omega_p` 需要根据该节点在全局物理域边界上的位置修正。

### AMR 处理

第一阶段建议仅支持单层网格：

$$
max\_level = 0 .
$$

如果必须支持 AMR，应避免粗细层重复计数。可选策略：

1. 只对 finest composite rho 过滤，再按 WarpX 现有同步路径回写。
2. 按层过滤，但对被细层覆盖的粗层节点及其控制体积使用 mask 排除。
3. 对 `rho_fp` 和 `rho_cp` 分别处理，并保证 reflux/sync 后的源项一致。

建议第一版在 `max_level > 0` 时直接报错，等单层控制算例验证通过后再扩展。

## 算法步骤

每次静电场求解前执行：

1. 保存过滤前的总电荷和可选诊断范数：

   $$
   Q_{\mathrm{before}} = \sum_p \rho_p \Omega_p .
   $$

2. 对每个有效三维节点计算 `(r_p,z_p)`、节点控制体积 `\Omega_p` 和权重 `W_{a p}`。
3. 累积 `Q_a` 与 `D_a`。
4. 对 `Q_a` 和 `D_a` 做 MPI 全局归约。
5. 计算 `\bar{\rho}_a = Q_a / D_a`。
6. 用同一组权重回填得到 `rho^f_p`，覆盖 `rho_fp` 有效区域。
7. 重新执行 rho 物理边界和 guard cell 更新。
8. 计算过滤后总电荷：

   $$
   Q_{\mathrm{after}} = \sum_p \rho^f_p \Omega_p .
   $$

9. 输出守恒误差：

   $$
   \epsilon_Q =
   \frac{|Q_{\mathrm{after}}-Q_{\mathrm{before}}|}
        {\max(|Q_{\mathrm{before}}|, Q_{\mathrm{scale}})} .
   $$

10. 输出过滤强度：

    $$
    R_\rho =
    \frac{
    \left(\sum_p (\rho_p-\rho^f_p)^2 \Omega_p\right)^{1/2}
    }{
    \left(\sum_p \rho_p^2 \Omega_p\right)^{1/2}
    } .
    $$

其中 `Q_scale` 可取

$$
Q_{\mathrm{scale}} = \sum_p |\rho_p| \Omega_p .
$$

避免净电荷接近零时相对误差失真。

## 诊断输出

建议至少输出：

```text
insert_ecdi_filter_charge_before
insert_ecdi_filter_charge_after
insert_ecdi_filter_charge_abs_sum
insert_ecdi_filter_charge_relative_error
insert_ecdi_filter_rho_l2_ratio
```

如果需要更强的物理对照，应增加周向模谱诊断。对固定 `(r_i,z_j)` bin，定义

$$
\rho_m(r_i,z_j)
=
\frac{
\sum_p \rho_p \Omega_p W_{i j,p} e^{-i m \theta_p}
}{
\sum_p \Omega_p W_{i j,p}
}.
$$

模能量可定义为

$$
P_m =
\sum_{i,j} |\rho_m(r_i,z_j)|^2 D_{i j}.
$$

对照时重点观察 `m=1,2,...` 是否在完整算例中增长，而在控制算例中被静电反馈抑制。

## 验证计划

### 单元级验证

1. 均匀密度：

   $$
   \rho(x,y,z)=\rho_0
   $$

   过滤后应满足 `rho_f = rho_0`，总电荷误差接近机器精度。

2. 纯 `m=1` 扰动：

   $$
   \rho(r,\theta,z)=\rho_0 + A(r,z)\cos\theta
   $$

   过滤后应只剩 `rho_0`。

3. 纯 `m=2` 扰动：

   $$
   \rho(r,\theta,z)=\rho_0 + A(r,z)\cos(2\theta)
   $$

   过滤后应只剩 `rho_0`。

4. 守恒检查：

   $$
   \left|\sum_p \rho^f_p \Omega_p - \sum_p \rho_p \Omega_p\right|
   $$

   应在归约和浮点误差范围内。

### 算例级验证

使用同一个 checkpoint 分叉运行：

```text
common startup -> checkpoint -> full 3D branch
                          -> m=0 charge-filtered branch
```

两个分支保持粒子数、随机数种子、碰撞模型、电离模型、边界条件、时间步长和输出频率一致。
比较：

```text
electron temperature
electron axial/radial/azimuthal energy
J_e dot E
discharge current
ion energy distribution
wall current
rho_m and phi_m spectra
```

控制算例结论应表述为：完整算例与 `m=0` 源项过滤算例之间的差异来自周向非轴对称
静电场反馈，包括 ECDI 及相关周向模的作用。

## 风险和限制

1. 过滤的是 `rho -> E` 的反馈通道，不是所有 ECDI 相关物理机制。
2. 如果三维边界条件或几何本身不是轴对称，即使 `rho_f` 是轴对称，三维 Poisson 解也可能
   保留边界诱导的非轴对称场。
3. 若净电荷很小，必须使用绝对电荷和 L2 范数共同评估过滤强度。
4. AMR 支持需要额外设计，第一版应限制为单层网格。
5. 过滤过强会同时移除真实周向输运和电子加热，控制算例结果不能直接解释为实验中
   “无 ECDI”的真实状态。

## 建议实施顺序

1. 增加参数读取和空 hook，默认关闭。
2. 实现单层三维分片常数投影和守恒回填。
3. 在 `LabFrameExplicitES::ComputeSpaceChargeField` 的 Poisson 前接入 hook。
4. 增加总电荷守恒和 `R_rho` 诊断。
5. 用解析密度场做离线或小网格验证。
6. 用同一 checkpoint 分叉做完整算例和控制算例对照。
7. 在确认需求后再扩展双线性权重、AMR 和周向模谱诊断。
