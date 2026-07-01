# zmin 混合边界谱 Schur 修正实现方案

## 1. 目标

在不修改 WarpX/AMReX nodal Poisson 求解器的前提下，在 `Source/Insert`
中独立实现 zmin 面上的非均匀第二类边界修正。

当前约束是：

- zmin 面上的阳极 Dirichlet 电势已经由 WarpX 原静电求解流程处理；
- AMReX nodal Poisson 求解器只提供齐次 Neumann 形式，不能直接传入非齐次
  Neumann 通量；
- 该修正只针对无 EB、规则矩形域、3D lab-frame electrostatic 求解路径；
- 大部分实现应放在 `Source/Insert` 内，减少对 WarpX 主流程的侵入。

因此实现采用后处理叠加：

\[
\phi_{\text{total}} = \phi_{\text{warpx}} + \phi_{\text{corr}},
\]

\[
\nabla^2 \phi_{\text{corr}} = 0.
\]

其中 `phi_warpx` 已经包含体电荷、阳极电势和 WarpX 已有边界条件。
`phi_corr` 只补足 zmin 非阳极区域上的目标 Neumann 通量。

## 2. 修正问题定义

设 zmin 面为

\[
\Gamma_0 = \{z=z_{\min}\}.
\]

该面分为：

- `Gamma_D`：阳极 Dirichlet 区域，由 WarpX 原求解器负责；
- `Gamma_N`：需要施加非均匀 Neumann 通量的区域。

修正场满足：

\[
\phi_{\text{corr}} = 0
\quad \text{on } \Gamma_D,
\]

并在 `Gamma_N` 上满足：

\[
\partial_n \phi_{\text{corr}}
=
g_{\text{target}} - \partial_n \phi_{\text{warpx}}.
\]

若物理模型以表面电荷给出：

\[
\sigma_s = -\epsilon_0 \partial_n \phi,
\]

则：

\[
g_{\text{target}} = -\frac{\sigma_s}{\epsilon_0}.
\]

zmin 面外法向为：

\[
n=-\hat{z},
\qquad
\partial_n = -\partial_z.
\]

所以实现中如果用正 z 方向差分估计 \(\partial_z\phi\)，需要显式处理符号：

\[
\partial_n \phi_{\text{warpx}}
\approx
-\frac{\phi_{\text{warpx}}(k_0+1)-\phi_{\text{warpx}}(k_0)}{\Delta z}.
\]

## 3. 谱 Schur 方程的 zmin 形式

把文档 `spectral_boundary_schur_mixed_bc.md` 中的 `x=0` 推导换成 zmin 面后，
横向方向为 `x,y`，传播方向为 `z`。

对 zmin 面上的修正电势迹

\[
u(x,y) = \phi_{\text{corr}}(x,y,z_{\min})
\]

做二维 sine 展开。对每个模式：

\[
k_{mn}
=
\sqrt{
\left(\frac{m\pi}{L_x}\right)^2
+
\left(\frac{n\pi}{L_y}\right)^2
},
\]

\[
\lambda_{mn} = k_{mn}\coth(k_{mn}L_z).
\]

Dirichlet-to-Neumann 算子为：

\[
\Lambda
=
\operatorname{IDST}
\operatorname{diag}(\lambda_{mn})
\operatorname{DST}.
\]

在 Neumann 区域上的 Schur 方程为：

\[
P_N\Lambda P_N^T u_N
=
g_N,
\]

其中：

\[
g_N
=
\left(
g_{\text{target}} - \partial_n \phi_{\text{warpx}}
\right)\big|_{\Gamma_N}.
\]

若以表面电荷输入：

\[
g_N =
-\frac{\sigma_s}{\epsilon_0}
- \partial_n \phi_{\text{warpx}}.
\]

阳极区域不进入未知量，扩展到全 zmin 面时填零：

\[
u = P_N^T u_N,
\qquad
u|_{\Gamma_D}=0.
\]

## 4. 推荐代码结构

新增一个独立模块：

```text
Source/Insert/Fields/SpectralBoundarySchur.cpp
Source/Insert/Fields/SpectralBoundarySchur.h
```

在 `Source/Insert/CMakeLists.txt` 中加入该 `.cpp`。

建议命名空间：

```cpp
namespace Insert {
namespace SpectralBoundarySchur {
    void Initialize();
    void ApplyZMinCorrection();
    void Finalize();
}
}
```

对外入口保持很窄：

```cpp
namespace Insert {
    void ApplyElectrostaticBoundaryCorrection();
}
```

第一版入口只负责从 zmin 面表面电荷密度求解边界未知量 `u_N`。后续实现体内
重建后，主流程再在 WarpX 原 Poisson 解和原电场计算之后调用该入口，作为
`Efield_fp += E_corr` 的后处理。

如果后续确实需要诊断 `phi_corr`，可以额外存储一个 Insert 私有 `MultiFab`。
第一版不生成 `phi_corr`，也不回写 `phi_fp`，避免影响 WarpX 后续 guard cell
与 Dirichlet 处理逻辑。

## 5. 数据和配置

建议新增输入开关：

```text
insert.schur_boundary.enabled = 0/1
insert.schur_boundary.face = zmin
insert.schur_boundary.max_iter = 100
insert.schur_boundary.rel_tol = 1e-10
insert.schur_boundary.abs_tol = 0.0
insert.schur_boundary.rebuild_volume_field = 0/1
```

Schur 算子使用全横向 sine 模态。每个方向的模态数固定为该方向参与 DST 的
节点个数：

```text
num_modes_x = nx_internal
num_modes_y = ny_internal
```

其中 `nx_internal` 和 `ny_internal` 是进入 sine transform 的节点数。如果 zmin
面数组包含横向 Dirichlet 端点，则端点不进入 DST，自由节点数为
`nx_internal = nx_face - 2`、`ny_internal = ny_face - 2`；如果调用方已经只传入
内部面节点，则 `nx_internal = nx_face`、`ny_internal = ny_face`。第一版不提供
运行时模态截断参数。
第一版不实现体内重建，`rebuild_volume_field` 固定为 `0`；后续实现
`phi_corr`/`E_corr` 后再启用该开关。

第一版只支持：

- `face = zmin`
- `WarpX_DIMS = 3`
- 无 EB
- 单层 level 0
- 规则矩形、非 AMR 修正

第一版接口先假设外部传入一个 zmin 面上的
`amrex::Gpu::DeviceVector<amrex::Real>`，其中 x 方向快变：

```text
sigma_s_device[i + nx_face * j] = sigma_s(i,j)
```

`nx_face` 和 `ny_face` 必须与 Schur 模块使用的 zmin 面节点布局一致。数组值为
已经沉积完成的表面电荷密度，单位按
\(\sigma_s = -\epsilon_0 \partial_n\phi\) 解释。

表面通量来源可以分阶段实现：

1. 常数或解析 parser：便于验证；
2. zmin 面 `DeviceVector`：便于接入外部提供的表面电荷密度；
3. WarpX 当前尚未实现壁面电荷沉积；待该能力实现后，Schur 模块不负责沉积，
   而是直接接收已经沉积好的 \(\sigma_s(x,y)\) 表面电荷密度。

## 6. 面 mask 构造

复用当前阳极判断逻辑，例如 `InsertBoundaryPhi.cpp` 中的
`IsHallAnodeRingNode` 所表达的 zmin 阳极区域。

构造二维面 mask。第一版直接使用一个二维 `bool` 数组表示自由 Neumann
节点，避免在 Schur 模块内部重复维护多套 mask：

```text
NeumannMask(i,j) = true   if interior zmin node belongs to Gamma_N
NeumannMask(i,j) = false  if node belongs to Gamma_D or transverse boundary endpoint
```

注意事项：

- sine transform 只作用在横向 Dirichlet 端点内部点；
- `x=xmin/xmax` 与 `y=ymin/ymax` 的端点如果按接地 Dirichlet 处理，不进入
  DST 自由度；
- 阳极节点对应 `NeumannMask = false`，并在扩展到全 zmin 面时保持
  `phi_corr = 0`；
- `NeumannMask` 不应覆盖阳极节点，否则会重复处理阳极边界；
- 如果后续调用方已经能提供 zmin 面 bool mask，则 Schur 模块直接使用该数组；
  否则由共享阳极几何函数在初始化阶段生成一次。

建议将阳极几何判断抽出为可复用函数，避免在
`InsertBoundaryPhi.cpp` 和 Schur 模块中复制半径、电压、中心等参数。

## 7. 每步执行流程

### 7.1 WarpX 原求解

保持当前流程：

```text
1. WarpX deposit rho
2. Insert::SetBoundaryPhi() 设置阳极电势
3. Insert::BuildPhiOversetMasks() 固定阳极 Dirichlet 节点
4. WarpX/AMReX Poisson 求解 phi_warpx
5. WarpX 根据 phi_warpx 计算 Efield_fp
```

这些步骤已经处理 zmin 阳极电势，不在 Schur 模块中重复。

### 7.2 构造 Schur RHS

在 zmin 的 `Gamma_N` 上计算：

```text
normal_grad_warpx = -(phi_warpx(i,j,k0+1) - phi_warpx(i,j,k0)) / dz
g_rhs = g_target(i,j) - normal_grad_warpx
```

若输入为表面电荷：

```text
g_target = -sigma_s(i,j) / epsilon0
```

其中 `sigma_s(i,j)` 来自 `ZMinWallChargeDeposit` 已累计的 zmin 面壁面电荷。

第一版实际接入接口按 x 快变的 device 数组读取：

```text
sigma_s(i,j) = sigma_s_device[i + nx_face * j]
```

只将 `Gamma_N` 上的 `g_rhs` 放入 CG 右端项。

### 7.3 CG 求解边界未知量

Schur 算子应用：

```text
apply_A(v_N):
    face = 0
    face[NeumannMask] = v_N

    face_hat = DST2(face)
    q_hat[m,n] = lambda[m,n] * face_hat[m,n]
    q_face = IDST2(q_hat)

    return q_face[NeumannMask]
```

求解：

```text
P_N Lambda P_N^T u_N = g_rhs
```

这里不乘 \(\epsilon_0\)，因为 RHS 已经是法向导数单位。如果直接用
\(\sigma_s\) 作为 RHS，则算子与 RHS 需要统一改写为
\(\epsilon_0\Lambda u = -\sigma_s - \epsilon_0\partial_n\phi_{\text{warpx}}\)。

这里使用全模态 Schur 算子，不做模态截断。这样 \(P_N\Lambda P_N^T\) 保持完整
边界 Schur 系统的对称正定结构。

### 7.4 重建修正场

第一版暂不实现体内重建，只求解并保留 zmin 边界未知量 `u_N`，可选输出
`u_face` 或 `u_hat` 供后续阶段使用。本节保留为后续实现方案。

把 `u_N` 扩展为全 zmin 面：

```text
u_face = 0
u_face[NeumannMask] = u_N
```

做 DST 得到 `u_hat[m,n]`。

对每个 z 层：

\[
R_{mn}(z)
=
\frac{\sinh(k_{mn}(L_z - (z-z_{\min})))}
{\sinh(k_{mn}L_z)}.
\]

后续体内重建采用分层模态截断和矩阵乘法实现，不再默认对每个 z 层做完整
`IDST2`。Schur 边界方程仍使用全模态；下面的截断只用于把已经求得的 `u_hat`
近似重建到体内。

对于 `Nz` 个 z 层，按距离 zmin 的归一化深度分层：

```text
0        <= z/Lz < 1/16 : Mx = nx_internal,     My = ny_internal
1/16     <= z/Lz < 1/8  : Mx = nx_internal / 2, My = ny_internal / 2
1/8      <= z/Lz < 1/4  : Mx = nx_internal / 4, My = ny_internal / 4
1/4      <= z/Lz < 1/2  : Mx = nx_internal / 8, My = ny_internal / 8
1/2      <= z/Lz <= 1   : Mx = nx_internal / 16, My = ny_internal / 16
```

实际实现中 `Mx` 和 `My` 至少为 1，并向下取整到合法模态数。分层边界处不做
谱滤波或过渡平滑，直接切换保留的模态数。

对每个 z 层，先形成保留模态上的谱系数：

```text
A[m,n,z] = R[m,n,z] * u_hat[m,n]
```

然后利用 sine 基函数的可分离性，用两次矩阵乘法重建物理空间电势：

```text
tmp[i,n,z]      = sum_m Sx[i,m] * A[m,n,z]
phi_corr[i,j,z] = sum_n tmp[i,n,z] * Sy[j,n]
```

其中：

```text
Sx[i,m] = sin(m * pi * x_i / Lx)
Sy[j,n] = sin(n * pi * y_j / Ly)
```

推荐实现为按 z 层或按分层区间的 batched GEMM/tiled matrix multiply。对于
`Nx = Ny = Nz = 256` 且使用上面的分层规则，矩阵乘法重建的数量级约为
`2.4e9` FLOP，显著低于直接 `ijk` 双重模态求和的约 `4.0e11` FLOP。

电场修正：

```text
E_corr = -grad(phi_corr)
Efield_fp += E_corr
```

后续实现可以先重建 `phi_corr` 到临时 nodal `MultiFab`，然后复用与
`phi_fp -> Efield_fp` 一致的有限差分方式生成 `E_corr`，这样更容易对齐 WarpX
场布局。优化版再考虑直接在谱空间重建 `E_corr`。

## 8. Transform 实现策略

优先选择已有依赖，避免新增外部库。

候选方案：

1. 使用 AMReX FFT 接口实现 DST 的奇延拓版本；
2. 若当前构建没有启用 `WarpX_FFT`，第一版可先实现 CPU-only 验证路径；
3. 后续再补 GPU/多 MPI rank 的并行面 transform。

第一版建议明确限制并在运行时检查：

- zmin 面是否完整落在单进程或已经 gather 到 IOProcessor；
- transform 大小是否等于内部节点数；
- 并行运行时是否需要 broadcast `u_N` 或 `E_corr`。

如果目标算例必须 MPI 并行，推荐实现面数据 gather/scatter：

```text
MultiFab zmin face -> host contiguous face array on IOProcessor
IOProcessor runs CG/DST
IOProcessor reconstructs needed correction
scatter correction back to Efield_fp owners
```

该方案实现简单，但体内重建和 scatter 成本较高。后续可替换为分布式
AMReX FFT。

## 9. 数值稳定和预处理

预处理阶段缓存：

```text
k[m,n]
lambda[m,n] = k[m,n] * coth(k[m,n] * Lz)
NeumannMask
```

后续体内重建阶段缓存：

```text
Sx[i,m] = sin(m * pi * x_i / Lx)
Sy[j,n] = sin(n * pi * y_j / Ly)
LayerSchedule(z) -> (Mx, My)
可选 R[m,n,k] 和 dRdz[m,n,k]
```

如果 `R[m,n,k]` 内存成本过高，可以在生成 `A[m,n,z]` 时按需计算。

计算 `coth` 时避免直接使用大参数 `sinh/cosh`：

```text
if a > threshold:
    coth(a) = 1
else:
    coth(a) = (1 + exp(-2a)) / (1 - exp(-2a))
```

第一版 CG 预处理可用对角近似。若收敛慢，再实现整面 DtN 逆作为近似预处理：

```text
M^{-1} r = IDST2( DST2(r_face) / lambda[m,n] )
```

mask 会导致模式耦合，所以该预处理只是近似。

实现上只需要保存一个 `NeumannMask` bool 数组。所有不属于 `Gamma_N` 的点，
包括阳极 Dirichlet 节点和横向端点，都通过 `NeumannMask = false` 隐式处理为
`phi_corr = 0`。

## 10. 与现有 Insert 代码的关系

建议调整：

- 将阳极几何配置从 `InsertBoundaryPhi.cpp` 的匿名命名空间迁出；
- 新增 `Boundary/HallAnodeGeometry.h/.cpp` 或等价小工具；
- `AnodeVoltage()`、`BuildPhiOversetMasks()`、Schur mask 构造共用同一判断函数；
- 新增 Schur 模块只依赖该几何工具和 WarpX 场数据。

不建议调整：

- 不修改 `PoissonSolver.H` 中 AMReX linop 构造；
- 不修改 AMReX `MLNodeTensorLaplacian`；
- 不在 `PoissonBoundaryHandler` 中伪造非齐次 Neumann；
- 不让 Schur 模块重新设置 zmin 阳极电势。

## 11. 验证计划

### 11.1 单独数学验证

构造解析解：

\[
\phi(x,y,z) =
\sin\frac{\pi x}{L_x}
\sin\frac{\pi y}{L_y}
\frac{\sinh(k(L_z-z))}
{\sinh(kL_z)}
\]

其中：

\[
k = \sqrt{(\pi/L_x)^2 + (\pi/L_y)^2}.
\]

在 zmin 给定对应 Neumann 通量，第一版检查 Schur 求得的 `u_face` 是否与解析
边界迹一致。体内 `phi_corr` 重建留到后续阶段验证。

### 11.2 WarpX 耦合验证

使用无粒子或固定 \(\rho=0\) 算例：

```text
phi_warpx = 阳极 Dirichlet 解
u_N = Schur 求得的边界未知量
```

检查：

- 阳极节点不进入 `u_N`，扩展到 `u_face` 时为 0；
- zmin 非阳极区域的 Schur RHS 和 `A u_N` 残差满足容差；
- 传入的 `sigma_s_device[i + nx_face * j]` 索引与面节点布局一致；
- 关闭 Schur 开关时结果与当前 WarpX 完全一致。

### 11.3 回归风险检查

需要覆盖：

- `insert.schur_boundary.enabled = 0` 的默认路径；
- `HALL3D` 下阳极 mask 与 overset mask 一致；
- 非 3D、EB、AMR、多 level 情况下给出明确报错或自动禁用；
- MPI 运行模式下的数据 gather/scatter 一致性。

## 12. 分阶段实施

### 阶段一：最小可验证版本

- 新建 Schur 模块和入口；
- 只支持 3D、level 0、无 EB；
- 输入为 x 快变的 zmin 面 `amrex::Gpu::DeviceVector<amrex::Real>` 表面电荷密度；
- Schur 入口只读取已累计的 zmin 壁面电荷作为 `sigma_s` 数据源；
- CG + CPU DST；在 `cpu_serial` 后端中，可先把 device RHS 拷贝到 host 后求解；
- 求解到边界未知量 `u_N` 为止，暂不重建 `phi_corr`，暂不叠加 `Efield_fp`；
- 每个方向的模态数固定为该方向参与 DST 的节点个数，不提供运行时截断；
- 写一个小规模输入算例做符号和边界未知量验证。

### 阶段二：接入实际表面模型

- 从外部壁面电荷沉积/表面模型读取已经沉积完成的 \(\sigma_s(x,y)\)；
- 直接使用 zmin 面 `NeumannMask` bool 数组；
- 增加诊断输出：`u_N`、`u_face`、`u_hat`、RHS 残差、CG 迭代数；
- 支持运行时可调容差和最大迭代数。

### 阶段三：性能和并行化

- 实现体内 `phi_corr`/`E_corr` 重建并叠加到 `Efield_fp`；
- 替换 CPU serial DST 为 AMReX FFT 或分布式 transform；
- 用分层模态截断和两步矩阵乘法/GEMM 重建 `phi_corr`；
- 分层边界直接切换模态数，不做谱滤波或过渡平滑；
- 预存 sine 基矩阵、分层表和必要的传播因子，或按需计算传播因子；
- 后续可避免完整体内 `phi_corr` 存储，直接生成 `E_corr`；
- 增加预处理器降低 CG 迭代数。

## 13. 主要风险

- 文档推导假设其余边界为接地 Dirichlet；若实际 WarpX 边界不是该形式，
  Schur 修正与原问题会有模型误差；
- WarpX 原解中的 \(\partial_n\phi_{\text{warpx}}\) 必须从同一网格布局上计算，
  否则 RHS 会有半格偏移或符号错误；
- zmin 面端点节点、阳极环边缘节点和 DST 内部自由度的对应关系必须统一；
- `sigma_s_device` 使用 x 快变布局，任何调用方传入前都必须保证和
  `NeumannMask(i,j)` 使用同一面索引；
- 分层重建阶段不做边界滤波，模态数硬切换可能在分层边界附近引入小的
  `phi_corr` 或 `E_corr` 层间不连续；
- 若并行时面数据跨多个 rank，CPU gather/scatter 版本会成为瓶颈；
- 如果回写 `phi_fp`，可能影响 WarpX 后续 guard cell 和诊断；第一版不回写
  `phi_fp`，也不叠加 `Efield_fp`。

## 14. 推荐接入点

第一版只求解边界未知量，可在 WarpX 原 Poisson 解之后、需要诊断或后续模块读取
`u_N` 时调用。后续实现 `E_corr` 重建后，推荐在
`LabFrameExplicitES::ComputeSpaceChargeField` 中形成如下流程：

```text
computePhi(...)
Insert::SetPhiGuards()
WarpX 原逻辑由 phi_fp 计算 Efield_fp
Insert::ApplyElectrostaticBoundaryCorrection()
```

当前 `LabFrameExplicitES.cpp` 中 `SetPhiGuards()` 已经紧跟 `computePhi(...)`
执行，随后才调用 `computeE(...)`。因此体内重建实现后，推荐在 `computeE(...)`
之后调用 Schur 修正入口，只叠加 `Efield_fp`，不改变 `phi_fp` 和 guard cell。
核心原则是：

- Schur 不参与 AMReX MLMG 求解；
- Schur 不覆盖阳极 Dirichlet 电势；
- 第一版 Schur 只求解边界未知量；
- 后续 Schur 只补充总电场，使 zmin 非阳极区域满足目标非齐次 Neumann 条件。
