# ECDI 电荷密度轴对称过滤算法规格（WarpX 实现参考）

> 本文档从 `ChargeTest` 离线验证代码提炼，给出可在 WarpX 中直接复现的离散算法、公式、约束与迁移映射。阅读前应先了解 `ECDI_SUPPRESSION_PLAN.md` 的物理背景。
>
> 当前实现范围：三维、单层网格（`max_level = 0`）、单 MPI rank。算法本身已通过独立离线验证；本文重点约束 WarpX 迁移时的离散索引、边界和接口语义。

---

## 1. 算法目标

在三维静电 Particle-in-Cell（PIC）模拟中，于每次 Poisson 求解前，将三维 nodal 电荷密度 `rho` 投影到轴对称的周向平均分量（`m=0`），再回填到三维网格，从而人为切断周向非轴对称电荷扰动对静电场的反馈。

---

## 2. 数据布局与符号约定

### 2.1 三维网格（Cartesian，Nodal）

| 符号 | 含义 | 备注 |
|------|------|------|
| `nx, ny, nz` | 各方向 cell 数 | — |
| `dx, dy, dz` | 均匀网格步长 | > 0 |
| `prob_lo[3]` | 物理域下角 | — |
| `prob_hi[3]` | 物理域上角 | — |
| `n_nodes_x = nx + 1` | x 方向节点数 | nodal 布局 |
| `n_nodes_y = ny + 1` | y 方向节点数 | nodal 布局 |
| `n_nodes_z = nz + 1` | z 方向节点数 | nodal 布局 |
| `N = n_nodes_x * n_nodes_y * n_nodes_z` | 总节点数 | `long` 型防溢出 |
| `p` | 三维节点展平索引 | `p = ((k * n_nodes_y + j) * n_nodes_x + i)`，使用 `long` |
| `(x_p, y_p, z_p)` | 节点坐标 | `x(i) = prob_lo[0] + i * dx` |
| `rho_p` | 节点电荷密度 | 输入 `rho_in[p]` |
| `Omega_p` | 节点控制体积 | 见第 3.1 节 |
| `r_p` | 节点到旋转轴距离 | `sqrt((x_p - x_c)^2 + (y_p - y_c)^2)` |

### 2.2 二维轴对称投影网格（径向 bin + z 层）

| 符号 | 含义 | 备注 |
|------|------|------|
| `nr` | 径向 bin 数 | >= 0 |
| `dr` | 径向步长 | > 0（若 `nr > 0`） |
| `(x_c, y_c)` | 旋转轴中心 | 由运行时参数指定 |
| `a = (i_r, k)` | 二维 bin 索引 | `i_r = 0 ... nr`，`k = 0 ... n_nodes_z-1` |
| `bar_rho_a` | 轴对称平均密度 | 投影结果，长度 `(nr+1) * n_nodes_z` |
| `W_{a,p}` | 从节点 `p` 到 bin `a` 的权重 | 满足 `sum_a W_{a,p} = 1` |

### 2.3 径向网格范围

径向二维网格必须覆盖三维横截面内的最远 nodal 节点。自动生成 `nr` 时：

```
r_domain_max = max over all nodal (i,j) of sqrt((x(i)-x_c)^2 + (y(j)-y_c)^2)
nr = ceil(r_domain_max / dr)
rmax = nr * dr
```

因此最远节点（通常是横截面边界上的角/棱节点）满足 `r <= rmax`，落在二维径向网格内部。`compute_weights` 中的 `r >= rmax` 分支仅作为浮点舍入和手动参数不匹配时的安全保护，不应成为自动参数下的常规路径。

---

## 3. 核心子算法

### 3.1 节点控制体积 `computeNodeControlVolumes`

控制体积必须与 nodal 离散一致，物理边界需截断：

```
dV = dx * dy * dz

for k = 0 ... n_nodes_z-1:
    fz = 1.0; if (k == 0 || k == nz) fz *= 0.5

    for j = 0 ... n_nodes_y-1:
        fy = 1.0; if (j == 0 || j == ny) fy *= 0.5

        for i = 0 ... n_nodes_x-1:
            fx = 1.0; if (i == 0 || i == nx) fx *= 0.5

            p  = ((k * n_nodes_y + j) * n_nodes_x + i)
            Omega_p = dV * fx * fy * fz
```

**边界规则**：
- 面边界节点：`0.5 * dV`
- 边边界节点：`0.25 * dV`
- 角边界节点：`0.125 * dV`
- 内部节点：`1.0 * dV`

### 3.2 线性径向权重 `compute_weights(r)`

使用一维线性 hat function（linear interpolation kernel），返回长度为 `nr+1` 的向量 `w`，满足 `sum(w) == 1`。

**算法**：

```
if (nr == 0):
    w[0] = 1.0
    return w

rmax = nr * dr

if (r >= rmax):
    w[nr] = 1.0
    return w

i_left  = floor(r / dr)   // clamp to [0, nr-1]
i_right = i_left + 1
r_left  = i_left * dr
r_right = i_right * dr
denom   = r_right - r_left   // = dr，除非退化

if (denom < 1e-15):
    w[i_left] = 1.0
    return w

w[i_right] = (r - r_left) / denom
w[i_left]  = (r_right - r) / denom
return w
```

**关键性质**：
- 对任意 `r >= 0`，至多两个相邻分量非零。
- `sum_{i_r} w[i_r] == 1`（严格成立）。
- `r >= rmax` 时边界归一到 `w[nr] = 1`，防止越界。
- **投影与回填必须使用同一套权重函数**；这是电荷守恒的充要条件。

---

## 4. 主过滤算法 `filterChargeDensity`

### 4.1 输入输出

```cpp
void filterChargeDensity(
    const Geometry3D& geom,
    const CylindricalBins& bins,
    const std::vector<double>& omega,   // Omega_p
    const std::vector<double>& rho_in,  // rho_p
    std::vector<double>& rho_out,       // rho^f_p（输出）
    FilterDiagnostics* diag = nullptr
);
```

### 4.2 步骤 1：可选前置诊断

若 `diag != nullptr`，遍历所有节点计算：

```
q_before = sum_p rho_p * Omega_p
q_scale  = sum_p |rho_p| * Omega_p
l2_before = sum_p rho_p^2 * Omega_p
```

- `q_scale` 用于净电荷接近零时的相对误差归一化，避免失真。

### 4.3 步骤 2：投影（3D → 2D）

**按 `z = k` 层独立处理**（z 方向完全重合，无需插值）。

对每一层 `k = 0 ... n_nodes_z - 1`：

```
初始化 Q[0...nr] = 0.0
初始化 D[0...nr] = 0.0

for j = 0 ... n_nodes_y - 1:
    y = geom.y(j) - bins.center[1]
    for i = 0 ... n_nodes_x - 1:
        x = geom.x(i) - bins.center[0]
        r = sqrt(x*x + y*y)

        weights = bins.compute_weights(r)   // w[0...nr]

        p       = ((k * n_nodes_y + j) * n_nodes_x + i)
        rho_val = rho_in[p]
        vol     = omega[p]

        for ir = 0 ... nr:
            Q[ir] += rho_val * vol * weights[ir]
            D[ir] += vol * weights[ir]

for ir = 0 ... nr:
    a = k * (nr + 1) + ir
    if (D[ir] > 0.0):
        bar_rho[a] = Q[ir] / D[ir]
    else:
        bar_rho[a] = 0.0
        empty_bins++
```

**要点**：
- `Q[ir]` 是加权电荷量之和；`D[ir]` 是加权控制体积之和。
- 必须使用 `D[ir]` 做归一，而不是解析环形体积。这自动兼容非圆柱域和边界裁剪。
- `D[ir] <= 0` 时置 `bar_rho = 0`（除零保护）。

### 4.4 步骤 3：回填（2D → 3D）

使用**与投影完全相同的权重函数** `compute_weights(r)`：

```
for k = 0 ... n_nodes_z - 1:
    for j = 0 ... n_nodes_y - 1:
        y = geom.y(j) - bins.center[1]
        for i = 0 ... n_nodes_x - 1:
            x = geom.x(i) - bins.center[0]
            r = sqrt(x*x + y*y)

            weights = bins.compute_weights(r)

            p   = ((k * n_nodes_y + j) * n_nodes_x + i)
            val = 0.0
            for ir = 0 ... nr:
                a   = k * (nr + 1) + ir
                val += bar_rho[a] * weights[ir]
            rho_out[p] = val
```

### 4.5 步骤 4：可选后置诊断

若 `diag != nullptr`：

```
q_after  = sum_p rho_out[p] * Omega_p
l2_diff  = sum_p (rho_in[p] - rho_out[p])^2 * Omega_p
max_abs  = max_p |rho_in[p] - rho_out[p]|
max_rel  = max_p |rho_in[p] - rho_out[p]| / max(|rho_in[p]|, 1e-30)

charge_relative_error = |q_after - q_before| / q_scale   (若 q_scale > 0)
rho_l2_ratio          = sqrt(l2_diff / l2_before)         (若 l2_before > 0)
```

---

## 5. 守恒证明

离散全局电荷守恒严格成立，前提是：

1. 投影和回填使用同一套权重 `W_{a,p}`。
2. `bar_rho_a = Q_a / D_a`，其中 `D_a = sum_p Omega_p W_{a,p}`。

**证明**：

```
Q_after = sum_p rho^f_p * Omega_p
        = sum_p Omega_p * sum_a bar_rho_a * W_{a,p}
        = sum_a bar_rho_a * sum_p Omega_p * W_{a,p}
        = sum_a bar_rho_a * D_a
        = sum_a Q_a
        = sum_a sum_p rho_p * Omega_p * W_{a,p}
        = sum_p rho_p * Omega_p * sum_a W_{a,p}
        = sum_p rho_p * Omega_p
        = Q_before
```

**WarpX 实现约束**：
- 回填时严禁使用另一套几何插值（如双线性/三线性插值与投影的权重不一致）。
- 节点控制体积 `Omega_p` 在投影和回填诊断中必须一致。

---

## 6. 数值安全与边界处理

| 场景 | 处理策略 |
|------|----------|
| `D_a <= 0`（空 bin） | `bar_rho_a = 0.0`，计数 `empty_bins` |
| `r >= rmax` | 权重归一到最外层 bin `w[nr] = 1.0` |
| `dr` 极小或退化 | 若 `denom < 1e-15`，`w[i_left] = 1.0` |
| 净电荷接近零 | 相对误差分母用 `q_scale = sum \|rho_p\| Omega_p` |
| 边界节点体积 | 面×0.5、边×0.25、角×0.125 |

---

## 7. WarpX / AMReX 迁移映射

以下将离线验证代码中的 `std::vector` 实现映射到 WarpX 的 `MultiFab` 体系。当前仅考虑单层网格和单 MPI rank；因此不需要跨 rank 的全局归约。

| 离线代码 | WarpX / AMReX 对应 | 备注 |
|----------|-------------------|------|
| `Geometry3D` | `warpx.Geom(lev)` / `amrex::Geometry` | 使用 `ProbLo()`、`ProbHi()`、`CellSize()`、`Domain()` |
| `std::vector<double> rho_in/out` | `ablastr::fields::MultiLevelScalarField rho_fp` | 单层时使用 `rho_fp[0]`，类型为 `amrex::MultiFab*` |
| `std::vector<double> omega` | 现场计算或预计算为临时 `MultiFab` | 控制体积只需几何信息，可每层算一次 |
| 展平索引 `p` | `amrex::MFIter` + `amrex::Box` | 遍历 `rho_fp[0]` 的 nodal valid/tile box |
| 节点坐标 `(x_p, y_p)` | `ProbLo + i * dx` | 不使用 `CellCenter`；节点坐标由 nodal index 直接计算 |
| `compute_weights(r)` | host/device 可调函数 | 建议 `AMREX_GPU_HOST_DEVICE` 内联；逻辑极简单，可直接展开 |
| `Q_a, D_a` 局部累加 | 每层 `amrex::Vector<amrex::Real>` | 大小 `(nr+1) * n_nodes_z` |
| MPI 全局归约 | 暂不需要 | 当前单 MPI rank；未来扩展 MPI 时再加入 `ParallelDescriptor::ReduceRealSum` |
| `bar_rho` | `amrex::Vector<amrex::Real>`（host）或临时 `MultiFab`（2D） | 数据量小，建议 host vector |
| 回填写入 `rho_out` | `amrex::MFIter` 遍历 + `fab(i,j,k,n)` 赋值 | 只写 valid region，不直接写 guard cell；第一版要求或假定 `nComp()==1` |
| guard cell 填充 | `rho_fp[lev]->FillBoundary(geom.periodicity())` 或与当前 rho 同步路径一致的本地填充 | `ApplyRhofieldBoundary` 不负责一般 guard cell 填充 |
| 物理边界处理 | `warpx.ApplyRhofieldBoundary(lev, rho_fp[lev], PatchType::fine)` | 仅处理反射/PEC/PECInsulator/PMC 等 rho 边界逻辑 |
| 诊断输出 | `amrex::Print()` 或 `WarpX::RecordWarning` | 建议每步或按 `filter_interval` 输出 |

### 7.1 建议接入点

在 `LabFrameExplicitES::ComputeSpaceChargeField` 中：

```cpp
const MultiLevelScalarField rho_fp = fields.get_mr_levels(FieldType::rho_fp, max_level);
const MultiLevelScalarField rho_cp = fields.get_mr_levels(
    FieldType::rho_cp, max_level, skip_lev0_coarse_patch);

mpc.DepositCharge(rho_fp, 0.0_rt);
if (mfl) {
    const int lev = 0;
    mfl->DepositCharge(fields, *rho_fp[lev], lev);
}

const Vector<std::unique_ptr<MultiFab>> rho_buf(num_levels);
warpx.SyncRho(rho_fp, rho_cp, amrex::GetVecOfPtrs(rho_buf));

for (int lev = 0; lev < num_levels; lev++) {
    warpx.ApplyRhofieldBoundary(lev, rho_fp[lev], PatchType::fine);
}

if (insert_ecdi_control_enabled) {
    Insert::FilterRhoForECDIControl(rho_fp, max_level);
    rho_fp[0]->FillBoundary(warpx.Geom(0).periodicity());
    warpx.ApplyRhofieldBoundary(0, rho_fp[0], PatchType::fine);
}

computePhi(rho_fp, phi_fp, ...)
computeE(Efield_fp, phi_fp, ...)
```

过滤算子只覆盖 `rho_fp[0]` 的有效节点区域。`FillBoundary` 用于刷新 guard cells；`ApplyRhofieldBoundary` 用于重新施加 WarpX 源码中已有的 rho 物理边界处理。

### 7.2 AMR 与 MPI 限制

- **第一版仅支持 `max_level = 0`**。
- **第一版仅支持单 MPI rank**，不做 `Q`/`D` 的跨 rank 归约。
- 若未来扩展 AMR，需处理粗细层重复计数：finest composite 过滤后同步，或对粗层被覆盖节点使用 mask 排除。
- 若未来扩展 MPI 或多 rank 分解，需增加跨 rank 全局归约；若 nodal box 在分块边界存在共享节点，还需使用唯一节点所有权策略，避免投影和诊断重复计数。

---

## 8. 运行时参数建议

```ini
insert.ecdi_control.enabled = 0          # 默认关闭
insert.ecdi_control.axis = z             # 旋转轴，初期仅支持 z
insert.ecdi_control.center = 0.0 0.0     # (x_c, y_c)
insert.ecdi_control.nr = -1              # 径向 bin 数；-1 表示自动
insert.ecdi_control.dr = -1.0            # 径向步长；-1 表示自动（建议 min(dx, dy)）
insert.ecdi_control.filter_interval = 1  # 每隔多少步过滤一次
insert.ecdi_control.diagnostics = 1      # 是否输出守恒误差与过滤强度
```

自动参数规则：
- 若 `dr < 0`，取 `dr = min(dx, dy)`。
- 若 `nr < 0`，先计算所有横截面 nodal 节点到过滤轴的最大距离 `r_domain_max`，再取 `nr = ceil(r_domain_max / dr)`。
- 若用户显式给定 `nr` 和 `dr`，必须检查 `nr * dr >= r_domain_max`；否则应报错，而不是静默把外侧节点压到最外层 bin。

---

## 9. 验证基准（ChargeTest 单元测试）

以下 9 项测试应在 WarpX 实现完成后复现，作为回归基准：

| 测试名 | 验证性质 | 阈值建议 |
|--------|----------|----------|
| `weight_normalization` | 对所有节点，`sum(weights) == 1` | `TOL = 1e-14` |
| `transpose_symmetry` | 投影权重与回填权重数值相同 | `TOL = 1e-14` |
| `uniform_field` | 常数场过滤后不变，守恒误差机器精度 | `TOL = 1e-14` |
| `m1_perturbation` | `rho = rho0 + A*cos(theta)` → `rho0` | `TOL = 1e-14` |
| `m2_perturbation` | `rho = rho0 + A*cos(2*theta)` → `rho0` | `TOL = 1e-14` |
| `charge_conservation` | 随机场过滤后全局电荷守恒 | `TOL = 1e-14` |
| `approximate_idempotence` | 轴对称场两次过滤的相对偏差 | `TOL_LOOSE = 1e-2` |
| `control_volume_sum` | `sum(Omega_p) == 物理域体积` | `TOL = 1e-14` |
| `r0_regularity` | 轴心邻域无显著伪影 | `TOL_LOOSE = 1e-2` |

---

## 10. 与 `ECDI_SUPPRESSION_PLAN.md` 的关系

- `ECDI_SUPPRESSION_PLAN.md`：物理计划、风险控制、算例级验证策略、实施顺序。
- **本文档**：离散算法精确规格、公式、代码映射、单元测试基准。

在 WarpX 中正式实现时，建议两文档对照阅读：
- 若需理解**为什么做**，参考 `ECDI_SUPPRESSION_PLAN.md`。
- 若需理解**怎么做**，参考本文档。
