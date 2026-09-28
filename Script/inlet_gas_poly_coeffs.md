# 入口气体拟合多项式系数记录

记录各入口工况下中性气体 (Xe) 出口剖面的 3 阶多项式拟合系数。

## 拟合形式

以孔内局部归一化半径 `u = r / hole_radius ∈ [0,1]` 为自变量：

```
f(u) = c0 + c1*u + c2*u^2 + c3*u^3
```

拟合的物理量：

| 前缀    | 物理量                  | 单位        |
|---------|-------------------------|-------------|
| `jf_`   | 数通量 flux             | 1/(m²·s)    |
| `tx_`   | x 方向平动温度 T_x      | K           |
| `ty_`   | y 方向平动温度 T_y      | K           |
| `tz_`   | z 方向平动温度 T_z      | K           |
| `mvz_`  | 平均轴向速度 mean_vz    | m/s         |

在输入文件中的用法（以 3d_hall 为例）：

- 位置采样：`P(r) ∝ r * flux(u)`（几何权重 r 显式写在 parser 表达式中）
- 速度：`vx/vy` 为高斯，`sigma = sqrt(kb * T(u) / m_xe)`；
  `vz` 为截断 vz > 0 的漂移高斯，`mean = mean_vz(u)`，`sigma = sqrt(kb * T_z(u) / m_xe)`

## 工况 1：壁面 600K（通道直径 1mm）

- **来源**：ZmaxRadialExitAnalyzer（3d_xe_dsmc 壁面 600K 工况）
- **通道直径**：1mm（`R = 0.5 mm`）
- **使用位置**：`Script/3d_xe_exit_inlet`（600K 拟合的采样验证算例）、
  `Script/3d_hall_dsmc_init`（`xe_neutral_inlet`，48 个直径 0.5mm 周向小孔）、
  `Script/3d_hall_dielectric_hole`（`xe_neutral_inlet`，48 个直径 0.5mm 周向小孔）

```
jf_c0 =  8.2652e23   jf_c1 = -2.2830e23   jf_c2 =  2.9381e23   jf_c3 = -4.5040e23
tx_c0 =  4.6232e2    tx_c1 =  8.8315e0    tx_c2 = -2.1431e0    tx_c3 =  6.5293e1
ty_c0 =  4.6121e2    ty_c1 =  1.6222e1    ty_c2 = -1.7777e1    ty_c3 =  7.4966e1
tz_c0 =  2.9017e2    tz_c1 = -2.6955e1    tz_c2 =  6.2493e1    tz_c3 = -5.0484e1
mvz_c0 = 3.1477e2    mvz_c1 = 2.5567e0    mvz_c2 = -2.3046e1   mvz_c3 = -4.8893e0
```

## 工况 2：壁面 400K（通道直径 1mm）

- **来源**：ZmaxRadialExitAnalyzer（3d_xe_dsmc 壁面 400K 工况，
  `Script/3d_xe_dsmc`，漫反射壁面温度 400K、镜面反射比例 0.05、注入温度 400K）
- **通道直径**：1mm（归一化半径 `u = r/R`，`R = 0.5 mm`）
- **拟合精度**：加权 RMS 相对误差 — 通量 1.01%，T_x/T_y/T_z ≤ 0.23%，mean_vz 0.08%；
  mean_vx/mean_vy 按对称性恒为 0（加权 RMS ≤ 0.06 m/s）
- **使用位置**：`Script/3d_hall`（`xe_neutral_inlet`）

```
jf_c0 =  1.1065e24   jf_c1 = -2.5115e23   jf_c2 =  1.4286e23   jf_c3 = -4.7771e23
tx_c0 =  3.0396e2    tx_c1 =  2.3711e1    tx_c2 = -4.6530e1    tx_c3 =  8.1279e1
ty_c0 =  3.0243e2    ty_c1 =  3.4419e1    ty_c2 = -6.5369e1    ty_c3 =  9.0793e1
tz_c0 =  1.9873e2    tz_c1 = -1.5765e1    tz_c2 =  3.2015e1    tz_c3 = -2.8509e1
mvz_c0 = 2.6902e2    mvz_c1 = -3.6125e0   mvz_c2 = -1.0281e1   mvz_c3 = -1.6341e1
```

## 工况 3：壁面/注入 800K（通道直径 1mm）

- **来源**：ZmaxRadialExitAnalyzer（3d_xe_dsmc 壁面 800K 工况，
  `Script/3d_xe_dsmc`，新版 `insert.analytic_walls` 框架，
  漫反射壁面温度 800K、镜面反射比例 0.05、注入温度 800K）
- **通道直径**：1mm（归一化半径 `u = r/R`，`R = 0.5 mm`）
- **拟合精度**：加权 RMS 相对误差 — 通量 1.14%，T_x/T_y/T_z ≤ 0.22%，mean_vz 0.13%；
  mean_vx/mean_vy 按对称性恒为 0（加权 RMS ≤ 0.08 m/s）
- **使用位置**：待接入（参考 `Script/3d_hall` 的 `xe_neutral_inlet` 用法）

```
jf_c0 =  8.9952e23   jf_c1 = -2.6127e23   jf_c2 =  3.7691e23   jf_c3 = -5.2273e23
tx_c0 =  6.1551e2    tx_c1 =  1.0779e1    tx_c2 =  7.3039e0    tx_c3 =  7.3112e1
ty_c0 =  6.1659e2    ty_c1 =  4.2385e0    ty_c2 =  2.2056e1    ty_c3 =  6.3314e1
tz_c0 =  3.8264e2    tz_c1 = -3.3719e1    tz_c2 =  8.0867e1    tz_c3 = -6.5401e1
mvz_c0 = 3.5923e2    mvz_c1 = 6.7369e0    mvz_c2 = -3.4595e1   mvz_c3 = 2.1311e0
```

## 工况 4：注入 800K / 壁面 600K（通道直径 0.5mm）

- **来源**：ZmaxRadialExitAnalyzer（`Script/3d_xe_dsmc` 缩小通道工况，
  新版 `insert.analytic_walls` 框架，注入温度 800K、漫反射壁面温度 600K、
  镜面反射比例 0.05、质量流率 0.6 mg/s 不变）
- **通道直径**：0.5mm（归一化半径 `u = r/R`，`R = 0.25 mm`）
- **拟合精度**：加权 RMS 相对误差 — 通量 1.44%，T_x/T_y/T_z ≤ 0.20%，mean_vz 0.10%；
  mean_vx/mean_vy 按对称性恒为 0（加权 RMS ≤ 0.08 m/s）
- **使用位置**：待接入（参考 `Script/3d_hall` 的 `xe_neutral_inlet` 用法）

```
jf_c0 =  3.0072e24   jf_c1 = -9.4486e23   jf_c2 =  1.1450e24   jf_c3 = -1.8228e24
tx_c0 =  4.5632e2    tx_c1 =  1.1027e1    tx_c2 =  3.0145e-1   tx_c3 =  6.2892e1
ty_c0 =  4.5546e2    ty_c1 =  1.8788e1    ty_c2 = -1.6058e1    ty_c3 =  7.2252e1
tz_c0 =  2.9410e2    tz_c1 = -2.2501e1    tz_c2 =  5.0391e1    tz_c3 = -4.6670e1
mvz_c0 = 3.2603e2    mvz_c1 = 5.5370e0    mvz_c2 = -3.6283e1   mvz_c3 = 2.4365e-1
```

## 工况 5：注入 800K / 壁面 600K，1/48 流量（通道直径 0.5mm）

- **来源**：ZmaxRadialExitAnalyzer（`Script/3d_xe_dsmc` 缩小通道工况，
  质量流率 0.6/48 = 0.0125 mg/s（对应 3d_hall 整机 0.6 mg/s 的 48 孔分流），
  `xe_weight = 1e5` 保持宏粒子数不变；注入温度 800K、漫反射壁面温度 600K、
  镜面反射比例 0.05）
- **通道直径**：0.5mm（归一化半径 `u = r/R`，`R = 0.25 mm`）
- **流态**：近自由分子流（Kn ~ 10–40），用于与工况 4（同几何同温度、
  48 倍流量、Kn ~ 0.1–0.5）对比气-气碰撞对出口剖面的影响
- **拟合精度**：加权 RMS 相对误差 — 通量 1.32%（cos 形更优，0.62%），
  T_x/T_y ≤ 0.61%，T_z 0.12%，mean_vz 0.18%；
  mean_vx/mean_vy 按对称性恒为 0（加权 RMS ≤ 0.14 m/s）
- **使用位置**：真实流量下的推荐系数（参考 `Script/3d_hall` 的
  `xe_neutral_inlet` 用法）

```
jf_c0 =  3.1093e22   jf_c1 = -1.2104e22   jf_c2 =  2.5988e22   jf_c3 = -2.3639e22
tx_c0 =  4.5345e2    tx_c1 = -8.9290e1    tx_c2 =  2.3191e2    tx_c3 = -1.3822e2
ty_c0 =  4.5155e2    ty_c1 = -7.6317e1    ty_c2 =  2.0671e2    ty_c3 = -1.2413e2
tz_c0 =  2.7916e2    tz_c1 = -1.3112e1    tz_c2 =  3.6252e1    tz_c3 = -3.0710e1
mvz_c0 = 2.8110e2    mvz_c1 = 1.8574e1    mvz_c2 = -4.9945e1   mvz_c3 = 3.0711e1
```

## 工况 6：注入 800K / 壁面 600K，1/48 流量，全漫反射（通道直径 0.5mm）

- **来源**：ZmaxRadialExitAnalyzer（`Script/3d_xe_dsmc`，同工况 5 但
  `p_specular = 0`，即壁面 100% 漫反射，用于评估镜面份额的敏感性上限）
- **通道直径**：0.5mm（归一化半径 `u = r/R`，`R = 0.25 mm`）
- **拟合精度**：加权 RMS 相对误差 — 通量 1.19%（cos 形更优，0.71%），
  T_x/T_y ≤ 0.57%，T_z 0.13%，mean_vz 0.21%；
  mean_vx/mean_vy 按对称性恒为 0（加权 RMS ≤ 0.2 m/s）
- **结论**：与工况 5（5% 镜面）几乎一致 — 温度/速度差 < 1%，
  仅出口通量低约 7%（透射率下降，反流增多）

```
jf_c0 =  2.8906e22   jf_c1 = -1.0278e22   jf_c2 =  2.1909e22   jf_c3 = -2.0619e22
tx_c0 =  4.4883e2    tx_c1 = -8.1788e1    tx_c2 =  2.1619e2    tx_c3 = -1.2381e2
ty_c0 =  4.4576e2    ty_c1 = -6.4611e1    ty_c2 =  1.9233e2    ty_c3 = -1.1460e2
tz_c0 =  2.7917e2    tz_c1 = -1.1341e1    tz_c2 =  3.6000e1    tz_c3 = -3.2665e1
mvz_c0 = 2.8057e2    mvz_c1 = 2.2278e1    mvz_c2 = -5.7676e1   mvz_c3 = 3.5114e1
```

<!-- 新增工况请复制上方模板，注明拟合来源与使用位置 -->

## 工况 7：注入 800K / 壁面 600K，1/48 流量，全漫反射，5 倍几何缩比（通道直径 0.5mm）

- **来源**：ZmaxRadialExitAnalyzer（`Script/3d_xe_dsmc` 缩比工况 `l_factor = 5`：
  几何/半径缩小 5 倍（R = 0.25mm → 50 um，20 径向档），`dt` 除以 5，
  质量流率 0.0125 mg/s 由注入源内部除以 `l_factor^2`，`xe_weight = 1e5/l_factor^3 = 800`，
  Xe-Xe 碰撞截面使用 `Xe_Xe_VHS_elastic_x5.dat`（截面 ×5）保持 Knudsen 数不变；
  注入温度 800K、壁面 600K、全漫反射（`p_specular` 未启用）；
  稳态累积 step 10000–19999，约 2.4e7 宏粒子）
- **通道直径**：0.5mm 物理尺寸（归一化半径 `u = r/R`，`R = 0.25 mm`；
  缩比仿真内 `R = 50 um`，归一化坐标相同）
- **拟合精度**：加权 RMS 相对误差 — 通量 1.21%（cos 形更优，0.70%），
  T_x 0.54%，T_y 0.50%，T_z 0.15%，mean_vz 0.20%；
  mean_vx/mean_vy 按对称性恒为 0（加权 RMS ≤ 0.17 m/s）
- **结论**：与工况 6（未缩比、同温度、同流量、全漫反射）几乎一致 —
  各系数差异 < 1%（如 jf_c0 2.8933e22 vs 2.8906e22），
  验证了几何缩比 + 截面 ×l_factor 的 Kn 相似性成立，缩比仿真可直接替代
  原几何算例计算进气分布
- **使用位置**：缩比正式仿真的推荐进气系数（全漫反射版；若需 5% 镜面
  反射版本，须在 `Script/3d_xe_dsmc` 启用 `p_specular = 0.05` 后重跑，
  对应未缩比结果为工况 5）

```
jf_c0 =  2.8933e22   jf_c1 = -1.0528e22   jf_c2 =  2.2350e22   jf_c3 = -2.0845e22
tx_c0 =  4.4684e2    tx_c1 = -6.8512e1    tx_c2 =  1.9290e2    tx_c3 = -1.1217e2
ty_c0 =  4.4720e2    ty_c1 = -7.0169e1    ty_c2 =  1.9720e2    ty_c3 = -1.1474e2
tz_c0 =  2.7951e2    tz_c1 = -1.2349e1    tz_c2 =  3.5171e1    tz_c3 = -3.0965e1
mvz_c0 = 2.8075e2    mvz_c1 = 2.2326e1    mvz_c2 = -5.9251e1   mvz_c3 = 3.6617e1
```
