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

## 工况 1：壁面 600K

- **来源**：ZmaxRadialExitAnalyzer（3d_xe_dsmc 壁面 600K 工况）
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

## 工况 2：壁面 400K

- **来源**：ZmaxRadialExitAnalyzer（3d_xe_dsmc 壁面 400K 工况，
  `Script/3d_xe_dsmc`，漫反射壁面温度 400K、镜面反射比例 0.05、注入温度 400K）
- **归一化半径**：`u = r/R`，`R = 0.5 mm`
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

## 工况 3：壁面/注入 800K

- **来源**：ZmaxRadialExitAnalyzer（3d_xe_dsmc 壁面 800K 工况，
  `Script/3d_xe_dsmc`，新版 `insert.analytic_walls` 框架，
  漫反射壁面温度 800K、镜面反射比例 0.05、注入温度 800K）
- **归一化半径**：`u = r/R`，`R = 0.5 mm`
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

<!-- 新增工况请复制上方模板，注明拟合来源与使用位置 -->
