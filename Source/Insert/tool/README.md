# Cross-Section Tools

独立编译的碰撞截面数据处理工具，用于为 WarpX `ScatteringProcess` 准备输入文件。

---

## 编译

```bash
cd Source/Insert/tool
make
```

产物输出到 `bin/` 目录：
- `bin/CrossSectionInterpolator`
- `bin/CrossSectionScaler`
- `bin/SpokeAzimuthalDensity`
- `bin/XeXeVHSGenerator`
- `bin/ZmaxRadialExitAnalyzer`

---

## 路径配置 `paths.conf`

工具运行时读取当前目录下的 `paths.conf`，将其内容作为**基础路径**。所有相对路径的参数都会自动拼接该基础路径。

默认内容：
```
../../../build/bin
```

**解析规则**：
`paths.conf` 必须与可执行文件放在**同一目录**下（即 `bin/paths.conf`）。`Makefile` 编译时会自动从源码目录同步复制到 `bin/`。

工具运行时从**可执行文件所在目录**读取 `paths.conf`，因此 `paths.conf` 中的相对路径始终相对于 `bin/` 目录解析。无论从哪个工作目录运行工具，数据路径都指向同一个位置。

**行为规则**：
- 存在 `paths.conf` 时，相对路径参数定向到配置指定的目录。
- 传入**绝对路径**时不受 `paths.conf` 影响。
- 没有 `paths.conf` 时，相对路径相对于**当前工作目录**。

> 建议在 `Source/Insert/tool/` 目录下运行工具，这样相对路径会自动映射到 `build/bin/`。

---

## 工具一：CrossSectionInterpolator

将能量非均匀分布的碰撞截面数据插值为等间距能量网格。

### 用法

```bash
CrossSectionInterpolator <input_file> <energy_step> [output_file]
```

| 参数 | 说明 |
|------|------|
| `input_file` | 输入截面文件（两列：`能量[eV]` `截面[m²]`） |
| `energy_step` | 目标能量步长（单位 eV） |
| `output_file` | **可选**。省略时默认输出为 `<input_file>_uniform.dat` |

### 示例

```bash
cd Source/Insert/tool
./bin/CrossSectionInterpolator my_xs.dat 0.5
```

- 读取 `build/bin/my_xs.dat`（通过 `paths.conf` 重定向）
- 输出 `build/bin/my_xs_uniform.dat`（默认输出文件名基于输入路径）

---

## 工具二：CrossSectionScaler

将碰撞截面的所有 `sigma` 值乘以一个常数因子。

### 用法

```bash
CrossSectionScaler <input_file> <output_file> <scale_factor>
```

| 参数 | 说明 |
|------|------|
| `input_file` | 输入截面文件 |
| `output_file` | 输出文件路径 |
| `scale_factor` | 缩放倍数（如 `2.0`、`0.5`） |

### 示例

```bash
cd Source/Insert/tool
./bin/CrossSectionScaler my_xs_uniform.dat my_xs_scaled.dat 2.0
```

- 读取 `build/bin/my_xs_uniform.dat`（通过 `paths.conf` 重定向）
- 所有截面值乘以 `2.0`
- 输出 `build/bin/my_xs_scaled.dat`

---

## 输入/输出文件格式

纯 ASCII 两列数据，与 WarpX `ScatteringProcess::readCrossSectionFile` 兼容：

```text
0.0    1.0e-20
1.0    2.5e-20
5.0    4.0e-20
```

- 每行两个值：能量（eV）、截面（m²）
- 空行自动忽略
- **不要添加 `#` 注释头**——WarpX 的读取逻辑不支持注释行

---

## 典型工作流

```bash
cd Source/Insert/tool

# Step 1: 将原始数据插值到均匀能量网格
./bin/CrossSectionInterpolator raw_xs.dat 1.0

# Step 2: 对结果进行缩放（例如调整反应率）
./bin/CrossSectionScaler raw_xs_uniform.dat final_xs.dat 0.8
```

最终 `build/bin/final_xs.dat` 可直接作为 WarpX 碰撞截面输入。

---

## 工具三：SpokeAzimuthalDensity

输出 spoke 工况下中性气体（`neutral_spoke`）和等离子体（`multi_spoke`）的**归一化周向密度**曲线（spoke 数 1–4），用于 Origin 绘图。

### 用法

```bash
SpokeAzimuthalDensity [output_file] [num_points]
```

| 参数 | 说明 |
|------|------|
| `output_file` | **可选**。输出文件路径，默认 `spoke_azimuthal_density.dat` |
| `num_points` | **可选**。周向采样点数，默认 3601（0.1° 步长） |

### 示例

```bash
cd Source/Insert/tool
./bin/SpokeAzimuthalDensity spoke_density.dat
```

### 输出格式

纯 ASCII 表格，首行为 `#` 注释的列名（Origin 导入时作为 Long Name）：

```text
# theta[deg]    neutral_n1    plasma_n1    neutral_n2    plasma_n2    ...
0.00000000e+00  2.97632914e-01  7.00092220e-01  ...
```

- 共 9 列：角度（度）+ 1–4 个 spoke 的中性/等离子体密度各 4 列
- 每条曲线按各自最大值归一化（峰值 = 1）
- 分布公式与 `Source/Insert/Injection/HallDistribution1D.cpp` 完全一致
- 形状参数取自当前 spoke 工况（`Script/3d_hall_spoke`）：中性耗尽宽度 0.30π、最小比 1/6、下降指数 4、reverse=1；等离子体峰 σ=π/8，相位相对中性最小值向峰方向偏移 15°
- **相位经过平移**：整体平移 180°，中性峰在 200°、等离子体峰在 161°，坐标取整便于绘图标注，仅为绘图方便，不改变形状
- 修改默认参数需编辑源文件顶部的常量并重新编译

---

## 工具四：XeXeVHSGenerator

用 VHS（Variable Hard Sphere）模型生成 Xe–Xe 中性原子弹性碰撞截面文件，供 WarpX DSMC 模块使用。物理规范见 `xe_coll.md`。

### 用法

```bash
XeXeVHSGenerator [output_file] [num_points] [E_min_eV] [E_max_eV]
```

| 参数 | 说明 |
|------|------|
| `output_file` | **可选**。输出文件路径，默认 `Xe_Xe_VHS_elastic.dat` |
| `num_points` | **可选**。对数间隔能量点数，默认 1000 |
| `E_min_eV` | **可选**。最小质心系能量（eV），默认 1e-6 |
| `E_max_eV` | **可选**。最大质心系能量（eV），默认 10 |

### 示例

```bash
cd Source/Insert/tool
./bin/XeXeVHSGenerator
```

- 输出 `build/bin/Xe_Xe_VHS_elastic.dat`（通过 `paths.conf` 重定向）

### 物理模型

- VHS 参数：T_ref = 273 K，d_ref = 5.74e-10 m，omega = 0.85，m_Xe = 131.293 u
- 截面公式：`sigma(E) = pi*d_ref^2 / Gamma(2.5-omega) * (E_ref/E)^(omega-1/2)`，其中 `E_ref = k_B*T_ref/e`
- 即 `sigma ∝ E^(-0.35)`
- 第一列为质心系相对能量（eV），对数间隔、严格递增，不含 E = 0

### 输出与自检

输出为纯两列 ASCII（`能量[eV]` `截面[m²]`），无注释行，与 WarpX `ScatteringProcess::readCrossSectionFile` 兼容。程序运行后打印若干测试能量点的截面值及 log-log 斜率（应约为 -0.35），并在截面非正或网格非严格递增/单调下降时报错退出。

WarpX 输入中使用方式：

```text
XeXe.type = dsmc
XeXe.species = neutral_Xe
XeXe.scattering_processes = elastic
XeXe.elastic_cross_section = Xe_Xe_VHS_elastic.dat
XeXe.elastic_scattering_angle_model = isotropic
```

---

## 工具五：ZmaxRadialExitAnalyzer

分析 WarpX `zmax_radial_exit_diag` 诊断的输出（`zmax_radial/density.dat` 与 `vdist.dat`），给出各半径档的平均速度、方向温度以及各速度分量与麦克斯韦分布的偏离程度。

### 用法

```bash
ZmaxRadialExitAnalyzer <zmax_radial_dir> [mass_kg] [output_file]
```

| 参数 | 说明 |
|------|------|
| `zmax_radial_dir` | 包含 `density.dat` 和 `vdist.dat` 的目录 |
| `mass_kg` | **可选**。粒子质量，默认 2.18017e-25（Xe 原子） |
| `output_file` | **可选**。默认 `<zmax_radial_dir>/radial_analysis.dat` |

### 示例

```bash
cd Source/Insert/tool
./bin/ZmaxRadialExitAnalyzer zmax_radial
```

### 输出列

`r_center_m  weight  num_macro  mean_vx  mean_vy  mean_vz  drift_speed  T_x  T_y  T_z  T_avg  TV_x  TV_y  TV_z  KL_x  KL_y  KL_z`

- 方向温度 `T_i = m·sigma_i^2 / k_B`，`T_avg` 为三分量平均
- `TV`（全变差距离）与 `KL`（KL 散度，nats）衡量该分量直方图与同均值、同方差麦克斯韦分布的差异，0 表示完全一致
- 直方图按 bin 中心取值；权重为累积宏粒子权重，归一化在内部完成

### 径向剖面拟合

工具还会对 `flux(r)`、`T_x/T_y/T_z(r)`、`mean_vz(r)` 做以 `num_macro` 为权重的最小二乘拟合（`u = r/R ∈ [0,1]`，`R` 取 `density.dat` 最后一档的 `r_hi`）：

- 多项式 `polyN`：`y = c0 + c1 u + … + cN u^N`（N = 0–3，正规方程 + Gauss-Jordan 部分主元消元）
- 形状族（仅通量）：`pow` = `A(1-u²)^p`，`cos` = `A·cos(πu/2)^p`；`p` 网格搜索（0.1–4，步长 0.02），`A` 对固定 `p` 解析求解

stdout 打印每个量全部模型的系数、加权 RMS 相对误差和最大相对误差，并标出最佳模型；`mean_vx`/`mean_vy` 按对称性应恒为 0，只报告加权绝对 RMS（m/s）。各量的测量值与最佳拟合值写入 `<zmax_radial_dir>/radial_fit.dat`（每量两列，便于绘图）。
