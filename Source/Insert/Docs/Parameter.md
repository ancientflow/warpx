# Insert 诊断与注入参数

本文档整理 `Source/Insert` 当前运行时诊断和 Hall 注入输入参数。参数名按
AMReX `ParmParse` 写法列出，例如 `my_constants.foo` 和 `insert.foo`。

## 运行时诊断

### 总体节奏

边界吸收粒子相关诊断共享同一个间隔：

```text
my_constants.hall_diag_interval = 10
```

默认值为 `10`，内部会限制为至少 `1`。以下路径均使用该间隔读取
`ParticleBoundaryBuffer`：

- `AnodeCurrentCalc`
- `ZMinWallChargeDeposit` 的电荷沉积部分
- `ThrustCalc`
- `BeamDivergenceCalc`
- `IEDFCalc`
- `ClearHallBoundaryParticleCache`
- `SecondaryEmission`
- `AnodeIonNeutralization`
- `NeutralAtomEBInteraction`

这意味着边界粒子缓存的读取和清理使用同一时步判断，避免诊断间隔不一致导致
粒子在被统计前被清掉。

### 开关与输出

| 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `my_constants.particle_number_diag` | `0` | 打印每个 species 当前粒子数。需要编译宏 `NUMP`。 |
| `my_constants.collision_record_diag` | `0` | 写出碰撞产生的电子和离子宏粒子数到 `collision_record.dat`。需要编译宏 `COLLISION_RECORD`。 |
| `my_constants.anode_current_diag` | `0` | 统计 zmin 和阳极环相关电流计数。需要 `HALL3D`。 |
| `my_constants.zmin_wall_charge_diag` | `0` | 沉积 zmin 非阳极环壁面电荷。需要 `HALL3D`。 |
| `my_constants.thrust_diag` | `0` | 统计出口离子轴向动量并计算推力。需要 `HALL3D`。 |
| `my_constants.beam_divergence_diag` | `0` | 统计出口离子束流发散角。需要 `HALL3D`。 |
| `my_constants.iedf_diag` | `0` | 统计出口离子能量分布函数。需要 `HALL3D`。 |
| `my_constants.clear_hall_boundary_particle_cache_diag` | `0` | 按 `hall_diag_interval` 清理边界粒子缓存。需要 `HALL3D`。 |

`clear_hall_boundary_particle_cache_diag` 应在所有读取边界缓存的诊断之后执行。
当前 `Insert::AfterDiagnostics()` 中的调用顺序已满足这一点。

### 中性原子 EB 壁面作用

```text
insert.neutral_atom_eb.enabled = 1
insert.neutral_atom_eb.species = xe_netural
insert.neutral_atom_eb.model = diffuse
xe_netural.save_particles_at_eb = 1
```

该功能默认关闭，触发间隔复用 `my_constants.hall_diag_interval`。当前仅建立
基于 `VariableCountCopyTransformBoundaryBuffer` 的 Decision、Transform、目标粒子分配
以及漫反射/镜面反射设备算子骨架，壁面作用逻辑尚未实现。`model` 可取
`diffuse` 或 `specular`，默认值为 `diffuse`。处理后不会单独清理中性原子
buffer；如需与电子、离子统一清理，应同时设置：

```text
my_constants.clear_hall_boundary_particle_cache_diag = 1
```

### 阳极电流

```text
my_constants.anode_current_diag = 1
my_constants.anode_current_path = "anode_current.dat"
```

`anode_current_path` 默认值为 `anode_current.dat`。输出列为：

```text
time    zmin_electron    zmin_ion    anode_electron    anode_electron_cut
```

统计对象：

- `electrons` 在 zlo 边界缓存中的权重。
- `xe_ions` 在 zlo 边界缓存中可回溯命中阳极环的权重。
- 电子当前位置或回溯到 zmin 后命中阳极环的权重。

阳极环几何来自 `ReadHallAnodeRingConfig()`：

```text
my_constants.L = <physical length>
my_constants.l_factor = <scale factor>
my_constants.voltage = <stored in config>
```

其中阳极环中心为 `(ProbLo(0) + L/l_factor/2, ProbLo(1) + L/l_factor/2)`，
半径范围固定为 `0.021/2/l_factor` 到 `0.031/2/l_factor`。

### zmin 壁面电荷

```text
my_constants.zmin_wall_charge_diag = 1
my_constants.zmin_wall_charge_dir = "zmin_wall_charge"
my_constants.zmin_wall_charge_write_interval = 100
```

参数：

| 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `my_constants.zmin_wall_charge_dir` | `"zmin_wall_charge"` | 写出目录。空字符串会回退到默认目录。 |
| `my_constants.zmin_wall_charge_write_interval` | `100` | 文件写出步号间隔。 |
| `my_constants.zmin_wall_charge_interval` | `100` | 兼容旧参数名；现在仅作为写出间隔读取。 |

壁面电荷沉积不使用写出间隔，而是使用 `my_constants.hall_diag_interval`。
写出文件名格式：

```text
<zmin_wall_charge_dir>/zmin_wall_charge_00000100.dat
```

文件内容是从运行开始累计到当前写出步的 zmin 壁面电荷密度
`sigma_s(x,y)`，单位为 `C/m^2`。沉积使用粒子的物理电荷符号：电子贡献为负，
离子贡献为正。数组按 y 节点为行、x 节点为列写出，x 是最快变化索引。

当该累计量用于 zmin 非齐次 Neumann 边界时，约定 zmin 面外法向
`n = -z`，并施加 `dphi/dn = sigma_s/epsilon0`，等价于
`dphi/dz = -sigma_s/epsilon0`。因此 `InsertBoundaryPhi` 的 guard cell 更新为
`phi(k0-1) = phi(k0+1) + 2*dz*sigma_s/epsilon0`。

### 推力

```text
my_constants.thrust_diag = 1
my_constants.thrust_diag_path = "thrust.dat"
```

`thrust_diag_path` 默认值为 `thrust.dat`。统计出口边界：

```text
xlo xhi ylo yhi zhi
```

输出列为：

```text
time    dt    ion_weight    thrust_N    axial_momentum_kg_m_per_s
```

`dt` 为本次诊断与上次诊断之间的物理时间；首次统计时如果无法由时间差得到，
使用 `hall_diag_interval * WarpX::getdt(0)` 估算。

### 束流发散

```text
my_constants.beam_divergence_diag = 1
my_constants.beam_divergence_path = "beam_divergence.dat"
```

`beam_divergence_path` 默认值为 `beam_divergence.dat`。输出列为：

```text
time    dt    ion_weight    divergence_angle_rad    divergence_angle_deg
axial_momentum_kg_m_per_s    transverse_momentum_kg_m_per_s
```

发散角按 `atan2(transverse_momentum, axial_momentum)` 计算。

### IEDF

```text
my_constants.iedf_diag = 1
my_constants.iedf_path = "iedf.dat"
my_constants.iedf_bins = 200
my_constants.iedf_min_eV = 0
my_constants.iedf_max_eV = 500
```

参数：

| 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `my_constants.iedf_path` | `"iedf.dat"` | 输出文件。 |
| `my_constants.iedf_bins` | `200` | 能量 bin 数，内部限制为至少 `1`。 |
| `my_constants.iedf_min_eV` | `0` | 能量下限。 |
| `my_constants.iedf_max_eV` | `500` | 能量上限；若不大于下限，会设为 `min + 1`。 |

输出列为：

```text
time    bin_center_eV    bin_lo_eV    bin_hi_eV    ion_weight
current_A    pdf_per_eV
```

每次诊断块写完后会额外写一个空行。

## Hall 注入

### 入口参数

当前 `HALL3D` / `HALL3D_INIT` 注入入口统一走 `HallInjector`。全局入口参数为：

```text
insert.initial_sources = source1 source2
insert.continuous_sources = source3 source4
my_constants.dt = 3e-12
```

- `insert.initial_sources` 在初始化阶段执行一次。
- `insert.continuous_sources` 在每个粒子注入步执行。
- 两个数组为空时，不会创建任何默认 source。
- 连续注入使用 `my_constants.dt` 作为传入 source 的时间步。若未设置，该值为 `0`，
  依赖 `dt` 的速率模型不会产生粒子。

### Source 基本结构

普通 source 使用 source 名作为参数前缀：

```text
<source>.rate = current
<source>.species = electrons
<source>.macro_weight = elec_weight

<source>.position.coordinate_system = cylindrical
<source>.position.r.distribution = area_uniform
...

<source>.velocity.coordinate_system = cartesian
<source>.velocity.vx.distribution = gaussian
...
```

多物种 source 使用同一组位置样本：

```text
<source>.shared_species = electrons xe_ions
```

每个物种的粒子权重读取优先级为：

```text
<source>.<species>.weight
<source>.macro_weight
1.0
```

每个物种的速度分布读取优先级为：

```text
<source>.<species>.velocity.*
<source>.velocity.*
```

可选 source 级参数：

| 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `<source>.x_offset` | `0.0` | 对采样位置 x 加偏移。 |
| `<source>.y_offset` | `0.0` | 对采样位置 y 加偏移。 |
| `<source>.batch_multiplier` | 未启用 | 批量注入阈值倍数；设置时必须 `>= 1.0`。 |

### 注入数量模型

`<source>.rate` 支持以下取值：

| rate | 必填参数 | 可选参数 | 宏粒子期望数 |
| --- | --- | --- | --- |
| `fixed_count` | `count` | 无 | `count` |
| `density_volume` | `density`, `volume`, `macro_weight` | 无 | `density * volume / macro_weight` |
| `current` | `current`, `macro_weight` | `l_factor = 1.0` | `current * dt / (l_factor^2 * q_e * macro_weight)` |
| `mass_flow` | `mass_flow`, `mass`, `macro_weight` | `l_factor = 1.0` | `mass_flow * dt / (l_factor^2 * mass * macro_weight)` |

`fixed_count` 和 `density_volume` 不使用 `dt`，适合一次性初始化。
`current` 和 `mass_flow` 用于连续注入。

### 坐标系统

位置分布默认使用 `cylindrical`，速度分布默认使用 `cartesian`。

| coordinate_system | position 轴名 | velocity 轴名 |
| --- | --- | --- |
| `cartesian` | `x`, `y`, `z` | `vx`, `vy`, `vz` |
| `cylindrical` | `r`, `theta`, `z` | `vr`, `vtheta`, `vz` |
| `local_normal` | 不支持 | `vnormal`, `vt1`, `vt2` |

`local_normal` 速度坐标会根据采样位置中的法向量转换到笛卡尔速度。
普通坐标分布产生的位置法向量默认为 `(0, 0, 1)`。

### 一维分布

每个坐标轴通过 `<prefix>.distribution` 指定一维分布：

| distribution | 必填参数 | 可选参数 | 说明 |
| --- | --- | --- | --- |
| `constant` | `value` | 无 | 固定值。 |
| `uniform` | `min`, `max` | 无 | 区间均匀分布，要求 `max > min`。 |
| `area_uniform` | `min`, `max` | 无 | 面积均匀半径分布，要求 `min >= 0` 且 `max > min`。 |
| `gaussian` | `sigma` | `mean = 0` | 高斯分布，要求 `sigma >= 0`。 |
| `positive_gaussian` | `sigma` | `mean = 0` | 重采样直到值非负。 |
| `single_spoke` | `sigma` 或 `spoke_sigma` | `center` / `spoke_center = pi`, `num_bins = 1024` | 单 spoke 角向分布。 |
| `multi_spoke` | `sigma` 或 `spoke_sigma` | `spoke_count = 1`, `spoke_phase = 0`, `num_bins = 1024` | 多 spoke 角向分布。 |
| `neutral_spoke` | `ion_width` | `min_ratio = 0.25`, `drop_exponent = 4.0`, `phase` / `phase_offset` / `spoke_phase = 0`, `reverse = 0`, `num_bins = 1024` | 中性气体 spoke 耗尽分布。 |
| `neutral_spoke_depletion` | 同 `neutral_spoke` | 同 `neutral_spoke` | `neutral_spoke` 的别名。 |
| `discrete` | `values` | `weights` | 离散采样；未给 `weights` 时等权。 |
| `tabulated` | `values`, `pdf` | `num_bins = 1024` | 表格 PDF，`values` 必须严格递增。 |
| `parser` | `min`, `max`, `function(<axis>)` 或 `function` | `num_bins = 1024` | 解析表达式 PDF。 |

### hole_array_plane 位置耦合分布

位置分布可以使用特殊耦合分布：

```text
<source>.position.coupled_distribution = hole_array_plane
<source>.hole_count = 48
<source>.hole_radius = <radius>
<source>.hole_ring_radius = <ring_radius>
<source>.z = 0.0
```

参数：

| 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `<source>.hole_count` | `48` | 孔数量，必须 `> 0`。 |
| `<source>.hole_radius` | 必填 | 单个孔半径，必须 `> 0`。 |
| `<source>.hole_ring_radius` | 若缺省则读取 `<source>.ring_radius` | 孔中心所在环半径，必须 `>= 0`。 |
| `<source>.z` | `0.0` | 注入平面 z 坐标。 |

使用该模式时，`<source>.position.*` 的普通三轴分布不会用于位置采样。
速度仍按普通 velocity 分布配置。

### 平均电离源

编译宏 `IONIZATION_SOURCE_INJECT` 打开时，`insert.initial_sources` 或
`insert.continuous_sources` 可列出特殊 source：

```text
insert.continuous_sources = average_ionization_source
```

该 source 不读取普通 `<source>.*` 参数，而是读取 `ionization_source_fab` 目录：

```text
ionization_source_fab/metadata.txt
ionization_source_fab/node_source_rate
ionization_source_fab/cell_cdf
```

`metadata.txt` 需要包含：

```text
center = <x_center> <y_center>
r_min = <value>
r_max = <value>
z_min = <value>
z_max = <value>
nr = <value>
nz = <value>
A_tot = <total_rate>
elec_weight = <weight>
```

运行时会要求 `metadata.txt` 中的 `elec_weight` 与当前
`my_constants.elec_weight` 一致。注入物种固定为 `electrons` 和 `xe_ions`，
权重均为 `my_constants.elec_weight`。

### 常用配置片段

连续阴极电子注入：

```text
insert.continuous_sources = cathode_electron

cathode_electron.species = electrons
cathode_electron.rate = current
cathode_electron.current = Ic
cathode_electron.macro_weight = elec_weight
cathode_electron.l_factor = l_factor
cathode_electron.x_offset = 0.025 / l_factor
cathode_electron.y_offset = 0.025 / l_factor

cathode_electron.position.coordinate_system = cylindrical
cathode_electron.position.r.distribution = area_uniform
cathode_electron.position.r.min = 0.016 / l_factor
cathode_electron.position.r.max = 0.020 / l_factor
cathode_electron.position.theta.distribution = uniform
cathode_electron.position.theta.min = 0.0
cathode_electron.position.theta.max = 2*pi
cathode_electron.position.z.distribution = uniform
cathode_electron.position.z.min = 0.042 / l_factor
cathode_electron.position.z.max = 0.046 / l_factor

cathode_electron.velocity.coordinate_system = cartesian
cathode_electron.velocity.vx.distribution = gaussian
cathode_electron.velocity.vx.sigma = 592982
cathode_electron.velocity.vy.distribution = gaussian
cathode_electron.velocity.vy.sigma = 592982
cathode_electron.velocity.vz.distribution = gaussian
cathode_electron.velocity.vz.sigma = 592982
```

中性气体质量流量注入：

```text
insert.continuous_sources = xe_neutral_inlet

xe_neutral_inlet.species = xe_netural
xe_neutral_inlet.rate = mass_flow
xe_neutral_inlet.mass_flow = m_dot
xe_neutral_inlet.mass = m_xe
xe_neutral_inlet.macro_weight = xe_weight
xe_neutral_inlet.l_factor = l_factor
xe_neutral_inlet.x_offset = 0.025 / l_factor
xe_neutral_inlet.y_offset = 0.025 / l_factor

xe_neutral_inlet.position.coupled_distribution = hole_array_plane
xe_neutral_inlet.hole_count = 48
xe_neutral_inlet.hole_radius = 0.001 / l_factor
xe_neutral_inlet.hole_ring_radius = (0.021 + 0.031) / 4 / l_factor
xe_neutral_inlet.z = 0.0

xe_neutral_inlet.velocity.coordinate_system = cartesian
xe_neutral_inlet.velocity.vx.distribution = gaussian
xe_neutral_inlet.velocity.vx.sigma = sqrt(kb * Tx / m_xe)
xe_neutral_inlet.velocity.vy.distribution = gaussian
xe_neutral_inlet.velocity.vy.sigma = sqrt(kb * Ty / m_xe)
xe_neutral_inlet.velocity.vz.distribution = positive_gaussian
xe_neutral_inlet.velocity.vz.mean = vz0
xe_neutral_inlet.velocity.vz.sigma = sqrt(kb * Tz / m_xe)
```

## 相关源文件

- `Source/Insert/Diagnostics/InsertRuntimeDiagnostics.cpp`
- `Source/Insert/Utils/InsertUtils.cpp`
- `Source/Insert/Injection/HallInjector.cpp`
- `Source/Insert/Injection/HallRateModel.cpp`
- `Source/Insert/Injection/HallCoordinateDistribution.cpp`
- `Source/Insert/Injection/HallDistribution1D.cpp`
- `Source/Insert/Injection/HallInjectionSource.cpp`
