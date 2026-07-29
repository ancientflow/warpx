# Hall 注入配置说明

本文档说明 `Source/Insert/Injection` 中当前 Hall 注入框架的输入文件配置方法。

## 总体行为

注入入口由 `insert` 命名空间下的两个数组控制：

```text
insert.initial_sources = source1 source2
insert.continuous_sources = source3 source4
```

- `insert.initial_sources`：初始化阶段执行一次，对应 `HallInjector::InitializePlasma()`。
- `insert.continuous_sources`：每个注入步执行，对应 `HallInjector::InjectParticles()`。
- 如果这两个数组没有显式写出，内部 source 列表为空，不会进行任何默认初始化或默认注入。
- source 名称没有固定集合。除特殊的 `average_ionization_source` 外，每个 source 都按 `<source_name>.*` 参数解析。

## Source 基本结构

每个普通 source 需要配置三部分：

```text
<source>.rate = ...
<source>.position....
<source>.species = ...
# 或
<source>.shared_species = ...
<source>.velocity....
# 或 <source>.<species>.velocity....
```

### 物种列表

单物种 source：

```text
my_source.species = electrons
```

多物种共享同一组位置样本：

```text
my_source.shared_species = electrons xe_ions
```

每个物种的宏粒子权重按如下优先级读取：

```text
my_source.electrons.weight = elec_weight
```

如果没有物种专属权重，则使用 source 级别权重：

```text
my_source.macro_weight = elec_weight
```

如果两者都没有，默认权重为 `1.0`。注意：部分 rate model 本身也需要 `macro_weight` 参数，见下一节。

## 注入数量模型

`<source>.rate` 支持以下取值。

### fixed_count

每次执行 source 时注入固定宏粒子数：

```text
my_source.rate = fixed_count
my_source.count = 1000
```

### density_volume

按密度、体积和宏粒子权重计算一次性宏粒子数：

```text
my_source.rate = density_volume
my_source.density = 1e18
my_source.volume = 1.0e-8
my_source.macro_weight = elec_weight
```

宏粒子数为：

```text
density * volume / macro_weight
```

### current

按电流和时间步计算宏电子数：

```text
my_source.rate = current
my_source.current = Ic
my_source.macro_weight = elec_weight
my_source.l_factor = l_factor
```

`l_factor` 可省略，默认 `1.0`。

### mass_flow

按质量流量计算宏粒子数：

```text
my_source.rate = mass_flow
my_source.mass_flow = m_dot
my_source.mass = 2.179e-25
my_source.macro_weight = xe_weight
my_source.l_factor = l_factor
```

`l_factor` 可省略，默认 `1.0`。

### batch_multiplier

连续注入时可以设置批量注入阈值：

```text
my_source.batch_multiplier = 100.0
```

设置后，source 会累积期望宏粒子数，只有累积超过批量阈值时才一次性注入一批。未设置时使用普通小数累积器，每步按整数部分注入。

## 位置分布

普通位置分布使用：

```text
my_source.position.coordinate_system = cylindrical
```

支持的 position 坐标系：

- `cartesian`：轴名为 `x y z`
- `cylindrical`：轴名为 `r theta z`

source 级别还可以加平移：

```text
my_source.x_offset = 0.025 / l_factor
my_source.y_offset = 0.025 / l_factor
```

### 笛卡尔位置示例

```text
my_source.position.coordinate_system = cartesian
my_source.position.x.distribution = constant
my_source.position.x.value = 0.01
my_source.position.y.distribution = constant
my_source.position.y.value = 0.02
my_source.position.z.distribution = uniform
my_source.position.z.min = 0.0
my_source.position.z.max = 0.01
```

### 柱坐标位置示例

```text
my_source.position.coordinate_system = cylindrical
my_source.position.r.distribution = area_uniform
my_source.position.r.min = 0.0105 / l_factor
my_source.position.r.max = 0.0155 / l_factor
my_source.position.theta.distribution = multi_spoke
my_source.position.theta.spoke_count = 4
my_source.position.theta.spoke_sigma = pi / 4
my_source.position.theta.spoke_phase = 0
my_source.position.z.distribution = uniform
my_source.position.z.min = 0.001 / l_factor
my_source.position.z.max = 0.004 / l_factor
```

## 一维分布类型

每个坐标轴使用 `<prefix>.distribution` 指定分布。

### constant

```text
axis.distribution = constant
axis.value = 0.0
```

### uniform

```text
axis.distribution = uniform
axis.min = 0.0
axis.max = 1.0
```

### area_uniform

按面积均匀采样半径，常用于柱坐标 `r`：

```text
axis.distribution = area_uniform
axis.min = 0.01
axis.max = 0.02
```

### gaussian

```text
axis.distribution = gaussian
axis.mean = 0.0
axis.sigma = 1.0
```

`mean` 可省略，默认 `0.0`。

### positive_gaussian

正值截断高斯，常用于入口 `vz`：

```text
axis.distribution = positive_gaussian
axis.mean = vz0
axis.sigma = sqrt(kb * Tz / mass)
```

### single_spoke

```text
axis.distribution = single_spoke
axis.spoke_center = pi
axis.spoke_sigma = pi / 4
axis.num_bins = 1024
```

也可以使用 `center` 和 `sigma` 参数名。`num_bins` 可省略，默认 `1024`。

### multi_spoke

```text
axis.distribution = multi_spoke
axis.spoke_count = 4
axis.spoke_sigma = pi / 8
axis.spoke_phase = 0
axis.num_bins = 1024
```

`spoke_phase` 可省略，默认 `0`。`num_bins` 可省略，默认 `1024`。

### neutral_spoke

中性气体初始化用的 spoke 耗尽分布。它用于构造一个稳定 spoke 扫过后的中性气体
角向分布：电离区内中性密度快速下降，电离区外用一条直线恢复到下一个电离前沿。

```text
axis.distribution = neutral_spoke
axis.ion_width = 0.30*pi
axis.min_ratio = 0.25
axis.drop_exponent = 4.0
axis.phase = 0.0
axis.reverse = 0
axis.num_bins = 1024
```

也可以使用名称 `neutral_spoke_depletion`。

令：

```text
phi = mod(theta - phase, 2*pi)              # reverse = 0
phi = mod(phase - theta, 2*pi)              # reverse = 1
n_min = min_ratio * n_max
```

电离区内：

```text
0 <= phi < ion_width
s = phi / ion_width
n(phi) = n_min + (n_max - n_min) * (1 - s)^drop_exponent
```

电离区外：

```text
ion_width <= phi < 2*pi
s = (phi - ion_width) / (2*pi - ion_width)
n(phi) = n_min + (n_max - n_min) * s
```

参数含义：

- `ion_width`：电离前沿到电离后沿的角宽，必须满足 `0 < ion_width < 2*pi`。
- `min_ratio`：电离后沿最低中性密度与最高中性密度之比，默认 `0.25`。
- `drop_exponent`：电离区内快速下降指数，默认 `4.0`。越大，下降越集中在电离前沿附近。
- `phase`：电离前沿相位，默认 `0.0`。也可写作 `phase_offset` 或 `spoke_phase`。
- `reverse`：相位取反开关，默认 `0`。也可写作 `reverse_phase` 或 `phase_reverse`。
- `num_bins`：数值 CDF 采样 bin 数，默认 `1024`。

### discrete

从离散值中采样：

```text
axis.distribution = discrete
axis.values = 0.0 1.0 2.0
axis.weights = 1.0 2.0 1.0
```

`weights` 可省略，此时所有值等权。

### tabulated

按表格 PDF 采样：

```text
axis.distribution = tabulated
axis.values = 0.0 0.5 1.0
axis.pdf = 0.0 1.0 0.0
axis.num_bins = 1024
```

### parser

使用 parser 表达式作为 PDF：

```text
axis.distribution = parser
axis.min = 0.0
axis.max = 2*pi
axis.function(theta) = 1.0 + 0.5*cos(theta)
axis.num_bins = 1024
```

变量名由轴名决定，例如 `theta` 轴使用 `function(theta)`。

## 速度分布

速度分布使用：

```text
my_source.velocity.coordinate_system = cartesian
```

支持的 velocity 坐标系：

- `cartesian`：轴名为 `vx vy vz`
- `cylindrical`：轴名为 `vr vtheta vz`
- `local_normal`：轴名为 `vnormal vt1 vt2`

如果某个物种配置了专属速度分布，则使用：

```text
my_source.electrons.velocity.coordinate_system = cartesian
```

否则使用 source 级别速度分布：

```text
my_source.velocity.coordinate_system = cartesian
```

### 各向同性热速度示例

```text
my_source.velocity.coordinate_system = cartesian
my_source.velocity.vx.distribution = gaussian
my_source.velocity.vx.sigma = 592982
my_source.velocity.vy.distribution = gaussian
my_source.velocity.vy.sigma = 592982
my_source.velocity.vz.distribution = gaussian
my_source.velocity.vz.sigma = 592982
```

### 中性气体入口速度示例

```text
xe_neutral_inlet.velocity.coordinate_system = cartesian
xe_neutral_inlet.velocity.vx.distribution = gaussian
xe_neutral_inlet.velocity.vx.sigma = sqrt(kb * Tx / mass)
xe_neutral_inlet.velocity.vy.distribution = gaussian
xe_neutral_inlet.velocity.vy.sigma = sqrt(kb * Ty / mass)
xe_neutral_inlet.velocity.vz.distribution = positive_gaussian
xe_neutral_inlet.velocity.vz.mean = vz0
xe_neutral_inlet.velocity.vz.sigma = sqrt(kb * Tz / mass)
```

## 圆锥面位置分布

`truncated_cone_surface` 将给定径向区间上的一维母线绕 `z` 轴旋转为
轴对称圆锥面。采样器只根据斜率、区间和参考点进行几何坐标运算：

```text
xe_neutral_inlet.position.coupled_distribution = truncated_cone_surface
xe_neutral_inlet.position.slope = eb_k
xe_neutral_inlet.position.r_min = inlet_r_min
xe_neutral_inlet.position.r_max = inlet_r_max
xe_neutral_inlet.position.r_reference = eb_a1
xe_neutral_inlet.position.z_reference = 0.0
xe_neutral_inlet.position.theta_min = 0.0
xe_neutral_inlet.position.theta_max = 2*pi
```

采样坐标满足
`z = z_reference + slope * (r - r_reference)`。`r_reference` 默认等于
`r_min`，`z_reference`、`theta_min` 分别默认为 `0`，`theta_max` 默认为
`2*pi`。径向采样概率正比于 `r`，已经包含圆锥面积元中的极坐标权重，
因此宏粒子不需要再乘径向权重。圆锥轴心仍由 source 的 `x_offset` 和
`y_offset` 指定。完整二维采样使用两个独立均匀随机数：第一个对 `r^2`
采样，第二个对 `theta` 采样；`z` 随后由上述母线方程确定，不能再同时
配置独立的 `position.z.distribution`。

## 孔阵列平面位置分布

孔阵列配置为：

```text
xe_neutral_inlet.position.coupled_distribution = hole_array_plane
xe_neutral_inlet.hole_count = 48
xe_neutral_inlet.hole_radius = 0.001 / l_factor
xe_neutral_inlet.hole_ring_radius = (0.021 + 0.031) / 4 / l_factor
xe_neutral_inlet.z = 0.0
```

`hole_ring_radius` 可用旧参数名 `ring_radius` 替代。使用 `hole_array_plane` 时，`position.coordinate_system` 和 `position.r/theta/z` 不参与采样。

## 特殊 source：average_ionization_source

当编译期开启 `IONIZATION_SOURCE_INJECT` 时，可以显式加入平均电离源：

```text
insert.continuous_sources = average_ionization_source
```

该 source 不按普通 `<source>.*` 解析位置和速度，而是从 `IonizationSourceSampler` 读取预计算电离源数据，并自动注入：

- `electrons`
- `xe_ions`

使用时要求：

```text
my_constants.elec_weight
```

与电离源数据中的 electron weight 一致。

## 示例：初始 spoke 等离子体

```text
insert.initial_sources = channel_plasma

channel_plasma.rate = density_volume
channel_plasma.density = 1e18
channel_plasma.volume = ((0.0155/l_factor)^2 - (0.0105/l_factor)^2) * pi * (0.003/l_factor)
channel_plasma.macro_weight = elec_weight
channel_plasma.shared_species = electrons xe_ions
channel_plasma.x_offset = 0.025 / l_factor
channel_plasma.y_offset = 0.025 / l_factor

channel_plasma.position.coordinate_system = cylindrical
channel_plasma.position.r.distribution = area_uniform
channel_plasma.position.r.min = 0.0105 / l_factor
channel_plasma.position.r.max = 0.0155 / l_factor
channel_plasma.position.theta.distribution = single_spoke
channel_plasma.position.theta.spoke_center = pi
channel_plasma.position.theta.spoke_sigma = pi / 4
channel_plasma.position.z.distribution = uniform
channel_plasma.position.z.min = 0.001 / l_factor
channel_plasma.position.z.max = 0.004 / l_factor

channel_plasma.electrons.weight = elec_weight
channel_plasma.electrons.velocity.coordinate_system = cartesian
channel_plasma.electrons.velocity.vx.distribution = gaussian
channel_plasma.electrons.velocity.vx.sigma = 592982
channel_plasma.electrons.velocity.vy.distribution = gaussian
channel_plasma.electrons.velocity.vy.sigma = 592982
channel_plasma.electrons.velocity.vz.distribution = gaussian
channel_plasma.electrons.velocity.vz.sigma = 592982

channel_plasma.xe_ions.weight = elec_weight
channel_plasma.xe_ions.velocity.coordinate_system = cartesian
channel_plasma.xe_ions.velocity.vx.distribution = gaussian
channel_plasma.xe_ions.velocity.vx.sigma = 1212.41
channel_plasma.xe_ions.velocity.vy.distribution = gaussian
channel_plasma.xe_ions.velocity.vy.sigma = 1212.41
channel_plasma.xe_ions.velocity.vz.distribution = gaussian
channel_plasma.xe_ions.velocity.vz.sigma = 1212.41
```

## 示例：阴极电子连续注入

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

## 示例：Xe 中性气体耗尽分布初始化

```text
insert.initial_sources = initial_xe_neutral

initial_xe_neutral.species = xe_netural
initial_xe_neutral.rate = density_volume
initial_xe_neutral.density = n_xe0
initial_xe_neutral.volume = volume_xe
initial_xe_neutral.macro_weight = xe_weight
initial_xe_neutral.x_offset = 0.025 / l_factor
initial_xe_neutral.y_offset = 0.025 / l_factor

initial_xe_neutral.position.coordinate_system = cylindrical
initial_xe_neutral.position.r.distribution = area_uniform
initial_xe_neutral.position.r.min = r_min
initial_xe_neutral.position.r.max = r_max
initial_xe_neutral.position.theta.distribution = neutral_spoke
initial_xe_neutral.position.theta.ion_width = 0.30*pi
initial_xe_neutral.position.theta.min_ratio = 0.25
initial_xe_neutral.position.theta.drop_exponent = 4.0
initial_xe_neutral.position.theta.phase = neutral_phase
initial_xe_neutral.position.theta.reverse = 0
initial_xe_neutral.position.z.distribution = uniform
initial_xe_neutral.position.z.min = z_min
initial_xe_neutral.position.z.max = z_max

initial_xe_neutral.velocity.coordinate_system = cartesian
initial_xe_neutral.velocity.vx.distribution = gaussian
initial_xe_neutral.velocity.vx.sigma = sqrt(kb * Tx / 2.179e-25)
initial_xe_neutral.velocity.vy.distribution = gaussian
initial_xe_neutral.velocity.vy.sigma = sqrt(kb * Ty / 2.179e-25)
initial_xe_neutral.velocity.vz.distribution = gaussian
initial_xe_neutral.velocity.vz.sigma = sqrt(kb * Tz / 2.179e-25)
```

## 当前限制

- 不再提供任何隐式默认 source。所有初始化和注入都必须通过输入文件显式列入 `insert.initial_sources` 或 `insert.continuous_sources`。
- 当前没有输入参数对应旧代码中的 `setTimeStepOverride(5.6e-10)`。如果需要复现该行为，需要新增 source 级别的 `dt_override` 参数解析。
- `hole_array_plane` 只支持位置采样；速度仍通过普通 velocity distribution 配置。
