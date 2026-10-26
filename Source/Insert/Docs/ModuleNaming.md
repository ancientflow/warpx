# 模块命名与适用范围

更新日期：2026-10-26。

通用代码按数值方法、数据职责和物理过程命名，不以 CHT、PHT 或霍尔推力器作为
模块分类。具体装置通过输入中的几何、物种、材料、边界及源配置表达。此次调整
修改命名、构建引用和文档，并为诊断参数提供兼容映射，不改变求解与采样算法。

## 本次命名迁移

| 原名称 | 当前名称与职责 |
|---|---|
| `cht` | `Insert::Poisson`：节点材料 Poisson 算子和边界辅助函数 |
| `CHT_MLNodeRobinLaplacian.*` | `MLNodeRobinLaplacian.*` |
| `CHT_MonopoleBoundary.H` / `CHT_NodeRobin_K.H` | `MonopoleBoundary.H` / `NodeRobin_K.H` |
| `HallElectrostaticMaterial` | `ElectrostaticMaterial`：介电系数、固定电势 mask 和指定电势 |
| `HallWallCharge` | `WallCharge`：解析壁面持久电荷与网格源 |
| `HallElectrostaticConfig.H` | `ElectrostaticMaterialConfig.H` |
| `WARPX_USE_HALL_ELECTROSTATIC_MATERIALS` | `WARPX_USE_ELECTROSTATIC_MATERIALS` |
| `HallInjector` / `HallInjectionSource` | `ParticleInjector` / `ParticleInjectionSource` |
| `HallRateModel` | `InjectionRateModel` |
| `HallCoordinateDistribution` / `HallCoordinateTransform` | `CoordinateDistribution` / `CoordinateTransform` |
| `HallDistribution1D` / `HallPositionSampler` | `Distribution1D` / `PositionSampler` |
| `Hall*` 分布、速率和坐标辅助类型 | 删除 `Hall` 前缀，保留分布或坐标的具体含义 |
| `InitializeHallInjection` / `InjectHallParticles` | `InitializeParticleInjection` / `InjectConfiguredParticles` |
| `ClearHallBoundaryParticleCache` | `ClearBoundaryParticleCache` |
| `HALL3D` / `HALL3D_INIT` | `WARPX_INSERT_3D` / `WARPX_INSERT_INITIALIZATION_3D` |
| `HallAnodeRingConfig` | `LegacyAnnularElectrodeConfig`：保留旧 zmin 环形电极尺寸约定的兼容助手 |
| `Examples/Tests/hall_injection` | `Examples/Tests/particle_injection`：通用坐标采样测试 |

对应文件名、头文件保护宏、包含路径、CMake 源列表、调用点及文档引用同步更新。
注入设计和配置文档分别改为 `PARTICLE_INJECTOR_DESIGN.md` 和
`PARTICLE_INJECTION_CONFIG.md`，取消“注入框架只服务霍尔推力器”的设计限制。
编译宏的启用条件和默认开关保持原有行为；内部 C++ 名称不保留旧命名别名。

两个运行时参数提供中性名称：

| 兼容旧名称 | 推荐名称 |
|---|---|
| `my_constants.hall_diag_interval` | `my_constants.boundary_particle_diag_interval` |
| `my_constants.clear_hall_boundary_particle_cache_diag` | `my_constants.clear_boundary_particle_cache_diag` |

`Insert::BackwardCompatibility()` 读取旧参数并补充新名称；同时指定新旧名称时必须
取相同值，否则明确报错。既有输入文件可以继续使用，不会静默丢失旧配置。

## 分类审查与仍然具体的模型

- `Fields/NodeRobin`、`Fields/ElectrostaticMaterial` 和 `Fields/WallCharge` 按数学
  模型与源数据分类。原有三维、静态网格、边界类型等实现限制仍然有效，重命名
  不意味着增加了维度、几何或边界条件支持。
- `Injection` 提供按配置组合的粒子源、坐标分布和速率模型，物种由源参数指定。
  为兼容现有输入，位置默认柱坐标、速度默认笛卡尔坐标；可以显式选择坐标系。
- `CeramicInteraction`、`StainlessSteelInteraction`、Xe 碰撞截面工具和 ECDI 模型
  有具体材料或物理过程含义，保留这些名称。上游混合 PIC/Ohm 求解中的 Hall 项
  是物理项，不是装置分类，也保持原名。
- 真实装置输入路径（如 `Script/3d_hall_*`）、文献和历史算例记录保留其名称。
  它们不定义通用材料或注入 API 的适用范围。
- `LegacyAnnularElectrodeConfig` 仍含旧环形电极的固定半径及坐标约定，不能当作
  任意几何生成器。一般电极使用材料隐式函数或配置的解析壁面。
- `InsertInjection.cpp` 中的旧 Xe/阴极/环形初始化函数，以及可选平均电离源中的
  既有物种和温度约定，仍属于具体源模型；本次没有重写其物理参数。
- `ThrustCalc`、`BeamDivergenceCalc`、`IEDFCalc` 是显式启用的出口诊断，当前
  仍固定读取 `xe_ions`，不能据其通用诊断目录推断已支持任意离子物种。
  后续推广需单独配置物种并核对质量、电荷和出口定义。

## 验证

- CLion 构建仓库内 `build` 的 `app_3d` 通过（CUDA / MPI / DOUBLE）。
- 重命名后的 `inputs_test_3d_particle_rotating_axis` 正常运行；原有
  `analysis_rotating_axis.py` 通过，未修改统计断言或容差。
  分析使用 `build/naming-validation-env` 中的 yt/scipy，输出见
  `build/naming-injection/analysis.log`。
- 两层材料 Dirichlet 和单层 monopole 直接运行均退出 0，固定节点误差为 0；
  日志分别为 `build/naming-dirichlet/run.log` 和 `build/naming-monopole/run.log`。
- 旧参数和新参数分别运行均通过，见 `build/naming-legacy-parameters/run.log`
  和 `build/naming-neutral-parameters/run.log`；新旧间隔参数冲突时按预期中止，
  见 `build/naming-conflicting-parameters/run.log`。
- 源码扫描未发现旧通用类型、命名空间、编译宏或包含路径残留；旧运行时参数
  仅作为兼容映射保留。`git diff --check` 通过。

未新增测试算法或 CTest 项；现有注入测试仅迁移目录和名称。当前构建仍关闭
`BUILD_TESTING`，上述验证使用已有输入直接运行二进制和分析脚本。
