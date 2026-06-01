# Insert Reorganization Plan

## 背景和目标

`Source/Insert` 当前承载的是长期私有功能。由于这些功能预计不会向 WarpX
官方提交 PR，整理目标不是满足上游合并标准，而是降低长期跟踪
`development` 分支时的冲突成本，并让私有物理模型更容易定位、验证和继续扩展。

核心原则：

- 官方文件只保留最薄的 hook 和分发逻辑。
- 私有物理模型、诊断、注入、局部背景密度和辅助模板尽量都放在 `Source/Insert`。
- 避免在 WarpX 通用头文件中加入只服务私有模型的 overload、宏和全局状态。
- 保留当前独特功能，但把功能边界整理清楚，便于长期维护差异。

## 当前主要问题

1. `WarpXSimulationFunction.h` 和 `WarpXInsertFunction.h` 过大，混合了注入、
   诊断、边界电势、背景密度、二次电子发射和测试代码。
2. `BackgroundMCCCollision.cpp` 中包含关键私有碰撞模型修改，导致私有模型和
   WarpX 官方 MCC 实现深度混合。
3. `FilterCopyTransform.H` 增加了私有耦合电离使用的 overload，扩大了官方
   粒子创建工具的维护冲突面。
4. `global_background_density` 当前定义在头文件中，长期存在重复定义和 include
   顺序风险。
5. 私有配置使用手写宏分散控制，部分宏存在拼写或作用域问题，例如
   `BENCHMAKR2D`、无条件 `WAVE1D`。
6. 若干官方核心文件有非功能性格式化 diff，会增加未来 rebase/merge 冲突。

## 目标目录结构

建议逐步整理为：

```text
Source/Insert/
  BackgroundCoupledDensity.cpp
  BackgroundCoupledDensity.h
  BackgroundMCCCoupled.H
  BackgroundMCCCoupled.cpp
  BackgroundMCCCoupledI.H
  FilterCopyTransformCoupled.H
  InsertBoundaryParticles.cpp
  InsertBoundaryParticles.h
  InsertConfig.H
  InsertInjection.cpp
  InsertInjection.h
  InsertRuntimeDiagnostics.cpp
  InsertRuntimeDiagnostics.h
  InsertBackgroundDensity.cpp
  InsertBackgroundDensity.h
  InsertHooks.cpp
  InsertHooks.h
  InsertState.cpp
  InsertState.h
  InsertBoundaryPhi.cpp
  InsertBoundaryPhi.h
```

说明：

- `InsertHooks.*`：对外暴露给 WarpX 官方文件调用的薄入口。
- `InsertState.*`：集中保存私有全局状态，例如背景密度对象。
- `BackgroundMCCCoupled.*`：承载从 `BackgroundMCCCollision.cpp` 移出的
  耦合 MCC 模型。
- `FilterCopyTransformCoupled.H`：承载当前加到官方 `FilterCopyTransform.H`
  中的私有 overload。
- `InsertConfig.H`：集中管理私有宏和功能开关。

## 阶段 1：低风险清理

目标：先修正明显问题，不改变大结构。

状态：已完成。

1. 已确认 `WarpXSimulationFunction.h` 中没有未使用的 `amrex::MLEBABecLap ml_eb`
   和对应新增 include。
2. 已修正 `BENCHMAKR2D` 为 `BENCHMARK_2D`。
3. 已删除 `WarpXInsertFunction.h` 中无条件 `#define WAVE1D`，只保留
   `WarpXSimulationConfig.h` 中按维度启用。
4. 已删除 `BeforeCollision/AfterCollision` 的 `if_split` 参数，背景密度更新和
   清理只保留原来的非 split 分支时序。
5. 已确认 `git diff --check` 无错误。
6. 已将 `global_background_density` 的定义从头文件移动到 `.cpp`，头文件只保留
   `extern` 或通过 `InsertState` 访问。

验收标准：

- 当前功能入口行为不变。
- 编译通过。
- `git diff --check` 无错误。

## 阶段 2：建立稳定 hook 层

目标：让官方核心文件只依赖一个轻量入口。

状态：已完成；由于不创建新文件，稳定 hook 层先复用现有
`WarpXInsert.h/cpp`。

1. 已将当前 `WarpXInsert.h/cpp` 中的入口函数迁移或重命名为：
   - `Insert::Initialize()`
   - `Insert::BeforeStep()`
   - `Insert::ParticleInjection()`
   - `Insert::SetBoundaryPhi()`
   - `Insert::SetPhiGuards()`
   - `Insert::BeforeCollision(int step)`
   - `Insert::AfterCollision(int step)`
   - `Insert::AfterDiagnostics()`
2. 已将 `WarpXEvolve.cpp`、`LabFrameExplicitES.cpp` 等官方文件改为只 include
   `Insert/WarpXInsert.h`。
3. 已确认阶段 2 中官方演化和场求解 hook 调用点不再直接依赖私有实现头。
   `MCC_EXCITATION`、`COLLISION_RECORD` 等编译期功能开关必须继续生效，因此
   相关实现文件和头文件中的 `WarpXFunctionConfig.h` 依赖暂时保留到阶段 4
   统一迁移。`MCC_DENSITY_MID` 已在阶段 4.1 删除，不再作为有效功能开关。

验收标准：

- 官方核心文件中的私有 diff 只剩 hook 调用。
- hook 函数内部可通过宏决定是否执行具体私有功能。

## 阶段 3：拆分大头文件

状态：已完成。

目标：将巨型 header 拆成按功能组织的 `.cpp/.h`。阶段 3 新文件统一使用
`Insert` 前缀，避免和 WarpX 官方已有或未来新增的通用模块名冲突。所有仅用于
人工检查、临时输出或一次性验证的检查函数全部抛弃，不再迁移到新文件。

1. 注入：
   - 已新建 `InsertInjection.*`。
   - 已迁移 `CathodeInjection3D`、`XeInjection`、`XeFastInjection`、`PlasmaInit`
     等粒子注入和初始化注入函数。
2. 诊断：
   - 已新建 `InsertRuntimeDiagnostics.*`。
   - 已迁移保留运行时仍需要的诊断，例如 `ParticleNumber`、`AnodeCurrentCalc`、
     `ShowAndWriteIonzationNum`。
   - 已抛弃检查类函数，例如 `PhiExamine`、`XeRhoExamine`、`PhiBCExamine`、
     `DataExamine` 及其入口，不再保留相关宏入口。
3. 边界粒子处理：
   - 已新建 `InsertBoundaryParticles.*`。
   - 已迁移二次电子发射相关逻辑，包括 `SecondaryEmissionFilter`、
     `SecondaryEmissionTransform`、`SecondaryEmission`。
   - 后续如有边界粒子反射、吸收缓存处理等私有逻辑，也归入该分类。
4. 边界电势处理：
   - 已新建 `InsertBoundaryPhi.*`。
   - 已迁移 `AnodeVoltage`、`DirichletPhiGuardSet`、`VoltageAdjustment`。
5. 背景密度 hook：
   - 已新建 `InsertBackgroundDensity.*`。
   - 已迁移 `global_background_density`、`GlobalBackgroundDensityInit`、
     `GlobalBackgroundDensityUpdate`、`GlobalBackgroundDensityClean`。
6. 碰撞：
   - 阶段 3 未迁移 `BackgroundMCCCollision.cpp` 中的耦合碰撞主体。
   - 已只做配合当前拆分所需的 include、声明和调用调整。
   - 实际耦合 MCC 模型迁移统一放到阶段 4。
7. 已只在需要模板或 GPU lambda 的位置保留局部实现。

验收标准：

- `WarpXSimulationFunction.h` 和 `WarpXInsertFunction.h` 可以删除，或缩减为
  临时兼容 include。
- 每个功能文件的 include 更小，编译依赖更清楚。
- 检查类函数和对应入口已删除，不再参与编译。
- 碰撞行为保持不变，耦合 MCC 主体仍留在原位置等待阶段 4。

## 阶段 4：新增 Insert 私有耦合 MCC 碰撞类

目标：不再把私有耦合 MCC 算法嵌入官方 `BackgroundMCCCollision`。允许在
`Source/Insert` 中复制官方 `BackgroundMCCCollision` 的必要代码，形成独立的
私有碰撞类，并通过新的碰撞类型接入 `CollisionHandler`。官方目录中的长期
diff 应缩小到类型分发和删除历史私有改动。

阶段 4 继续保留阶段 4.1 的结论：`MCC_DENSITY_MID` 中点密度路径已删除，私有
耦合 MCC 只保留形函数插值读取局部背景密度、按沉积阶数更新 `ground_rho` 的
路径。

### 4.1 删除中点密度路径

状态：已完成。

已移除 `MCC_DENSITY_MID` 相关实现和配置入口：

- 已删除 `WarpXFunctionConfig.h` 中的 `MCC_DENSITY_MID` 开关。
- 已删除 `BackgroundCoupledDensity.h/cpp` 中按 `MCC_DENSITY_MID` 分支维护
  `std::unique_ptr<amrex::MultiFab>` 的路径，只保留
  `amrex::Vector<amrex::MultiFab>`。
- 已删除 `BackgroundMCCCollision.cpp` 和后续迁移目标中的
  `#ifdef/#ifndef MCC_DENSITY_MID` 分支，只保留形函数插值读取局部背景密度、
  并按沉积阶数更新密度的路径。
- 已删除 `FilterCopyTransform.H` 或私有迁移头中由 `MCC_DENSITY_MID` 控制的
  中点读取逻辑。

删除后已执行构建验证，确认当前耦合碰撞只剩单一路径，再开始移动算法。

### 4.2 新建 Insert 私有碰撞类

状态：已完成。

在 `Source/Insert` 中新建独立碰撞类，例如：

```text
Source/Insert/CoupledBackgroundMCCCollision.H
Source/Insert/CoupledBackgroundMCCCollision.cpp
Source/Insert/CoupledBackgroundMCCCollisionI.H
```

类定义为：

```cpp
class CoupledBackgroundMCCCollision final : public CollisionBase
{
public:
    explicit CoupledBackgroundMCCCollision(std::string const& collision_name);

    void doCollisions(
        amrex::Real cur_time,
        amrex::Real dt,
        MultiParticleContainer* mypc) override;
};
```

注意：虽然文件位于 `Source/Insert`，碰撞类本身不放入 `namespace Insert`，而是
保持和官方 `CollisionBase`、`BackgroundMCCCollision` 一致的全局命名空间。
这样 `CollisionHandler` 可以用和官方碰撞类相同的方式统一分发。

该类可以复制官方 `BackgroundMCCCollision` 中必要的构造、过程读取、概率初始化、
普通散射和电离逻辑。复制是有意为之：用 Insert 内部维护一份私有碰撞实现，换取
官方 MCC 文件长期低冲突。

私有类内部直接保存和处理：

- `ground_rho_index`
- `excitation_product`
- `have_excitation`
- `excitation_rho_index`
- `sigma_max`
- `ioni_sigma_max`
- 背景粒子容器和 `BackgroundCoupledDensity` 访问
- `elec_weight`、`ncell`、`sim_L`、`inv_gap`
- 背景粒子扣减、删除、激发态替换、局部密度更新

验收标准：

- `CoupledBackgroundMCCCollision` 已在 `Source/Insert` 中定义。
- `CoupledBackgroundMCCCollision` 使用全局命名空间，和官方碰撞类保持一致。
- 官方 `BackgroundMCCCollision.H/.cpp` 暂时可以不改。
- 复制代码允许暂时包含较大实现，后续再按功能拆分。

### 4.3 迁移耦合散射和耦合电离主体

把当前官方 `BackgroundMCCCollision.cpp` 中的私有耦合主体复制/迁移到
`CoupledBackgroundMCCCollision`：

- `doBackgroundCollisionsWithinTileCouple<depos_order>`
- `doBackgroundIonizationCouple<depos_order>`
- `MCC_EXCITATION` 下的背景粒子替换和激发态粒子生成
- `MCC_DELETE` 下的背景粒子权重扣减和删除
- 局部背景密度插值和 `ground_rho` 扣减更新

迁移后 Insert 私有类可以直接调用自身成员，不需要通过官方
`BackgroundMCCCollision` 传入 `ctx` 或官方参数。

验收标准：

- 私有耦合散射和耦合电离代码已位于 Insert 私有类。
- `MCC_EXCITATION` 仍能完整关闭激发态相关逻辑。
- `MCC_DENSITY`、`MCC_DELETE`、`COLLISION_RECORD` 等功能开关仍保持关闭/启用作用。

### 4.4 移动私有 FilterCopyTransform overload

将当前添加到官方 `FilterCopyTransform.H` 的私有 overload 移动到：

```text
Source/Insert/FilterCopyTransformCoupled.H
```

Insert 私有碰撞类 include 这个私有头。官方 `FilterCopyTransform.H` 恢复为接近
上游状态。

验收标准：

- 官方 `FilterCopyTransform.H` 不再依赖 Insert 配置。
- 私有 overload 只被 Insert 私有碰撞类使用。
- 私有 overload 不再包含 `MCC_DENSITY_MID` 分支。

### 4.5 接入新的碰撞类型

在官方 `CollisionHandler.cpp` 中增加最小分发：

```cpp
#include "Insert/CoupledBackgroundMCCCollision.H"

...
else if (type == "insert_background_mcc") {
    allcollisions[i] =
        std::make_unique<CoupledBackgroundMCCCollision>(collision_names[i]);
}
```

输入文件中私有耦合 MCC 使用：

```text
<collision_name>.type = insert_background_mcc
```

官方 `background_mcc` 保持官方语义，不再承载私有耦合模型。

验收标准：

- `CollisionHandler.cpp` 只新增 include 和一个类型分发分支。
- 私有输入文件已改用 `insert_background_mcc`。
- 官方 `background_mcc` 路径仍可独立编译和运行。

### 4.6 清理官方 BackgroundMCCCollision 的历史私有改动

在 Insert 私有类可用后，清理官方文件：

- 删除 `BackgroundMCCCollision.H` 中的私有成员：
  - `m_ground_rho_index`
  - `m_excitation_product`
  - `m_have_excitation`
  - `m_excitation_rho_index`
  - `m_sigma_max`
  - `m_ioni_sigma_max`
- 删除 `BackgroundMCCCollision.cpp` 中对 `global_background_density`、
  `BackgroundCoupledDensity`、`WarpXFunc.h` 和私有宏的直接依赖。
- 删除官方 `BackgroundMCCCollision.cpp` 中的耦合散射、电离和激发态替换逻辑。
- 保留官方普通 MCC 行为，尽量恢复接近上游。

验收标准：

- 耦合 MCC 的核心算法代码不再位于 `BackgroundMCCCollision.cpp`。
- `BackgroundMCCCollision` 不再有非官方耦合参数成员。
- 官方 MCC 文件仍可清楚表示官方 `background_mcc` 行为。
- 私有功能主要集中在 `Source/Insert`。
- `git diff --check` 和 `cmake --build build -j 8` 通过。

## 阶段 5：减少官方文件格式化 diff

目标：降低长期 merge/rebase 冲突。

1. 还原 `AddParticles.cpp` 中与功能无关的格式化改动。
2. 还原 include 顺序中非必要调整。
3. 对 `BackgroundMCCCollision.cpp` 只保留与私有模型接入有关的最小 diff。
4. 对 `LabFrameExplicitES.cpp`、`WarpXEvolve.cpp`、`PhysicalParticleContainer.cpp`
   保持小范围 hook 调用，不做邻近代码重排。

验收标准：

- `git diff origin/development...HEAD --stat` 中官方核心文件的 diff 明显缩小。
- 私有功能主要集中在 `Source/Insert`。

## 阶段 6：验证

每个阶段至少执行：

```bash
git diff --check
cmake --build build -j 8
```

涉及 Python bindings 时执行：

```bash
cmake --build build --target pip_install_nodeps -j 8
```

涉及碰撞模型时优先运行私有输入文件对应的小规模算例。若存在 checksum 差异，只把
物理量、粒子数、背景密度扣减和日志输出作为主要判断依据，不以官方 checksum 为目标。

## 推荐实施顺序

1. 阶段 1：先修正低风险问题。
2. 阶段 2：建立稳定 hook 层。
3. 阶段 4.4：先移动私有 `FilterCopyTransform` overload。
4. 阶段 4.1 到 4.3：移动耦合 MCC 核心模型。
5. 阶段 3：拆分注入、诊断、边界和二次电子发射。
6. 阶段 5：清理官方文件非功能性 diff。
7. 阶段 6：按功能分批验证。

## 风险点

- 耦合 MCC 模板函数包含 GPU lambda，移动时应优先保留 header-only 或 `.I`
  include 形式，避免 CUDA/HIP/SYCL 编译和链接问题。
- 不应为了移动代码而把 `BackgroundMCCCollision` 的大量私有成员改成 public。
  优先使用 context 结构显式传参。
- 背景粒子扣减依赖 cell bin、粒子权重和局部密度场同步，移动过程中应保持原有
  操作顺序。
- `PUSH_GAP`、碰撞 supercycle/subcycle、背景密度更新清理三者存在时序耦合，
  每次改动都要用小算例确认粒子数和背景密度扣减行为。
