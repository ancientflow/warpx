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
