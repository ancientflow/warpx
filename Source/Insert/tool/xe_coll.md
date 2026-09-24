为 WarpX 的 DSMC 模块生成 Xe–Xe 中性原子弹性碰撞截面文件。

要求只计算 VHS 碰撞截面，不实现 DSMC 碰撞算法。

## 1. VHS 参数

采用 Xe–Xe：

$$
m_{\rm Xe}=131.293\,u
$$

其中

$$
u=1.66053906660\times10^{-27}\ {\rm kg}.
$$

因此：

$$
m_{\rm Xe}\approx2.18017\times10^{-25}\ {\rm kg}.
$$

VHS 参数：

$$
T_{\rm ref}=273\ {\rm K}
$$

$$
d_{\rm ref}=5.74\times10^{-10}\ {\rm m}
$$

$$
\omega=0.85
$$

Xe–Xe 约化质量：

$$
\mu=\frac{m_{\rm Xe}}{2}.
$$

---

## 2. WarpX 输入能量定义

WarpX DSMC 截面文件第一列要求使用质心系相对动能：

$$
E_{\rm cm}
=
\frac12\mu g^2
$$

其中：

$$
g=|\mathbf v_1-\mathbf v_2|.
$$

由于文件中的能量单位为 eV：

$$
E_{\rm cm,J}=E_{\rm cm,eV}e
$$

其中：

$$
e=1.602176634\times10^{-19}\ {\rm J/eV}.
$$

因此：

$$
g=
\sqrt{
\frac{2E_{\rm cm,eV}e}{\mu}
}.
$$

---

## 3. VHS 截面

使用：

$$
\boxed{
\sigma_{\rm VHS}(g)
=
\frac{\pi d_{\rm ref}^2}
{\Gamma(2.5-\omega)}
\left[
\frac{2k_BT_{\rm ref}}
{\mu g^2}
\right]^{\omega-\frac12}
}
$$

利用：

$$
E_{\rm cm}=\frac12\mu g^2
$$

可直接写成能量形式：

$$
\boxed{
\sigma_{\rm VHS}(E_{\rm cm})
=
\frac{\pi d_{\rm ref}^2}
{\Gamma(2.5-\omega)}
\left(
\frac{k_BT_{\rm ref}}
{E_{\rm cm}}
\right)^{\omega-\frac12}
}
$$

注意这里：

- $k_BT_{\rm ref}$ 和 $E_{\rm cm}$ 必须使用相同能量单位；
- 推荐全部转换为 eV。

定义：

$$
E_{\rm ref}
=
\frac{k_BT_{\rm ref}}{e}.
$$

因此最终应直接使用：

$$
\boxed{
\sigma_{\rm VHS}(E_{\rm eV})
=
\frac{\pi d_{\rm ref}^2}
{\Gamma(2.5-\omega)}
\left(
\frac{E_{\rm ref,eV}}{E_{\rm eV}}
\right)^{\omega-\frac12}
}
$$

对于 Xe：

$$
\omega-\frac12=0.35
$$

所以：

$$
\boxed{
\sigma(E)\propto E^{-0.35}
}
$$

---

## 4. 输出文件格式

生成例如：

```text
Xe_Xe_VHS_elastic.dat
```

文件必须只有两列数值：

```text
E_eV sigma_m2
```

例如：

```text
1.0e-5  ...
1.1e-5  ...
1.2e-5  ...
...
1.0e+0  ...
```

第一列必须严格递增。

不要加入单位字符串。

注释行如果没有必要则不写，以避免解析兼容性问题。

---

## 5. 建议能量范围

针对室温至气体分配器喷孔中的中性 Xe，主要相对速度约为几十到上千 m/s，对应的质心系能量很低。

建议先覆盖：

$$
10^{-6}\ {\rm eV}
\le E_{\rm cm}
\le
10\ {\rm eV}.
$$

使用对数间隔。

例如生成：

```text
1000
```

个能量点。

这样既能覆盖低温 Xe–Xe 碰撞，也给非平衡高速尾部留出充分范围。

不要包含：

$$
E=0
$$

因为 VHS 模型在 $E\rightarrow0$ 时：

$$
\sigma(E)\rightarrow\infty.
$$

---

## 6. 数值检查

程序生成数据后输出若干测试值，例如：

```text
E = 1e-4 eV
E = 1e-3 eV
E = 1e-2 eV
E = 1e-1 eV
E = 1 eV
```

检查：

1. 所有截面均为正；
2. 能量越高，截面单调下降；
3. log-log 坐标下斜率应约为：

$$
-0.35.
$$

也就是：

$$
\frac{d\ln\sigma}{d\ln E}\approx-0.35.
$$

---

## 7. WarpX 输入

截面文件用于一个 DSMC elastic scattering process，例如：

```text
collisions.collision_names = XeXe

XeXe.type = dsmc
XeXe.species = neutral_Xe
XeXe.scattering_processes = elastic

XeXe.elastic_cross_section = Xe_Xe_VHS_elastic.dat
XeXe.elastic_scattering_angle_model = isotropic
```

具体 species 名称按照现有输入文件调整。

核心任务只是生成：

```text
E_cm[eV]    sigma[m^2]
```

这两列数据。