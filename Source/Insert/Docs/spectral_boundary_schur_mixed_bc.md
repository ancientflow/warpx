# 规则矩形域混合边界 Laplace 子问题的谱边界 Schur 方法

## 1. 问题背景

考虑规则矩形区域

\[
\Omega=(0,L_x)\times(0,L_y)\times(0,L_z).
\]

在静电问题中，总电势可按线性叠加分解为

\[
\phi=\phi_\rho+\phi_\sigma,
\]

其中 \(\phi_\rho\) 由体电荷密度 \(\rho\) 驱动，\(\phi_\sigma\) 用来表示边界表面电荷 \(\sigma_s\) 产生的附加场。本文只讨论 \(\phi_\sigma\) 对应的 Laplace 子问题：

\[
\nabla^2\phi=0 \quad \text{in } \Omega.
\]

设 \(x=0\) 面为混合边界：

\[
\Gamma_0=\{x=0,\ 0<y<L_y,\ 0<z<L_z\}
\]

并分为两个不重叠区域

\[
\Gamma_0=\Gamma_D\cup\Gamma_N,\qquad \Gamma_D\cap\Gamma_N=\emptyset.
\]

边界条件为

\[
\phi=0 \quad \text{on } \Gamma_D,
\]

\[
-\epsilon_0\frac{\partial \phi}{\partial n}=\sigma_s
\quad \text{on } \Gamma_N,
\]

其余五个面均为接地 Dirichlet：

\[
\phi=0
\quad \text{on }
x=L_x,\ y=0,\ y=L_y,\ z=0,\ z=L_z.
\]

这里 \(n\) 是从计算域 \(\Omega\) 指向外部的单位法向。在 \(x=0\) 面上，

\[
n=-\hat{x},\qquad \frac{\partial}{\partial n}=-\frac{\partial}{\partial x}.
\]

目标是在不修改三维 Poisson/MLMG 求解器非齐次 Neumann 支持的前提下，用一个只定义在 \(x=0\) 面上的二维边界问题表示 \(\phi_\sigma\)。

## 2. 基本思想

谱边界 Schur 方法的核心是：

1. 将整个 \(x=0\) 面上的边界电势迹记为

   \[
   u(y,z)=\phi(0,y,z).
   \]

2. 给定 \(u\)，可以在矩形域内解析地解出满足其余五面接地的 Laplace 解。

3. 因此可以定义一个 Dirichlet-to-Neumann 算子 \(\Lambda\)，把边界电势 \(u\) 映射为同一面上的外法向导数：

   \[
   \Lambda u=\left.\frac{\partial \phi}{\partial n}\right|_{x=0}.
   \]

4. 混合边界条件变成只在 \(x=0\) 面上的二维边界方程：

   \[
   u=0 \quad \text{on } \Gamma_D,
   \]

   \[
   -\epsilon_0\Lambda u=\sigma_s \quad \text{on } \Gamma_N.
   \]

这样，三维体内自由度被解析消元，只剩 Neumann 区域上的未知边界电势。

## 3. Dirichlet-to-Neumann 算子的连续推导

给定整个 \(x=0\) 面的 Dirichlet 数据 \(u(y,z)\)，考虑辅助问题：

\[
\nabla^2\phi=0,
\]

\[
\phi(0,y,z)=u(y,z),
\]

\[
\phi=0
\quad \text{on }
x=L_x,\ y=0,\ y=L_y,\ z=0,\ z=L_z.
\]

由于 \(y,z\) 方向均为接地 Dirichlet，使用双重 sine 展开：

\[
u(y,z)=
\sum_{m=1}^{\infty}\sum_{n=1}^{\infty}
u_{mn}
\sin\frac{m\pi y}{L_y}
\sin\frac{n\pi z}{L_z}.
\]

令

\[
k_{mn}=
\sqrt{
\left(\frac{m\pi}{L_y}\right)^2+
\left(\frac{n\pi}{L_z}\right)^2
}.
\]

对每个 \((m,n)\) 模式，体内解满足

\[
\frac{d^2 X_{mn}}{dx^2}-k_{mn}^2 X_{mn}=0,
\]

并满足

\[
X_{mn}(0)=1,\qquad X_{mn}(L_x)=0.
\]

因此

\[
X_{mn}(x)=
\frac{\sinh(k_{mn}(L_x-x))}
{\sinh(k_{mn}L_x)}.
\]

体内解为

\[
\phi(x,y,z)=
\sum_{m,n}
u_{mn}
\frac{\sinh(k_{mn}(L_x-x))}
{\sinh(k_{mn}L_x)}
\sin\frac{m\pi y}{L_y}
\sin\frac{n\pi z}{L_z}.
\]

对 \(x\) 求导：

\[
\left.\frac{\partial\phi}{\partial x}\right|_{x=0}
=
\sum_{m,n}
\left[
-k_{mn}\coth(k_{mn}L_x)u_{mn}
\right]
\sin\frac{m\pi y}{L_y}
\sin\frac{n\pi z}{L_z}.
\]

在 \(x=0\) 面上外法向导数为

\[
\frac{\partial\phi}{\partial n}
=-\frac{\partial\phi}{\partial x}.
\]

所以

\[
\Lambda u
=
\left.\frac{\partial\phi}{\partial n}\right|_{x=0}
=
\sum_{m,n}
\lambda_{mn}u_{mn}
\sin\frac{m\pi y}{L_y}
\sin\frac{n\pi z}{L_z},
\]

其中

\[
\lambda_{mn}=k_{mn}\coth(k_{mn}L_x).
\]

因此 \(\Lambda\) 在 sine 谱空间中是对角算子。

## 4. 混合边界 Schur 方程

定义 \(P_N\) 为从整个 \(x=0\) 面提取 Neumann 区域 \(\Gamma_N\) 上自由度的限制算子。其转置 \(P_N^T\) 将 Neumann 区域上的未知量放回整个面，并在 Dirichlet 区域填零。

令

\[
u_N=u|_{\Gamma_N}.
\]

由于 \(\Gamma_D\) 接地，整个面上的电势可写为

\[
u=P_N^T u_N.
\]

在 Neumann 区域要求

\[
-\epsilon_0\Lambda u=\sigma_s.
\]

限制到 \(\Gamma_N\) 得到边界 Schur 方程：

\[
-\epsilon_0 P_N\Lambda P_N^T u_N=\sigma_N,
\]

其中

\[
\sigma_N=\sigma_s|_{\Gamma_N}.
\]

也可写为正定形式

\[
A u_N=b,
\]

\[
A=\epsilon_0 P_N\Lambda P_N^T,
\]

\[
b=-\sigma_N.
\]

如果直接以外法向导数

\[
g_N=\left.\frac{\partial \phi}{\partial n}\right|_{\Gamma_N}
\]

作为输入，则边界方程为

\[
P_N\Lambda P_N^T u_N=g_N.
\]

若采用表面电荷约定

\[
\sigma_s=-\epsilon_0\frac{\partial\phi}{\partial n},
\]

则

\[
\epsilon_0 P_N\Lambda P_N^T u_N=-\sigma_N.
\]

只要 \(\Gamma_D\) 非空，且 Neumann 区域没有完全脱离接地区域，离散后的 \(A\) 通常为对称正定矩阵，可使用 CG 求解。

## 5. 算子应用流程

边界矩阵 \(A\) 不需要显式组装。对任意 Neumann 区未知向量 \(v_N\)，计算

\[
w_N=A v_N
\]

的流程如下。

1. 扩展到整个 \(x=0\) 面：

   \[
   v=P_N^T v_N.
   \]

   即 Neumann 区域填 \(v_N\)，Dirichlet 区域填 0。

2. 对 \(v(y,z)\) 做二维 sine transform，得到 \(v_{mn}\)。

3. 在谱空间逐模态相乘：

   \[
   \widehat{q}_{mn}=\epsilon_0\lambda_{mn}v_{mn}.
   \]

4. 做二维 inverse sine transform，得到

   \[
   q(y,z)=\epsilon_0\Lambda v.
   \]

5. 限制到 Neumann 区域：

   \[
   w_N=P_N q.
   \]

在伪代码中：

```text
function apply_A(v_N):
    v_face = 0 on whole xmin face
    v_face[NeumannMask] = v_N

    v_hat = DST2(v_face)
    q_hat[m,n] = eps0 * lambda[m,n] * v_hat[m,n]
    q_face = IDST2(q_hat)

    return q_face[NeumannMask]
```

然后用 CG 求解

```text
A u_N = -sigma_N
```

即可得到混合边界下 \(x=0\) 面 Neumann 区域的未知电势。

## 6. 体内电势重建

求得 \(u_N\) 后，将其扩展为整个面上的 \(u=P_N^T u_N\)，做 sine transform 得到 \(u_{mn}\)。

体内任意 \(x_i\) 处：

\[
\phi(x_i,y,z)=
\sum_{m,n}
u_{mn}
R_{mn}(x_i)
\sin\frac{m\pi y}{L_y}
\sin\frac{n\pi z}{L_z},
\]

其中传播因子为

\[
R_{mn}(x)=
\frac{\sinh(k_{mn}(L_x-x))}
{\sinh(k_{mn}L_x)}.
\]

数值上可对每个 \(x_i\) 层做：

```text
phi_hat[m,n] = R[m,n,x_i] * u_hat[m,n]
phi[x_i,:,:] = IDST2(phi_hat)
```

## 7. 电场重建

电场为

\[
\mathbf{E}=-\nabla\phi.
\]

### 7.1 法向分量

\[
\frac{\partial \phi}{\partial x}
=
\sum_{m,n}
u_{mn}
R'_{mn}(x)
\sin\frac{m\pi y}{L_y}
\sin\frac{n\pi z}{L_z},
\]

其中

\[
R'_{mn}(x)=
-k_{mn}
\frac{\cosh(k_{mn}(L_x-x))}
{\sinh(k_{mn}L_x)}.
\]

因此

\[
E_x=-\frac{\partial\phi}{\partial x}
=
\sum_{m,n}
u_{mn}
k_{mn}
\frac{\cosh(k_{mn}(L_x-x))}
{\sinh(k_{mn}L_x)}
\sin\frac{m\pi y}{L_y}
\sin\frac{n\pi z}{L_z}.
\]

在 \(x=0\) 面：

\[
E_x(0)=\Lambda u.
\]

因为 \(n=-\hat{x}\)，法向电场为

\[
E_n=\mathbf{E}\cdot n=-E_x.
\]

表面电荷条件

\[
\sigma_s=\epsilon_0 E_n
\]

等价于

\[
\sigma_s=-\epsilon_0\Lambda u.
\]

### 7.2 切向分量

\[
\frac{\partial \phi}{\partial y}
=
\sum_{m,n}
u_{mn}R_{mn}(x)
\frac{m\pi}{L_y}
\cos\frac{m\pi y}{L_y}
\sin\frac{n\pi z}{L_z},
\]

\[
\frac{\partial \phi}{\partial z}
=
\sum_{m,n}
u_{mn}R_{mn}(x)
\frac{n\pi}{L_z}
\sin\frac{m\pi y}{L_y}
\cos\frac{n\pi z}{L_z}.
\]

因此切向场可用 sine/cosine 混合变换，或在物理空间中对重建后的 \(\phi\) 做有限差分。若目标是与 WarpX 现有场布局一致，通常后者实现更直接。

## 8. 离散实现注意事项

### 8.1 网格位置

若 \(y,z\) 方向的接地边界位于网格端点，sine transform 作用于内部点。对 \(N_y,N_z\) 个内部自由度，常见离散点为

\[
y_j=j\Delta y,\quad j=1,\ldots,N_y,
\]

\[
z_k=k\Delta z,\quad k=1,\ldots,N_z.
\]

端点 \(y=0,L_y,z=0,L_z\) 为 Dirichlet 0，不作为 transform 自由度。

WarpX/AMReX 的 nodal 数据包括端点节点；若直接耦合，需要明确：

- 哪些节点属于 sine transform 的内部自由度；
- 端点接地节点是否固定为 0；
- \(\Gamma_D\) 与 \(\Gamma_N\) mask 是否只作用于内部面节点。

### 8.2 连续符号与离散符号

本文使用外法向导数：

\[
\Lambda u=\partial_n\phi.
\]

在 \(x=0\) 面，

\[
\partial_n=-\partial_x.
\]

若实现中使用 \(d\phi/dx\)，则 DtN 符号相反：

\[
\left.\frac{\partial\phi}{\partial x}\right|_{x=0}
=-\Lambda u.
\]

表面电荷条件若写作

\[
\sigma_s=-\epsilon_0\partial_n\phi,
\]

则边界方程为

\[
\epsilon_0 P_N\Lambda P_N^T u_N=-\sigma_N.
\]

### 8.3 稳定计算 \(\lambda_{mn}\)

\[
\lambda_{mn}=k_{mn}\coth(k_{mn}L_x).
\]

当 \(k_{mn}L_x\) 很大时，\(\coth(k_{mn}L_x)\approx 1\)。不要直接用容易溢出的 \(\sinh,\cosh\) 计算。可用

\[
\coth a=\frac{1+e^{-2a}}{1-e^{-2a}}.
\]

当 \(a\) 很大时直接取 \(\coth a=1\)。

### 8.4 体内传播因子

传播因子

\[
R_{mn}(x)=
\frac{\sinh(k_{mn}(L_x-x))}
{\sinh(k_{mn}L_x)}
\]

也建议用指数形式：

\[
R_{mn}(x)
=
e^{-k_{mn}x}
\frac{1-e^{-2k_{mn}(L_x-x)}}{1-e^{-2k_{mn}L_x}}
\]

或等价稳定形式，避免 \(k_{mn}L_x\) 大时溢出。

### 8.5 预处理器

边界 Schur 系统的矩阵为

\[
A=\epsilon_0 P_N\Lambda P_N^T.
\]

若 \(\Gamma_N\) 是整个 \(x=0\) 面，\(P_N\) 为恒等，\(A\) 在 sine 空间对角，可直接求解。

若 \(\Gamma_N\) 是面内一部分，mask 会导致模式耦合，需要 CG/GMRES。常用预处理：

- 用整面 DtN 的逆作为近似预处理；
- 用对角近似；
- 对固定 patch mask 预构造低秩或稀疏近似预条件器。

## 9. 复杂度与存储

令

\[
N_f=N_yN_z
\]

为 \(x=0\) 面上的自由度数量。

一次 Schur 算子应用需要一次 2D DST 和一次 inverse DST：

\[
O(N_f\log N_f).
\]

若 CG 迭代次数为 \(N_{\text{iter}}\)，边界求解成本为

\[
O(N_{\text{iter}}N_f\log N_f).
\]

体内重建若对每个 \(x\) 层做 inverse DST：

\[
O(N_xN_f\log N_f).
\]

对 \(256^3\) 网格：

\[
N_f=256^2=65536.
\]

建议预存：

- \(k_{mn}\)，大小 \(O(N_f)\)；
- \(\lambda_{mn}\)，大小 \(O(N_f)\)；
- transform plan；
- mask；
- 可选的 \(R_{mn}(x_i)\)，大小 \(O(N_xN_f)\)。

不建议预存每个表面模式对应的完整三维 Green 响应场，因为存储量为

\[
O(N_f\cdot N_xN_yN_z),
\]

对 \(256^3\) 已达到 TB 级别。

## 10. 与三维多重网格的比较

三维多重网格直接在体内求解，理论复杂度约为

\[
O(N_xN_yN_z).
\]

谱边界 Schur 方法把混合边界 Laplace 子问题降到二维边界方程：

\[
O(N_{\text{iter}}N_yN_z\log(N_yN_z))
\]

加上必要时的体内重建成本。

因此它适合以下情况：

- 几何为规则矩形；
- 其余五面为接地 Dirichlet；
- 混合边界 patch 固定；
- 表面电荷 \(\sigma_s\) 随时间变化，但边界类型和几何不变；
- 需要把表面电荷贡献作为附加场叠加到 WarpX 原 Poisson 解上。

它不适合以下情况：

- 复杂 EB 几何；
- 非矩形边界；
- 变系数介电常数；
- AMR patch 上需要严格局部守恒；
- 混合边界区域剧烈变化并导致预处理失效。

## 11. 推荐数值流程

预处理阶段：

```text
1. 读取几何 Lx, Ly, Lz 和网格尺寸。
2. 构造 xmin 面上的 Dirichlet/Neumann mask。
3. 建立 2D DST/IDST plan。
4. 预计算 k[m,n]。
5. 预计算 lambda[m,n] = k[m,n] coth(k[m,n] Lx)。
6. 可选：预计算 R[m,n,i] 和 dRdx[m,n,i]。
7. 构造 Schur 算子 apply_A。
8. 可选：构造预处理器。
```

每个时间步：

```text
1. 从表面充电模型得到 sigma_s(y,z)。
2. 在 Neumann mask 上取 RHS: b = -sigma_s。
3. 用 CG 解 A u_N = b。
4. 将 u_N 扩展为全 xmin 面电势 u。
5. 用 DST 得到 u_hat。
6. 重建 phi_sigma 或 E_sigma。
7. 将 E_sigma 叠加到 WarpX 原有电场。
```

若只需要场而不需要电势，可直接重建 \(\mathbf{E}_\sigma\)，避免存储完整 \(\phi_\sigma\)。

## 12. 小结

谱边界 Schur 方法利用规则矩形域的可分离结构，将三维 Laplace 混合边界问题压缩为 \(x=0\) 面上的二维边界问题。其数学核心是 Dirichlet-to-Neumann 算子

\[
\Lambda=
\operatorname{IDST}
\operatorname{diag}\left(k_{mn}\coth(k_{mn}L_x)\right)
\operatorname{DST}.
\]

混合边界方程为

\[
\epsilon_0 P_N\Lambda P_N^T u_N=-\sigma_N.
\]

该方法避免了直接修改 WarpX/AMReX nodal Poisson 求解器以支持非齐次 Neumann 边界，同时能正确体现接地 Dirichlet patch 的感应电荷效应。对于规则矩形、固定混合边界区域、表面电荷随时间变化的问题，它是一个存储可控、可预处理、并且与原有 WarpX 静电解可线性叠加的方案。
