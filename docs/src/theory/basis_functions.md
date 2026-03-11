# 基函数 (Basis Functions)

矩量法将连续的积分方程转化为离散的矩阵方程，其核心在于将未知电流密度展开为一组基函数的线性组合。

## 1. RWG 基函数 (Rao-Wilton-Glisson)

RWG 基函数是定义在平面三角形网格上的矢量基函数，广泛用于表面积分方程 (EFIE, MFIE, CFIE)。

### 1.1 定义
对于第 ``n`` 条公共边，连接两个三角形 ``T_n^+`` 和 ``T_n^-``。
```math
\mathbf{f}_n(\mathbf{r}) = \begin{cases}
\frac{l_n}{2A_n^+} (\mathbf{r} - \mathbf{v}_n^+), & \mathbf{r} \in T_n^+ \\
\frac{l_n}{2A_n^-} (\mathbf{v}_n^- - \mathbf{r}), & \mathbf{r} \in T_n^- \\
0, & \text{otherwise}
\end{cases}
```
其中 ``l_n`` 是边长，``A_n^\pm`` 是面积，``\mathbf{v}_n^\pm`` 是相对顶点。

### 1.2 与 EMSuite 实现一致的统一写法

教材中常把 RWG 写成正半基函数与负半基函数两段。EMSuite 的实现则把“正负号”单独存到 support-local 的符号 ``s_{n,i} \in \{+1,-1\}`` 中，并把公共边长始终存为正数 ``l_n > 0``。于是可以把两段公式统一写成
```math
\mathbf{f}_n(\mathbf{r}) = \begin{cases}
\dfrac{\tilde l_{n,i}}{2 A_{n,i}} \left( \mathbf{r} - \mathbf{v}_{n,i}^{\mathrm{opp}} \right), & \mathbf{r} \in T_{n,i},\; i \in \{1,2\} \\
0, & \text{otherwise}
\end{cases}
```
其中
```math
	ilde l_{n,i} = s_{n,i} \, l_n,
\qquad
s_{n,1} = +1,
\qquad
s_{n,2} = -1.
```

这种写法与传统的 ``T_n^+ / T_n^-`` 二分写法完全等价，但在推导实现公式时更贴近仓库内部数据结构：

*   `edge_length` 始终是无符号公共边长 ``l_n``。
*   `support[i]` 给出第 ``i`` 个支撑三角形 ``T_{n,i}``。
*   `local_edge_idx[i]` 用来确定该三角形内与公共边相对的顶点 ``\mathbf{v}_{n,i}^{\mathrm{opp}}``。
*   `signs[i]` 提供 ``s_{n,i}``，用于统一处理正半/负半的方向差异。

因此，在实现相关的公式中，把符号“藏”进带符号边长 ``\tilde l_{n,i}`` 往往比显式写正负半基函数更直接。

### 1.3 性质
1.  **法向连续性**：流过公共边的电流法向分量连续，无电荷积累。
2.  **散度有限**：
    ```math
    \nabla \cdot \mathbf{f}_n(\mathbf{r}) = \begin{cases}
    \frac{l_n}{A_n^+}, & \mathbf{r} \in T_n^+ \\
    -\frac{l_n}{A_n^-}, & \mathbf{r} \in T_n^-
    \end{cases}
    ```
    每个三角形面上的总电荷为零（整体电中性）。

对于统一写法，上式也可以写成
```math
\nabla \cdot \mathbf{f}_n(\mathbf{r}) = \frac{\tilde l_{n,i}}{A_{n,i}},
\qquad \mathbf{r} \in T_{n,i},\; i \in \{1,2\}.
```
这样在 EFIE 的标量势项和电荷项推导中，符号会自然跟随 ``\tilde l_{n,i}`` 进入局部系数。

### 1.4 装配时的支撑配对

一个 RWG-RWG 矩阵元不是“边对边”的单个积分，而是四个支撑子三角形配对的求和：
```math
Z_{mn} = \sum_{i=1}^2 \sum_{j=1}^2 Z_{mn}^{(i,j)},
```
其中每个子项都只在 ``T_{m,i} \times T_{n,j}`` 上积分。若采用上面的带符号边长记号，则每个子项都可以直接复用同一套几何公式，而不再额外区分正、负半基函数。

## 2. SWG 基函数 (Schaubert-Wilton-Glisson)

SWG 基函数是定义在四面体网格上的矢量基函数，用于体积分方程 (VIE) 或穿透性介质问题。

### 2.1 定义
对于第 ``n`` 个公共面，连接两个四面体 ``T_n^+`` 和 ``T_n^-``。
```math
\mathbf{f}_n(\mathbf{r}) = \begin{cases}
\frac{A_n}{3V_n^+} (\mathbf{r} - \mathbf{v}_n^+), & \mathbf{r} \in T_n^+ \\
\frac{A_n}{3V_n^-} (\mathbf{v}_n^- - \mathbf{r}), & \mathbf{r} \in T_n^- \\
0, & \text{otherwise}
\end{cases}
```
其中 ``A_n`` 是公共面面积，``V_n^\pm`` 是体积。

### 2.2 性质
类似于 RWG，SWG 保证了穿过公共面的通量连续性，且散度在四面体内为常数。

在实现层面，SWG 也常采用“公共面面积为正、方向由局部符号或顶点顺序决定”的存储方式。因此，上面对 RWG 的“几何尺度与方向分离”思路同样适用于 SWG、PWC 和混合面-体耦合推导。

## 3. 分片常数基函数 (Piecewise Constant Basis Functions)

分片常数 (PWC) 基函数通常用于展开标量物理量（如电荷密度、标量电位）或作为测试函数。

### 3.1 定义
在每个单元 ``T_n`` 上定义脉冲函数：
```math
P_n(\mathbf{r}) = \begin{cases}
1, & \mathbf{r} \in T_n \\
0, & \text{otherwise}
\end{cases}
```

### 3.2 应用
*   **电荷密度展开**：在求解混合位积分方程时，有时独立展开电荷密度 ``\rho(\mathbf{r}) = \sum q_n P_n(\mathbf{r})``。
*   **点匹配法 (Point Matching)**：作为测试函数时，相当于在单元中心进行狄拉克 ``\delta`` 采样。

## 4. 屋顶基函数 (Rooftop Basis Functions)

屋顶基函数定义在矩形（2D）或六面体（3D）网格上，是 RWG/SWG 在结构化网格上的对应形式。

### 4.1 定义 (3D 六面体)
对于连接两个六面体单元 ``V_n^+`` 和 ``V_n^-`` 的公共面（假设法向为 ``x`` 方向）：
```math
\mathbf{f}_n(\mathbf{r}) = \hat{x} \Lambda(x) \Pi(y) \Pi(z)
```
其中 ``\Lambda(x)`` 是三角形函数（屋顶形状），``\Pi(y), \Pi(z)`` 是矩形脉冲函数。
```math
\Lambda(x) = \begin{cases}
\frac{x - x_{i-1}}{\Delta x}, & x_{i-1} \le x \le x_i \\
\frac{x_{i+1} - x}{\Delta x}, & x_i \le x \le x_{i+1}
\end{cases}
```

### 4.2 性质
*   **分量解耦**：通常 ``\mathbf{f}_x, \mathbf{f}_y, \mathbf{f}_z`` 分别定义在不同的边/面上。
*   **散度有限**：类似于 RWG，保证法向连续。
*   **应用**：常用于 FDTD 算法或基于体素的矩量法 (Voxel-based MoM)。

## 5. BC 基函数 (Buffa-Christiansen)

BC 基函数是 RWG 基函数的对偶基函数，定义在重心对偶网格（Barycentric Dual Mesh）上。

### 5.1 应用场景
在 MFIE 的离散中，如果测试函数和基函数都使用 RWG，会导致矩阵条件数不佳。采用 BC 基函数作为测试函数（RWG 作为基函数），可以显著改善 MFIE 的精度和稳定性，这被称为 **Galerkin MFIE** 的修正形式。

### 5.2 构造
BC 基函数是 RWG 基函数的线性组合，但在几何上旋转了 90 度，模拟磁流的流动特性。

## 6. 高阶基函数 (Higher-Order Basis Functions)

为了提高精度或减少未知量数目，可在每个单元上使用高阶多项式基函数。

### 6.1 Graglia-Wilton-Peterson (GWP) 基函数
定义在弯曲三角形上，能够更好地拟合曲面几何。
阶数 ``p`` 决定了基函数的多项式次数。
*   ``p=1``：类似 RWG。
*   ``p=2``：包含边上的高阶变化和面内的涡旋模式。

### 6.2 优势
*   **几何拟合**：配合曲面网格，减少几何离散误差。
*   **收敛速度**：具有指数收敛特性（``hp``-refinement）。
