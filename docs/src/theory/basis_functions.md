# 基函数 (Basis Functions)

矩量法将连续的积分方程转化为离散的矩阵方程，其核心在于将未知电流密度展开为一组基函数的线性组合。

## 1. RWG 基函数 (Rao-Wilton-Glisson)

RWG 基函数是定义在平面三角形网格上的矢量基函数，广泛用于表面积分方程 (EFIE, MFIE, CFIE)。

### 1.1 定义
对于第 \$n\$ 条公共边，连接两个三角形 \$T_n^+\$ 和 \$T_n^-\$。
\$\$
\mathbf{f}_n(\mathbf{r}) = \begin{cases}
\frac{l_n}{2A_n^+} (\mathbf{r} - \mathbf{v}_n^+), & \mathbf{r} \in T_n^+ \\
\frac{l_n}{2A_n^-} (\mathbf{v}_n^- - \mathbf{r}), & \mathbf{r} \in T_n^- \\
0, & \text{otherwise}
\end{cases}
\$\$
其中 \$l_n\$ 是边长，\$A_n^\pm\$ 是面积，\$\mathbf{v}_n^\pm\$ 是相对顶点。

### 1.2 性质
1.  **法向连续性**：流过公共边的电流法向分量连续，无电荷积累。
2.  **散度有限**：
    \$\$
    \nabla \cdot \mathbf{f}_n(\mathbf{r}) = \begin{cases}
    \frac{l_n}{A_n^+}, & \mathbf{r} \in T_n^+ \\
    -\frac{l_n}{A_n^-}, & \mathbf{r} \in T_n^-
    \end{cases}
    \$\$
    每个三角形面上的总电荷为零（整体电中性）。

## 2. SWG 基函数 (Schaubert-Wilton-Glisson)

SWG 基函数是定义在四面体网格上的矢量基函数，用于体积分方程 (VIE) 或穿透性介质问题。

### 2.1 定义
对于第 \$n\$ 个公共面，连接两个四面体 \$T_n^+\$ 和 \$T_n^-\$。
\$\$
\mathbf{f}_n(\mathbf{r}) = \begin{cases}
\frac{A_n}{3V_n^+} (\mathbf{r} - \mathbf{v}_n^+), & \mathbf{r} \in T_n^+ \\
\frac{A_n}{3V_n^-} (\mathbf{v}_n^- - \mathbf{r}), & \mathbf{r} \in T_n^- \\
0, & \text{otherwise}
\end{cases}
\$\$
其中 \$A_n\$ 是公共面面积，\$V_n^\pm\$ 是体积。

### 2.2 性质
类似于 RWG，SWG 保证了穿过公共面的通量连续性，且散度在四面体内为常数。

## 3. 分片常数基函数 (Piecewise Constant Basis Functions)

分片常数 (PWC) 基函数通常用于展开标量物理量（如电荷密度、标量电位）或作为测试函数。

### 3.1 定义
在每个单元 \$T_n\$ 上定义脉冲函数：
\$\$
P_n(\mathbf{r}) = \begin{cases}
1, & \mathbf{r} \in T_n \\
0, & \text{otherwise}
\end{cases}
\$\$

### 3.2 应用
*   **电荷密度展开**：在求解混合位积分方程时，有时独立展开电荷密度 \$\rho(\mathbf{r}) = \sum q_n P_n(\mathbf{r})\$。
*   **点匹配法 (Point Matching)**：作为测试函数时，相当于在单元中心进行狄拉克 \$\delta\$ 采样。

## 4. 屋顶基函数 (Rooftop Basis Functions)

屋顶基函数定义在矩形（2D）或六面体（3D）网格上，是 RWG/SWG 在结构化网格上的对应形式。

### 4.1 定义 (3D 六面体)
对于连接两个六面体单元 \$V_n^+\$ 和 \$V_n^-\$ 的公共面（假设法向为 \$x\$ 方向）：
\$\$
\mathbf{f}_n(\mathbf{r}) = \hat{x} \Lambda(x) \Pi(y) \Pi(z)
\$\$
其中 \$\Lambda(x)\$ 是三角形函数（屋顶形状），\$\Pi(y), \Pi(z)\$ 是矩形脉冲函数。
\$\$
\Lambda(x) = \begin{cases}
\frac{x - x_{i-1}}{\Delta x}, & x_{i-1} \le x \le x_i \\
\frac{x_{i+1} - x}{\Delta x}, & x_i \le x \le x_{i+1}
\end{cases}
\$\$

### 4.2 性质
*   **分量解耦**：通常 \$\mathbf{f}_x, \mathbf{f}_y, \mathbf{f}_z\$ 分别定义在不同的边/面上。
*   **散度有限**：类似于 RWG，保证法向连续。
*   **应用**：常用于 FDTD 算法或基于体素的矩量法 (Voxel-based MoM)。

## 5. BC 基函数 (Buffa-Christiansen)

BC 基函数是 RWG 基函数的对偶基函数，定义在重心对偶网格（Barycentric Dual Mesh）上。

### 3.1 应用场景
在 MFIE 的离散中，如果测试函数和基函数都使用 RWG，会导致矩阵条件数不佳。采用 BC 基函数作为测试函数（RWG 作为基函数），可以显著改善 MFIE 的精度和稳定性，这被称为 **Galerkin MFIE** 的修正形式。

### 3.2 构造
BC 基函数是 RWG 基函数的线性组合，但在几何上旋转了 90 度，模拟磁流的流动特性。

## 4. 高阶基函数 (Higher-Order Basis Functions)

为了提高精度或减少未知量数目，可在每个单元上使用高阶多项式基函数。

### 4.1 Graglia-Wilton-Peterson (GWP) 基函数
定义在弯曲三角形上，能够更好地拟合曲面几何。
阶数 \$p\$ 决定了基函数的多项式次数。
*   \$p=1\$：类似 RWG。
*   \$p=2\$：包含边上的高阶变化和面内的涡旋模式。

### 4.2 优势
*   **几何拟合**：配合曲面网格，减少几何离散误差。
*   **收敛速度**：具有指数收敛特性（\$hp\$-refinement）。
