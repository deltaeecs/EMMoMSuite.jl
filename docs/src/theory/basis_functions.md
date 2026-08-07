# 基函数

矩量法把连续积分方程离散成线性代数系统，而离散质量在很大程度上取决于基函数与测试函数的选择。EMMoMSuite 当前覆盖的核心基函数包括表面 RWG、体 SWG 以及若干分片常数和结构化网格基函数。

## 1. RWG 基函数

RWG 基函数定义在共享公共边的两个三角形上，是表面积分方程中最常见的电流展开基。

### 1.1 经典定义

设第 n 条公共边连接两个三角形 $T_n^+$ 与 $T_n^-$，则 RWG 基函数写为

$$
\mathbf{f}_n(\mathbf{r}) = \begin{cases}
\dfrac{l_n}{2 A_n^+} (\mathbf{r} - \mathbf{v}_n^+), & \mathbf{r} \in T_n^+ \\
\dfrac{l_n}{2 A_n^-} (\mathbf{v}_n^- - \mathbf{r}), & \mathbf{r} \in T_n^- \\
0, & \text{otherwise}
\end{cases}
$$

其中 $l_n$ 为公共边长，$A_n^\pm$ 为两个支撑三角形的面积，$\mathbf{v}_n^\pm$ 为对应的对顶点。

### 1.2 与 EMMoMSuite 实现一致的统一写法

EMMoMSuite 在内部把几何尺度和方向符号拆开存储。对第 n 个 RWG，可以把两个支撑三角形统一记为 $T_{n,i}$，其中 $i \in \{1,2\}$，并定义带符号边长

$$
	ilde l_{n,i} = s_{n,i} \, l_n,
\qquad
s_{n,1} = +1,
\qquad
s_{n,2} = -1.
$$

于是 RWG 基函数可统一写为

$$
\mathbf{f}_n(\mathbf{r}) = \begin{cases}
\dfrac{\tilde l_{n,i}}{2 A_{n,i}} \left( \mathbf{r} - \mathbf{v}_{n,i}^{\mathrm{opp}} \right), & \mathbf{r} \in T_{n,i},\; i \in \{1,2\} \\
0, & \text{otherwise}
\end{cases}
$$

这种写法与经典正半边和负半边的表示完全等价，但更接近实现中的数据布局：

- `edge_length` 存储无符号公共边长 $l_n$。
- `support[i]` 指向第 i 个支撑三角形。
- `local_edge_idx[i]` 用于恢复对顶点 $\mathbf{v}_{n,i}^{\mathrm{opp}}$。
- `signs[i]` 提供局部方向符号 $s_{n,i}$。

### 1.3 关键性质

1. 法向分量在公共边上连续，因此适合表示守恒电流。
2. 散度在每个支撑三角形上为常数，有利于 EFIE 中标量势项与电荷项的离散。

经典形式下，散度为

$$
\nabla \cdot \mathbf{f}_n(\mathbf{r}) = \begin{cases}
\dfrac{l_n}{A_n^+}, & \mathbf{r} \in T_n^+ \\
-\dfrac{l_n}{A_n^-}, & \mathbf{r} \in T_n^-
\end{cases}
$$

统一写法下则变成

$$
\nabla \cdot \mathbf{f}_n(\mathbf{r}) = \dfrac{\tilde l_{n,i}}{A_{n,i}},
\qquad
\mathbf{r} \in T_{n,i},\; i \in \{1,2\}.
$$

### 1.4 装配视角

一个 RWG-RWG 矩阵元通常不是单个面元积分，而是四个支撑子三角形配对项的求和：

$$
Z_{mn} = \sum_{i=1}^{2} \sum_{j=1}^{2} Z_{mn}^{(i,j)}.
$$

因此实现时最关键的是把几何、符号和局部积分核组织成可复用的统一公式。

## 2. SWG 基函数

SWG 基函数定义在共享公共面的两个四面体上，常用于体积分方程和穿透介质问题。

### 2.1 定义

设第 n 个公共面连接四面体 $V_n^+$ 与 $V_n^-$，则对应基函数可写为

$$
\mathbf{f}_n(\mathbf{r}) = \begin{cases}
\dfrac{A_n}{3 V_n^+} (\mathbf{r} - \mathbf{v}_n^+), & \mathbf{r} \in V_n^+ \\
\dfrac{A_n}{3 V_n^-} (\mathbf{v}_n^- - \mathbf{r}), & \mathbf{r} \in V_n^- \\
0, & \text{otherwise}
\end{cases}
$$

其中 $A_n$ 为公共面面积，$V_n^\pm$ 为四面体体积。

### 2.2 性质

SWG 与 RWG 的思想一致，核心特征是跨公共面的通量连续，且散度在单元内为常数。对 VEFIE 或面体混合耦合问题，这种局部守恒结构很重要。

## 3. 分片常数基函数

分片常数基函数通常用于标量量的展开，例如电荷密度、标量势或辅助测试空间。

### 3.1 定义

在单元 $T_n$ 上定义

$$
P_n(\mathbf{r}) = \begin{cases}
1, & \mathbf{r} \in T_n \\
0, & \text{otherwise}
\end{cases}
$$

### 3.2 应用

- 用于独立展开电荷密度等标量未知量。
- 作为点匹配或单元平均测试的基础空间。

## 4. 屋顶基函数

屋顶基函数定义在结构化网格上，可以看作 RWG 和 SWG 在规则单元上的对应形式，常见于矩形面元或六面体体元。

### 4.1 一个典型形式

若公共面法向取 x 方向，则可写为

$$
\mathbf{f}_n(\mathbf{r}) = \hat{\mathbf{x}} \, \Lambda(x) \, \Pi(y) \, \Pi(z)
$$

其中 $\Lambda(x)$ 是沿法向的屋顶函数，$\Pi(y)$ 和 $\Pi(z)$ 是矩形脉冲函数。

### 4.2 适用场景

- 规则几何上的结构化离散。
- 体素化 MoM 或其他基于正交网格的离散方案。

## 5. BC 基函数

BC 基函数是 RWG 的对偶测试空间，常用于改进 MFIE 或 Calderon 相关离散的稳定性。

### 5.1 作用

若 MFIE 中基函数和测试函数都选 RWG，矩阵条件数与精度往往不理想。使用 BC 作为测试空间、RWG 作为展开空间，可以显著改善离散品质。

### 5.2 直观理解

BC 基函数可看作在重心对偶网格上的旋转型流函数，其几何行为更接近磁流测试空间。

## 6. 高阶基函数

高阶基函数通过在单元内引入更高次多项式，提高几何拟合与场分布逼近能力。

### 6.1 典型形式

Graglia-Wilton-Peterson 等高阶基函数常定义在曲面三角形上，通过阶数 $p$ 控制多项式复杂度。

- $p = 1$ 时可退化到类似 RWG 的一阶表示。
- $p = 2$ 及以上时可以表达更复杂的边沿和面内变化模式。

### 6.2 优点

- 在光滑几何上可减少几何离散误差。
- 在相同精度目标下，可能降低总未知量规模。
- 适合与曲面单元和 hp 细化策略配合使用。
