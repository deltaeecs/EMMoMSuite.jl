# 基函数

矩量法把连续积分方程离散成线性代数系统，离散质量在很大程度上取决于基函数与
测试函数的选择。EMMoMSuite 覆盖的核心基函数包括表面 RWG、体 SWG、屋顶 RBF、
分片常数 PWC（四面体/六面体）以及若干结构化网格基函数。本章公式均来自
论文第 2 章"常见基函数"一节。

## 1. 统一矢量基函数记号

RWG、SWG 与 RBF（屋顶）基函数具有统一形式（论文式 (2-34)~(2-36)）：

$$
\bm{f}_n(\bm{r}) =
\begin{cases}
\dfrac{a_n^{\pm}}{C_\Omega J^{\pm}}\, \bm{\rho}, & \bm{r} \in \Omega^{\pm} \\
\bm{0}, & \bm{r} \notin \Omega^{\pm}
\end{cases}
$$

其中：

- $\Omega^+$、$\Omega^-$ 为共享公共边/公共面的"正"、"负"网格单元；
- $a_n$ 为公共边长或公共面积，$a_n^{\pm} = \pm a_n$ 把方向符号转移到系数上；
- $J^{\pm}$ 为单元雅可比（三角形的面积、四面体/六面体的体积）；
- $\bm{\rho} = \bm{r} - \bm{r}_0^{\pm}$，$\bm{r}_0^{\pm}$ 为公共边/面的对侧
  自由顶点（六面体为场点映射到公共面对面上的点）；
- 常数 $C_\Omega$ 使得散度为常数（论文式 (2-36)）：

$$
\nabla \cdot \bm{f}_n(\bm{r}) = \frac{a_n^{\pm}}{J^{\pm}}
$$

因此 $C_\Omega$ 对三角形网格取 2、四面体取 3、其它网格一般取 1。
这种"符号提前折叠"的写法使同一基函数在正、负单元上有相同的表达形式，
矩阵填充时无需在循环内判断正负，与 EMMoMSuite 中 `edge_length`（无符号长度）、
`support`（支撑单元）、`local_edge_idx`（恢复对顶点）和 `signs`（局部符号）
的数据布局一致。

## 2. RWG 基函数

RWG 基函数定义在一对共享公共边的三角形 $T_n^+$、$T_n^-$ 上，经典形式为：

$$
\bm{f}_n(\bm{r}) =
\begin{cases}
\dfrac{l_n}{2A_n^+}(\bm{r} - \bm{v}_n^+), & \bm{r} \in T_n^+ \\
\dfrac{l_n}{2A_n^-}(\bm{v}_n^- - \bm{r}), & \bm{r} \in T_n^- \\
\bm{0}, & \text{otherwise}
\end{cases}
$$

其中 $l_n$ 为公共边长，$A_n^{\pm}$ 为支撑三角形面积，$\bm{v}_n^{\pm}$ 为对顶点。

**关键性质**：

1. 法向分量在公共边上连续，适合表示守恒电流；
2. 散度在每个支撑三角形上为常数，利于 EFIE 中矢量势项与标量势项离散。

### RWG-RWG 阻抗矩阵元

将 RWG 基函数代入 EFIE 的 $\bm{Z}^{SS}$ 定义，得到（论文式 (2-39)）：

$$
Z^{SS}_{mn} = {\rm j}k\eta \sum_{t \in \{+,-\}} \sum_{s \in \{+,-\}}
\left\{ \frac{l_m^t l_n^s}{4A_m^t A_n^s}
\int_{S_m^t} \int_{S_n^s}
\left[ \bm{\rho}_m(\bm{r}) \cdot \bm{\rho}_n(\bm{r}') - \frac{4}{k^2} \right]
G(R)\, dS' dS \right\}
$$

其中 $l_m^{\pm} = \pm l_m$、$l_n^{\pm} = \pm l_n$ 为带符号边长。该矩阵元具有
对称性 $Z^{SS}_{mn} = Z^{SS}_{nm}$，装配时可只计算一次填充两个位置。

## 3. SWG 基函数

SWG 基函数定义在一对共享公共面的四面体上，用于体积分方程与穿透介质问题：

$$
\bm{f}_n(\bm{r}) =
\begin{cases}
\dfrac{A_n}{3V_n^+}(\bm{r} - \bm{v}_n^+), & \bm{r} \in V_n^+ \\
\dfrac{A_n}{3V_n^-}(\bm{v}_n^- - \bm{r}), & \bm{r} \in V_n^- \\
\bm{0}, & \text{otherwise}
\end{cases}
$$

其中 $A_n$ 为公共面面积，$V_n^{\pm}$ 为四面体体积。SWG 与 RWG 思想一致：
跨公共面的通量连续、单元内散度为常数，对 VEFIE 或面体混合耦合问题保持局部守恒。

### SWG/RBF 的 $\bm{Z}^{VV}$ 矩阵元

将 SWG 或 RBF 基函数代入体 EFIE 得到六项积分（论文式 (2-40)）：

$$
\begin{aligned}
Z^{VV}_{mn} = & \sum_{t,s \in \{+,-\}} \Bigg\{
\frac{1}{{\rm j}\omega C_\Omega^2} \frac{a_m^t}{J_m^t} \frac{a_n^s}{J_n^s} \frac{1}{\varepsilon_n^s}
\int_{V_m^t} \bm{\rho}_m \cdot \bm{\rho}_n\, dV \\
& + \frac{{\rm j}\eta k}{C_\Omega^2} \frac{a_m^t}{J_m^t} \frac{a_n^s}{J_n^s} \kappa_n^s
\int_{V_m^t} \int_{V_n^s} \bm{\rho}_m \cdot \bm{\rho}_n\, G(R)\, dV' dV \\
& - \frac{{\rm j}\eta}{k} \frac{a_m^t}{J_m^t} \frac{a_n^s}{J_n^s} \kappa_n^s
\int_{V_m^t} \int_{V_n^s} G(R)\, dV' dV \\
& + \frac{(s)\,{\rm j}\eta}{k} \frac{a_m^t}{J_m^t} \kappa_n^s
\int_{V_m^t} \int_{\Gamma_n} G(R)\, dS' dV \\
& + \frac{(t)\,{\rm j}\eta}{k} \frac{a_n^s}{J_n^s} \kappa_n^s
\int_{\Gamma_m} \int_{V_n^s} G(R)\, dV' dS \\
& - \frac{(t)(s)\,{\rm j}\eta}{k} \kappa_n^s
\int_{\Gamma_m} \int_{\Gamma_n} G(R)\, dS' dS \Bigg\}
\end{aligned}
$$

其中 $\Gamma$ 为基函数的公共面，$\kappa_n = \varepsilon_n / \varepsilon_0$ 为相对
介电常数（介质对比度），$(t)$、$(s)$ 为正负号本身。前两项来自体 EFIE 的质量项
与 $\mathcal{L}$ 算子展开，后三项仅在测试/源函数为半基函数（位于介质边界）时
才需要计算。

## 4. 分片常数基函数 (PWC)

PWC 基函数对网格类型无限制，可处理四面体、六面体乃至混合剖分，因此更为普适。
其定义（论文式 (2-38)）为：

$$
\bm{f}_n(\bm{r}) = \bm{\hat{e}}_i, \qquad \bm{r} \in \Omega_n
$$

其中 $\bm{\hat{e}}_i$ 为 $x$、$y$ 或 $z$ 方向的单位向量，即每个体单元内
有 3 个朝向的 PWC 未知量。

### PWC 的 $\bm{Z}^{VV}$ 矩阵元

$$
Z^{VV}_{mn} = \int_{V_m} \frac{1}{{\rm j}\omega\varepsilon_n} \bm{\hat{e}}_m \cdot \bm{\hat{e}}_n\, dV
+ {\rm j}k\eta\, \bm{\hat{e}}_m \cdot \left[ \kappa_n
\int_{V_m} \int_{V_n} \left( \overline{\bm{I}} + \frac{1}{k^2}\nabla\nabla \right) G(R)\, dV' dV
\right] \cdot \bm{\hat{e}}_n
$$

其中的二重积分可作并矢展开以方便处理奇异性（论文式 (2-41)）：

$$
\left(\overline{\bm{I}} + \frac{1}{k^2}\nabla\nabla\right) G =
\left(\overline{\bm{I}} - \bm{\hat{R}}\bm{\hat{R}}\right) G
- \frac{\overline{\bm{I}} - 3\bm{\hat{R}}\bm{\hat{R}}}{k^2 R}
\left({\rm j}k + \frac{1}{R}\right) G(R)
$$

## 5. 屋顶基函数 (RBF)

屋顶基函数定义在一对共享公共面的六面体上，是 RWG/SWG 在规则六面体单元上的
对应形式。若公共面法向取 $x$ 方向，可写为：

$$
\bm{f}_n(\bm{r}) = \bm{\hat{x}}\, \Lambda(x)\, \Pi(y)\, \Pi(z)
$$

其中 $\Lambda(x)$ 为沿法向的屋顶函数，$\Pi(y)$、$\Pi(z)$ 为矩形脉冲函数。
其矩阵元公式与 SWG 共用第 3 节的六项统一形式（$C_\Omega = 1$）。

## 6. BC 基函数与高阶基函数

- **BC 基函数**是 RWG 的对偶测试空间（重心对偶网格上的旋转型流函数），
  用于改善 MFIE 与 Calderón 型离散的条件数与稳定性。
- **高阶基函数**（如 Graglia-Wilton-Peterson）在单元内引入更高次多项式，
  $p=1$ 时可退化到类似 RWG 的一阶表示，$p \geq 2$ 可表达更复杂的边沿与
  面内变化模式，在光滑几何上可降低总未知量规模。
