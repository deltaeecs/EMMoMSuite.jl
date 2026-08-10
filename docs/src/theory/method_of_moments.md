# 矩量法实现

本章聚焦矩量法（MoM）的数学原理与代码离散实现：伽辽金离散、各类积分方程的
矩阵元推导、高斯数值积分以及论文提出的"奇异值提取 + 递推公式"整体奇异性
处理方案。公式均来自论文第 2 章。

## 1. 伽辽金离散 (Galerkin Discretization)

设连续算子方程为（论文式 (2-16)）：

$$
L(\bm{f}) = \bm{g}
$$

将未知函数展开为基函数线性组合 $\bm{f} \approx \sum_{n=1}^{N} a_n \bm{f}_n$，
并用测试函数 $\bm{w}_m$ 取内积，得到线性代数系统（论文式 (2-17)~(2-19)）：

$$
\sum_{n=1}^{N} \langle \bm{w}_m, L(\bm{f}_n) \rangle\, a_n = \langle \bm{w}_m, \bm{g} \rangle
$$

即 $\bm{Z}\bm{I} = \bm{V}$，其中：

$$
Z_{mn} = \int \bm{w}_m \cdot L(\bm{f}_n)\, d\Omega, \qquad
V_m = \int \bm{w}_m \cdot \bm{g}\, d\Omega
$$

伽辽金法取 $\bm{w}_m = \bm{f}_m$，矩阵对称性最好。

## 2. 矩阵元素计算

### 2.1 EFIE 矩阵（RWG-RWG）

表面 EFIE 的阻抗矩阵元（论文式 (2-39)）：

$$
Z^{SS}_{mn} = {\rm j}k\eta \sum_{t,s \in \{+,-\}}
\frac{l_m^t l_n^s}{4A_m^t A_n^s}
\int_{S_m^t} \int_{S_n^s}
\left[ \bm{\rho}_m(\bm{r}) \cdot \bm{\rho}_n(\bm{r}') - \frac{4}{k^2} \right] G(R)\, dS' dS
$$

代码实现中（`calc_interaction!` / `calc_self_interaction!` /
`calc_near_interaction!`）把双面积分展开为 4 个支撑子三角形配对，先累加
$(\bm{\rho}_m \cdot \bm{\rho}_n - 4/k^2) G$ 的加权和，最后统一乘
$l_m l_n / (A_m A_n)$ 与全局因子 ${\rm j}k\eta/(16\pi)$。

### 2.2 MFIE 矩阵

MFIE 离散形式包含质量矩阵项与 $\mathcal{K}$ 算子主值项：

$$
Z^{SS}_{mn} = \eta \left[
\frac{1}{2}\int_{S_m} \bm{f}_m \cdot \bm{f}_n\, dS
- \int_{S_m} \bm{f}_m \cdot \left( \hat{\bm{n}} \times \int_{S_n}
\bm{f}_n \times \nabla' G\, dS' \right) dS
\right]
$$

其中 $\frac{1}{2}$ 项来自 $\mathcal{K}$ 算子在光滑表面上的主值跳跃项。
代码中 `calc_self_term!` 计算质量项（$\eta/(8A)$ 缩放），`calc_k_term_fast!`
计算 $\mathcal{K}$ 项，其核为 $(\bm{\rho}_m\cdot\bm{R})(\hat{\bm{n}}_t\cdot\bm{\rho}_n) - (\hat{\bm{n}}_t\cdot\bm{R})(\bm{\rho}_m\cdot\bm{\rho}_n)$
乘 $({\rm j}k + 1/R) e^{-{\rm j}kR}/R^2$。

### 2.3 VIE 矩阵（SWG/RBF 与 PWC）

- SWG/RBF 展开体等效电流：六项积分公式见 `basis_functions.md` 第 3 节
  （论文式 (2-40)）。
- PWC 展开：质量项 + 并矢格林函数积分（论文式 (2-42)）。
- 面-体耦合项 $Z^{SV}$、$Z^{VS}$（RWG-SWG 见论文式 (2-43)~(2-44)，
  RWG-PWC 见 (2-45)~(2-46)），配合 $\kappa_n = \varepsilon_n/\varepsilon_0$
  的比例系数填入面体混合矩阵：

$$
\begin{bmatrix}
Z_{SS} & Z_{SV} \\
Z_{VS} & Z_{VV}
\end{bmatrix}
\begin{bmatrix} I_S \\ I_V \end{bmatrix}
=
\begin{bmatrix} V_S \\ V_V \end{bmatrix}
$$

### 2.4 数值积分

所有矩阵元内的积分使用高斯求积转换为加权求和（论文式 (2-47)）：

$$
\int_\Omega \bm{f}(\bm{r})\, d\Omega = J_\Omega \sum_{i=1}^{N_G} w_i \bm{f}(\bm{r}_i)
$$

其中 $J_\Omega$ 为网格单元雅可比（面积/体积），$N_G$ 为求积点数。
EMMoMSuite 对三角形远场用 4 点规则、近邻/自项用 7 点规则；四面体默认 5 点；
六面体 8 点。三角形/四面体/六面体的求积点与权重见
`Geometry/GaussQuadrature.jl`。

## 3. 奇异性处理 (Singularity Treatment)

当二重积分的积分区域重合（$R = 0$）或非常接近时，格林函数奇异或剧烈变化，
纯数值求积不准确。论文采用**奇异值提取法**（论文式 (2-48)~(2-57)），
将 $G(R)$ 在 $R = 0$ 附近泰勒展开：

$$
G(R) = \sum_{n=0}^{\infty} \frac{(-{\rm j}k)^n}{4\pi\, n!} R^{n-1}
= \sum_{n=0}^{\infty} C_G^n R^{n-1}
$$

对面积分/体积分中出现的两种基本积分：

$$
\mathcal{I}_{G\Omega} = \int_\Omega G(R)\, d\Omega', \qquad
\overline{\bm{\mathcal{I}}}_{G\Omega}^m = \int_\Omega \frac{\bm{R}}{R^m} G(R)\, d\Omega'
$$

交换求和与积分次序，并定义 $\mathcal{I}_{R\Omega}^n = \int_\Omega R^n d\Omega'$、
$\overline{\bm{\mathcal{I}}}_{R\Omega}^n = \int_\Omega \bm{R} R^n d\Omega'$。
对**面网格**，利用递推公式把面积分化到边界线上的线积分（论文式 (2-52)~(2-55)）：

$$
(n+2)\mathcal{I}_{RS,i}^n = n\, d^2 \mathcal{I}_{RS}^{n-2} + \sum_i P_{0i} \mathcal{I}_{Rl,i}^n
$$

$$
(n+1)\mathcal{I}_{Rl,i}^n = l_i^+ R_{i+}^n - l_i^- R_{i-}^n + n R_{0i}^2 \mathcal{I}_{Rl,i}^{n-2}
$$

$$
\overline{\bm{\mathcal{I}}}_{RS}^n = -\frac{1}{n+2} \sum_i \bm{\hat{u}}_i \mathcal{I}_{Rl,i}^{n+2}
+ d\, \bm{\hat{n}}\, \mathcal{I}_{RS}^n
$$

初值为（论文式 (2-56)）：

$$
\mathcal{I}_{Rl,i}^{-1} = \ln\frac{R_{i+} + l_{i+}}{R_{i-} + l_{i-}}, \qquad
\mathcal{I}_{Rl,i}^{0} = l_i, \qquad
\mathcal{I}_{RS,i}^{0} = A_S
$$

$$
\mathcal{I}_{RS,i}^{-1} = \sum_j \left(P_{0i}\mathcal{I}_{Rl,i}^{-1} - |d_i| \beta_i\right), \qquad
\mathcal{I}_{RS,i}^{-3} = \frac{1}{|d_i|} \sum_j \beta_i
$$

其中 $d$ 为场点到投影点的有向距离，$P_{0i}$ 为投影点到边的有向距离，
$l_i^{\pm}$、$R_{i\pm}$、$R_{0i}$、$\beta_i$ 的几何含义见论文图 2-9。

对**体网格**，利用 $\nabla' R^n = -n\bm{R}R^{n-2}$ 与
$\nabla'\cdot(\bm{R}R^n) = -(n+3)R^n$，把体积分转换到包围体网格的面上的积分
（论文式 (2-57)~(2-58)）：

$$
\mathcal{I}_{RV}^n = -\frac{1}{n+3} \sum_j d_j \mathcal{I}_{RS,j}^n, \qquad
\overline{\bm{\mathcal{I}}}_{RV}^n = -\frac{1}{n+2} \sum_j \hat{\bm{n}}_j \mathcal{I}_{RS,j}^{n+2}
$$

该方案把体积分奇异性转换到面上、面上的奇异性转换到线上，计算线积分后
用递推公式得到所有高阶项。论文采用 14 阶展开（展开阶数在 8 阶以上即可在
0.2 波长内获得高于 $10^{-5}$ 的精度），并在 0.2 波长范围内只处理重合与相邻
网格的奇异性。`Singularities.jl` 中的 `faceSingularityIgIvecg`、
`volumeSingularityIgIvecg`、`singularF1/F21/F22` 即该方案的实现。
