# 快速算法

传统矩量法在显式装配与矩阵向量乘阶段需要 $O(N^2)$ 的时间与存储。对电大尺寸
问题，EMMoMSuite 通过多层快速多极子算法（MLFMA）把远场相互作用的复杂度
降低到接近 $O(N \log N)$，并采用论文第 4 章提出的 Lebedev 求积与球面非规则
矢量插值进一步降低常数因子。

## 1. 多层快速多极子算法

### 1.1 球面波加法定理与转移函数

MLFMA 的核心是对自由空间标量格林函数做球面波展开（论文式 (2-38)~(2-41)）：

$$
\frac{e^{-{\rm j}k|\bm{R} + \bm{d}|}}{4\pi|\bm{R} + \bm{d}|}
= -\frac{{\rm j}k}{4\pi} \sum_{l=0}^{\infty} (-{\rm j})^{2l} (2l+1)
j_l(kd)\, h_l^{(2)}(kR)\, P_l(\hat{\bm{d}} \cdot \hat{\bm{R}})
$$

其中 $j_l$ 为第一类球贝塞尔函数，$h_l^{(2)}$ 为第二类球汉克尔函数，
$P_l$ 为勒让德多项式。对 $(-{\rm j})^l j_l(kd) P_l$ 做平面波展开：

$$
(-{\rm j})^l j_l(kd) P_l(\hat{\bm{d}} \cdot \hat{\bm{R}})
= \frac{1}{4\pi} \int_\Omega e^{-{\rm j}k\hat{\bm{k}} \cdot \bm{d}}
P_l(\hat{\bm{k}} \cdot \hat{\bm{R}})\, d^2\hat{\bm{k}}
$$

于是带偏置的格林函数可写成单位球面上的积分，其中**转移函数**为
（论文式 (2-41)）：

$$
T_\tau(k, \hat{\bm{k}}, \bm{R}) =
\frac{-{\rm j}k}{(4\pi)^2} \sum_{l=0}^{\tau} (-{\rm j})^l (2l+1)
h_l^{(2)}(kR)\, P_l(\hat{\bm{k}} \cdot \hat{\bm{R}})
$$

当源点、场点各自在盒子中心附近变化时转移函数保持不变，据此可加速矩阵向量乘。
球面积分近似为加权求和（论文式 (2-43)）：

$$
\int_\Omega f(\hat{\bm{k}})\, d^2\hat{\bm{k}} \approx \sum_{p=1}^{N_p} W_p f(\hat{\bm{k}}_p), \qquad
\sum_{p=1}^{N_p} W_p = 4\pi
$$

![快速多极子算法示意图：r_a 与 r_b 分别为包围场、源两点的盒子中心（论文图 2-11）](figures/fma.png)

### 1.2 远场矩阵向量乘积

以面 EFIE 为例，远场部分（论文式 (2-45)~(2-49)）：

$$
b_m = {\rm j}k\eta \int_\Omega
\bm{\mathcal{R}}_S(\bm{f}_m^S, \hat{\bm{k}}) \cdot
T_\tau(k, \hat{\bm{k}}, \bm{R}_{ba})
\sum_{n=1}^{N} x_n \bm{\mathcal{F}}_S(\bm{f}_n^S, \hat{\bm{k}})\, d^2\hat{\bm{k}}
$$

其中辐射/配置积分算子为：

$$
\bm{\mathcal{R}}_S(\bm{f}_m^S, \hat{\bm{k}}) = \int_{S_m}
\left(\overline{\bm{I}} - \hat{\bm{k}}\hat{\bm{k}}\right) \cdot \bm{f}_m^S(\bm{r})\,
e^{-{\rm j}k\hat{\bm{k}}\cdot(\bm{r} - \bm{r}_a)}\, dS
$$

$$
\bm{\mathcal{F}}_S(\bm{f}_n^S, \hat{\bm{k}}) = \int_{S_n}
\left(\overline{\bm{I}} - \hat{\bm{k}}\hat{\bm{k}}\right) \cdot \bm{f}_n^S(\bm{r}')\,
e^{-{\rm j}k\hat{\bm{k}}\cdot(\bm{r}_b - \bm{r}')}\, dS'
$$

在球坐标下 $\overline{\bm{I}} - \hat{\bm{k}}\hat{\bm{k}} = \mathrm{diag}(0,1,1)$，
即只保留 $\theta$、$\phi$ 切向分量。多层形式把 $\bm{r} - \bm{r}'$ 逐层拆分，
聚合（upward）、转移（translation）、配置（downward）三阶段独立复用。

![多层快速多极子算法示意图：虚线围成的小盒子为子层盒子，实线大方框为其父层盒子（论文图 2-13）](figures/mlfma.png)

### 1.3 截断项经验公式

转移函数截断项数按盒子尺寸选取（论文式 (2-42)）：

$$
\tau(l) \approx 1.73\, k a_l + 2.16\, d_0^{2/3} (k a_l)^{1/3}, \qquad d_0 = 3
$$

其中 $a_l$ 为第 $l$ 层盒子边长，$d_0$ 为精度位数。代码中
`truncation_kernel(rel_l)` 传入 $rel_l = a_l/\lambda$，等价地写为
$2\pi\sqrt{3}\, rel_l + 2.16\, d_0^{2/3} (2\pi\, rel_l)^{1/3}$
（注意 $1.73 \cdot 2\pi \approx 2\pi\sqrt{3}$；实现中的精度参数取 9，得到的
截断更保守）。

### 1.4 层间插值与反插值

各层转移函数阶数不同，所需球面采样率也不同，因此聚合/配置时需要在层间
匹配采样点（论文式 (4-2)~(4-3)）：

$$
\bm{\mathcal{F}}(\hat{\bm{k}}^{l-1}) = \bm{\Gamma}^{l-1,l} \bm{\mathcal{F}}(\hat{\bm{k}}^{l})
$$

其中 $\bm{\Gamma}^{l-1,l}$ 为第 $l$ 层到父层第 $l-1$ 层的插值矩阵。
**反插值不是下采样**，而是由"同层两种采样率下转移结果相等"推导得到的
伴随关系：转移并乘上求积权重后，左乘插值矩阵的转置即可完成向子层的配置
（论文式 (4-8)~(4-13)）：

$$
\left(\bm{\Gamma}^{l-1,l}\right)^T \bm{\mathcal{T}}(\hat{\bm{k}}^{l-1})
= \bm{\mathcal{T}}(\hat{\bm{k}}^{l})
$$

### 1.5 算法流程

1. 八叉树划分：叶层盒子边长一般取 $0.20\lambda \sim 0.25\lambda$，盒子边长
   $a_0 = 2^i a_f$ 逐层倍增（`OctreeBuilder.build_octree`）。

   ![八叉树的二维邻盒子、次邻盒子分布示意图（论文图 2-14）](figures/octree2D.png)

2. 计算近场矩阵（CSR/CSC，叶层 + 邻盒子），远场封装为算子。
3. 远场算子：聚合（叶层辐射积分 + 逐层插值/相移）、转移（同层远亲盒子）、
   配置（反插值/相移 + 叶层测试）。
4. 完整矩阵向量乘积 $Zx = Z_{near} x + Z_{far}(x)$。

## 2. Lebedev 求积与球面非规则矢量插值

### 2.1 为什么用 Lebedev

球面高斯求积把球面展开为 $\theta \in [0,\pi]$、$\phi \in [0,2\pi]$ 的矩形区域，
$\theta$ 方向取 $\tau+1$ 个 Gauss-Legendre 点、$\phi$ 方向均匀取
$2(\tau+1)$ 个点，共 $2(\tau+1)^2$ 个点（求积效率仅 $2/3$，两极附近冗余）。
Lebedev 求积把单位球面按直角坐标划分为 8 个旋转对称区域并求解非线性方程组，
得到更均匀的点分布，求积效率接近 1，点数约为球面高斯的 $2/3$。

![典型八叉树层截断项下的球面高斯求积采样点分布（论文图 4-3）](figures/gq_points_sets.png)

![典型八叉树层截断项下的 Lebedev 球面采样点分布（论文图 4-4）](figures/lq_points_sets.png)

### 2.2 球面高斯求积与拉格朗日局部插值

球面高斯求积的采样点落在 $\theta$、$\phi$ 的规则格点上，两层之间可以方便地
做拉格朗日局部插值：以待插值点附近 $N_{ip} \times N_{ip}$ 个采样点为支撑构造
插值多项式（论文式 (4-4)~(4-7)），插值矩阵与函数值无关、可在预处理阶段计算，
且保持稀疏。

![截断项 tau 分别为 5 和 7 的两层间球面高斯求积点分布（论文图 4-6）](figures/laginterp_k2f.png)

单步插值每行有 $N_{ip}^2$ 个非零元，完成一次插值的计算量为
$N_p^{l-1} N_{ip}^2$；分步插值先在 $\theta$ 方向一维插值得到中间态、再沿
$\phi$ 方向插值，计算量降为 $\left(N_p^{l-1} + N_{p,\theta}^{l-1}\right) N_{ip}$，
且中间态可被多个待插值点复用。反插值使用插值矩阵的转置（分步插值则两步转置
后调换顺序，论文 4.1.4 节）。

![利用 4x4 插值窗口的拉格朗日插值策略：(a) 单步插值；(b) 分步插值（论文图 4-7）](figures/laginterp_1and2step.png)

### 2.3 球面非规则矢量插值

Lebedev 点不在规则网格上，且 $\theta$、$\phi$ 展开极不均匀，不能按标量分量
分别插值。论文规定：极点处（$\theta = 0, \pi$）切向矢量取直角坐标 $x, y$
（$\theta=\pi$ 取 $x, -y$）分量，其余位置取 $\theta$、$\phi$ 分量；把两分量
拼接为列向量后，层间插值写成 $2N_p^{l-1} \times 2N_p^l$ 的分块形式
（论文式 (4-20)）：

$$
\begin{bmatrix}
\bm{\mathcal{F}}_\theta(\hat{\bm{k}}^{l-1}) \\
\bm{\mathcal{F}}_\phi(\hat{\bm{k}}^{l-1})
\end{bmatrix}
=
\begin{bmatrix}
\bm{\Gamma}_{\theta\theta}^{l-1,l} & \bm{\Gamma}_{\theta\phi}^{l-1,l} \\
\bm{\Gamma}_{\phi\theta}^{l-1,l} & \bm{\Gamma}_{\phi\phi}^{l-1,l}
\end{bmatrix}
\begin{bmatrix}
\bm{\mathcal{F}}_\theta(\hat{\bm{k}}^{l}) \\
\bm{\mathcal{F}}_\phi(\hat{\bm{k}}^{l})
\end{bmatrix}
$$

其中交叉块 $\bm{\Gamma}_{\theta\phi}$、$\bm{\Gamma}_{\phi\theta}$ 连接了两个
分量的信息，使插值成为矢量插值而非分量独立插值。反插值仍为插值矩阵的转置。

### 2.4 插值矩阵计算

**辐射函数数据集**。各类基函数的辐射积分可归纳为通用形式（论文式 (4-21)）：

$$
\bm{\mathcal{F}}(\hat{\bm{k}}) = C_f \left(\overline{\bm{I}} - \hat{\bm{k}}\hat{\bm{k}}\right)
\cdot \bm{\hat{\rho}}\, e^{-{\rm j}k\hat{\bm{k}} \cdot (\bm{r}_b - \bm{r}')}
$$

随机生成 $\bm{\hat{\rho}}$ 与尺度在盒子大小的 $\bm{r}_b - \bm{r}'$，在同一辐射源
的两层采样点上批量计算（实部、虚部分别作为数据列，数据集翻倍），得到
$\mathbb{F}(\hat{\bm{k}}^{l-1})$ 与 $\mathbb{F}(\hat{\bm{k}}^{l})$。

**矩阵初始化（稀疏模式）**。对每个父层待插值点，选距离最近的 $N_k$ 个子层
采样点作为插值点（4 个子块同步标记），其余位置为零，保证矩阵高度稀疏。
采样点按直角坐标 $x$、$y$、$z$ 值排序（满足 $\hat{\bm{k}}_i = -\hat{\bm{k}}_{N_p^l+1-i}$），
不同阶数的 Lebedev 点之间共享 14 个固定点（6 个轴点 $\pm\hat{x}, \pm\hat{y}, \pm\hat{z}$
与 8 个立方体角点 $(\pm\sqrt{3}/3, \pm\sqrt{3}/3, \pm\sqrt{3}/3)$），
这些点不需要插值，对应行只有 1 个非零元"1"。

![Lebedev 求积层间采样点分布与层间插值矩阵示意图：(a) 截断项 5 与 7 的采样点分布；(b) 不同插值点数下的稀疏模式，图中更稀疏的行对应 14 个层间共享点（论文图 4-11）](figures/lebedev_k2f_pattern.png)

**伪逆法（逐行）**。完整求解 $\bm{\Gamma}^{l-1,l} = \mathbb{F}(\hat{\bm{k}}^{l-1})\mathbb{F}^{\dagger}(\hat{\bm{k}}^{l})$
得到的是稠密矩阵。论文改为逐行计算
（论文式 (4-24)~(4-26)）：对第 $p$ 行，提取非零元 $\bm{\gamma}_p$ 与列索引
集合 $C_p$（$2N_k$ 个，两个子矩阵的行），求解

$$
\bm{\gamma}_p = \mathbb{F}_p(\hat{\bm{k}}^{l-1})\,
\mathbb{F}_{C_p}^{\dagger}(\hat{\bm{k}}^{l})
$$

行满秩条件由 $N_d > 2N_{pl}$ 降为 $N_d > 2N_k$，数据集规模只需匹配插值点数。
EMMoMSuite 中 `pinv2interpW.jl` 即按此实现（先反距离权重初始化稀疏模式，
再逐行 `pinv` 求解）；同时提供 `SHInterp.jl` 的球谐精确/局部/混合解析权重
作为确定性的替代路径。

**神经网络训练**（论文式 (4-27)~(4-29)）。把插值问题化为损失函数
$\mathrm{Loss}(\bm{\Gamma}) = \|\bm{\Gamma}\bm{\mathcal{F}}(\hat{\bm{k}}^{l}) - \bm{\mathcal{F}}(\hat{\bm{k}}^{l-1})\|^{\mathcal{P}}$
的最小化，用梯度下降
（论文采用 Nesterov 加速法）更新稀疏插值矩阵。训练收敛结果与伪逆法基本一致，
实际计算 CPU 资源充足时推荐伪逆法。

### 2.5 精度与效率

- 插值误差定义（论文式 (4-30)）：
  $\epsilon_i = |\tilde{\bm{\mathcal{F}}} - \bm{\mathcal{F}}| / \max|\bm{\mathcal{F}}|$。
- 拉格朗日插值 $4\times4$、$6\times6$ 窗口的平均误差约为 $6\times10^{-4}$、
  $6\times10^{-5}$；Lebedev 矢量插值用 8 个、18 个插值点达到同等精度。

![Lebedev 矢量插值在典型八叉树层下的平均插值误差 epsilon_i 随插值点数 N_k 的变化（论文图 4-13）](figures/lebedev_interp_error.png)

![Lebedev 矢量插值几个典型算例中的远场矩阵向量乘积相对误差分布（论文图 4-14）](figures/lebedev_mv_error.png)

- 相同精度下单点计算量与拉格朗日单步插值一致；Lebedev 采样点约为球面高斯的
  $2/3$，总体时间、空间复杂度降低约 $1/3$。论文算例中叶层插值点数取 9、
  其余层取 8。

### 2.6 高阶层的实现策略

Lebedev 数据集最高支持 131 阶多项式（论文给出 131 阶的数据集）。当某层的
多项式阶数 $p = 2\tau+1$ 超过数据集上限时，`LebedevSortedPoints.high_order_nodes`
采用 Fibonacci 准均匀格点（等权重 $4\pi/n$）作为替代，保证任意高阶都能构造
求积点；`LVI.jl` 的 `levelIntegralInfoCal` 会对此发出警告。

## 3. 低频 MLFMA 与 ACA

标准平面波形式的 MLFMA 在低频下会出现 low-frequency breakdown（展开形式随
盒尺寸与波数缩小而病态）。常见解决路线包括低频稳定的多极子展开形式、
归一化平面波展开或低频重标定。

自适应交叉近似（ACA）是核无关的代数压缩方法，对远场块 $\bm{Z}_{block}$ 构造
$\bm{Z}_{block} \approx \bm{U}\bm{V}^{T}$（$\bm{U} \in \mathbb{C}^{M\times r}$，
$\bm{V} \in \mathbb{C}^{N\times r}$，$r \ll M, N$），不依赖核函数的解析加法定理，
对复杂介质核更灵活，但常数因子与稳定性需要单独评估。

## 4. ACA 实现（`FastAlgorithms.ACA`）

EMMoMSuite 已实现部分主元 ACA（Gibson《Method of Moments》Ch9 Algorithm 6）：

- `aca(getrow, getcol, m, n; tol, maxrank, recompress)`：按行/列采样构造
  $\bm{Z}_{block} \approx \bm{U}\bm{V}^{T}$（**转置约定，无共轭**，适配复数对称
  阻抗矩阵）；收敛判据为 $\|\bm{u}_k\|\|\bm{v}_k\| \le \mathrm{tol}\,\|\tilde{\bm{Z}}\|_F$，
  $\|\tilde{\bm{Z}}\|_F^2$ 用递推估计；含零行/零块早期终止。
- `recompress!`：QR/SVD 再压缩（$\tau_{SVD} \approx 10\,\tau_{ACA}$），典型再压缩
  20–30%。
- `BlockEvaluator` / `eval_block`：按全局基函数索引求值 `Z[rows, cols]`，支持
  EFIE/MFIE/CFIE 的 RWG 块，供 ACA 按行/列稀疏采样，避免装配整块矩阵。
- `ACAOperator <: AbstractIntegralOperator`：复用 MLFMA 八叉树聚类与近场稀疏
  装配；叶层非邻盒子对按 ACA 压缩为低秩块，实现 `mul!` 与 GMRES 直接兼容。
  对称矩阵用转置语义应用下三角块（`V*(Uᵀ*x)`），避免重复压缩。

## 5. MLACA 实现（`FastAlgorithms.ACA.MLACAOperator`）

H-矩阵风格多层递归块压缩（Gibson Ch10 思想，基于八叉树多层结构）：

- 可容许对（非邻盒子）→ ACA 压缩整个子树块（可跨越多个叶层盒子）；
- 近邻对 → 下钻子盒子对；叶层近邻对由近场稀疏矩阵覆盖；
- 对角块 → 递归到子层；叶层自/邻对在近场中。

该结构保证每个叶层基函数对恰好被近场或某个低秩块覆盖一次，无重叠、无漏算。
当前支持复数对称算子（EFIE）；非对称问题使用 `ACAOperator(symmetric=false)`。

## 6. 参数配置与实测结果（2026-08-11）

- `nInterp`（层间插值点数）与 `precision_digits`（截断公式精度参数 $d_0$，默认
  9.0）不再硬编码，可通过 `MLFMAOperator`/`ACAOperator`/`MLACAOperator` 构造参数
  配置，默认值保持现状。
- 实测（EFIE 球体，本机 Julia 1.12，1 线程）：

| 用例 | 方法 | MatVec 相对误差 | 解相对误差 | 压缩率 | 说明 |
|------|------|-----------------|------------|--------|------|
| N=792, 叶层 0.25λ | MLFMA | 3.4e-4 | 1.9e-3 | 75.2%* | 近场稀疏 + 远场隐式 |
| N=792, 叶层 0.25λ | ACA | 5.6e-5 | 1.9e-4 | 40.3% | 1252 个低秩块 |
| N=252, 叶层 0.125λ | ACA | 2.8e-6 | 2.1e-5 | 17.7% | 13304 块（单层） |
| N=252, 叶层 0.125λ | MLACA | 1.6e-5 | 3.5e-5 | 22.3% | 3652 块（多层，块数少 73%） |
| N=252, 叶层 0.125λ | MLFMA | 2.1 | 1.3 | — | 小叶层失效（breakdown） |

`*` MLFMA 的"压缩率"仅计近场稀疏存储（远场不显式存储），与 ACA 的显式低秩存储
口径不同，不可直接比较。MLACA 在四层八叉树上以更少块数实现更高压缩率，且
N=252/0.125λ 用例显示 MLFMA 在该尺度的 breakdown——这正是文档 §3 所述
低频/小叶层失效，ACA/MLACA 作为核无关路径形成互补。

预条件实测（N=252 ACA 算子，GMRES，`abstol=1e-6`）：Identity 63 迭代，
ILU(0.01) 20 迭代，SPAI 19 迭代，BlockJacobi 400 迭代未收敛（叶层块过小）。
ILU/SPAI 可直接用 `ILUPreconditioner(op)` / `SPAIPreconditioner(op)` 便捷构造。

## 7. 非对称 MLACA 与 PMCHW 多基函数支持

- **非对称 MLACA**：`MLACAOperator(symmetric=false)` 对每个盒子对的两个方向
  分别压缩（适合 CFIE 等非对称算子）；对称时仍用转置语义避免重复压缩。
- **PMCHW 多基函数块求值**：`PMCHWBlockEvaluator` 按 2N 系统全局索引求值
  远场块（J/M 双通道），四个子块 EJ/HM 用 L 算子、EM/HJ 用 ±K^PMCHW；
  `ACAOperator(pmchw, basis, ...)` 与 `MLACAOperator(pmchw, basis, ...)`
  直接支持 PMCHW（2N×2N，非对称双向压缩）。注意每个 L 算子必须独立缓冲
  （`efie_interaction!` 末尾会整体乘 factor，同一缓冲连续调用会重复缩放）。
- PMCHW 该实现本机 cond≈4.3e6（N=150 球体），解误差受条件数放大；算子级
  门控以 MatVec 误差与 GMRES 算子残差为准。

## 8. 直接块 LU（多 RHS）

`block_lu(op)`（`FastAlgorithms.BlockLUModule`）对 ACA/MLACA 算子做叶层分块
直接 LU（Gibson Ch9 Algorithm 7）：
- 对角块 `A_bb = Z_bb − Σ_{p<b} L_bp U_pb`，用**无主元 LU**
  （`lu(A, NoPivot())`，块公式要求 `A_bb = L_bb U_bb` 精确成立）；
- 离对角块 `L_sb = (Z_sb − Σ L_sp U_pb) U_bb⁻¹`、`U_bs = L_bb⁻¹ (Z_bs − Σ L_bp U_ps)`，
  可用 ACA 再压缩控制存储；
- `block_lu_solve(F, B)` 前代 + 回代支持多 RHS；`F \ b` 支持单 RHS。

## 9. 更大规模实测（2026-08-11，本地 1 线程，稠密参照）

| 用例 | 方法 | N | MatVec 误差 | 解误差 | 压缩率 | 说明 |
|------|------|---|------------|--------|--------|------|
| EFIE 球 | MLFMA | 1734 | 6.2e-5 | 5.9e-4 | 75.0%* | 76 迭代 |
| EFIE 球 | ACA | 1734 | 2.8e-5 | 3.3e-4 | 57.0% | 求解 0.19s |
| EFIE 球 | MLFMA | 2280 | 4.1e-5 | 7.1e-4 | 74.8%* | 500 迭代 82.2s |
| EFIE 球 | ACA/MLACA | 2280 | 2.5e-5 | 5.0e-4 | 60.7% | 500 迭代 1.8s |
| CFIE 球 | ACA/MLACA | 792 | 1.9e-5 | 2.0e-5 | -3.6% | 7 迭代（cond≈22） |
| PMCHW 球 | ACA/MLACA | 600 | 1.2e-4 | 1.2 | 4.8% | cond≈4.3e6 限制 |
| 低频 EFIE | ACA/MLACA | 792 | 5.1e-6 | 2.9e-3 | 42.5% | 30 MHz，0.03λ 叶层 |
| EFIE 球 | ACA/MLACA | 11352 | 5.1e-6 | 7.2 | 71.5% | GMRES+ILU 500 迭代未收敛（EFIE 稠密网格预条件挑战） |
| CFIE 球 | ACA/MLACA | 11352 | 2.2e-6 | 2.6e-5 | 67.8% | 100 迭代收敛（cond≈1.2e3） |

`*` MLFMA 压缩率仅计近场稀疏存储（口径不同）。全部用例 `finite/NaN` 检查通过。
N=2280 时 ACA/MLACA 单次 MatVec 远快于 MLFMA（同 500 迭代下求解 1.8s vs 82.2s），
且压缩率随 N 增大而提升；小 N（792）CFIE 非对称双向压缩开销超过收益（负压缩率），
属正常现象。N=11352（约 1.1 万未知量）本地实测：EFIE 压缩率 71.5%、MatVec 误差
5.1e-6，但 ILU 预条件 GMRES 500 迭代不收敛（残差停滞 ~2.5），属 EFIE 稠密网格
预条件挑战（算子本身精确）；良态 CFIE 同规模 100 迭代收敛、解误差 2.6e-5、
压缩率 67.8%。该大规模基准脚本
（`benchmark/run_large_fast_solvers_benchmark.jl`）仅供本地手动运行，未接入 CI。

## 10. H 矩阵 H-LU 与 H2 扩展

### 10.1 实现

- `FastAlgorithms.ACA.HMatrixModule`：`HMatrixNode`（`:dense`/`:lowrank`/`:split`）
  与 `hmatrix_from_mlaca`——从 MLACA/ACA 算子的分层低秩结构重建显式 H 矩阵树
  （对称算子反方向用转置因子 `V*(Uᵀ)`；PMCHW 2N 自动展开 J/M 双通道）。
- `h_lu!`：标准分层块 LU（对角先更新 `A_bb = Z_bb − Σ L_bp U_pb` 再递归分解；
  离对角 `L_sb = (Z_sb − Σ L_sp U_pb) U_bb⁻¹`、`U_bs = L_bb⁻¹ (Z_bs − Σ L_bp U_ps)`；
  递归右除 U / 左除 L）。
- `h_lu_solve`：前代 + 回代多 RHS 直接求解（全局索引输入/输出）。
- 默认 `recompress=false` 精确分解；`recompress=true` 用 ACA 截断离对角因子块
  （误差校验 ≤10·tol，否则回退稠密）。

### 10.2 实测结果（2026-08-11，本地 1 线程，`benchmark/benchmark_hl_lu.jl`）

| 用例 | N | 方法 | 因子化 | 求解 | 残差 | 压缩率 |
|------|---|------|--------|------|------|--------|
| EFIE（MLACA 多层树） | 150 | H-LU 精确 | 0.37s | 0.008s | 9.4e-16 | 0% |
| EFIE | 792 | BlockLU | 2.06s | 0.008s | 1.5e-15 | 0% |
| EFIE | 792 | H-LU 精确 | 0.50s | 0.014s | 1.7e-15 | 0% |
| EFIE | 792 | H-LU 再压缩 | 0.76s | 0.014s | 7.8e-4 | **12.7%**（944 低秩块） |
| CFIE（非对称） | 150 | H-LU 精确 | 0.14s | 0.003s | 4.5e-16 | 0% |
| PMCHW（2N） | 300 | H-LU 精确 | 0.21s | 0.005s | 1.1e-15 | 0% |
| PMCHW | 300 | H-LU 再压缩 | 0.23s | 0.005s | 0.50 | 6.2%（cond≈4e6 不稳定） |
| 低频 EFIE | 150 | H-LU 精确 | 0.06s | 0.003s | 4.0e-16 | 0% |

结论：H-LU 精确分解比叶层 BlockLU 快 4-14 倍、残差 ~1e-15；再压缩在良态系统上
以 ≤10·tol 受控误差换取存储（N=792 12.7%），病态系统（PMCHW）不稳定，默认关闭。

### 10.3 H2 与 H2-LU（设计，下一阶段）

- H2 矩阵：嵌套基（cluster bases U_t/V_t + 父子转移矩阵），MatVec 用嵌套基加速；
- **H2-LU 为 opt-in 实验特性**：按 H2Mat4Ham 经验，仅当残差/解误差/秩/内存对照
  H-LU 无回归时才启用；每个实验保留稠密参照门控。
