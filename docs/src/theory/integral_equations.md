# 积分方程 (Integral Equations)

本章基于等效原理推导 EMMoMSuite 使用的各类积分方程：PEC 的 S-EFIE / S-MFIE /
CFIE，介质的 PMCHWT，非均匀介质的 V-EFIE / V-MFIE，以及金属-介质混合的
VSIE。公式均来自论文第 2 章。

## 1. 理想导体 (PEC) 表面积分方程

金属在微波频段近似为理想电导体。PEC 表面 $S$ 上的边界条件为
（论文式 (2-25)）：

$$
\left. \left[\bm{E}^i + \bm{E}^s\right] \right|_t = 0, \qquad
\bm{J}_S^m = 0, \qquad
\bm{J}_S = \hat{\bm{n}} \times \left(\bm{H}^i + \bm{H}^s\right)
$$

即切向电场为零、面磁流为零。

### 1.1 面电场积分方程 (S-EFIE)

用面等效源替换金属目标并强制切向电场边界条件，得到 S-EFIE（论文式 (2-29)）：

$$
\left. \bm{E}^i(\bm{r}) \right|_t = -\left. \eta \mathcal{L}\left[\bm{J}_S(\bm{r}')\right] \right|_t, \qquad \bm{r} \in S
$$

其显式形式为第一类 Fredholm 积分方程：

$$
\hat{\bm{n}} \times \bm{E}^i(\bm{r}) = \hat{\bm{n}} \times \left[
{\rm j}\omega\mu \int_S \bm{J}_S G\, dS'
- \frac{1}{{\rm j}\omega\varepsilon} \nabla \int_S \nabla' \cdot \bm{J}_S\, G\, dS'
\right]
$$

**特点**：
- 适用于开域（薄板、线天线）与闭域问题。
- 属于第一类 Fredholm 积分方程，离散矩阵条件数较差。
- 在闭合目标的谐振频率处解不唯一。

### 1.2 面磁场积分方程 (S-MFIE)

强制磁场边界条件得到 S-MFIE（论文式 (2-30)）：

$$
\hat{\bm{n}} \times \bm{H}^i(\bm{r}) = \bm{J}_S + \hat{\bm{n}} \times \mathcal{K}\left[\bm{J}_S(\bm{r}')\right], \qquad \bm{r} \in S
$$

当观测点趋近光滑表面时，$\mathcal{K}$ 算子的主值积分产生 $+\frac{1}{2}\bm{J}_S$ 项：

$$
\frac{1}{2}\bm{J}_S(\bm{r}) - \hat{\bm{n}} \times \int_{PV} \bm{J}_S(\bm{r}') \times \nabla' G(\bm{r}, \bm{r}')\, dS' = \hat{\bm{n}} \times \bm{H}^i(\bm{r})
$$

**特点**：
- 仅适用于闭合曲面。
- 属于第二类 Fredholm 积分方程，条件数较好。
- 在内谐振频率处解不唯一。

### 1.3 混合场积分方程 (CFIE)

为消除内谐振问题并改善收敛，取 EFIE 与 MFIE 的线性组合（论文式 (2-31)）：

$$
\text{CFIE} = \alpha\, \text{EFIE} + (1-\alpha)\, \eta\, \text{MFIE}
$$

其中 $\alpha \in [0, 1]$ 为组合系数（论文采用经验值 0.6，EMMoMSuite 默认 0.5），
$\eta$ 为波阻抗，用于平衡两个方程的数值量级。CFIE 在所有频率下解唯一，
条件数优于 EFIE。

## 2. 介质体积分方程 (VIE) 与 PMCHWT

### 2.1 体电场/磁场积分方程 (V-EFIE / V-MFIE)

对非均匀介质采用体等效源后，空间总场为入射场与散射场之和，得到
V-EFIE 与 V-MFIE（论文式 (2-32)~(2-33)）：

$$
\bm{E}^i(\bm{r}) = \bm{E}(\bm{r}) - \eta \mathcal{L}\left[\bm{J}_V\right] - \mathcal{K}\left[\bm{J}_V^m\right], \qquad \bm{r} \in V
$$

$$
\bm{H}^i(\bm{r}) = \bm{H}(\bm{r}) - \frac{1}{\eta} \mathcal{L}\left[\bm{J}_V^m\right] + \mathcal{K}\left[\bm{J}_V\right], \qquad \bm{r} \in V
$$

当 $\mu = \mu_0$ 时 $\bm{J}_V^m = 0$，只剩 V-EFIE。

### 2.2 金属-介质混合 (VSIE)

金属表面用面等效源、介质区域用体等效源，得到 VS-EFIE / VS-MFIE
（论文式 (2-34)~(2-35)，无磁体流时为 (2-36)~(2-37)）：

$$
\left. \bm{E}^i \right|_t = -\left.\left\{ \eta \mathcal{L}[\bm{J}_S] + \eta \mathcal{L}[\bm{J}_V] + \mathcal{K}[\bm{J}_V^m] \right\}\right|_t, \qquad \bm{r} \in S
$$

$$
\bm{E}^i(\bm{r}) = \bm{E}(\bm{r}) - \eta \mathcal{L}[\bm{J}_S] - \eta \mathcal{L}[\bm{J}_V] - \mathcal{K}[\bm{J}_V^m], \qquad \bm{r} \in V
$$

### 2.3 PMCHWT 方程

对均匀介质体（内外区域分别为 $(\varepsilon_1,\mu_1)$ 与 $(\varepsilon_2,\mu_2)$），
未知量为表面等效电流 $\bm{J}$ 与等效磁流 $\bm{M}$。PMCHWT 利用切向场连续条件，
把内、外区域算子叠加为 $2N \times 2N$ 的分块系统：

$$
\begin{bmatrix}
\mathcal{L}_1 + \mathcal{L}_2 & \mathcal{K}_1 + \mathcal{K}_2 \\
-\mathcal{K}_1 - \mathcal{K}_2 & \frac{1}{\eta_1}\mathcal{L}_1 + \frac{1}{\eta_2}\mathcal{L}_2
\end{bmatrix}
\begin{bmatrix} \bm{J} \\ \bm{M} \end{bmatrix}
=
\begin{bmatrix} -\bm{E}^{i}|_t \\ -\bm{H}^{i}|_t \end{bmatrix}
$$

矩阵显式区分 EJ/EM/HJ/HM 四个块。EMMoMSuite 中 `PMCHWBlockOperators` 即按
该块结构实现；`assemble_K_offdiag` 提供不含质量矩阵主值项的 $\mathcal{K}$ 块
（用于 EM 与 HJ 块）。
