# 积分方程 (Integral Equations)

本章详细推导用于求解电磁散射与辐射问题的各类积分方程。

## 1. 理想导体 (PEC) 表面积分方程

考虑位于自由空间背景中的理想导体，表面为 ``S``。入射波为 ``(\mathbf{E}^{inc}, \mathbf{H}^{inc})``。

### 1.1 电场积分方程 (EFIE)

EFIE 基于 PEC 表面的切向电场边界条件：
```math
\hat{n} \times (\mathbf{E}^{inc} + \mathbf{E}^{scat}) = 0, \quad \mathbf{r} \in S
```
散射场由表面感应电流 ``\mathbf{J}`` 产生：
散射场由表面感应电流 ``\mathbf{J}`` 产生：
```math
\mathbf{E}^{scat}(\mathbf{J}) = -j\omega\mu \int_S \mathbf{J}(\mathbf{r}') G(\mathbf{r}, \mathbf{r}') dS' + \frac{1}{j\omega\epsilon} \nabla \int_S \nabla' \cdot \mathbf{J}(\mathbf{r}') G(\mathbf{r}, \mathbf{r}') dS'
```
代入边界条件得到 EFIE：
```math
\hat{n} \times \mathbf{E}^{inc}(\mathbf{r}) = \hat{n} \times \left[ j\omega\mu \int_S \mathbf{J} G dS' - \frac{1}{j\omega\epsilon} \nabla \int_S \nabla' \cdot \mathbf{J} G dS' \right]
```
**特点**：
*   适用于开域（如薄板）和闭域问题。
*   属于第一类 Fredholm 积分方程，条件数较差。
*   在内谐振频率处解不唯一。

### 1.2 磁场积分方程 (MFIE)

MFIE 基于 PEC 表面的切向磁场边界条件：
```math
\hat{n} \times \mathbf{H}^{tot} = \mathbf{J} \implies \mathbf{J} = \hat{n} \times (\mathbf{H}^{inc} + \mathbf{H}^{scat})
```
散射磁场为：
```math
\mathbf{H}^{scat}(\mathbf{J}) = \int_S \mathbf{J}(\mathbf{r}') \times \nabla' G(\mathbf{r}, \mathbf{r}') dS'
```
当观测点 ``\mathbf{r}`` 趋近于光滑表面时，积分主值产生奇异项 ``\frac{1}{2}\mathbf{J}(\mathbf{r})``。
```math
\frac{1}{2}\mathbf{J}(\mathbf{r}) - \hat{n} \times \int_{PV} \mathbf{J}(\mathbf{r}') \times \nabla' G(\mathbf{r}, \mathbf{r}') dS' = \hat{n} \times \mathbf{H}^{inc}(\mathbf{r})
```
**特点**：
*   仅适用于闭合曲面。
*   属于第二类 Fredholm 积分方程，条件数较好。
*   在内谐振频率处解不唯一。

### 1.3 混合场积分方程 (CFIE)

为了消除内谐振问题，采用 EFIE 和 MFIE 的线性组合：
```math
	ext{CFIE} = \alpha \cdot \text{EFIE} + (1-\alpha) \cdot \eta \cdot \text{MFIE}
```
其中 ``\eta = \sqrt{\mu/\epsilon}`` 为波阻抗，``\alpha`` 为加权系数（通常取 0.2 ~ 0.8）。
CFIE 在所有频率下均有唯一解，且条件数优于 EFIE。

## 2. 介质体积分方程 (Dielectric Surface Integral Equations)

考虑均匀介质体，外部区域 ``E_1, H_1 (\epsilon_1, \mu_1)``，内部区域 ``E_2, H_2 (\epsilon_2, \mu_2)``。边界为 ``S``。
未知量为表面的等效电流 ``\mathbf{J}`` 和等效磁流 ``\mathbf{M}``。

### 2.1 PMCHWT 方程

Poggio-Miller-Chang-Harrington-Wu-Tsai (PMCHWT) 方程是最常用的介质表面积分方程。
利用切向场连续性：
```math
\hat{n} \times (\mathbf{E}_1^{inc} + \mathbf{E}_1^{scat}) = \hat{n} \times \mathbf{E}_2^{tot}
```
```math
\hat{n} \times (\mathbf{H}_1^{inc} + \mathbf{H}_1^{scat}) = \hat{n} \times \mathbf{H}_2^{tot}
```
在 PMCHWT 形式中，我们将内外部的积分算子组合：
```math
\left[ \hat{n} \times \mathbf{E}_1^{scat}(\mathbf{J}, \mathbf{M}) + \hat{n} \times \mathbf{E}_2^{scat}(-\mathbf{J}, -\mathbf{M}) \right]_{tan} = -\hat{n} \times \mathbf{E}^{inc}
```
```math
\left[ \hat{n} \times \mathbf{H}_1^{scat}(\mathbf{J}, \mathbf{M}) + \hat{n} \times \mathbf{H}_2^{scat}(-\mathbf{J}, -\mathbf{M}) \right]_{tan} = -\hat{n} \times \mathbf{H}^{inc}
```
这构成了一个 ``2N \times 2N`` 的矩阵方程系统。

## 3. 体积分方程 (Volume Integral Equation, VIE)

对于非均匀介质，需采用体积分方程。
定义体等效极化电流：
```math
\mathbf{J}_{vol} = j\omega (\epsilon(\mathbf{r}) - \epsilon_0) \mathbf{E}(\mathbf{r})
```
总电场为入射场与体电流产生的散射场之和：
```math
\mathbf{E}(\mathbf{r}) = \mathbf{E}^{inc}(\mathbf{r}) + \int_V \overline{\mathbf{G}}_e(\mathbf{r}, \mathbf{r}') \cdot \mathbf{J}_{vol}(\mathbf{r}') dV'
```
将 ``\mathbf{J}_{vol}`` 代入，得到关于总场 ``\mathbf{E}`` 的积分方程：
```math
\mathbf{E}(\mathbf{r}) - \int_V \overline{\mathbf{G}}_e(\mathbf{r}, \mathbf{r}') \cdot [j\omega (\epsilon(\mathbf{r}') - \epsilon_0) \mathbf{E}(\mathbf{r}')] dV' = \mathbf{E}^{inc}(\mathbf{r})
```
通常使用 SWG 基函数对电通量密度 ``\mathbf{D}`` 进行离散求解。
