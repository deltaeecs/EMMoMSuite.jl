# 基础电磁理论 (Fundamental Electromagnetics)

本章介绍计算电磁学所需的基础理论，包括麦克斯韦方程组、势函数理论、格林函数以及等效原理。

## 1. 麦克斯韦方程组 (Maxwell's Equations)

在各向同性、线性、均匀的介质区域 ``V`` 中，时谐电磁场（时间因子 ``e^{j\omega t}``）满足麦克斯韦方程组：

```math
\nabla \times \mathbf{E} = -j\omega\mu \mathbf{H} - \mathbf{M}
```
```math
\nabla \times \mathbf{H} = j\omega\epsilon \mathbf{E} + \mathbf{J}
```
```math
\nabla \cdot \mathbf{E} = \frac{\rho_e}{\epsilon}
```
```math
\nabla \cdot \mathbf{H} = \frac{\rho_m}{\mu}
```

其中：
*   ``\mathbf{E}, \mathbf{H}`` 分别为电场强度和磁场强度。
*   ``\mathbf{J}, \mathbf{M}`` 分别为外加电流密度和磁流密度（磁流为数学上的等效引入）。
*   ``\rho_e, \rho_m`` 分别为电荷密度和磁荷密度。
*   ``\epsilon, \mu`` 分别为介质的介电常数和磁导率。

电流与电荷满足连续性方程：
```math
\nabla \cdot \mathbf{J} = -j\omega\rho_e, \quad \nabla \cdot \mathbf{M} = -j\omega\rho_m
```

## 2. 势函数 (Potentials)

为了求解麦克斯韦方程组，通常引入磁矢势 ``\mathbf{A}`` 和电矢势 ``\mathbf{F}``（或标量势 ``\Phi_e, \Phi_m``）。

### 2.1 洛伦兹规范 (Lorenz Gauge)

```math
\nabla \cdot \mathbf{A} = -j\omega\mu\epsilon \Phi_e
```
```math
\nabla \cdot \mathbf{F} = -j\omega\mu\epsilon \Phi_m
```

### 2.2 亥姆霍兹方程 (Helmholtz Equations)

势函数满足非齐次矢量亥姆霍兹方程：
```math
\nabla^2 \mathbf{A} + k^2 \mathbf{A} = -\mu \mathbf{J}
```
```math
\nabla^2 \mathbf{F} + k^2 \mathbf{F} = -\epsilon \mathbf{M}
```
其中 ``k = \omega\sqrt{\mu\epsilon}`` 为波数。

### 2.3 场的势函数表示

```math
\mathbf{E} = -j\omega\mathbf{A} - \nabla\Phi_e - \frac{1}{\epsilon} \nabla \times \mathbf{F}
```
```math
\mathbf{H} = -j\omega\mathbf{F} - \nabla\Phi_m + \frac{1}{\mu} \nabla \times \mathbf{A}
```

利用洛伦兹规范，可消去标量势：
```math
\mathbf{E} = -j\omega\mathbf{A} + \frac{1}{j\omega\mu\epsilon} \nabla(\nabla \cdot \mathbf{A}) - \frac{1}{\epsilon} \nabla \times \mathbf{F}
```
```math
\mathbf{H} = -j\omega\mathbf{F} + \frac{1}{j\omega\mu\epsilon} \nabla(\nabla \cdot \mathbf{F}) + \frac{1}{\mu} \nabla \times \mathbf{A}
```

## 3. 格林函数 (Green's Functions)

### 3.1 标量格林函数 (Scalar Green's Function)

三维自由空间中标量亥姆霍兹方程 ``(\nabla^2 + k^2)G = -\delta(\mathbf{r}-\mathbf{r}')`` 的解为：
```math
G(\mathbf{r}, \mathbf{r}') = \frac{e^{-jkR}}{4\pi R}, \quad R = |\mathbf{r} - \mathbf{r}'|
```

势函数的积分解为：
```math
\mathbf{A}(\mathbf{r}) = \mu \int_V \mathbf{J}(\mathbf{r}') G(\mathbf{r}, \mathbf{r}') dV'
```
```math
\mathbf{F}(\mathbf{r}) = \epsilon \int_V \mathbf{M}(\mathbf{r}') G(\mathbf{r}, \mathbf{r}') dV'
```

### 3.2 并矢格林函数 (Dyadic Green's Functions)

为了直接表示场与源的关系，引入并矢格林函数 ``\overline{\mathbf{G}}_e``。
```math
\mathbf{E}(\mathbf{r}) = -j\omega\mu \int_V \overline{\mathbf{G}}_e(\mathbf{r}, \mathbf{r}') \cdot \mathbf{J}(\mathbf{r}') dV'
```
自由空间电并矢格林函数为：
```math
\overline{\mathbf{G}}_e(\mathbf{r}, \mathbf{r}') = \left( \overline{\mathbf{I}} + \frac{1}{k^2} \nabla\nabla \right) G(\mathbf{r}, \mathbf{r}')
```
其中 ``\overline{\mathbf{I}}`` 是单位并矢。这对应于 ``\mathbf{E} = -j\omega (\mathbf{A} + \frac{1}{k^2}\nabla\nabla\cdot\mathbf{A})`` 的形式。

## 4. 等效原理 (Equivalence Principle)

### 4.1 惠更斯原理 (Huygens' Principle)

场源分布在区域 ``V`` 内，边界为 ``S``。外部区域的场可以通过边界 ``S`` 上的等效电流 ``\mathbf{J}_s`` 和等效磁流 ``\mathbf{M}_s`` 唯一确定。

### 4.2 Love 等效原理 (Surface Equivalence Principle)

在边界 ``S`` 上建立等效面电流和面磁流：
```math
\mathbf{J}_s = \hat{n} \times \mathbf{H}
```
```math
\mathbf{M}_s = \mathbf{E} \times \hat{n}
```
其中 ``\hat{n}`` 为指向外部区域的单位法向量。
这些等效源在外部区域产生与原始源相同的场，在内部区域产生零场（零场定理）。

### 4.3 体等效原理 (Volume Equivalence Principle)

对于非均匀介质体 ``V``，其介电常数为 ``\epsilon(\mathbf{r})``，磁导率为 ``\mu(\mathbf{r})``。我们可以将其替换为背景介质（通常为真空 ``\epsilon_0, \mu_0``），并引入体等效极化电流 ``\mathbf{J}_{eq}`` 和磁流 ``\mathbf{M}_{eq}`` 来模拟介质的存在。

根据麦克斯韦方程：
```math
\nabla \times \mathbf{H} = j\omega\epsilon \mathbf{E} = j\omega\epsilon_0 \mathbf{E} + j\omega(\epsilon - \epsilon_0)\mathbf{E}
```
```math
\nabla \times \mathbf{E} = -j\omega\mu \mathbf{H} = -j\omega\mu_0 \mathbf{H} - j\omega(\mu - \mu_0)\mathbf{H}
```

定义体等效电流：
```math
\mathbf{J}_{eq}(\mathbf{r}) = j\omega(\epsilon(\mathbf{r}) - \epsilon_0) \mathbf{E}(\mathbf{r})
```
定义体等效磁流：
```math
\mathbf{M}_{eq}(\mathbf{r}) = j\omega(\mu(\mathbf{r}) - \mu_0) \mathbf{H}(\mathbf{r})
```

此时，总场 ``\mathbf{E}, \mathbf{H}`` 可以看作是由入射场（源在外部）和体等效源产生的散射场之和：
```math
\mathbf{E} = \mathbf{E}^{inc} + \mathbf{E}^{scat}(\mathbf{J}_{eq}, \mathbf{M}_{eq})
```
```math
\mathbf{H} = \mathbf{H}^{inc} + \mathbf{H}^{scat}(\mathbf{J}_{eq}, \mathbf{M}_{eq})
```
这就是体积分方程 (VIE) 的物理基础。

### 4.4 表面积分方程的基础

利用等效原理，我们可以将散射体替换为背景介质，并在其表面分布等效电流 ``\mathbf{J}_s``（对于 PEC，``\mathbf{M}_s=0``）。通过强制满足边界条件（如切向电场为零），建立积分方程求解 ``\mathbf{J}_s``。
