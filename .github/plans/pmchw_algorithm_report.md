# PMCHW (Poggio-Miller-Chang-Harrington-Wu-Tsai) 算法完整理论报告

> **目的**：精确提炼 PMCHW 表面积分方程的完整理论和 EMSuite 实现流程，覆盖从物理推导到 MLFMA 加速的全链路，精度足以指导 Phase 15 `PMCHWMLFMAOperator` 开发。  
> **参考**：EMSuite 代码 + CEM 知识库 (Gibson MoM Ch.3 §3.6, Ch.11 Algorithm 14)  
> **时间约定**：$e^{-j\omega t}$ ⟹ Green 函数 $G = \frac{e^{-jkR}}{4\pi R}$, Hankel 函数用 $h_l^{(2)}$  
> **注意**：Legacy (`MoM_Kernels`) 不包含 PMCHW 实现，PMCHW 完全在 EMSuite 中开发。

---

## 总体数据流

```
Nastran .nas 文件
  ↓ MeshProcess
[node(3,N_node), triangles(3,N_tri)]
  ↓ RWG 基函数构建
[TriangleInfo 数组, RWG 基函数数组 (N个)]
  ↓ PMCHW(freq, εᵣ, μᵣ) 算子构造
[k₀, η₀, k₁, η₁]
  ↓ 矩阵装配 / MLFMA 构造
[Z(2N×2N) 或 PMCHWMLFMAOperator]
  ↓ 激励向量 V(2N)
[V_E(1:N), V_H(N+1:2N)]
  ↓ 求解 Z·I = V
[I_2N = [I_J; I_M]]
  ↓ 双流 RCS
[σ(θ,ϕ)]
```

---

## 0. 物理基础：面等效原理与 PMCHW 推导

### 0.1 面等效原理 (Surface Equivalence Principle)

对于均匀介质体（外部区域 $R_0$，内部区域 $R_1$），在介质界面 $S$ 上引入等效面电流 $\mathbf{J}$ 和等效面磁流 $\mathbf{M}$：

$$\mathbf{J}(\mathbf{r}) = \hat{n} \times [\mathbf{H}_0(\mathbf{r}) - \mathbf{H}_1(\mathbf{r})]$$
$$\mathbf{M}(\mathbf{r}) = -\hat{n} \times [\mathbf{E}_0(\mathbf{r}) - \mathbf{E}_1(\mathbf{r})]$$

其中 $\hat{n}$ 为界面法向（从 $R_1$ 指向 $R_0$）。

利用消光定理 (Extinction Theorem)，$R_0$ 和 $R_1$ 的场可分别用各自区域的 Green 函数表示。关键的线性相关性：

$$\mathbf{J}_0(\mathbf{r}) = -\mathbf{J}_1(\mathbf{r}), \quad \mathbf{M}_0(\mathbf{r}) = -\mathbf{M}_1(\mathbf{r})$$

即两侧的等效电流方向相反，非独立未知量。

### 0.2 L 算子和 K 算子定义

**L 算子** (EFIE 型)：
$$(\mathbf{L}\mathbf{X})(\mathbf{r}) = \left(1 + \frac{1}{k^2}\nabla\nabla\cdot\right) \int_S G(\mathbf{r}, \mathbf{r}') \mathbf{X}(\mathbf{r}') d\mathbf{r}'$$

**K 算子** (MFIE 型)：
$$(\mathbf{K}\mathbf{X})(\mathbf{r}) = \nabla \times \int_S G(\mathbf{r}, \mathbf{r}') \mathbf{X}(\mathbf{r}') d\mathbf{r}'$$

用算子表示区域 $R_l$ 中的总场 (Gibson Eq.3.75-3.76)：

$$\mathbf{E}(\mathbf{r}) = -j\omega\mu_l (\mathbf{L}\mathbf{J})(\mathbf{r}) - (\mathbf{K}\mathbf{M})(\mathbf{r})$$
$$\mathbf{H}(\mathbf{r}) = -j\omega\varepsilon_l (\mathbf{L}\mathbf{M})(\mathbf{r}) + (\mathbf{K}\mathbf{J})(\mathbf{r})$$

### 0.3 边界积分方程

在界面 $S$ 上应用边界条件，得到 EFIE 和 MFIE（Gibson §3.6.2, Eq.3.178-3.179）：

**EFIE (电场边界)**：
$$\frac{1}{2}\hat{n} \times \mathbf{M}_l(\mathbf{r}) = -j\omega\mu_l (\mathbf{L}\mathbf{J}_l)(\mathbf{r}) - (\mathbf{K}\mathbf{M}_l)(\mathbf{r}) - \hat{n} \times \mathbf{E}^i_l(\mathbf{r})$$

**MFIE (磁场边界)**：
$$\frac{1}{2}\hat{n} \times \mathbf{J}_l(\mathbf{r}) = -j\omega\varepsilon_l (\mathbf{L}\mathbf{M}_l)(\mathbf{r}) + (\mathbf{K}\mathbf{J}_l)(\mathbf{r}) - \hat{n} \times \mathbf{H}^i_l(\mathbf{r})$$

### 0.4 PMCHW 方程的形成

PMCHW 方法将介质界面两侧的 EFIE/MFIE 相加（Gibson §3.6.3.1）：

1. 利用 $\mathbf{J}_0 = -\mathbf{J}_1$, $\mathbf{M}_0 = -\mathbf{M}_1$ 合并列（消除线性相关未知量）
2. 将两侧同类方程的行相加（消除超定性，得到方阵）

相加后的 E 方程（外部 EFIE + 内部 EFIE）：
$$[j\omega\mu_0 (\mathbf{L}_0\mathbf{J}) + j\omega\mu_1 (\mathbf{L}_1\mathbf{J})] + [(\mathbf{K}_0\mathbf{M}) + (\mathbf{K}_1\mathbf{M})] = \mathbf{E}^i$$

相加后的 H 方程（外部 MFIE + 内部 MFIE）：
$$-[(\mathbf{K}_0\mathbf{J}) + (\mathbf{K}_1\mathbf{J})] + [j\omega\varepsilon_0 (\mathbf{L}_0\mathbf{M}) + j\omega\varepsilon_1 (\mathbf{L}_1\mathbf{M})] = \mathbf{H}^i$$

> **关键**：PMCHW 的 EFIE/MFIE 相加消去了 $\frac{1}{2}$ 对角项——两侧的 $\pm\frac{1}{2}$ 对角项符号相反，相加后为零。因此 PMCHW 的 K 算子是**纯主值积分 (PV, Principal Value)**，无对角项 ($\frac{1}{2}$ mass matrix)。

### 0.5 PMCHW 与 MFIE/CFIE 的 K 算子区别

| 特性 | MFIE K | PMCHW K |
|------|--------|---------|
| 对角项 $\frac{1}{2}\hat{n}\times\mathbf{X}$ | ✅ 有，贡献 mass matrix | ❌ 无，相加消去 |
| 测试函数 | $\hat{n} \times \mathbf{f}_m$ (旋转后的) | $\mathbf{f}_m$ (直接) |
| 积分核 | $(\hat{n} \times \mathbf{f}_m) \cdot (\nabla G \times \mathbf{f}_n)$ | $\mathbf{f}_m \cdot (\nabla G \times \mathbf{f}_n)$ |
| η₀ 预乘 | ✅ CFIE 中 $\eta_0$ 预乘 | ❌ 无 |

---

## 1. 材料参数与波数 (`PMCHW` 结构体)

### 1.1 构造函数

```julia
PMCHW(freq::FT, eps_r, mu_r = 1.0)
```

### 1.2 参数计算

| 参数 | 公式 | 说明 |
|------|------|------|
| $k_0$ | $2\pi f / c_0$ | 自由空间波数（实数） |
| $\eta_0$ | $\sqrt{\mu_0/\varepsilon_0} \approx 376.73\;\Omega$ | 自由空间波阻抗（实数） |
| $k_1$ | $k_0 \sqrt{\varepsilon_r \mu_r}$ | 介质波数（复数，含损耗） |
| $\eta_1$ | $\eta_0 \sqrt{\mu_r / \varepsilon_r}$ | 介质波阻抗（复数） |

### 1.3 物理常数

| 常数 | 值 |
|------|------|
| $c_0$ | $299\,792\,458.0$ m/s |
| $\mu_0$ | $4\pi \times 10^{-7}$ H/m |
| $\varepsilon_0$ | $8.854\ldots \times 10^{-12}$ F/m |
| $\eta_0$ | $376.730\ldots\;\Omega$ |

### 1.4 有损介质处理

对于 Green 函数核中的波数，使用 $k_1$ 的实部近似（忽略指数衰减）：

```julia
k1_r = FT(_k_real(k1_c))       # Re(k₁)，用于 exp(-jkR)/R 计算
eta1_r = FT(abs(real(eta1_c)))  # |Re(η₁)|，仅用于 EFIE 记录字段
```

**精度说明**：对无损介质（$\varepsilon_r$ 为实数），结果精确。对有损介质，Green 函数核使用实数 $k$ 是近似。但前系数 (factor) 仍使用完整复数 $k_1, \eta_1$。

---

## 2. 2N×2N 阻抗矩阵结构

### 2.1 矩阵分块定义

$$Z = \begin{bmatrix} Z^{EJ} & Z^{EM} \\ Z^{HJ} & Z^{HM} \end{bmatrix} \in \mathbb{C}^{2N \times 2N}$$

$$Z \begin{pmatrix} \mathbf{I}_J \\ \mathbf{I}_M \end{pmatrix} = \begin{pmatrix} \mathbf{V}_E \\ \mathbf{V}_H \end{pmatrix}$$

### 2.2 各子块公式

#### Z^EJ [1:N, 1:N] — L 算子 (EFIE 型)

$$Z^{EJ}_{mn} = \sum_{l \in \{0,1\}} \frac{jk_l \eta_l}{16\pi} \cdot l_m \cdot l_n \int\!\!\!\int \left[(\boldsymbol{\rho}_m \cdot \boldsymbol{\rho}_n) - \frac{4}{k_l^2}\right] \frac{e^{-jk_l R}}{R}\; dA_t\, dA_s$$

**代码因子**：`factor = im * k_c * eta_c / (16π)`, 即 $jk\eta / (16\pi)$

**物理含义**：$Z^{EJ} = L(k_0, \eta_0) + L(k_1, \eta_1)$，等同于两个 EFIE 的叠加。矢量势项 $\propto j\omega\mu$，标量势项 $\propto j/(\omega\varepsilon)$。

**积分核展开**：
$$\int\!\!\!\int [\boldsymbol{\rho}_m \cdot \boldsymbol{\rho}_n] \cdot G \; dA_t dA_s \quad \text{(矢量势项)}$$
$$-\frac{4}{k^2} \int\!\!\!\int G \; dA_t dA_s \quad \text{(标量势项, } \nabla\cdot\mathbf{f}_m = \pm l_m/A\text{)}$$

#### Z^HM [N+1:2N, N+1:2N] — L 算子 (对偶 EFIE 型, η⁻¹)

$$Z^{HM}_{mn} = \sum_{l \in \{0,1\}} \frac{jk_l}{\eta_l \cdot 16\pi} \cdot l_m \cdot l_n \int\!\!\!\int \left[(\boldsymbol{\rho}_m \cdot \boldsymbol{\rho}_n) - \frac{4}{k_l^2}\right] \frac{e^{-jk_l R}}{R}\; dA_t\, dA_s$$

**代码因子**：`factor = im * k_c / (eta_c * 16π)`, 即 $jk / (\eta \cdot 16\pi)$

**物理含义**：$Z^{HM} = L_e(k_0, \eta_0) + L_e(k_1, \eta_1)$，"反转 η" 的 EFIE。矢量势项 $\propto j\omega\varepsilon$，标量势项 $\propto j/(\omega\mu)$。

> **Z^EJ 与 Z^HM 的对偶关系**：$Z^{EJ}$ 使用 $jk\eta$ (= $j\omega\mu$), $Z^{HM}$ 使用 $jk/\eta$ (= $j\omega\varepsilon$)。两者结构完全相同，仅前系数中 $\eta \leftrightarrow 1/\eta$。

#### Z^EM [1:N, N+1:2N] — K 算子 (PMCHW 专用)

$$Z^{EM}_{mn} = \sum_{l \in \{0,1\}} \frac{1}{16\pi} \cdot l_m \cdot l_n \int\!\!\!\int \mathbf{f}_m(\mathbf{r}) \cdot [\nabla G_l(\mathbf{r}, \mathbf{r}') \times \mathbf{f}_n(\mathbf{r}')] \; dA_t \, dA_s$$

**积分核推导**：

Green 函数梯度：
$$\nabla G = -\left(jk + \frac{1}{R}\right) \frac{e^{-jkR}}{R} \hat{R} = -\text{temp} \cdot \mathbf{R}_\text{vec}$$

其中 $\mathbf{R}_\text{vec} = \mathbf{r} - \mathbf{r}'$, $R = |\mathbf{R}_\text{vec}|$, 

$$\text{temp} = \left(jk + \frac{1}{R}\right) \frac{e^{-jkR}}{R^2} \cdot w_i \cdot w_j$$

代入：
$$\mathbf{f}_m \cdot (\nabla G \times \mathbf{f}_n) = -\text{temp} \cdot \boldsymbol{\rho}_m \cdot (\mathbf{R}_\text{vec} \times \boldsymbol{\rho}_n)$$

**最终 kernel**：`Z_local[m,n] += -dot(rho_m, cross(rvec, rho_n)) * temp`

**缩放因子**：`Z_local[m,n] *= lm * ln * η/(16π)`

> **注意**：`calc_k_pmchw_term!` 中使用 `eta_div_16pi = mfie.eta / (16π)`。当 `mfie.eta = 1`（通过 `MFIE(0, k, 1, gq)` 构造），实际因子为 $1/(16\pi)$。

#### Z^HJ [N+1:2N, 1:N] — 结构不变量

$$Z^{HJ} = -Z^{EM} \quad \text{（精确成立，代码中直接取反）}$$

**代码**：`Z[N+1:2N, 1:N] .= -K0`

### 2.3 矩阵结构总结表

| 子块 | 位置 | 算子 | 前系数 | 物理 | 对角 |
|------|------|------|--------|------|------|
| $Z^{EJ}$ | [1:N, 1:N] | $L_0 + L_1$ | $jk\eta/(16\pi)$ | $j\omega\mu$ | 有 (奇异项) |
| $Z^{EM}$ | [1:N, N+1:2N] | $K_0^{PV} + K_1^{PV}$ | $1/(16\pi)$ | — | 无 (PV) |
| $Z^{HJ}$ | [N+1:2N, 1:N] | $-K_0^{PV} - K_1^{PV}$ | $-1/(16\pi)$ | — | 无 (PV) |
| $Z^{HM}$ | [N+1:2N, N+1:2N] | $L_{e,0} + L_{e,1}$ | $jk/(\eta \cdot 16\pi)$ | $j\omega\varepsilon$ | 有 (奇异项) |

### 2.4 对比 EFIE 单块矩阵

| 特性 | EFIE | PMCHW |
|------|------|-------|
| 矩阵大小 | $N \times N$ | $2N \times 2N$ |
| 未知数 | J (N个) | J + M (2N个) |
| 方程 | 1个 (EFIE) | 2个 (E方程 + H方程) |
| L 算子 | 1个区域 ($k_0, \eta_0$) | 2个区域叠加 |
| K 算子 | 无 | 有 (PV, 无对角) |
| 矩阵对称性 | 对称 | $Z^{HJ} = -Z^{EM}$ (反对称对) |

---

## 3. 矩阵装配算法

### 3.1 Z^EJ 和 Z^HM 的装配

Z^EJ 和 Z^HM 复用 EFIE 的 `assemble_impedance_matrix`，仅修改前系数 (factor)：

```julia
# Z^EJ: factor = jk₀η₀/(16π) 和 jk₁η₁/(16π)
efie_ej0 = _l_block_operator(k0, eta0, k0_c, eta0_c, :EJ)
efie_ej1 = _l_block_operator(k1_r, eta1_r, k1_c, eta1_c, :EJ)

# Z^HM: factor = jk₀/(η₀·16π) 和 jk₁/(η₁·16π)  
efie_hm0 = _l_block_operator(k0, eta0, k0_c, eta0_c, :HM)
efie_hm1 = _l_block_operator(k1_r, eta1_r, k1_c, eta1_c, :HM)
```

`_l_block_operator` 中的 `mode` 参数决定 factor 计算：
- `:EJ` → `factor = jk · η / (16π)` → $j\omega\mu$ 比例
- `:HM` → `factor = jk / (η · 16π)` → $j\omega\varepsilon$ 比例

### 3.2 Z^EM 的装配 (`assemble_K_pmchw_offdiag`)

```
1. 构造 MFIE 壳体对象: MFIE(0, k, 1, gq)  [eta=1 → 因子 1/(16π)]
2. 预计算每个三角形的高斯点: quad_points[t] = 4点高斯求积
3. 三角形对循环:
   ├─ 跳过自身 (t_test == t_src): K 算子自对角为零 (PMCHW 无质量矩阵)
   └─ 非对角对: calc_k_pmchw_term!(Z_local, ...)
4. 全局装配: assemble_generic(..., symmetric=false)
```

### 3.3 K^PMCHW 积分核的计算细节

对于每个三角形对 (test × source)，3×3 局部子矩阵：

```
for j = 1:n_pts (源求积点)
    rho_n[1:3] = r_src[j] - v_src[1:3]
    for i = 1:n_pts (测试求积点)
        rvec = r_test[i] - r_src[j]
        R = |rvec|
        if R < 1e-12: skip (奇异点)
        
        temp = (jk + 1/R) × exp(-jkR)/R² × w_i × w_j
        rho_m[1:3] = r_test[i] - v_test[1:3]
        
        for n in 1:3, m in 1:3:
            Z_local[m,n] += -dot(rho_m[m], cross(rvec, rho_n[n])) × temp
```

缩放：`Z_local[m,n] *= lm × ln × η/(16π)` → 当 η=1 时为 `lm × ln / (16π)`

### 3.4 奇异积分处理

**L 算子** (Z^EJ, Z^HM)：

自对角元素 ($m = n$, 同三角形) 使用**解析奇异积分**：

$$G_\text{int} = (G - 1/R) + 1/R$$

- $(G - 1/R)$ 部分：光滑，数值积分
- $1/R$ 部分：解析公式
  - $F_1(a,b,c)$：$\iint 1/R \; dS dS'$  
  - $F_{21}$, $F_{22}$：$\iint \rho \cdot \rho' / R \; dS dS'$

**K 算子** (Z^EM, Z^HJ)：

自对角 $(t_\text{test} = t_\text{src})$：直接跳过。物理原因——PMCHW 的 K 算子的自对角 mass matrix 项 $\pm\frac{1}{2}$ 在两侧相加时精确抵消。

---

## 4. 激励向量

### 4.1 平面波激励

$$\mathbf{V} = \begin{pmatrix} \mathbf{V}_E \\ \mathbf{V}_H \end{pmatrix} \in \mathbb{C}^{2N}$$

#### V_E [1:N]：EFIE 标准激励

$$V_{E,m} = \int_S \mathbf{f}_m(\mathbf{r}) \cdot \mathbf{E}^{inc}(\mathbf{r}) \; dS$$

其中 $\mathbf{E}^{inc}(\mathbf{r}) = \mathbf{E}_0 \exp(jk_0 \hat{k} \cdot \mathbf{r})$（平面波）。

复用 EFIE 的激励计算：`excitation_vector(EFIE(freq), source, basis)`。

#### V_H [N+1:2N]：PMCHW 磁场激励

$$V_{H,m} = \int_S \mathbf{f}_m(\mathbf{r}) \cdot \mathbf{H}^{inc}(\mathbf{r}) \; dS$$

其中 $\mathbf{H}^{inc}(\mathbf{r}) = \frac{\hat{k} \times \mathbf{E}^{inc}(\mathbf{r})}{\eta_0}$。

> **与 MFIE 激励的关键区别**：
> - MFIE：$\eta_0 \int \mathbf{f} \cdot (\hat{n} \times \mathbf{H}^{inc}) \; dS$ （含面法向叉积 + η₀ 因子）
> - PMCHW V_H：$\int \mathbf{f} \cdot \mathbf{H}^{inc} \; dS$ （直接点积，无 $\hat{n}$×，无 η₀ 因子）

### 4.2 Delta-Gap 激励

$$V_E[\text{idx}] = V_\text{gap} \cdot l_\text{idx}, \qquad V_H = \mathbf{0}$$

缝隙电压仅激励 E 方程（前 N 行），H 方程的激励为零。M 电流完全由 K 算子耦合提供。

---

## 5. 求解

$$Z \cdot \begin{pmatrix} \mathbf{I}_J \\ \mathbf{I}_M \end{pmatrix} = \begin{pmatrix} \mathbf{V}_E \\ \mathbf{V}_H \end{pmatrix}$$

- **直接法**：LU 分解，适用于小系统
- **迭代法**：GMRES + MLFMA（4遍远场加速），适用于大系统
- **预条件子**：可用 ILU 或块 Jacobi（需覆盖 2N×2N 结构）

---

## 6. 双流 RCS 后处理

### 6.1 辐射积分

$$\mathbf{N}(\theta, \phi) = \int_S \mathbf{J}(\mathbf{r}') e^{jk_0 \hat{r} \cdot \mathbf{r}'} \; dS' \quad \text{（电流辐射积分）}$$
$$\mathbf{L}(\theta, \phi) = \int_S \mathbf{M}(\mathbf{r}') e^{jk_0 \hat{r} \cdot \mathbf{r}'} \; dS' \quad \text{（磁流辐射积分）}$$

其中 $\mathbf{J} = \sum_n I_{J,n} \mathbf{f}_n$, $\mathbf{M} = \sum_n I_{M,n} \mathbf{f}_n$。

实现：两次调用 `radiation_integral_rwg`，分别传入 $I_J$ 和 $I_M$。3点高斯求积。

### 6.2 散射远场合成

代码中的远场幅度模式（省略全局因子 $C$）：

$$E_\theta \equiv \eta_0 N_\theta + L_\phi$$
$$E_\phi \equiv \eta_0 N_\phi - L_\theta$$

完整远场：$\mathbf{E}^{\text{far}} = \frac{-jk_0 e^{-jk_0 r}}{4\pi r} (-E_\theta \hat{\theta} - E_\phi \hat{\phi})$，全局负号在 RCS 取模平方时消去。

### 6.3 RCS 计算

$$\sigma(\theta, \phi) = \frac{k_0^2}{4\pi} (|E_\theta|^2 + |E_\phi|^2) \quad \text{[m², 线性标度]}$$

$$\text{RCS}_\text{dB} = 10 \log_{10}(\sigma) \quad \text{[dBsm]}$$

> **与 EFIE RCS 的区别**：
> - EFIE：$\sigma = \frac{(k_0 \eta_0)^2}{4\pi} |\mathbf{N}|^2$（仅 J 贡献）
> - PMCHW：$\sigma = \frac{k_0^2}{4\pi} (|\eta_0 N_\theta + L_\phi|^2 + |\eta_0 N_\phi - L_\theta|^2)$（J + M 交叉耦合）

### 6.4 远场合成公式推导

从 Maxwell 方程出发，远区散射场（正比于 $e^{-jk_0 r}/r$）：

$$\mathbf{E}_\text{scat}^{\text{far}}(\mathbf{r}) = \frac{-jk_0 e^{-jk_0 r}}{4\pi r} \left[\eta_0 (\hat{r} \times (\hat{r} \times \mathbf{N})) + (\hat{r} \times \mathbf{L})\right]$$

球坐标正交关系：
- $\hat{r} \times \hat{\theta} = \hat{\phi}$
- $\hat{r} \times \hat{\phi} = -\hat{\theta}$

**第一项** $\hat{r} \times (\hat{r} \times \mathbf{N})$：

$$\hat{r} \times (\hat{r} \times \mathbf{N}) = \hat{r}(\hat{r}\cdot\mathbf{N}) - \mathbf{N} = -\mathbf{N}_\perp$$

远场中 $\mathbf{N}$ 无径向分量，故 $\hat{r} \times (\hat{r} \times \mathbf{N}) = -\mathbf{N}$，即 $\theta$ 分量 $= -N_\theta$, $\phi$ 分量 $= -N_\phi$。

**第二项** $\hat{r} \times \mathbf{L}$：

$$\hat{r} \times \mathbf{L} = \hat{r} \times (L_\theta \hat{\theta} + L_\phi \hat{\phi}) = L_\theta \hat{\phi} - L_\phi \hat{\theta}$$

即 $\theta$ 分量 $= -L_\phi$, $\phi$ 分量 $= +L_\theta$。

**合并** (乘以 $C = \frac{-jk_0 e^{-jk_0 r}}{4\pi r}$)：

$$E_\theta^{\text{full}} = C \cdot [\eta_0 \cdot (-N_\theta) + (-L_\phi)] = -C \cdot [\eta_0 N_\theta + L_\phi]$$
$$E_\phi^{\text{full}} = C \cdot [\eta_0 \cdot (-N_\phi) + L_\theta] = -C \cdot [\eta_0 N_\phi - L_\theta]$$

代码中定义 $E_\theta \equiv \eta_0 N_\theta + L_\phi$, $E_\phi \equiv \eta_0 N_\phi - L_\theta$（省略全局负号和 $C$ 因子），用于 RCS 的幅度计算。因为 $|E^{\text{full}}|^2 = |C|^2 |E^{\text{code}}|^2$，全局负号在取模平方时消失。

---

## 7. MLFMA 加速 (4遍远场算法)

### 7.1 总体架构

PMCHW 的 MLFMA 加速需要处理 4 个子块的远场相互作用，分解为 **4 遍** MLFMA：

| 遍 | 源电流 | 八叉树 | 写入目标 |
|----|--------|--------|----------|
| 1 | $\mathbf{I}_J$ (x[1:N]) | octree0 (k₀) | $y[1:N] += Z^{EJ}_\text{far}$, $y[N+1:2N] += Z^{HJ}_\text{far}$ |
| 2 | $\mathbf{I}_J$ (x[1:N]) | octree1 (k₁) | 同上（累加） |
| 3 | $\mathbf{I}_M$ (x[N+1:2N]) | octree0 (k₀) | $y[1:N] += Z^{EM}_\text{far}$, $y[N+1:2N] += Z^{HM}_\text{far}$ |
| 4 | $\mathbf{I}_M$ (x[N+1:2N]) | octree1 (k₁) | 同上（累加） |

每遍的流程：聚合(叶层) → 上推(插值) → 转移(Translation) → 下推(反插值) → 解聚(叶层)。

### 7.2 双八叉树

由于内/外区域波数不同（$k_0 \neq k_1$），必须使用**两套独立八叉树**：

| 八叉树 | 波数 | 波长 | 影响 |
|--------|------|------|------|
| octree0 | $k_0$ | $\lambda_0 = 2\pi/k_0$ | 叶盒边长、截断项数、极点采样 |
| octree1 | $k_1$ | $\lambda_1 = 2\pi/\text{Re}(k_1)$ | 同上 |

> **不能共享八叉树**：波数影响 (1) 叶层盒子边长 (以波长为参考), (2) 各层截断项数 $L$, (3) Lebedev/Gauss-Legendre 极点采样数, (4) 相移因子, (5) 转移函数。

### 7.3 叶层聚合 (`aggregate_leaf_pmchw!`)

PMCHW 的叶层聚合函数与 EFIE 的**完全一致**（辐射函数只依赖电流和基函数几何，与算子类型无关）：

$$\text{aggS}[p, \text{pol}, i_\text{cube}] \mathrel{+}= \sum_{n \in \text{cube}} I_n \sum_{gi} (\hat{e}_\text{pol} \cdot \boldsymbol{\rho}_{n,gi}) \cdot w_{gi} \cdot \frac{l_n}{2} \cdot s_n \cdot e^{+jk\hat{r}_p \cdot \mathbf{r}'_\text{local}}$$

其中：
- $\mathbf{r}'_\text{local} = \mathbf{r}_{gi} - \mathbf{c}_\text{cube}$ (相对于盒子中心)
- $s_n$ = RWG 符号 (±1)
- $\text{pol} \in \{1, 2\}$ 对应 $\hat{\theta}, \hat{\phi}$

遍 1-2 传入 `x[1:N]` (J 系数), 遍 3-4 传入 `x[N+1:2N]` (M 系数)。

> **与 EFIE aggSBF 的差异**：无差异。聚合函数不区分 J 或 M，只需把正确的系数向量传入。

### 7.4 上推、转移、下推

这三步与标准 EFIE MLFMA **完全相同**：

1. **上推** (`aggregate_upward!`)：从叶层到根层，对每层做 (a) 反插值（子层极点 → 父层极点）+ (b) 相移（子盒子中心 → 父盒子中心）
2. **转移** (`translate!`)：对每层的远亲盒子对，乘以转移函数 $\alpha_\text{Trans}[p] = \frac{-jk}{16\pi^2} W_p T_L(k, \hat{k}_p, \mathbf{D})$
3. **下推** (`disaggregate_downward!`)：从根层到叶层，对每层做 (a) 相移（父盒子中心 → 子盒子中心）+ (b) 插值（父层极点 → 子层极点）

> **关键**：遍 1-2 使用 octree0 (k₀) 或 octree1 (k₁)，其 Translation 函数用的波数不同，但算法结构相同。

> **叶层 Translation 跳过**：当前实现跳过叶层 (最细层) 的 Translation（`levelID == nLevels && continue`），叶层远亲对由 4-box 近场直接积分覆盖。这与 EFIE MLFMA 的策略一致。

### 7.5 叶层解聚——J 通道 (`disaggregate_leaf_pmchw_j!`)

当源电流为 J 时，远场 $\text{disaggG}$ 同时贡献给 $Z^{EJ}$（L 算子）和 $Z^{HJ}$（−K 算子）两行：

```
对每个叶层盒子内的基函数 bfID:
    (te_L, te_K) = _receive_terms(bf, basis, field, r0, k, η, poles)
    
    y[bfID]     += te_L × factor_EJ    # E-行: Z^EJ (L 算子)
    y[bfID + N] += te_K × factor_HJ    # H-行: Z^HJ (−K 算子)
```

因子：
- `factor_EJ = jkη / (4π)` — L 算子 MLFMA 因子 (直接法 $jk\eta/(16\pi)$ 的 4 倍)
- `factor_HJ = jk / (4π)` — −K 算子 MLFMA 因子 (**正号！** 因为 $Z^{HJ} = -K$，K 核带负号，两个负号抵消)

### 7.6 叶层解聚——M 通道 (`disaggregate_leaf_pmchw_m!`)

当源电流为 M 时，远场同时贡献给 $Z^{EM}$（K 算子）和 $Z^{HM}$（L 算子）两行：

```
对每个叶层盒子内的基函数 bfID:
    (te_L, te_K) = _receive_terms(bf, basis, field, r0, k, η, poles)
    
    y[bfID]     += te_K × factor_EM    # E-行: Z^EM (K 算子)
    y[bfID + N] += te_L × factor_HM    # H-行: Z^HM (L 算子)
```

因子：
- `factor_EM = -jk / (4π)` — K 算子 MLFMA 因子 (**负号！** K 核自身带负号)
- `factor_HM = jk / (η · 4π)` — L^M 算子 MLFMA 因子 (直接法 $jk/(\eta \cdot 16\pi)$ 的 4 倍)

### 7.7 接收函数 `_receive_terms`

同时计算 L 型和 K 型两种接收模式：

```
te_L = 0    # L-型接收：rho · E_inc
te_K = 0    # K-型接收：rho · (r̂ × E_inc)

for 每个三角形支撑, 每个求积点 r:
    r_local = r - cube_center
    w_f = sign × edge_length/2 × gq_weight
    
    for 每个极点 p:
        (Eθ, Eϕ) = disaggG[p, 1:2]
        phase = exp(-jk · r̂_p · r_local)
        
        # L-型：标准 EFIE 接收
        E_inc = (Eθ · θ̂_p + Eϕ · φ̂_p) × phase
        te_L += dot(rho, E_inc) × w_f
        
        # K-型：叉积接收（用于 K 算子）
        # r̂ × (Eθ θ̂ + Eϕ φ̂) = Eθ (r̂×θ̂) + Eϕ (r̂×φ̂) = Eθ φ̂ - Eϕ θ̂
        r̂×E = (Eθ · φ̂_p - Eϕ · θ̂_p) × phase
        te_K += dot(rho, r̂×E) × w_f

return (te_L, te_K)
```

### 7.8 MLFMA 系数链验证

#### Z^EJ 远场系数链

$$Z^{EJ,\text{far}}_{mn} = \underbrace{\text{te\_L}}_{\text{L-接收}} \cdot \underbrace{\frac{jk\eta}{4\pi}}_{\text{factor\_EJ}}$$

其中 te_L 包含：
- 解聚场 (`disaggG`)
- 转移因子 `αTrans = -jk/(16π²) · W_p · T_L`
- 聚合场 (`aggS`)

展开：

$$Z^{EJ,\text{far}} = \frac{l_m}{2} \cdot \frac{l_n}{2} \cdot \frac{-jk}{16\pi^2} \cdot W_p \cdot T_L \cdot \frac{jk\eta}{4\pi}$$

$$= \frac{l_m l_n}{4} \cdot \frac{-jk}{16\pi^2} \cdot jk\eta \cdot \frac{1}{4\pi} \cdot [(4\pi)^2 \cdot G_\text{int}]$$

利用 Addition Theorem $\sum W_p T_L e^{\pm jk\cdots} = \frac{(4\pi)^2}{-jk} G$：

$$= \frac{l_m l_n}{4} \cdot jk\eta \cdot G_\text{int} / (4\pi) = \frac{jk\eta}{16\pi} \cdot l_m l_n \cdot G_\text{int}$$

$$= Z^{EJ,\text{direct}} \quad \checkmark$$

#### Z^HM 远场系数链

同理，将 $\eta \to 1/\eta$：

$$Z^{HM,\text{far}} = \frac{jk}{\eta \cdot 16\pi} \cdot l_m l_n \cdot G_\text{int} = Z^{HM,\text{direct}} \quad \checkmark$$

#### Z^EM 远场系数链 (K 算子)

$$Z^{EM,\text{far}}_{mn} = \underbrace{\text{te\_K}}_{\text{K-接收}} \cdot \underbrace{\frac{-jk}{4\pi}}_{\text{factor\_EM}}$$

K 算子的直接法：`Z_local[m,n] *= lm * ln / (16π)` 并且核中含 $-(jk + 1/R) e^{-jkR}/R^2$ 梯度因子。

在远场 ($R \gg 1/k$)，Green 函数梯度的主导项：$\nabla G \approx -jk \hat{R} G$。

因此 K 算子的远场形式可通过 Addition Theorem 的梯度形式获得。远场系数链：

$$Z^{EM,\text{far}} = \frac{l_m}{2} \cdot \frac{l_n}{2} \cdot \frac{-jk}{16\pi^2} \cdot W_p \cdot T_L \cdot \frac{-jk}{4\pi}$$

利用 $\sum W_p T_L e^{\pm} = \frac{(4\pi)^2}{-jk} G$，括号中 $\frac{-jk}{16\pi^2} \cdot \frac{(4\pi)^2}{-jk} = 1$，再乘 $\frac{-jk}{4\pi}$：

$$= \frac{l_m l_n}{4} \cdot \frac{-jk}{4\pi} \cdot G = \frac{-jk}{16\pi} \cdot l_m l_n \cdot G$$

K 算子直接法远场极限：$\frac{l_m l_n}{16\pi} \cdot (-jk) G \cdot (\hat{R} \times \boldsymbol{\rho}_n) \cdot \boldsymbol{\rho}_m$

系数匹配 $\checkmark$

#### Z^HJ 系数链

$Z^{HJ} = -K$，factor_HJ = $+jk/(4\pi)$（正号），te_K 包含 `rho · (r̂ × E)`：

$$Z^{HJ,\text{far}} = \text{te\_K} \cdot \frac{jk}{4\pi} = -Z^{EM,\text{far}} \quad \checkmark$$

### 7.9 MLFMA 因子总结

| 子块 | 直接法因子 | MLFMA 因子 | 接收类型 | 比值 |
|------|-----------|-----------|---------|------|
| $Z^{EJ}$ | $jk\eta/(16\pi)$ | $jk\eta/(4\pi)$ | L-型 | ×4 |
| $Z^{HM}$ | $jk/(\eta \cdot 16\pi)$ | $jk/(\eta \cdot 4\pi)$ | L-型 | ×4 |
| $Z^{EM}$ | $1/(16\pi)$ | $-jk/(4\pi)$ | K-型 | — |
| $Z^{HJ}$ | $-1/(16\pi)$ | $+jk/(4\pi)$ | K-型 | — |

> **×4 因子来源**：聚合使用 `edge_length/2` 归一化，解聚也使用 `edge_length/2` 归一化，产生 $l_m l_n / 4$ ≠ 直接法的 $l_m l_n$。因此 MLFMA 因子 = 直接法因子 × 4。但 K 算子的 MLFMA 因子形式不同，因为 K 核是 $\nabla G$ 而非 $G$，MLFMA 转移函数的梯度由 `r̂ ×` 操作隐式实现。

---

## 8. 近场矩阵 (`assemble_near_field_pmchw`)

### 8.1 策略

当前实现采用**全矩阵 + 稀疏化**策略（非最高效，但正确）：

1. 装配完整 2N×2N 密集矩阵 `Z_full`
2. 根据 octree0 的叶层盒子邻近关系，确定近场对
3. 提取近场对的矩阵元素，构造 CSC 稀疏矩阵

### 8.2 近邻搜索范围

使用 **4-box** 搜索（±4 个 cube 宽度），覆盖 $9^3 = 729$ 个候选偏移方向。

> **与 EFIE 近场的区别**：
> - EFIE 近场：仅 octree0 的 N×N 邻近对
> - PMCHW 近场：octree0 的 N×N 邻近对 ×4（EJ, EM, HJ, HM 四个子块均需包含）

### 8.3 近场对注册

```
for each (i, j) in near basis function pairs:
    push!(near_pairs, (i,     j    ))   # Z^EJ 块
    push!(near_pairs, (i,     j + N))   # Z^EM 块
    push!(near_pairs, (i + N, j    ))   # Z^HJ 块
    push!(near_pairs, (i + N, j + N))   # Z^HM 块
```

### 8.4 待优化项

- **近场 octree 选择**：当前仅用 octree0 的邻近关系。理论上 octree1 可能有不同的邻近对（因波长不同）。如 $k_1 > k_0$（即 $\varepsilon_r > 1$），octree1 的叶盒更小，远亲关系可能不同。目前这是一个已知的近似。
- **内存效率**：先装配全矩阵再稀疏化浪费内存。可改为直接按邻近对计算。

---

## 9. mul! 完整流程

```julia
function mul!(y, A::PMCHWMLFMAOperator, x)
    fill!(y, 0)
    N = num_basis(A.basis)
    
    # ① 近场：2N×2N 稀疏矩阵乘
    mul!(y, A.Z_near, x)
    
    # ② 远场：4遍 MLFMA
    y_far = zeros(2N)
    
    # 遍1: J×k₀
    clear_agg!(octree0)
    aggregate_leaf_pmchw!(octree0, basis, x, sorted_ids0, 1:N, k0)
    up_translate_down!(octree0)
    disaggregate_leaf_pmchw_j!(octree0, ..., :k0)
    
    # 遍2: J×k₁
    clear_agg!(octree1)
    aggregate_leaf_pmchw!(octree1, basis, x, sorted_ids1, 1:N, k1)
    up_translate_down!(octree1)
    disaggregate_leaf_pmchw_j!(octree1, ..., :k1)
    
    # 遍3: M×k₀
    clear_agg!(octree0)
    aggregate_leaf_pmchw!(octree0, basis, x, sorted_ids0, (N+1):(2N), k0)
    up_translate_down!(octree0)
    disaggregate_leaf_pmchw_m!(octree0, ..., :k0)
    
    # 遍4: M×k₁
    clear_agg!(octree1)
    aggregate_leaf_pmchw!(octree1, basis, x, sorted_ids1, (N+1):(2N), k1)
    up_translate_down!(octree1)
    disaggregate_leaf_pmchw_m!(octree1, ..., :k1)
    
    y .+= y_far
end
```

### 9.1 清零时机

每遍开始前必须清零 aggS 和 disaggG：

| 步骤 | 清零目标 | 时机 |
|------|---------|------|
| `clear_agg!(octree)` | 所有层的 `aggS` 和 `disaggG` | 每遍开始前 |
| `translate!` 内部 | 当层 `disaggG` | 转移前 (由 `translate!` 执行) |

> **与 EFIE 的区别**：EFIE 只需 1 遍，PMCHW 需要 4 遍。每遍必须独立清零，否则上一遍的残留数据会污染结果。

---

## 10. EMSuite PMCHW 与 EFIE 的代码复用关系

| 组件 | EFIE | PMCHW | 复用程度 |
|------|------|-------|---------|
| L 算子 (矩阵装配) | `assemble_impedance_matrix(EFIE, ...)` | 直接调用，不同 factor | 100% 复用 |
| K 算子 (矩阵装配) | 无 | `calc_k_pmchw_term!` (新增) | 新增 |
| 奇异积分 | `Singularity.jl` | L 复用；K 跳过自对角 | 部分复用 |
| 激励 V_E | `excitation_vector(EFIE, ...)` | 直接调用 | 100% 复用 |
| 激励 V_H | 无 | `_pmchw_excitation_H` (新增) | 新增 |
| RCS | `radarCrossSection(RWG, ...)` | 新增双流版本 | 新增 |
| 辐射积分 | `radiation_integral_rwg` | 直接调用 (分别对 J 和 M) | 100% 复用 |
| MLFMA 聚合 | `aggregate_leaf!` | `aggregate_leaf_pmchw!` (结构相同) | 接口相同 |
| MLFMA 上推/转移/下推 | `aggregate_upward!`, `translate!`, `disaggregate_downward!` | 直接调用 | 100% 复用 |
| MLFMA 解聚 | `disaggregate_leaf!` | `disaggregate_leaf_pmchw_j/m!` (新增) | 结构新增 |

---

## 11. 关键公式速查表

### 11.1 子块因子表

| 子块 | 直接法因子 | 含义 | 对应 Gibson |
|------|-----------|------|------------|
| $Z^{EJ}$ | $\frac{jk\eta}{16\pi} = \frac{j\omega\mu}{16\pi}$ | 矢量势 $\propto j\omega\mu$, 标量势 $\propto \frac{j}{\omega\varepsilon}$ | Eq.3.185 |
| $Z^{HM}$ | $\frac{jk}{\eta \cdot 16\pi} = \frac{j\omega\varepsilon}{16\pi}$ | 矢量势 $\propto j\omega\varepsilon$, 标量势 $\propto \frac{j}{\omega\mu}$ | Eq.3.190 (对偶) |
| $Z^{EM}$ | $\frac{1}{16\pi}$ | K 算子 (PV) | Eq.3.187 |
| $Z^{HJ}$ | $-\frac{1}{16\pi}$ | $-K$ (结构不变量) | Eq.3.188 |

### 11.2 MLFMA 解聚因子表

| 子块 | MLFMA 因子 | 接收类型 | 符号 |
|------|-----------|---------|------|
| EJ (L) | $jk\eta/(4\pi)$ | te_L | + |
| HJ (−K) | $jk/(4\pi)$ | te_K | + |
| EM (K) | $-jk/(4\pi)$ | te_K | − |
| HM (L^M) | $jk/(\eta \cdot 4\pi)$ | te_L | + |

### 11.3 散射截面公式

$$\sigma = \frac{k_0^2}{4\pi}(|E_\theta|^2 + |E_\phi|^2), \quad E_\theta = \eta_0 N_\theta + L_\phi, \quad E_\phi = \eta_0 N_\phi - L_\theta$$

---

## 附录 A：EMSuite 关键文件索引

| 文件 | 功能 |
|------|------|
| `src/IntegralEquations/PMCHW.jl` | PMCHW 结构体 + L/K 算子 + 2N×2N 矩阵装配 |
| `src/IntegralEquations/Excitation.jl` | PMCHW 平面波/Delta-Gap 激励向量 |
| `src/IntegralEquations/EFIE.jl` | L 算子复用 (Z^EJ, Z^HM) |
| `src/IntegralEquations/Impedance.jl` | 通用矩阵装配框架 |
| `src/IntegralEquations/Singularities.jl` | 奇异积分 (F₁, F₂₁, F₂₂) |
| `src/FastAlgorithms/MLFMA/PMCHWMLFMAOperator.jl` | PMCHW MLFMA 线性算子 (mul!, 4遍) |
| `src/FastAlgorithms/MLFMA/Aggregation.jl` | 上推聚合 |
| `src/FastAlgorithms/MLFMA/Translation.jl` | 转移函数 |
| `src/FastAlgorithms/MLFMA/Disaggregation.jl` | 下推解聚 |
| `src/FastAlgorithms/MLFMA/OctreeBuilder.jl` | 八叉树构建 |
| `src/PostProcessing/RCS.jl` | 双流 RCS 后处理 |
| `src/PostProcessing/RadiationIntegral.jl` | 辐射积分 |

## 附录 B：与 Legacy MLFMA 报告的关系

| MLFMA 报告章节 | PMCHW 对应 | 差异 |
|---------|---------|------|
| M1/M2 (网格/RWG) | 完全复用 | 无 |
| §0 (Addition Theorem) | 完全复用 | 无 |
| §1.1 (八叉树构建) | 2个八叉树 (k₀, k₁) | 参数不同 |
| §1.2 (aggSBF/disaggSBF) | 不预计算 (实时计算) | 结构相同 |
| §1.3 (相移因子) | 完全复用 | 无 |
| §2 (mul!) | 4遍并行 | 结构新增 |
| §3 (系数链验证) | 4个子块分别验证 | 新增 K 算子链 |
| §5 (工程技巧) | 清零策略更复杂 | 4遍清零 |
| §6 (disaggG 清零) | 每遍必须独立清零 | 更严格 |
| §8 (求解器) | GMRES + 2N×2N 算子 | 规模翻倍 |
| §9 (RCS) | 双流合成 | 新增 L 积分 |

## 附录 C：PMCHW 特有的工程注意事项

### C.1 矩阵条件数

PMCHW 矩阵的条件数通常比 EFIE 差，因为：
- Z^EJ 和 Z^HM 的物理量纲不同（$\eta$ vs $1/\eta$）
- K 算子无对角项，缺乏正则化

可能需要适当的预条件子（如块对角预条件）来加速收敛。

### C.2 symmetry 利用

$Z^{HJ} = -Z^{EM}$ 是精确的结构不变量。在存储和计算中可利用这一关系：
- 直接法：只计算 K 一次，Z^HJ 取反
- MLFMA：遍 1/2 的 K 接收用正号，遍 3/4 的 K 接收用负号

### C.3 有损介质注意事项

- $k_1$ 为复数时，$e^{-jk_1 R}$ 包含衰减项 $e^{-\text{Im}(k_1) R}$
- **计算**（相移因子、转移函数、前系数 factor）均使用完整复数 $k_1, \eta_1$，保证相位和衰减正确
- **八叉树几何**（叶盒边长、截断项数、极点采样）使用 $\text{Re}(k_1)$ 对应的波长 $\lambda_1 = 2\pi / \text{Re}(k_1)$。对高损耗介质（大的 $\text{Im}(k_1)$），八叉树参数可能不再最优，但不影响正确性
- 代码: `λ1 = 2π / real(pmchw.k1)`, `JK = im * k` (k 取完整复数)
