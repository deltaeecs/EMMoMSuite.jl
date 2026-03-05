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

---

## 附录 D：模块验证方案

> **目的**：为每个 PMCHW 模块提供独立、可执行的验证方案。每个验证项包含：目标、方法、预期结果（含容差）、参考测试编号。
>
> **测试夹具约定**：除非另行说明，默认使用 $r = 0.5$ m 球体, $f = 300$ MHz, $\varepsilon_r = 4.0$, $\mu_r = 1.0$, `lat_divs=4, lon_divs=6` (N=54)。

### D.1 材料参数与波数验证 (§1)

**验证目标**：确认 PMCHW 结构体从 $(f, \varepsilon_r, \mu_r)$ 正确计算 $k_0, \eta_0, k_1, \eta_1$。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D1.1 | $k_0$ 计算 | `pmchw.k0 == 2π·f/c₀` | 精确匹配 | `< eps()` |
| D1.2 | $\eta_0$ 计算 | `pmchw.eta0 == √(μ₀/ε₀)` | $\approx 376.73\;\Omega$ | `< eps()` |
| D1.3 | $k_1$ 无损 | `pmchw.k1 == k₀√(ε_r·μ_r)` | $k_1 = 2k_0$ (当 $\varepsilon_r=4$) | `< eps()` |
| D1.4 | $\eta_1$ 无损 | `pmchw.eta1 == η₀√(μ_r/ε_r)` | $\eta_1 = η_0/2$ (当 $\varepsilon_r=4$) | `< eps()` |
| D1.5 | $k_1$ 有损 | `PMCHW(f, 4.0-0.1im)` | `Im(k₁) > 0` (衰减) | 定性 |
| D1.6 | $\eta_1$ 有损 | 同上 | `Im(η₁) ≠ 0` (复数) | 定性 |
| D1.7 | PEC 极限 | `PMCHW(f, 1e6)` | $\eta_1 \to 0$, $k_1 \to \infty$ | 趋势正确 |

**已有测试覆盖**：部分 (22.2 含 PEC 极限和有损测试，但未逐项验证 $k_1, \eta_1$)

```julia
# 推荐验证代码
pmchw = PMCHW(300e6, 4.0)
@test pmchw.k0 ≈ 2π * 300e6 / 299792458.0
@test pmchw.eta0 ≈ √(4π*1e-7 / (1/(299792458.0^2 * 4π*1e-7)))
@test pmchw.k1 ≈ pmchw.k0 * √(4.0)
@test pmchw.eta1 ≈ pmchw.eta0 / √(4.0)
```

---

### D.2 2N×2N 矩阵结构不变量验证 (§2)

**验证目标**：确认 4 个子块满足理论预测的对称性和代数关系。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D2.1 | $Z^{EJ}$ 对称 | `norm(Z_EJ - Z_EJ') / norm(Z_EJ)` | 复对称 (非 Hermitian) | `< 1e-8` |
| D2.2 | $Z^{HM}$ 对称 | `norm(Z_HM - Z_HM') / norm(Z_HM)` | 复对称 | `< 1e-8` |
| D2.3 | $Z^{HJ} = -Z^{EM}$ | `norm(Z_HJ + Z_EM) / norm(Z_EM)` | 精确为零 | `< 1e-10` |
| D2.4 | 范数比 | `norm(Z_EJ) / norm(Z_HM)` | $\approx \eta_0^2 \approx 141866$ | 2 个数量级以内 |
| D2.5 | 所有子块非零 | `norm(Z_block) > 0` | 4 块均非零 | 定性 |
| D2.6 | $Z^{EJ}$ 非 Hermitian | `Z_EJ ≠ Z_EJ'` (一般) | 复对称 ≠ Hermitian | 定性 |

**提取子块**：
```julia
Z = assemble_impedance_matrix(pmchw, basis)
N = num_basis(basis)
Z_EJ = Z[1:N, 1:N]
Z_EM = Z[1:N, N+1:2N]
Z_HJ = Z[N+1:2N, 1:N]
Z_HM = Z[N+1:2N, N+1:2N]
```

**已有测试覆盖**：✅ 全覆盖 (test 22.2)

---

### D.3 L 算子 (EFIE 核) 因子验证 (§3.1)

**验证目标**：确认 `_l_block_operator` 生成的 EFIE 对象 factor 值正确。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D3.1 | EJ factor (区域0) | 检查 `efie_ej0.factor` | $jk_0\eta_0/(16\pi)$ | `< eps()` |
| D3.2 | EJ factor (区域1) | 检查 `efie_ej1.factor` | $jk_1\eta_1/(16\pi)$ | `< eps()` |
| D3.3 | HM factor (区域0) | 检查 `efie_hm0.factor` | $jk_0/(\eta_0 \cdot 16\pi)$ | `< eps()` |
| D3.4 | HM factor (区域1) | 检查 `efie_hm1.factor` | $jk_1/(\eta_1 \cdot 16\pi)$ | `< eps()` |
| D3.5 | 单区域比值 | $Z^{EJ}_\text{region0} / Z^{HM}_\text{region0}$ | $= \eta_0^2$ (精确) | `< 1e-12` |

**验证方法（D3.5 详述）**：构造 $\varepsilon_r = 1, \mu_r = 1$ 的 PMCHW，此时 $k_0 = k_1$, $\eta_0 = \eta_1$，故 $Z^{EJ} = 2L(k_0,\eta_0)$, $Z^{HM} = 2L_e(k_0,\eta_0)$，逐元素比值精确为 $\eta_0^2$：

```julia
pmchw_vac = PMCHW(300e6, 1.0, 1.0)  # ε_r=1
Z_vac = assemble_impedance_matrix(pmchw_vac, basis)
# 逐元素比值 (跳过零元素)
ratio = Z_vac[1:N,1:N] ./ Z_vac[N+1:2N, N+1:2N]
@test all(r -> abs(r - pmchw_vac.eta0^2) < 1e-10, ratio[abs.(Z_vac[N+1:2N,N+1:2N]) .> 1e-20])
```

**已有测试覆盖**：⚠️ 仅范数比 (D2.4), 未逐元素验证 (D3.5)

---

### D.4 K 算子核验证 (§3.2–3.3)

**验证目标**：确认 `calc_k_pmchw_term!` 和 `assemble_K_pmchw_offdiag` 的正确性。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D4.1 | 输出形状 | `size(K) == (N, N)` | N×N | 精确 |
| D4.2 | 非零 | `norm(K) > 0` | 非零 | 定性 |
| D4.3 | 复数值 | `eltype(K) <: Complex` | 复数 | 精确 |
| D4.4 | 自对角为零 | 对任意 RWG $n$，其自对角元素应为零或极小 | $\|K[n,n]\| < \epsilon \cdot \|K\|_\infty$ | `< 1e-6` |
| D4.5 | 非对称 | `norm(K - K') / norm(K)` | 非对称 ($K \neq K^T$) | `> 1e-3` |
| D4.6 | 双区域叠加 | $K_\text{total} = K(k_0) + K(k_1)$ | 分别装配后相加 = 直接装配 | `< 1e-10` |
| D4.7 | η=1 因子 | 检查 K 算子使用 `MFIE(0, k, 1, gq)` | factor 为 $1/(16\pi)$，不含 $\eta$ | 代码审查 |

**D4.4 自对角验证的物理意义**：PMCHW 的 K 算子没有 mass matrix ($\frac{1}{2}$ 对角项被两侧相加消去)。因此自对角元素仅来自 PV 积分，且代码中 `t_test == t_src` 时跳过，所以 K 矩阵的"块对角"（共享三角形的基函数对）应显著小于非对角元素。但 $K[n,n]$ 不一定精确为零，因为同一 RWG 基函数的两个支撑三角形 $T^+, T^-$ 之间的 K 相互作用不为零。

**D4.6 验证代码**：
```julia
# 分别装配 K(k₀) 和 K(k₁)
K0 = assemble_K_pmchw_offdiag(pmchw.k0, basis)  
K1 = assemble_K_pmchw_offdiag(pmchw.k1, basis)
K_sum = K0 + K1

# 直接从 PMCHW 装配获取
Z = assemble_impedance_matrix(pmchw, basis)
K_direct = Z[1:N, N+1:2N]  # Z^EM = K₀ + K₁

@test norm(K_sum - K_direct) / norm(K_direct) < 1e-10
```

**已有测试覆盖**：✅ D4.1–D4.3 (test 22.1), ❌ D4.4–D4.6

---

### D.5 奇异积分验证 (§3.4)

**验证目标**：确认 L 算子的解析奇异项 ($F_1, F_{21}, F_{22}$) 和 K 算子的自对角跳过行为。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D5.1 | L 算子自对角非零 | $Z^{EJ}[n,n] \neq 0$ | 自对角元素是虚部主导（$jk\eta \cdot F_1$） | 定性 |
| D5.2 | L 算子自对角主导 | $\|Z^{EJ}_\text{diag}\|_\infty / \|Z^{EJ}_\text{offdiag}\|_\infty$ | $> 1$（自对角不弱于非对角） | 定性 |
| D5.3 | K 算子跳过验证 | 将 `t_test == t_src` 跳过条件移除 | 自对角贡献应极小（PV 积分几乎为零） | `< 1e-4 × ‖K‖` |

**已有测试覆盖**：❌

---

### D.6 激励向量验证 (§4)

**验证目标**：确认平面波和 Delta-Gap 激励计算正确。

#### D.6.1 平面波激励

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D6.1 | V 向量长度 | `length(V) == 2N` | 2N | 精确 |
| D6.2 | V_E 非零 | `norm(V[1:N]) > 0` | 非零 | 定性 |
| D6.3 | V_H 非零 | `norm(V[N+1:2N]) > 0` | 非零 | 定性 |
| D6.4 | 幅度比 | `norm(V[N+1:2N]) / norm(V[1:N])` | $\approx 1/\eta_0 \approx 2.65 \times 10^{-3}$ | 1 个数量级 |
| D6.5 | V_E 与 EFIE 一致 | `V[1:N] == excitation_vector(EFIE(freq), pw, basis)` | 精确相等 | `< eps()` |

**D6.4 的物理推导**：$V_H = \int \mathbf{f} \cdot \mathbf{H}^{inc}$, $V_E = \int \mathbf{f} \cdot \mathbf{E}^{inc}$。平面波中 $|\mathbf{H}^{inc}| = |\mathbf{E}^{inc}|/\eta_0$，故 $|V_H|/|V_E| \approx 1/\eta_0$。

#### D.6.2 Delta-Gap 激励

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D6.6 | V_H 全零 | `V[N+1:2N] == zeros(N)` | 精确零 | `< eps()` |
| D6.7 | 馈电元素 | `V[feed_idx] == V_gap × l_edge` | 精确 | `< eps()` |
| D6.8 | 非馈电元素 | `V[i] == 0` for `i ≠ feed_idx` | 精确零 | `< eps()` |

**已有测试覆盖**：✅ D6.1–D6.4 (test 22.3), ✅ D6.6–D6.8 (test 15.1–15.2)

---

### D.7 双流 RCS 验证 (§6)

**验证目标**：确认从双电流 $(\mathbf{I}_J, \mathbf{I}_M)$ 正确合成远场和 RCS。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D7.1 | RCS 非负 | `all(σ .>= 0)` | 物理约束 | 精确 |
| D7.2 | factor 正确 | 检查代码中 `k₀²/(4π)` | ≠ EFIE 的 $(k_0\eta_0)^2/(4\pi)$ | 代码审查 |
| D7.3 | M=0 退化 | 令 `I_M = zeros(N)`，比较 PMCHW RCS 与 EFIE RCS | 相等 ($\sigma_\text{PMCHW} = \sigma_\text{EFIE}$) | `< 1e-10` |
| D7.4 | J=0 退化 | 令 `I_J = zeros(N)`，仅 M 贡献 | RCS 仅含 $\|L\|^2$ 项 | 定性 |
| D7.5 | 远场分量合成 | 手动计算 $E_\theta = \eta_0 N_\theta + L_\phi$ | 与代码输出一致 | `< 1e-10` |
| D7.6 | PEC 极限一致性 | $\varepsilon_r \to \infty$ 时 PMCHW RCS → 某基准 | 求解收敛，RCS 有限 | 定性 |

**D7.3 详述——PMCHW vs EFIE RCS factor 差异**：
- EFIE RCS: $\sigma_\text{EFIE} = \frac{(k_0\eta_0)^2}{4\pi} |N|^2$
- PMCHW RCS: $\sigma_\text{PMCHW} = \frac{k_0^2}{4\pi} (|\eta_0 N_\theta + L_\phi|^2 + |\eta_0 N_\phi - L_\theta|^2)$

当 $\mathbf{I}_M = 0$ (即 $L=0$) 时：$\sigma_\text{PMCHW} = \frac{k_0^2 \eta_0^2}{4\pi} |N|^2 = \sigma_\text{EFIE}$。两者**在此特殊情况下相等**。

```julia
# D7.5 验证代码
N_θϕ = radiation_integral_rwg(r_info, basis, I_J)
L_θϕ = radiation_integral_rwg(r_info, basis, I_M)

E_θ_manual = pmchw.eta0 * N_θϕ[1] + L_θϕ[2]
E_ϕ_manual = pmchw.eta0 * N_θϕ[2] - L_θϕ[1]
σ_manual = pmchw.k0^2 / (4π) * (abs2(E_θ_manual) + abs2(E_ϕ_manual))

@test σ_manual ≈ σ_code
```

**已有测试覆盖**：❌ (RCS 完全未测试)

---

### D.8 MLFMA 聚合验证 (§7.3)

**验证目标**：确认叶层聚合 (`aggregate_leaf_pmchw!`) 正确计算辐射函数。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D8.1 | aggS 非零 | 单位电流输入后 `norm(aggS) > 0` | 非零 | 定性 |
| D8.2 | aggS 线性 | `aggS(2x) == 2·aggS(x)` | 精确线性 | `< 1e-12` |
| D8.3 | J/M 无差别 | 同样系数传入 1:N 或 N+1:2N | `aggS` 完全相同 | `< eps()` |
| D8.4 | 与 EFIE 一致 | 对比 EFIE `aggregate_leaf!` 输出 | 相同辐射函数 | `< 1e-12` |

**D8.3 的物理意义**：聚合函数只看 RWG 基函数几何和输入系数向量，不区分电流/磁流类型。传入同样的系数，无论从 x[1:N] 还是 x[N+1:2N] 提取，聚合结果必须相同。

**已有测试覆盖**：❌ (聚合未隔离测试)

---

### D.9 MLFMA 接收函数验证 (§7.7)

**验证目标**：确认 `_receive_terms` 正确产生 L 型和 K 型两种接收模式。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D9.1 | te_L 非零 | 给定非零 disaggG | `te_L ≠ 0` | 定性 |
| D9.2 | te_K 非零 | 给定非零 disaggG | `te_K ≠ 0` | 定性 |
| D9.3 | te_L/te_K 独立 | 在一般情况下 | `te_L ≠ te_K` | 定性 |
| D9.4 | K 型叉积关系 | 检查 $\hat{r} \times (\hat{\theta}, \hat{\phi})$ 展开 | $E_\theta \hat{\phi} - E_\phi \hat{\theta}$ | 代码审查 |
| D9.5 | 与 EFIE disagg 比对 | EFIE 解聚仅用 te_L | PMCHW te_L 与 EFIE 一致 | `< 1e-12` |

**D9.5 验证方法**：在相同八叉树和相同输入下，EFIE 的 `disaggregate_leaf!` 产生的结果应等于 PMCHW `_receive_terms` 返回的 `te_L` 乘以对应因子：

```julia
# 概念性伪代码 (需在测试中注入中间变量)
te_L_pmchw, te_K_pmchw = _receive_terms(bf, basis, field, r0, k, η, poles)
te_efie = disaggregate_leaf_efie_equivalent(...)  # EFIE 版本
@test te_L_pmchw ≈ te_efie
```

**已有测试覆盖**：❌ (接收函数未隔离测试)

---

### D.10 MLFMA 解聚因子验证 (§7.5–7.6)

**验证目标**：确认 J 通道和 M 通道解聚使用正确的因子和符号。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D10.1 | factor_EJ 值 | 数值检查 | $jk_l\eta_l/(4\pi)$ | `< eps()` |
| D10.2 | factor_HJ 值 | 数值检查 | $+jk_l/(4\pi)$（正号） | `< eps()` |
| D10.3 | factor_EM 值 | 数值检查 | $-jk_l/(4\pi)$（负号） | `< eps()` |
| D10.4 | factor_HM 值 | 数值检查 | $jk_l/(\eta_l \cdot 4\pi)$ | `< eps()` |
| D10.5 | J 通道写入位置 | disagg_j 写入 y[bfID] 和 y[bfID+N] | 同时写 E 行和 H 行 | 代码审查 |
| D10.6 | M 通道写入位置 | disagg_m 写入 y[bfID] 和 y[bfID+N] | 同时写 E 行和 H 行 | 代码审查 |
| D10.7 | HJ/EM 符号关系 | `factor_HJ / factor_EM` | $= -1$ | `< eps()` |

**D10.7 推导**：$Z^{HJ} = -Z^{EM}$。J 通道的 H 行和 M 通道的 E 行都使用 `te_K` 接收，因此符号关系由 factor 携带：`factor_HJ × te_K = -factor_EM × te_K`，即 `factor_HJ = -factor_EM`。代入验证：$+jk/(4\pi) = -(-jk/(4\pi))$ ✓。注意这两个因子都不含 $\eta$，因此比值为 $-1$ 而非 $-\eta$。

```julia
# 因子数值验证
k = pmchw.k0; η = pmchw.eta0
@test factor_EJ ≈ im * k * η / (4π)
@test factor_HJ ≈ im * k / (4π)           # 正号
@test factor_EM ≈ -im * k / (4π)           # 负号
@test factor_HM ≈ im * k / (η * 4π)
@test factor_HJ ≈ -factor_EM               # 反号关系
@test factor_EJ / factor_HM ≈ η^2          # 与直接法比值一致
```

**已有测试覆盖**：❌ (因子未单独验证)

---

### D.11 MLFMA 系数链端到端验证 (§7.8)

**验证目标**：通过 MLFMA matvec 与 Direct matvec 的对比，端到端验证所有系数正确。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D11.1 | matvec 整体 | `\|y_{mlfma} - y_{direct}\| / \|y_{direct}\|` | 远场截断误差 | `< 10%` |
| D11.2 | E 行分量 | `\|y[1:N]_{mlfma} - y[1:N]_{direct}\| / \|y[1:N]_{direct}\|` | 分量级误差 | `< 15%` |
| D11.3 | H 行分量 | `\|y[N+1:2N]_{mlfma} - y[N+1:2N]_{direct}\|` | 分量级误差 | `< 15%` |
| D11.4 | 随机输入 | 多个随机 x 向量 | 方差小于均值 | 统计 |

**D11.1 为何允许 10% 误差**：MLFMA 使用有限阶多级展开近似远场 Green 函数，对于小网格 (N=54) 远/近交界处截断误差较大。更大的网格 ($N > 500$) 误差通常 $< 3\%$。

```julia
Z_direct = assemble_impedance_matrix(pmchw, basis)
x = randn(ComplexF64, 2N)
y_direct = Z_direct * x

y_mlfma = similar(y_direct)
mul!(y_mlfma, op, x)

rel_err = norm(y_mlfma - y_direct) / norm(y_direct)
@test rel_err < 0.10
```

**已有测试覆盖**：✅ D11.1 (test 15.11)

---

### D.12 近场矩阵验证 (§8)

**验证目标**：确认稀疏近场矩阵正确提取了密集矩阵的邻近元素。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D12.1 | 稀疏性 | `nnz(Z_near) < (2N)^2` | 真稀疏 | 严格小于 |
| D12.2 | 4 子块非零 | 检查 4 个 N×N 子块均有非零元素 | 全部非零 | 定性 |
| D12.3 | HJ = -EM | `Z_near[N+1:2N, 1:N] + Z_near[1:N, N+1:2N]` | 零矩阵 | `< 1e-8 × ‖Z_near‖` |
| D12.4 | 对角元素匹配 | `diag(Z_near) ≈ diag(Z_full)` | 对角元素必在近场中 | `< 1e-12` |
| D12.5 | 近/远互补 | `Z_near + Z_far ≈ Z_full` | 分解完整 | `< 1e-10` |

**D12.5 验证方法**：Z_far 可通过 `Z_full - Z_near` 获得。但更直接的方式是检查近场矩阵的非零位置覆盖了所有自身和邻居基函数对：

```julia
Z_full = assemble_impedance_matrix(pmchw, basis)
Z_near = op.Z_near  # PMCHWMLFMAOperator 的近场稀疏矩阵

# 对角匹配
@test norm(diag(Z_near) - diag(Z_full)) / norm(diag(Z_full)) < 1e-12

# 稀疏性
@test nnz(Z_near) < (2N)^2
```

**已有测试覆盖**：✅ D12.1–D12.3 (test 15.8), ⚠️ D12.4–D12.5 (未单独验证)

---

### D.13 mul! 流程验证 (§9)

**验证目标**：确认 `mul!(y, A, x)` 正确执行近场+4遍远场，结果与直接法一致。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D13.1 | 零输入 | `mul!(y, A, zeros(2N))` | `y == zeros(2N)` | `< eps()` |
| D13.2 | 线性 | `mul!(y, A, α·x) ≈ α·mul!(y, A, x)` | 精确线性 | `< 1e-10` |
| D13.3 | 输出长度 | `length(y) == 2N` | 2N | 精确 |
| D13.4 | E/H 行均非零 | `norm(y[1:N]) > 0`, `norm(y[N+1:2N]) > 0` | 双行块均有贡献 | 定性 |
| D13.5 | 4 遍清零 | 连续两次 mul! 结果相同 | 无残留 | `< 1e-12` |
| D13.6 | GMRES 收敛 | 使用 MLFMA 算子求解 PMCHW | GMRES 迭代 < 100 | 收敛 |

**D13.5 验证代码**：
```julia
x = randn(ComplexF64, 2N)
y1 = zeros(ComplexF64, 2N)
y2 = zeros(ComplexF64, 2N)
mul!(y1, op, x)
mul!(y2, op, x)
@test norm(y1 - y2) / norm(y1) < 1e-12
```

**已有测试覆盖**：✅ D13.1–D13.3 (test 15.11), ⚠️ D13.4–D13.6 (部分覆盖)

---

### D.14 端到端求解验证

**验证目标**：确认 PMCHW 从网格到 RCS 的完整流程正确。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D14.1 | 直接法收敛 | LU 求解 2N×2N 系统 | 残差 $\|Zx-b\|/\|b\| < 10^{-10}$ | `< 1e-10` |
| D14.2 | I_J 非零 | `norm(I[1:N]) > 0` | 非零电流 | 定性 |
| D14.3 | I_M 非零 | `norm(I[N+1:2N]) > 0` | 非零磁流 | 定性 |
| D14.4 | I_M/I_J 比值 | `norm(I_M) / norm(I_J)` | $\approx 1/\eta_0$（量级） | 1–2 个数量级 |
| D14.5 | MLFMA Z_in | `Z_in_mlfma ≈ Z_in_direct` | 输入阻抗一致 | `< 5%` |
| D14.6 | RCS 物理合理 | 介质球 RCS $< $ PEC 球 RCS | 通过率 | 定性 |

**D14.4 的物理推导**：从矩阵结构 $Z^{EJ} I_J + Z^{EM} I_M = V_E$ 和 $\|Z^{EJ}\| \sim \eta_0^2 \|Z^{HM}\|$，可推出 $\|I_M\| / \|I_J\| \sim 1/\eta_0$。

**已有测试覆盖**：✅ D14.1–D14.3 (test 22.2 end-to-end), ✅ D14.5 (test 15.11 B2)

---

### D.15 MLFMA 4 遍隔离验证 (§7.1)

**验证目标**：确认每一遍 MLFMA 写入正确的 y 子块。

| # | 验证项 | 方法 | 预期 | 容差 |
|---|--------|------|------|------|
| D15.1 | 遍 1 (J×k₀) | 仅执行遍 1，检查 y[1:N] 和 y[N+1:2N] | 均非零 (L→E行, -K→H行) | 定性 |
| D15.2 | 遍 3 (M×k₀) | 仅执行遍 3，检查 y[1:N] 和 y[N+1:2N] | 均非零 (K→E行, L→H行) | 定性 |
| D15.3 | 遍 1+2 累加 | 遍 1 + 遍 2 = J 列的完整远场 | 对比直接法 $Z_\text{far}[:, 1:N] \cdot x_J$ | `< 15%` |
| D15.4 | 4 遍总和 | 全部 4 遍累加 | $\approx Z_\text{full} \cdot x$ (含近场) | `< 10%` |

**验证方法（需要代码修改以隔离单遍）**：

方案 A (推荐)：构造特殊输入使部分遍贡献为零
```julia
# 仅测试 J 列：令 I_M = 0
x_J_only = [randn(ComplexF64, N); zeros(ComplexF64, N)]
y_mlfma = zeros(ComplexF64, 2N); mul!(y_mlfma, op, x_J_only)
y_direct = Z_full * x_J_only
@test norm(y_mlfma - y_direct) / norm(y_direct) < 0.15

# 仅测试 M 列：令 I_J = 0  
x_M_only = [zeros(ComplexF64, N); randn(ComplexF64, N)]
y_mlfma = zeros(ComplexF64, 2N); mul!(y_mlfma, op, x_M_only)
y_direct = Z_full * x_M_only
@test norm(y_mlfma - y_direct) / norm(y_direct) < 0.15
```

**已有测试覆盖**：❌ (逐遍验证未实现)

---

### D.16 验证覆盖矩阵

| 模块 | 报告章节 | 验证项 | 已有测试 | 建议优先级 |
|------|---------|--------|---------|-----------|
| 材料参数 | §1 | D1.1–D1.7 | ⚠️ 部分 | 低 |
| 矩阵不变量 | §2 | D2.1–D2.6 | ✅ 全覆盖 | — |
| L 算子因子 | §3.1 | D3.1–D3.5 | ⚠️ 仅比值 | 中 |
| K 算子核 | §3.2 | D4.1–D4.7 | ⚠️ 部分 | 中 |
| 奇异积分 | §3.4 | D5.1–D5.3 | ❌ | 低 |
| 平面波激励 | §4.1 | D6.1–D6.5 | ✅ 全覆盖 | — |
| DeltaGap 激励 | §4.2 | D6.6–D6.8 | ✅ 全覆盖 | — |
| 双流 RCS | §6 | D7.1–D7.6 | ❌ | **高** |
| MLFMA 聚合 | §7.3 | D8.1–D8.4 | ❌ | 中 |
| MLFMA 接收 | §7.7 | D9.1–D9.5 | ❌ | 中 |
| MLFMA 解聚因子 | §7.5–7.6 | D10.1–D10.7 | ❌ | **高** |
| MLFMA 系数链 | §7.8 | D11.1–D11.4 | ✅ 整体 | — |
| 近场矩阵 | §8 | D12.1–D12.5 | ✅ 大部分 | 低 |
| mul! 流程 | §9 | D13.1–D13.6 | ✅ 大部分 | 低 |
| 端到端求解 | §5 | D14.1–D14.6 | ✅ 大部分 | — |
| 4 遍隔离 | §7.1 | D15.1–D15.4 | ❌ | **高** |

> **建议实现优先级**：D7 (RCS) > D10 (解聚因子) > D15 (4 遍隔离) > D4 (K 核隔离) > D3 (L 因子) > D8 (聚合) > D9 (接收) > D1 (材料) > D5 (奇异积分)
