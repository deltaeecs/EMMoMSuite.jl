# Legacy MLFMA 算法完整流程报告

> **目的**：精确提炼 Legacy (`MoM_Kernels`) 的 MLFMA 算法全流程,对比 EMSuite 实现, 定位系数差异。  
> **参考**：Legacy 代码 + CEM 知识库 (Gibson MoM Ch.21, Ergul MLFMA Ch.2, formulas.md §6)  
> **时间约定**：$e^{-j\omega t}$ ⟹ Green 函数 $G = \frac{e^{-jkR}}{4\pi R}$, Hankel 函数用 $h_l^{(2)}$, $\mathrm{JK\_0} = jk = j\omega/c$

---

## 总体数据流

```
Nastran .nas 文件
  ↓ MeshProcess (Pass1: 计数, Pass2: 解析)
[node(3,N_node), triangles(3,N_tri)]
  ↓ MeshAndBFs / RWG 构造
[TriangleInfo 数组, RWG 基函数数组]
  ↓ getOctreeAndReOrderBFs!
[八叉树(levels 1..nLevels), 重排后的基函数编号]
  ↓ calZnearCSC
[近场稀疏矩阵 Znear (CSC格式)]
  ↓ MLFMAIterator 构造
[MLFMA 算子 Zopt (含 aggSBF, disaggSBF, Znear)]
  ↓ Solvers.solve (GMRES/BiCGSTAB)
[电流系数 I(N_bf)]
  ↓ FarField / RCS
[E(θ,φ), σ(θ,φ)]
```

---

## M1. 网格读取 (`MeshProcess.jl`)

### 格式: Nastran 固定格式 `.nas`

两遍扫描:
1. **Pass 1**: 统计 GRID, CTRIA3, CTETRA, CHEXA 记录数 → 预分配数组
2. **Pass 2**: 解析节点坐标和单元连接关系
   - 旧节点编号 → 新编号 (1:nodenum)，通过 `nodeO2LID` 字典映射

### 输出:
```julia
struct MeshNodeTriTetraHexa{IT, FT}
    node[3, nodenum]          # 节点坐标
    triangles[3, trinum]      # 三角形连接关系
    tetrahedras[4, tetranum]  # 四面体连接关系
    hexahedras[8, hexanum]    # 六面体连接关系
end
```

---

## M2. RWG 基函数构建 (`RWG.jl`)

### 构建算法 (`rwgbfConstructerTrianglesInfoModifiers!`)

1. 创建边池: 每个三角形的3条边 → $3 \times N_\text{tri}$ 条候选
2. 按边的两个节点编号排序 → 相同边相邻
3. 分配 RWG 编号:
   - 出现2次的边 = **全基函数** (full BF), 跨两个三角形
   - 出现1次的边 = **半基函数** (half BF, 边界边)
4. 填充结构体: 边长、中心坐标、三角形引用

### RWG 数据结构:
```julia
struct RWG{IT, FT}
    isbd::Bool                # 是否为边界半基函数
    bfID::IT                  # 重编号后的ID
    edgel::FT                 # 边长 (含符号 ±)
    inGeo::MVector{2, IT}     # [tri₊, tri₋]
    inGeoID::MVector{2, IT}   # 在各三角形中的局部边序号(1:3)
    center::MVec3D{FT}        # 边中点 (用于八叉树分配)
end
```

### RWG 数学定义:

$$\mathbf{f}_n(\mathbf{r}) = \begin{cases} +\frac{l_n}{2A_+} \boldsymbol{\rho}_+ & \mathbf{r} \in T_+ \\ -\frac{l_n}{2A_-} \boldsymbol{\rho}_- & \mathbf{r} \in T_- \end{cases}$$

其中 $\boldsymbol{\rho}_\pm = \mathbf{r} - \mathbf{v}_\text{opp}^\pm$ (求积点到对顶点矢量).

散度: $\nabla \cdot \mathbf{f}_n = \pm \frac{l_n}{A_\pm}$

### 三角形几何量:
```julia
vertices[3]          # 3个顶点
edges[3]             # e₁, e₂, e₃ = 边向量
edgel[3]             # 带符号的边长
normal               # 法向量 = (e₁ × e₂) / |e₁ × e₂|
area                 # 面积 = ½|e₁ × e₂|
ρ[3,3]               # ρᵢ = 从第i条边到对顶点的向量
```

---

## 0. 理论基础：Addition Theorem

### 0.1 Green 函数多极展开

自由空间标量 Green 函数可分解为：

$$G(\mathbf{r}, \mathbf{r}') = \frac{e^{-jkR}}{4\pi R} \approx \frac{-jk}{(4\pi)^2} \int d^2\hat{k}\; e^{-j\mathbf{k}\cdot\mathbf{d}_1}\; T_L(k, \hat{k}, \mathbf{D})\; e^{+j\mathbf{k}\cdot\mathbf{d}_2}$$

其中：
- $\mathbf{r} = \mathbf{r}_a + \mathbf{d}_1$ (场点 = 场盒中心 + 局部偏移), $\mathbf{r}' = \mathbf{r}_b + \mathbf{d}_2$ (源点 = 源盒中心 + 局部偏移)
- $\mathbf{D} = \mathbf{r}_a - \mathbf{r}_b$ （盒子中心间向量，从源盒指向场盒）
- $T_L(k, \hat{k}, \mathbf{D}) = \sum_{l=0}^{L} (-j)^l (2l+1)\; h_l^{(2)}(k|\mathbf{D}|)\; P_l(\hat{k}\cdot\hat{D})$ （转移函数）
- $\frac{-jk}{(4\pi)^2}$ 是 **理论前系数**

### 0.2 球面积分的数值离散

球面积分 $\int d^2\hat{k}$ 用 $N_\text{poles}$ 个采样点离散化：

$$\int d^2\hat{k}\; f(\hat{k}) \approx \sum_{p=1}^{N_\text{poles}} W_p\; f(\hat{k}_p)$$

其中 $\sum_p W_p = 4\pi$（单位球面积）。

采样方案：
- $\theta$: Gauss-Legendre 求积, $N_\theta = L+1$ 个点
- $\phi$: 均匀采样, $N_\phi = 2(L+1)$ 个点
- 总极点数 $N_\text{poles} = N_\theta \times N_\phi = 2(L+1)^2$
- $W_p = W_\theta \cdot W_\phi$ （权重乘积）

### 0.3 截断项数

$$L = \lfloor 2\pi \frac{d}{\lambda}\sqrt{3} + 2.16 \cdot N_\text{digits}^{2/3} \cdot \left(2\pi\frac{d}{\lambda}\right)^{1/3} \rfloor$$

Legacy 默认 $N_\text{digits} = 3$。

---

## 1. 预处理阶段

### 1.1 八叉树构建 (`OctreeInfo`)

```
输入: 基函数中心坐标 bfCenters(3×N), 叶层盒子边长 leafCubeEdgel
输出: nLevels 层的八叉树

步骤:
1. 计算包围盒: CubeEdgel = max(Δx, Δy, Δz) + (√2 - 1) × leafEdgel (padding)
2. 计算层数: nLevels = ⌈log₂(CubeEdgel / leafEdgel)⌉
3. 计算大盒子: BigCubeEdgel = leafEdgel × 2^nLevels
4. 构建叶层: 按空间坐标将基函数分配到叶层盒子
5. 逐层向上构建: 每 2×2×2 个子盒子合并为一个父盒子
6. 重排基函数编号: 使同一盒子内的基函数编号连续 (关键优化!)
7. 计算邻近关系: 对每层计算 neighbors (≤27) 和 farneighbors (≤316)
```

**叶层盒子边长规则**：
- 三角形面网格：$d_\text{leaf} = 0.23\lambda$ (Legacy 默认)
- 四面体体网格：$d_\text{leaf} = 1.75 \times \bar{l}_\text{edge}$

**远亲定义**：在 $7^3$ 范围内偏移 $|i|>1 \lor |j|>1 \lor |k|>1$ 的盒子，最多 $7^3 - 3^3 = 316$ 个方向。远亲盒子必须由父盒子的邻盒子的子盒子中扣除自身邻盒子产生。

### 1.2 预计算辐射/配置函数 (`aggSBF`, `disaggSBF`)

Legacy **预计算**每个基函数在每个极点方向上的辐射模式，存为 `aggSBF(nPoles, 2, nBF)` 和 `disaggSBF(nPoles, 2, nBF)`。

#### RWG (EFIE) 辐射函数

$$\text{aggSBF}[p, \text{pol}, n] = \sum_{gi} (\hat{e}_\text{pol} \cdot \boldsymbol{\rho}_{n,gi}) \cdot w_{gi} \cdot \frac{l_n}{2} \cdot e^{+jk\hat{r}_p \cdot (\mathbf{r}_{gi} - \mathbf{c})}$$

其中：
- $p$ = 极点索引, $\text{pol} \in \{θ, ϕ\}$
- $\boldsymbol{\rho}_{n,gi} = \mathbf{r}_{gi} - \mathbf{v}_\text{opp}$ (自由顶点到求积点矢量)
- $l_n$ = RWG 边长（正值，符号通过 $\boldsymbol{\rho}$ 隐式编码）
- $w_{gi}/2$ = 三角形高斯求积权重 / 2
- $\mathbf{c}$ = 所在盒子中心
- $\hat{r}_p, \hat{\theta}_p, \hat{\phi}_p$ = 第 $p$ 个极点的方向矢量

> **注意**：Legacy 中 `ln = tri.edgel[ni]` 是**正值边长**（不含符号）。RWG 的 $\pm$ 方向通过 `ρ = rgs - tri.vertices[:,ni]` 隐式编码——每个三角形对应的自由顶点不同，自然产生正负方向。

#### EFIE 配置函数

$$\text{disaggSBF}[p, \text{pol}, n] = \overline{\text{aggSBF}[p, \text{pol}, n]}$$

即**取共轭**。物理含义：辐射用 $e^{+jk\hat{r}\cdot\mathbf{r}}$，配置用 $e^{-jk\hat{r}\cdot\mathbf{r}}$。

#### MFIE 配置函数（不同！）

MFIE 时 `disaggSBF` ≠ `conj(aggSBF)`：

$$\text{disaggSBF}_\text{MFIE}[p, \text{pol}, n] = \sum_{gi} (\hat{e}_\text{pol} \cdot (\boldsymbol{\rho} \times \hat{n} \times \hat{r}_p)) \cdot w_{gi} \cdot \frac{l_n}{2} \cdot e^{-jk\hat{r}_p \cdot (\mathbf{r}_{gi} - \mathbf{c})}$$

#### CFIE 配置函数  

$$\text{disaggSBF}_\text{CFIE} = \alpha \cdot \boldsymbol{\rho} + (1-\alpha) \cdot (\boldsymbol{\rho} \times \hat{n} \times \hat{r})$$

### 1.3 预计算相移因子 (`phaseShift2Kids`, `phaseShiftFromKids`)

对每层 $l \in [2, n_\text{Levels}-1]$，子盒子中心相对父盒子中心有 8 个可能偏移：

$$\Delta\mathbf{r}_\text{kid} = \frac{d_\text{parent}}{4} \cdot [(\pm 1, \pm 1, \pm 1)]^T = \frac{d_\text{child}}{2} \cdot [(\pm 1, \pm 1, \pm 1)]^T$$

$$\text{phaseShift2Kids}[p, i_\text{kid}] = e^{-jk\hat{r}_p \cdot \Delta\mathbf{r}_\text{kid}}$$

$$\text{phaseShiftFromKids}[p, i_\text{kid}] = \overline{\text{phaseShift2Kids}[p, i_\text{kid}]}$$

**工程技巧**：利用对称性，只计算4个偏移（$z<0$），另外4个（$z>0$）取共轭。

### 1.4 预计算转移因子 $\alpha_\text{Trans}$

对每层 $l \in [2, n_\text{Levels}]$ 的最多 316 个远亲方向：

$$\boxed{\alpha_\text{Trans}[p, i_\text{far}] = \frac{-jk}{(4\pi)^2} \cdot W_p \cdot \sum_{l=0}^{L} (-j)^l (2l+1)\; h_l^{(2)}(kR_{ab})\; P_l(\hat{r}_p \cdot \hat{R}_{ab})}$$

Legacy 代码：
```julia
mjKdiv16π² = -Params.JK_0 / (4π)^2      # = -jk / 16π²
αTrans[iPole, iFarNei] = αTransTemp * mjKdiv16π² * Wθϕs[iPole]
```

> **关键**：前系数 = $\frac{-jk}{(4\pi)^2} = \frac{-jk}{16\pi^2}$，**而非 $\frac{-jk}{4\pi}$**。$W_p$ 包含在 $\alpha_\text{Trans}$ 中。

### 1.5 预计算插值矩阵

采用两步 Lagrange 插值（先 $\phi$ 后 $\theta$）：

$$\mathbf{F}_\text{fine} = \mathbf{M}_\theta \cdot (\mathbf{M}_\phi \cdot \mathbf{F}_\text{coarse})$$

反插值（anterpolation - transposed interpolation）：

$$\mathbf{F}_\text{coarse} = \mathbf{M}_\phi^T \cdot (\mathbf{M}_\theta^T \cdot \mathbf{F}_\text{fine})$$

插值阶数默认 6 (即使用 6 个 Lagrange 基函数)。

### 1.6 近场矩阵组装

对叶层每个盒子及其邻盒子，直接计算 EFIE 矩阵元，组装成 CSC 稀疏矩阵 $\mathbf{Z}_\text{near}$。

EFIE 矩阵元（远场型，$R \geq R_\text{sglr}$）：

$$Z_{mn} = \frac{jk\eta}{16\pi} \cdot l_m \cdot l_n \int_{A_t} \int_{A_s} \left[(\boldsymbol{\rho}_m \cdot \boldsymbol{\rho}_n) - \frac{4}{k^2}\right] G(R)\; dA_t\; dA_s$$

其中 $G$ 已包含 $\frac{1}{4\pi R}$，$\frac{4}{k^2}$ 来自散度项 $(\nabla\cdot\mathbf{f}_m)(\nabla\cdot\mathbf{f}_n)$。

---

## 2. 运行时矩阵向量乘 (`mul!`)

### 2.0 总体流程

```
y = Z_near * x                           ← 近场
calZfarI!(Zopt, x):                       ← 远场
  ① aggOnBF!(leafLevel, aggSBF, x)       ← 叶层聚合
  ② agg2Level2!(levels, nLevels)          ← 聚合到第2层
  ③ transOnLevels!(levels, nLevels)       ← 各层转移
  ④ disagg2LeafLevel!(levels, nLevels)    ← 解聚到叶层
  ⑤ disaggOnBF!(leafLevel, disaggSBF, ZI) ← 叶层解聚
y += ZI                                   ← 合并
```

### 2.1 第①步：叶层聚合 (`aggOnBF!`)

$$\text{aggS}[p, \text{pol}, c] = \sum_{n \in \text{cube}(c)} I_n \cdot \text{aggSBF}[p, \text{pol}, n]$$

```
aggS .= 0    ← 每次mul!先清零
aggS[:,:,iCube] += IVec[n] * aggSBF[:,:,n]   ← 对基函数n累加
```

**物理含义**：将电流系数 $I_n$ 乘以预计算的辐射模式，得到盒子 $c$ 的辐射场在极点空间的表示。

### 2.2 第②步：聚合到第2层 (`agg2Level2!`)

从叶层逐层向上聚合到第2层（第1层是大盒子根，无意义）。

对每一层 $l$ 从 $(n_\text{Levels}-1)$ 到 $2$：

$$\text{aggS}_\text{parent}[p, \text{pol}, c] = \sum_{\text{kids}} \text{phaseShiftFromKids}[p, k_\text{in8}] \cdot [\mathbf{M}_\theta(\mathbf{M}_\phi \cdot \text{aggS}_\text{child}[:,:,\text{kid}])]_{[p, \text{pol}]}$$

```
tAggS .= 0                                   ← 父层清零
aggSInterped = θCSC * (ϕCSC * kAggS[:,:,kid]) ← 先ϕ后θ插值
tAggS[:,:,iCube] += phaseShiftFromKids[:,kIn8] .* aggSInterped  ← 相移+累加
```

**工程技巧**：
- 父层 aggS 先清零再累加 → 避免旧数据污染
- 插值顺序固定为 先ϕ后θ（与矩阵结构对应）
- 线程数 > 盒子数时用 BLAS 多线程，否则用 `@threads` 并行

### 2.3 第③步：各层转移 (`transOnLevels!`)

对第2层到叶层逐层执行：

$$\text{disaggG}[p, \text{pol}, c] = \sum_{f \in \text{farneighbors}(c)} \alpha_\text{Trans}[p, \text{idx}(c \leftarrow f)] \cdot \text{aggS}[p, \text{pol}, f]$$

```
disaggG .= 0                                     ← 每层先清零
disaggG[:,:,iCube] += αTrans[:,i1d] .* aggS[:,:,iFarNei]  ← 对远亲累加
```

**其中 `i1d`** 由远亲的相对3D偏移查 `αTransIndex` 表得到。

### 2.4 第④步：解聚到叶层 (`disagg2LeafLevel!`)

从第2层逐层向下到叶层：

$$\text{disaggG}_\text{child}[:,:,\text{kid}] \mathrel{+}= \mathbf{M}_\phi^T \cdot [\mathbf{M}_\theta^T \cdot (\text{phaseShift2Kids}[:,k_\text{in8}] \cdot \text{disaggG}_\text{parent}[:,:,c])]$$

```
disGshifted = phaseShift2Kids[:,kIn8] .* tCubeDisaggG  ← 相移
disGInterped = ϕCSCT * (θCSCT * disGshifted)          ← 先θ后ϕ反插值
kDisAggG[:,:,kCubeID] += disGInterped                  ← 累加到子层
```

> **注意**：反插值顺序是 先θ后ϕ（与正向插值的先ϕ后θ相反）。

> **关键**：子层 `kDisAggG` 在第③步 `transOnLevel!` 中**已经清零**，所以此处 `+=` 不会累积旧数据。

### 2.5 第⑤步：叶层解聚到基函数 (`disaggOnBF!`)

$$Z_\text{far}[n] = jk\eta \cdot \sum_{p=1}^{N_\text{poles}} \left[\text{disaggSBF}[p, 1, n] \cdot \text{disaggG}[p, 1, c] + \text{disaggSBF}[p, 2, n] \cdot \text{disaggG}[p, 2, c]\right]$$

```
ZI[n] = 0
for idx in 1:nPoles
    ZInTemp += disaggSBF[idx, 1, n] * disaggGCube[idx, 1] 
             + disaggSBF[idx, 2, n] * disaggGCube[idx, 2]
end
ZInTemp *= JK_0η     # = jkη ← 关键全局系数!
ZI[n] += ZInTemp
```

### 2.6 合并

$$y_n = \underbrace{(Z_\text{near} \cdot x)_n}_{\text{近场}} + \underbrace{Z_\text{far}[n]}_{\text{匹配 disaggOnBF! 的输出}}$$

---

## 3. 系数链完整验证

### 3.1 EFIE 直接法系数

$$Z_{mn}^\text{EFIE} = \underbrace{\frac{jk\eta}{16\pi}}_{\text{JKηdiv16π}} \cdot l_m \cdot l_n \int\!\!\!\int \left[(\boldsymbol{\rho}_m \cdot \boldsymbol{\rho}_n) - \frac{4}{k^2}\right] \frac{e^{-jkR}}{4\pi R}\; dA_t\, dA_s$$

注意 $G$ 已含 $1/(4\pi R)$，所以总系数 = $\frac{jk\eta}{16\pi} \cdot l_m \cdot l_n$。

### 3.2 MLFMA 远场系数链

$$Z_{mn}^{\text{far}} = \sum_p \left[\underbrace{\text{disaggSBF}_m[p]}_{\text{配置}}\right] \cdot \underbrace{\alpha_\text{Trans}[p]}_{\text{转移}} \cdot \left[\underbrace{\text{aggSBF}_n[p] \cdot I_n}_{\text{辐射}}\right] \cdot \underbrace{jk\eta}_{\text{全局因子}}$$

展开各项：

| 项 | 表达式 | 包含的因子 |
|---|---|---|
| `aggSBF` | $\int (\hat{e} \cdot \boldsymbol{\rho}_n)\; w\; \frac{l_n}{2}\; e^{+jk\hat{r}_p \cdot \mathbf{r}'_\text{local}}\; dA_s$ | $\frac{l_n}{2}$, $w$, $e^{+jk}$ |
| `disaggSBF` | $\overline{\text{aggSBF}_m}$ | $\frac{l_m}{2}$, $w$, $e^{-jk}$ |
| `αTrans` | $\frac{-jk}{16\pi^2} W_p \sum_l (...)$ | $\frac{-jk}{16\pi^2}$, $W_p$ |
| 全局因子 | $jk\eta$ | $jk\eta$ |

合并（忽略相位和基函数积分细节）：

$$Z_{mn}^{\text{far}} \propto \frac{l_m}{2} \cdot \frac{l_n}{2} \cdot \frac{-jk}{16\pi^2} \cdot jk\eta \cdot \underbrace{\sum_p W_p\; T_L\; e^{\pm jk\cdots}}_{\approx \int d^2\hat{k}\; T_L\; e^{\pm jk\cdots}}$$

利用 Addition Theorem：

$$\sum_p W_p\; T_L(p)\; e^{-j\mathbf{k}_p\cdot\mathbf{d}_1}\; e^{+j\mathbf{k}_p\cdot\mathbf{d}_2} \approx \int d^2\hat{k}\; T_L\; e^{-j\mathbf{k}\cdot\mathbf{d}_1}\; e^{+j\mathbf{k}\cdot\mathbf{d}_2} = \frac{(4\pi)^2}{-jk} \cdot G(\mathbf{r}, \mathbf{r}')$$

代入：

$$Z_{mn}^{\text{far}} = \frac{l_m}{2} \cdot \frac{l_n}{2} \cdot \frac{-jk}{16\pi^2} \cdot jk\eta \cdot \frac{16\pi^2}{-jk} \cdot G = \frac{l_m \cdot l_n}{4} \cdot jk\eta \cdot G$$

而直接法：

$$Z_{mn}^{\text{direct}} = \frac{jk\eta}{16\pi} \cdot l_m \cdot l_n \cdot G_{\text{int}} = \frac{l_m \cdot l_n}{4} \cdot jk\eta \cdot \frac{G_{\text{int}}}{4\pi}$$

因 $G = \frac{e^{-jkR}}{4\pi R}$ 已含 $\frac{1}{4\pi}$：

$$Z_{mn}^{\text{far}} = \frac{l_m \cdot l_n}{4} \cdot jk\eta \cdot \frac{e^{-jkR}}{4\pi R} = \frac{jk\eta}{16\pi} \cdot l_m \cdot l_n \cdot \frac{e^{-jkR}}{R}$$

直接法 $G_\text{int} = \frac{e^{-jkR}}{R}$（不含 $4\pi$，因为 $\frac{1}{4\pi}$ 在前系数里）。

$$\Rightarrow Z_{mn}^{\text{far}} = Z_{mn}^{\text{direct}} \quad \checkmark$$

---

## 4. EMSuite 与 Legacy 的系数分拆对比

### 4.1 转移因子 (`αTrans`)

| | Legacy | EMSuite | 比值 |
|---|---|---|---|
| `const_factor` | $\frac{-jk}{(4\pi)^2} = \frac{-jk}{16\pi^2}$ | $\frac{-jk}{4\pi}$ | EMSuite = $4\pi \times$ Legacy |

### 4.2 叶层解聚全局因子

| | Legacy | EMSuite | 比值 |
|---|---|---|---|
| 全局因子 | $jk\eta$ | $\frac{jk\eta}{4\pi}$ (= $4 \times \frac{jk\eta}{16\pi}$) | Legacy = $4\pi \times$ EMSuite |

### 4.3 总系数验证

- Legacy 总 = $\frac{-jk}{16\pi^2} \cdot jk\eta = \frac{k^2\eta}{16\pi^2}$ ✓
- EMSuite 总 = $\frac{-jk}{4\pi} \cdot \frac{jk\eta}{4\pi} = \frac{k^2\eta}{16\pi^2}$ ✓

**总系数一致**，但分拆方式不同。Legacy 把 $4\pi$ 因子放在 Translation 中（从 $\frac{1}{4\pi}$ 变成 $\frac{1}{16\pi^2}$），全局因子用完整的 $jk\eta$。

### 4.4 已确认完全一致的部分

| 组件 | Legacy | EMSuite | 状态 |
|---|---|---|---|
| aggSBF 辐射函数 | $\hat{e}\cdot\boldsymbol{\rho} \cdot w \cdot \frac{l}{2} \cdot e^{+jk}$ | 同左 | ✅ 一致 |
| disaggSBF (EFIE) | $\overline{\text{aggSBF}}$ | 同左 (实时计算 $e^{-jk}$) | ✅ 一致 |
| 相移因子 | $e^{-jk\hat{r}\cdot\Delta\mathbf{r}} / e^{+jk}$ | 同左 | ✅ 一致 |
| 远亲定义 | $7^3 - 3^3 = 316$ | 同左 | ✅ 一致 |
| 截断项数公式 | 同 §0.3 | 同左 | ✅ 一致 |
| 插值方式 | 两步 Lagrange (先ϕ后θ) | 同左 | ✅ 一致 |
| 反插值方式 | 转置 (先θ后ϕ) | 同左 | ✅ 一致 |
| Hankel 函数 | $h_l^{(2)}$ (第二类) | 同左 | ✅ 一致 |
| 总系数 | $k^2\eta / 16\pi^2$ | 同左 | ✅ 一致 |

---

## 5. Legacy 的关键工程技巧

### 5.1 预计算 aggSBF/disaggSBF (空间换时间)

Legacy 预计算 `aggSBF(nPoles, 2, nBF)` 数组，在每次 `mul!` 中只需做一次矩阵内积 `aggS += I[n] * aggSBF[:,:,n]`。

- **优点**：每次迭代只需 $O(N \cdot N_\text{poles})$ 乘法
- **代价**：需存储 $N_\text{poles} \times 2 \times N$ 复数矩阵
- **EMSuite 差异**：EMSuite 在每次 mul! 中重新计算辐射积分(不预计算 aggSBF)

### 5.2 对称相移因子 (8→4)

$z$ 轴以上和以下的子盒子偏移互为负向，其相移因子为共轭关系。只需计算4个，另外4个取共轭。

### 5.3 BLAS 线程 vs `@threads` 自适应

当盒子数 < 线程数时使用 BLAS 多线程做矩阵乘；否则用 `@threads` 并行对盒子循环。

### 5.4 提取公共指数项

在聚合/解聚内循环中，`exp(jk r̂·r_local) * w * l/2` 被提取为公共变量 `expWlntemp`，减少重复计算。

### 5.5 基函数全局重排

构建八叉树时将基函数按空间位置重排 (`kidsSorted`)，使同一盒子内的基函数编号连续 → `bfInterval` 为 `UnitRange`，便于向量化操作。

### 5.6 CSC 稀疏格式近场矩阵

近场矩阵按 CSC 格式存储，利用邻盒对的固定稀疏模式，预计算 `colPtr` 和 `rowIndices`。

### 5.7 分线程预分配临时内存

多线程循环中的临时数组按线程 id 索引 (`view(buffer, :, :, tid)`)，避免动态分配和竞争条件。

---

## 6. Legacy 的 `disaggG` 清零时机（关键！）

### 6.1 内存预分配

Legacy `memoryAllocationOnLevels!` 在八叉树构建完毕后，**对所有层** (2 到 nLevels) 预分配 `aggS` 和 `disaggG`，初始化为零。

### 6.2 每次 mul! 的清零顺序

```
① aggOnBF!:        aggS .= 0              (叶层 aggS 清零)
② agg2HighLevel!:  tAggS .= 0             (每个父层 aggS 在聚合前清零)
③ transOnLevel!:   disaggG .= 0           (每层 disaggG 在转移前清零)
④ disagg2KidLevel!: 不清零 kDisAggG       (子层 disaggG 已在③中被清零，此处只累加)
⑤ disaggOnBF!:     ZI 通过 setzero=true   (远场输出向量清零)
```

**关键点**：步骤③ `transOnLevels!` 对**所有层**（2 到 nLevels）的 `disaggG` 都清零了，包括叶层。然后步骤④ `disagg2KidLevel!` 对子层（叶层）做 `+=` 是安全的，因为叶层的 `disaggG` 已在该次 `transOnLevel!` 中被清零。

### 6.3 EMSuite 的清零时机

EMSuite `translate!` 也有 `fill!(level.disaggG, 0)` → 然后 `disaggregate_downward!` 的 `childDisaggG .+= childVal` 也是安全的。

**但有一个潜在问题**：EMSuite 的 `disaggregate_downward!` 中：
```julia
if !isdefined(childLevel, :disaggG)
    childLevel.disaggG = zeros(CT, ...)
end
# → 第二次 mul! 调用时，childLevel.disaggG 已存在，不会被清零！
# 但这不影响，因为 translate! 已经 fill! 过了。
```

实际上这不是问题，因为 `translate!` 对每层（包括叶层）都做了 `fill!(..., 0)`。只要 `translate!` 被正确执行，`disaggregate_downward!` 的累加就是安全的。

---

## 7. EMSuite 当前问题定位

### 7.1 已验证的事实

| 测试 | 结果 |
|---|---|
| nLevels=2 (叶层=根层) 精度 | 0.5%~9% ✅ |
| nLevels=3 精度 | 1005% ✗ |
| nLevels=4 精度 | 2918% ✗ |
| EFIE aggS == PMCHW aggS | ✅ 完全一致 |
| EFIE disaggG == PMCHW disaggG | ✅ 完全一致 |
| Σ_p αTrans[p] / 4πG | = 1.0000 ✅ (Translation 公式正确) |
| Σ_p W_p | = 4π ≈ 12.5664 ✅ (权重归一化正确) |
| 总系数 | Legacy == EMSuite ✅ |

### 7.2 推断

- Translation 公式和系数本身是正确的
- 总系数分拆虽然不同但等价
- aggSBF / disaggSBF 逻辑一致
- **问题极大概率在 aggregation upward 或 disaggregation downward 的 插值/反插值矩阵计算 或 相移因子计算 中**
- nLevels=2 正确表明叶层聚合、Translation、叶层解聚三步在单层情况下工作正常
- nLevels≥3 错误表明 upward/downward pass 的 插值或相移 存在数值错误

### 7.3 下一步排查方向

1. **对比 EMSuite 与 Legacy 的插值矩阵** (`θCSC`, `ϕCSC`) 的维度和数值
2. **对比相移因子** `phaseShiftFromKids` 的维度和数值
3. **验证 upward pass**：检查 `aggregate_upward!` 后 level 2 的 `aggS` 是否与 Legacy 一致
4. **验证 downward pass**：检查 `disaggregate_downward!` 后 level 3 的 `disaggG` 是否正确传递

---

## 附录 A：Legacy 关键常数表

| 符号 | Legacy 代码 | 定义 |
|---|---|---|
| $k_0$ | `Params.K_0` | $\omega/c_0$ |
| $jk_0$ | `Params.JK_0` | $jk_0$ |
| $\eta_0$ | `η_0` | $\sqrt{\mu_0/\varepsilon_0} \approx 376.73\;\Omega$ |
| $jk\eta$ | `Params.JKη_0 = JK_0 * η_0` | $jk_0\eta_0$ |
| $jk\eta/(16\pi)$ | `Params.JKηdiv16π` | EFIE 矩阵元前系数 |
| $4/k^2$ | `Params.C4divk²` | 散度项系数 |
| $-jk/(4\pi)^2$ | `mjKdiv16π²` | 转移因子前系数 |
| $R_\text{sglr}$ | `Params.Rsglr = 0.15λ` | 奇异性提取阈值 |

## 附录 B：Legacy 关键文件索引

| 文件 | 功能 |
|---|---|
| `MoM_Basics/src/ParametersSet.jl` | 物理/频率常数定义 |
| `MoM_Kernels/src/MLFMA/IterateOnOctree.jl` | mul! 远场核心：aggOnBF!, transOnLevel!, disaggOnBF! |
| `MoM_Kernels/src/MLFMA/AggOnBF/AggEFIE.jl` | EFIE aggSBF/disaggSBF 预计算 |
| `MoM_Kernels/src/MLFMA/AggOnBF/AggMFIE.jl` | MFIE aggSBF/disaggSBF 预计算 |
| `MoM_Kernels/src/MLFMA/AggOnBF/AggCFIE.jl` | CFIE aggSBF/disaggSBF 预计算 |
| `MoM_Kernels/src/MLFMA/PhaseShiftAndTransFactors.jl` | 相移因子 + αTrans 预计算 |
| `MoM_Kernels/src/MLFMA/LevelInfo.jl` | Level/Cube 数据结构 |
| `MoM_Kernels/src/MLFMA/OctreeInfo.jl` | 八叉树构建 |
| `MoM_Kernels/src/MLFMA/IntegralInterpolationInfo.jl` | 截断项数 + 采样点 |
| `MoM_Kernels/src/MLFMA/LagrangeInterpolation.jl` | Lagrange 插值矩阵计算 |
| `MoM_Kernels/src/MLFMA/MLFMAIterators.jl` | MLFMAIterator 构造 + mul! 封装 |
| `MoM_Kernels/src/Zmat/EFIE/EFIERWGTri.jl` | EFIE RWG 矩阵元直接计算 |
| `MoM_Kernels/src/Solvers.jl` | GMRES/BiCGSTAB 迭代求解器 |
| `MoM_Kernels/src/FastAlgorithm.jl` | MLFMA 管线构建入口 |
| `MoM_Basics/src/MeshProcess/MeshProcess.jl` | Nastran .nas 网格读取 |
| `MoM_Basics/src/BasisFunctions/RWG.jl` | RWG 基函数构造 |

---

## 8. 迭代求解器 (`Solvers.jl`)

### 8.1 求解器选择

```julia
solve(A::LinearMapType, b::Vector; solverT=:gmres, rtol=1e-3, 
      restart=200, Pl=Identity(), Pr=Identity())
```

| 求解器 | 选择符号 | 参数 |
|---|---|---|
| GMRES | `:gmres` (默认) | restart=200 |
| BiCGSTAB | `:bicgstab` | — |

### 8.2 收敛准则

$$\frac{\|\mathbf{r}_k\|}{\|\mathbf{r}_0\|} \leq \text{rtol}$$

默认 `rtol = 1e-3`。支持左/右预条件子 (`Pl`, `Pr`)。

### 8.3 输出

每次迭代输出残差范数到收敛历史文件。

---

## 9. 后处理：远场与 RCS

### 9.1 远场电场

$$\mathbf{E}_\text{far}(\hat{r}) = \frac{-jk_0\eta_0}{4\pi} \sum_{n} I_n \int_{S_n} \mathbf{J}(\mathbf{r}') e^{+jk_0\hat{r}\cdot\mathbf{r}'}\; dS'$$

其中 $\mathbf{J}(\mathbf{r}') = \sum_n I_n \mathbf{f}_n(\mathbf{r}')$ 是面电流密度。

远场积分的核 $e^{+jk\hat{r}\cdot\mathbf{r}'}$ 与 MLFMA 聚合中的 $e^{+jk\hat{r}\cdot\mathbf{r}'}$ 形式一致。

### 9.2 RCS (雷达散射截面)

$$\sigma(\theta, \phi) = \frac{(k_0\eta_0)^2}{4\pi} |\mathbf{N}(\theta, \phi)|^2$$

$$\text{RCS}_\text{dB} = 10 \log_{10} \sigma$$

其中 $\mathbf{N}(\theta, \phi)$ 是远场辐射积分（不含 $-jk\eta/(4\pi)$ 前系数）。

> **注意**：RCS 用功率维度 $10\log_{10}$，远场用场维度 $20\log_{10}$。

---

## 10. MLFMA 管线构建入口 (`FastAlgorithm.jl`)

### 10.1 构建顺序

```julia
function getImpedanceOpt(geosInfo, bfsInfo)
    # Step 1: 八叉树构建 + 基函数重排
    getOctreeAndReOrderBFs!(geosInfo, bfsInfo)
    #  ├─ 提取基函数中心 bfcenters[3, nbf]
    #  ├─ 计算包围盒、nLevels
    #  ├─ 递归构建层级 (nLevels → 1)
    #  ├─ 基函数重编号 (按空间位置)
    #  ├─ 预计算插值矩阵 (θ, ϕ, 6阶Lagrange)
    #  ├─ 预计算相移因子
    #  └─ 预计算转移因子

    # Step 2: 近场矩阵组装
    Znear = calZnearCSC(leafLevel, geosInfo, bfsInfo)
    #  └─ CSC 稀疏矩阵 [nbf × nbf]

    # Step 3: 构造 MLFMA 迭代算子
    Zopt = MLFMAIterator(Znear, octree, geosInfo, bfsInfo)
    #  ├─ 预计算 aggSBF[nPoles, 2, nbf]
    #  ├─ 预计算 disaggSBF[nPoles, 2, nbf]
    #  └─ 预分配 ZI[nbf] 工作空间
end
```

### 10.2 与 EMSuite 的结构差异

| 步骤 | Legacy | EMSuite |
|---|---|---|
| aggSBF/disaggSBF | 预计算并存储 | 每次 mul! 重新计算 |
| 内存分配 | 统一在构建时完成 | 惰性分配 (`isdefined` 检查) |
| 基函数重排 | 全局重排 (`kidsSorted`) | `sorted_ids` 映射 |
| 近场矩阵 | 在 `FastAlgorithm.jl` 中组装 | 在 `MLFMABuilder.jl` 中组装 |
