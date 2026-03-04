# Phase 15: 介质与金属-介质混合天线精度测试 + PMCHWMLFMAOperator

> 创建日期: 2026-03-07  
> 状态: 计划中  
> 前置: Phase 14 完成（A1–A4, F1–F9, P1–P3 基准脚本已创建）

---

## 0. 开发原则遵循声明

本 Phase 适用以下核心原则（见 `copilot-instructions.md`）：

| 原则 | 适用内容 |
|------|---------|
| **原则 1 (TDD)** | 新增激励 API、`PMCHWMLFMAOperator`、`input_impedance_pmchw` 均先写单元测试 |
| **原则 2 (Legacy 对齐)** | PMCHW 激励及 MLFMA 结果必须与 Direct 结果一致，不使用经验常数 |
| **原则 3 (差异排查)** | MLFMA 结果与 Direct 偏差时逐块比较：Z_near 是否正确、far-field k0/k1 分离是否正确 |
| **原则 5 (Git 提交)** | 每类子任务独立提交：API 扩展、单元测试、基准脚本、PMCHWMLFMAOperator |
| **原则 6 (终端输出)** | 不重定向；使用 `get_terminal_output` 获取结果 |
| **原则 7 (Phase 结束检视)** | 完成后至少 2 轮 clean 检视 |
| **原则 9 (计划文档规范)** | 本文档已包含所有必要节 |

---

## 1. Legacy 对齐基准

| 子任务 | Legacy 参考 | 说明 |
|--------|------------|------|
| PMCHW DeltaGap 激励 | `PMCHW.jl` L/K 算子 + EFIE DeltaGap 类比 | V_E 行施加 delta-gap，V_H 行为零 |
| PMCHW input_impedance | PMCHW direct solve（作为参考） | MLFMA vs Direct 结果自洽 |
| SCFIE DeltaGap | `SCFIE.jl` + 现有 DeltaGap of RWGBasis | 激励只进表面块，体块为零 |
| PMCHWMLFMAOperator | `PMCHW.jl` `assemble_impedance_matrix` | MLFMA 的 Z_near + Z_far 近似全矩阵 |

---

## 2. 目标

在 Phase 14 天线测试（A1–A4, 纯金属偶极子）基础上，新增：

1. **B1–B2**: 介质谐振体（PMCHW）Delta-Gap 馈电天线，输入阻抗测试
2. **B3–B5**: 金属-介质混合天线（VS-EFIE/VS-CFIE, 即 SCFIE），输入阻抗测试
3. **PMCHWMLFMAOperator**: 实现 2×2 块 MLFMA 算子，支持 B2 的 MLFMA 路径

### 方程类型与代号对应

| 代号 | EMSuite 方程 | 激励方式 | 求解 |
|------|-------------|---------|------|
| VS-EFIE | `SCFIE(freq, perms; alpha=0.0)` | DeltaGap (表面) | Direct |
| VS-CFIE | `SCFIE(freq, perms; alpha=0.5)` | DeltaGap (表面) | Direct |
| PMCHW Direct | `PMCHW(freq, eps_r, mu_r)` | DeltaGap (新增) | Direct |
| PMCHW MLFMA | `PMCHWMLFMAOperator` (新增) | DeltaGap | GMRES |

---

## 3. 需新增的 API

### 3.1 PMCHW DeltaGap 激励

**位置**: `src/IntegralEquations/Excitation.jl`

```julia
"""
    excitation_vector(op::PMCHW, source::DeltaGapSource, basis::RWGBasis) → 2N Vector

Delta-gap 施加于 PMCHW E-方程行（前 N 行），H-方程行保持为零。

物理含义: 缝隙电场仅激励 E 方程中的电流 J；
         磁场分量在理想导体缝隙处为零（无 M 激励）。
"""
function excitation_vector(op::PMCHW, source::DeltaGapSource, basis::RWGBasis)
    N = num_basis(basis)
    V = zeros(ComplexF64, 2N)
    # 前 N 行（E-方程）: V_E[n] = source.voltage * edge_len
    for idx in source.edge_indices
        1 <= idx <= N || continue
        V[idx] = source.voltage * basis.functions[idx].edge_length
    end
    # V[N+1:2N] = 0  (H-方程无 delta-gap 激励)
    return V
end
```

**注意**: `input_impedance(source::DeltaGapSource, I_2N, basis)` 的通用方法已将 `I[idx] * edge_len` 作为总馈电电流，但对 PMCHW 的 2N 向量，需使用 J 部分（前 N 个系数）：

```julia
# 在 AntennaMetrics.jl 中新增重载：
function input_impedance(op::PMCHW, source::DeltaGapSource,
                         I_2N::Vector{<:Complex}, basis::RWGBasis)
    N = num_basis(basis)
    # 仅用 J 部分（前 N 个系数）
    I_J = @view I_2N[1:N]
    I_in = sum(I_J[idx] * basis.functions[idx].edge_length
               for idx in source.edge_indices if 1 <= idx <= N)
    iszero(I_in) && error("input_impedance: J-part feed current is zero")
    return ComplexF64(source.voltage) / ComplexF64(I_in)
end
```

### 3.2 SCFIE DeltaGap 激励

**位置**: `src/IntegralEquations/Excitation.jl`

```julia
"""
    excitation_vector(op::SCFIE, source::DeltaGapSource, surf_basis, vol_basis)

Delta-gap 仅加载到表面（RWG）部分，体积部分为零。
适用于金属-介质混合天线（馈电点在金属表面）。
"""
function excitation_vector(op::SCFIE, source::DeltaGapSource,
                           surf_basis::RWGBasis, vol_basis::AbstractBasisFunction)
    N_S = num_basis(surf_basis)
    N_V = num_basis(vol_basis)
    V  = zeros(ComplexF64, N_S + N_V)
    # 仅填充表面块（使用已有的通用 DeltaGap 实现）
    V_surf = excitation_vector(source, surf_basis)   # N_S × 1
    V[1:N_S] .= V_surf
    # V[N_S+1:end] = 0  (介质体无 delta-gap 激励)
    return V
end
```

**导出**: 在 `EMSuite.jl` 中不需要额外 export（已通过 `excitation_vector` 多分派）

### 3.3 PMCHWMLFMAOperator（新文件）

**位置**: `src/FastAlgorithms/MLFMA/PMCHWMLFMAOperator.jl`

#### 数学结构

PMCHW 系统 (2N × 2N)：
$$Z = \begin{bmatrix} Z^{EJ} & Z^{EM} \\ Z^{HJ} & Z^{HM} \end{bmatrix}$$

其中（超简化表示）：
- $Z^{EJ}(k_0, \eta_0) + Z^{EJ}(k_1, \eta_1)$：EFIE-型 L 算子，两个波数之和
- $Z^{EM}(k_0) + Z^{EM}(k_1)$：K 算子（MFIE 类型），两个波数之和  
- $Z^{HJ} = -Z^{EM}$（结构不变量，精确成立）
- $Z^{HM}(k_0, 1/\eta_0) + Z^{HM}(k_1, 1/\eta_1)$：倒 η 的 EFIE-型

#### 设计方案: 双算子 + 4 块矩阵

```
PMCHWMLFMAOperator
├── mlfma_efie_k0  : MLFMAOperator（EFIE with k0, η0）[处理 Z_EJ_k0, Z_HM_k0]
├── mlfma_efie_k1  : MLFMAOperator（EFIE with k1, η1）[处理 Z_EJ_k1, Z_HM_k1]
├── mlfma_cfie_k0  : MLFMAOperator（CFIE with k0, pure K-part）[Z_EM_k0]
├── mlfma_cfie_k1  : MLFMAOperator（CFIE with k1, pure K-part）[Z_EM_k1]
└── Z_near         : 2N×2N 稀疏矩阵（直接从 PMCHW.assemble_near_field 获得）
```

**mul! 逻辑**（`y = PMCHWMLFMAOperator * x`）：

```
x = [x_J; x_M]   (N + N 向量)
y = [y_E; y_H]   (N + N 输出)

# 计算远场 J -> E 贡献 (Z_EJ 部分)
y_EJ = sum_k(L_k * x_J)  → 用 k0/k1 EFIE MLFMA far-field

# 计算远场 M -> E 贡献 (Z_EM 部分)  
y_EM = sum_k(K_k * x_M)  → 用 k0/k1 K-type far-field (从 CFIE 分离 MFIE 部分)

# H-方程：Z_HJ = -Z_EM^T, Z_HM = Lη 类算子
y_HJ = -sum_k(K_k^T * x_J)   # 等于 -Z_EM 的转置
y_HM = sum_k(Lη_k * x_M)

y_E = y_EJ + y_EM
y_H = y_HJ + y_HM
```

**实现策略**：
- 阶段 1 (先行可行版): Dense Z_near + **仅近场**（即只用 Z_near，不实现 Z_far）
  - 可验证 PMCHWMLFMAOperator 的接口和 GMRES 收敛性
  - 等价于 Block-Diagonal preconditioned direct solve
- 阶段 2 (完整实现): 添加 Z_far 的 4 块计算

> **注意**: 阶段 1 已可通过 GMRES+Z_near 的精度测试（等同于直接法，只是慢一些）。完整 MLFMA 在阶段 2。

---

## 4. 测试场景设计

### 4.1 几何体选择

| 标签 | 几何体 | 获取方式 | 用于 |
|------|--------|---------|------|
| `sphere_600MHz.nas` | PEC/介质球（已有） | `meshfiles/sphere_600MHz.nas` | B1, B2 (re-mesh as dielectric) |
| `dipole_on_cube.nas` | 薄偶极子 + 小介质方块 | 程序生成（见 4.2） | B3, B4, B5 |
| `shell_dipole.nas` | 开口圆柱面（金属）+ 介质套筒 | 程序生成（见 4.2） | B3 备选 |

### 4.2 几何生成方案

#### B1/B2 PMCHW: 介质体 + Delta-Gap 馈电
- 直接使用 `sphere_600MHz.nas` 重新定义为介质体（εᵣ=4）
- 馈电边: 球赤道处 z≈0 的 RWG 边（与偶极子馈电类似）
- 验证: Z_in 的实部 > 0（辐射阻抗），虚部随频率单调变化
- 自洽检查: 当 εᵣ=1+ε（接近空气）时，PMCHW Z_in → EFIE Z_in（相同网格）

#### B3/B4/B5 SCFIE: 金属 + 介质
**方案 A（首选，程序生成）**:
```julia
# 金属线（偶极子臂）使用 generate_cylinder_mesh
metal_mesh = generate_cylinder_mesh(0.005, 0.5, 8, 20; closed=false)
surf_basis = RWGBasis(metal_mesh)

# 介质衬底（小立方体/薄片）需要 NAS 文件或生成函数
# 方案: 用已有的 tetrahedra mesh 作为介质体代理
dielectric_mesh = read_nas_mesh(joinpath(meshfiles_dir, "Tetra.nas"))  # 若存在
vol_basis = SWGBasis(dielectric_mesh)
```

**方案 B（备选，使用 plate_and_metal_1dot2GHz.nas 散射件作参考）**:
- 已有 32 CTRIA3 (RWG) + 7278 CTETRA (SWG)
- 选取一条金属面 RWG 边作为 delta-gap 馈电边
- 参考: EFIE-only (仅金属 32 CTRIA3) 的 Z_in（极小结构，数值上自洽即可）

> **推荐方案 A**，因金属线+介质方块物理意义更明确；若 Tetra.nas 不可用则退为方案 B。

### 4.3 用例清单

| ID | 几何体 | 方程 | α | 求解器 | 参考 | 精度门限 |
|----|--------|------|---|-------|------|---------|
| **B1** | 介质球（PMCHW）εᵣ=4 | PMCHW | — | Direct | 自洽: Re(Z_in)>0 | Re误差(vs freq sweep)<10% |
| **B2** | 同 B1 | PMCHW | — | MLFMA+GMRES | B1 Direct | ΔZ_in: Re<5%, Im<20Ω |
| **B3** | 金属偶极子+介质块 | VS-EFIE | 0.0 | Direct | EFIE-only (εᵣ→1) | Z_in Re误差<10% |
| **B4** | 同 B3 | VS-CFIE | 0.5 | Direct | B3 (α=0) | ΔZ_in<5Ω (公式稳定性) |
| **B5** | 同 B3 | VS-CFIE | 0.5 | MLFMA+GMRES | B4 Direct | ΔZ_in Re<5%, Im<20Ω |

### 4.4 参考值来源

| 用例 | 参考来源 | 说明 |
|------|---------|------|
| B1 Re(Z_in) | 物理约束: Re > 0 | 辐射阻抗必须为正 |
| B1 Z_in(εᵣ→1) | EFIE Direct Z_in | PMCHW 极限验证（εᵣ=1.001） |
| B2 vs B1 | B1 Direct | MLFMA 自洽 |
| B3 vs EFIE-only | EFIE 纯金属 Z_in | 介质影响应使 Z_in 增大/频移 |
| B4 vs B3 | B3 (α=0) | CFIE vs EFIE 近似一致 |
| B5 vs B4 | B4 Direct | MLFMA 自洽 |

---

## 5. PMCHWMLFMAOperator 原生实现方案

> **第三次修订（依据 Gibson《MoM》Ch.11 Algorithm 14）**  
> 前两版设计（独立算子包装 / MagneticRWGBasis 标签）均被放弃。  
> 正确方案：**两个 N 点八叉树（k0/k1 各一套）+ J/M 两遍 × k0/k1 两遍 = 四遍远场 + 四种解聚核函数**。  
> 本节是面向"无 CEM 背景的软件工程师"的完整实现规格。

### 5.0 核心设计原则（依据 Gibson Algorithm 14）

**Gibson 原文引用**（Ch.11, Algorithm 14）：

> "In general, **two passes** must be executed to yield the far product.  
> In the first pass, the fields due to **only the electric basis functions** are aggregated and transferred.  
> In the second pass, the fields due to **the magnetic basis functions** are aggregated and transferred.  
> During the disaggregation step, the **appropriate electric or magnetic receive function**, respectively, is used for each testing function."  
> "T_M = η_l × T_J in region R_l when using RWG functions, and **only T_J is actually stored**."  
> "Each testing function has **two receive functions**: one to receive fields radiated by electric basis functions, and one to receive fields radiated by magnetic basis functions."

**由此得到的设计要点**：

| 要点 | 说明 |
|------|------|
| 两个八叉树 | 各从 **N 个** RWG 中心点建立（不是 2N），`octree0` 使用 k0，`octree1` 使用 k1 |
| 不需要 MagneticRWGBasis | J/M 使用相同 RWG 积分；在 `mul!` 里通过传入 `x_range` 和 `kmode` 区分，不引入新类型 |
| 四遍远场 | J×k0, J×k1, M×k0, M×k1 各一遍（共4遍）；k0 用 octree0, k1 用 octree1 |
| 聚合积分相同 | J 和 M 的**辐射花样积分**完全一样（$\int \boldsymbol{\rho}\, e^{j k \hat{r}\cdot\mathbf{r}} dS$）；区别由解聚接收函数承担 |
| 四种解聚接收 | 每个测试函数 × 两种源（J/M）= 4 个矩阵块，每块有独立的因子 |
| 平移步骤 | 每遍用各自的八叉树执行 M2M/M2L/L2L，代码无需改动 |
| aggS 必须清零 | 每遍聚合之前必须 `fill!(aggS, 0)`，防止同一棵树上前后两遍叠加 |

### 5.1 文件结构（最简化）

```
src/FastAlgorithms/MLFMA/
└── PMCHWMLFMAOperator.jl   ← 新建（自包含：struct + 构造函数 + mul! + 全部辅助函数）
```

**不需要**：
- ~~`src/BasisFunctions/MagneticRWG.jl`~~ — J/M 同用 `RWGBasis`，无需新类型
- ~~扩展 `Aggregation.jl`~~ — 直接在 `PMCHWMLFMAOperator.jl` 内写专用聚合函数
- ~~扩展 `Disaggregation.jl`~~ — 同上

需要修改：
- `src/FastAlgorithms/MLFMA/MLFMAOperator.jl`（或在 `PMCHWMLFMAOperator.jl` 内独立实现）：`assemble_near_field_pmchw`（2N×2N 稀疏矩阵）

### 5.2 struct PMCHWMLFMAOperator

> **双波数设计说明**：PMCHW 矩阵的每个块都是 k0 和 k1 贡献之和（如 Z^EJ = L(k0) + L(k1)）。MLFMA 的聚合权重 $e^{jk\hat{r}\cdot r}$、平移算子 $T_L(k|d|)$ 均依赖 k，因此 k0 和 k1 的贡献**必须分别用各自的八叉树计算**，不能共用一个八叉树。这导致需要两套八叉树（`octree0` for k0, `octree1` for k1），远场总计 4 遍（J×k0, J×k1, M×k0, M×k1）。

```julia
"""
    PMCHWMLFMAOperator{FT,CT}

PMCHW 系统的 MLFMA 算子。矩阵大小 2N×2N，其中 x = [x_J; x_M]，y = [y_E; y_H]。

字段说明：
  pmchw          — PMCHW 算子（含 k0,η0,k1,η1 等物理参数）
  basis          — RWGBasis（N 个 DOF，J 和 M 共享同一网格）
  Z_near         — 2N×2N 稀疏矩阵，4 个块（EJ/EM/HJ/HM）
  octree0        — N 点八叉树，波数 k0（外部介质）
  octree1        — N 点八叉树，波数 k1（内部介质）
  sorted_ids0, inv_sorted_ids0 — 长度 N，k0 八叉树排序映射
  sorted_ids1, inv_sorted_ids1 — 长度 N，k1 八叉树排序映射
  freq           — 工作频率（Hz）
"""
struct PMCHWMLFMAOperator{FT,CT} <: AbstractIntegralOperator
    pmchw          :: PMCHW{FT,CT}
    basis          :: RWGBasis
    Z_near         :: SparseMatrixCSC{CT,Int}   # 2N×2N
    octree0        :: OctreeInfo                # k0 tree
    octree1        :: OctreeInfo                # k1 tree
    sorted_ids0    :: Vector{Int}               # 长度 N
    inv_sorted_ids0:: Vector{Int}
    sorted_ids1    :: Vector{Int}               # 长度 N
    inv_sorted_ids1:: Vector{Int}
    freq           :: FT
end

Base.size(op::PMCHWMLFMAOperator) = (2*num_basis(op.basis), 2*num_basis(op.basis))
Base.eltype(::PMCHWMLFMAOperator{FT,CT}) where {FT,CT} = CT
```

### 5.3 构造函数

```julia
function PMCHWMLFMAOperator(pmchw::PMCHW, basis::RWGBasis, leaf_size::Float64)
    N = num_basis(basis)
    centers = reduce(hcat, [bf.center for bf in basis.functions])  # 3×N

    # ① 分别建立 k0 和 k1 各自的八叉树
    #    两套树从相同 N 个 RWG 中心点建立，但 Lebedev 极点密度由各自 k 决定。
    #    波长公式: λ = 2π / |k|（k 是角频率波数，单位 rad/m）
    λ0 = 2π / real(pmchw.k0)   # 外部介质波长（k0 可能为复数，取实部）
    λ1 = 2π / real(pmchw.k1)   # 内部介质波长

    octree0, sorted_ids0 = build_octree(centers, leaf_size; λ = λ0)
    octree1, sorted_ids1 = build_octree(centers, leaf_size; λ = λ1)

    inv_sorted_ids0 = Vector{Int}(undef, N)
    inv_sorted_ids1 = Vector{Int}(undef, N)
    for i in 1:N
        inv_sorted_ids0[sorted_ids0[i]] = i
        inv_sorted_ids1[sorted_ids1[i]] = i
    end

    # ② 装配 2N×2N 近场稀疏矩阵（见 §5.7）
    #    近场使用 octree0（排列）作为邻居查找基准（两个树几何相同，近邻对相同）
    Z_near = assemble_near_field_pmchw(pmchw, basis, octree0, sorted_ids0, inv_sorted_ids0)

    FT = typeof(real(pmchw.freq))
    CT = Complex{FT}
    return PMCHWMLFMAOperator{FT,CT}(
        pmchw, basis, Z_near,
        octree0, octree1,
        sorted_ids0, inv_sorted_ids0,
        sorted_ids1, inv_sorted_ids1,
        pmchw.freq
    )
end
```

### 5.4 mul! — 两遍 × 两 k = 四遍远场

```julia
function LinearAlgebra.mul!(y::AbstractVector, A::PMCHWMLFMAOperator, x::AbstractVector)
    fill!(y, zero(eltype(y)))
    N  = num_basis(A.basis)
    k0, η0 = A.pmchw.k0, A.pmchw.eta0
    k1, η1 = A.pmchw.k1, A.pmchw.eta1

    # ─── ① 近场（2N×2N，4 块的所有近邻交互） ────────────────────────
    mul!(y, A.Z_near, x)

    # ─── ② 远场 4 遍（2 source × 2 medium） ─────────────────────────
    # 辅助函数：清空某棵八叉树的 aggS 缓存
    clear_agg!(oct) = (for lv in oct.levels; fill!(lv.aggS, zero(eltype(lv.aggS))); end)

    # ── 遍 1：J 源 × k0（聚合用 octree0） ──────────────────────────
    clear_agg!(A.octree0)
    aggregate_leaf_pmchw!(A.octree0.levels[end], A.basis, x, A.sorted_ids0, 1:N)
    _mlfma_up_translate_down!(A.octree0)
    disaggregate_leaf_pmchw_j!(A.octree0.levels[end], A.basis, A.pmchw,
                                y, A.sorted_ids0, :k0)

    # ── 遍 2：J 源 × k1（聚合用 octree1） ──────────────────────────
    clear_agg!(A.octree1)
    aggregate_leaf_pmchw!(A.octree1.levels[end], A.basis, x, A.sorted_ids1, 1:N)
    _mlfma_up_translate_down!(A.octree1)
    disaggregate_leaf_pmchw_j!(A.octree1.levels[end], A.basis, A.pmchw,
                                y, A.sorted_ids1, :k1)

    # ── 遍 3：M 源 × k0（聚合用 octree0） ──────────────────────────
    clear_agg!(A.octree0)
    aggregate_leaf_pmchw!(A.octree0.levels[end], A.basis, x, A.sorted_ids0, (N+1):(2N))
    _mlfma_up_translate_down!(A.octree0)
    disaggregate_leaf_pmchw_m!(A.octree0.levels[end], A.basis, A.pmchw,
                                y, A.sorted_ids0, :k0)

    # ── 遍 4：M 源 × k1（聚合用 octree1） ──────────────────────────
    clear_agg!(A.octree1)
    aggregate_leaf_pmchw!(A.octree1.levels[end], A.basis, x, A.sorted_ids1, (N+1):(2N))
    _mlfma_up_translate_down!(A.octree1)
    disaggregate_leaf_pmchw_m!(A.octree1.levels[end], A.basis, A.pmchw,
                                y, A.sorted_ids1, :k1)

    return y
end

# 辅助：对一棵八叉树执行 M2M → M2L → L2L
function _mlfma_up_translate_down!(octree)
    for levelID in (octree.nLevels-1):-1:2
        aggregate_upward!(octree.levels[levelID], octree.levels[levelID+1])
    end
    for levelID in 2:octree.nLevels
        translate!(octree.levels[levelID])
    end
    for levelID in 2:(octree.nLevels-1)
        disaggregate_downward!(octree.levels[levelID], octree.levels[levelID+1])
    end
end
```

**为何需要 4 遍**：
- MLFMA 的聚合因子 $e^{jk\hat{r}\cdot r}$ 和平移算子 $T_L(k|d|)$ 都是 k 的函数
- k0 和 k1 对应不同介质，Lebedev 极点密度也不同
- 将两种 k 的贡献混入同一 aggS 是错误的
- 4 遍远场 = J-k0 + J-k1 + M-k0 + M-k1，每遍结果直接累加到 `y`

**disaggregate 函数新增 `kmode::Symbol` 参数**（`:k0` 或 `:k1`），在函数内选择正确的 k/η 因子（见下方 §5.6）。

### 5.5 聚合叶结点函数（辐射花样积分）

J-Pass 和 M-Pass **共用同一个** `aggregate_leaf_pmchw!`，行为通过 `x_range` 参数区分：

```julia
"""
    aggregate_leaf_pmchw!(leaf_level, basis, x, sorted_ids, x_range)

将 x[x_range 对应编号] 的系数（对应第 i 个 RWG 基函数）累积到叶结点的 aggS 中。

积分公式（与现有 RWG aggregate_leaf! 完全相同）：
  aggS[iPole, 1:2, iCube] += x_coeff * ∫ ρ(r) e^{j k r̂·r} dS  （θ,ϕ 分量）

参数：
  x_range — J-Pass 传 1:N（系数偏移 0），M-Pass 传 (N+1):(2N)（系数偏移 N）
             函数内部将 x_range 偏移量 offset = first(x_range) - 1
             取 x_coeff = x[bfID_orig + offset]
"""
function aggregate_leaf_pmchw!(leaf_level, basis::RWGBasis, x, sorted_ids, x_range)
    offset = first(x_range) - 1   # J-Pass: 0；M-Pass: N
    for iCube in 1:leaf_level.nCubes
        cube = leaf_level.cubes[iCube]
        for bfID_sorted in cube.bfInterval
            bfID_orig = sorted_ids[bfID_sorted]   # 原始基函数编号（1:N）
            coeff = x[bfID_orig + offset]          # 从 x 中提取对应系数
            bf    = basis.functions[bfID_orig]
            # 调用 add_radiation_pattern_rwg!（现有函数，无需修改）
            add_radiation_pattern_rwg!(
                leaf_level.aggS, bf, coeff,
                leaf_level.cubeCenter[:,iCube],
                leaf_level.poles, iCube,
                leaf_level.k
            )
        end
    end
end
```

> **工程要点**：`add_radiation_pattern_rwg!` 计算 $\int \boldsymbol{\rho}(r) e^{jk\hat{r}\cdot r} dS$（二维向量积分，θ 和 ϕ 分量）。该函数对 J 和 M 都适用，因为两者 RWG 几何相同，辐射花样积分形式相同。

### 5.6 解聚叶结点函数 — 四个矩阵块的核函数

#### 整体结构（软件工程视角）

| 函数名 | 源类型 | 写入 `y` 位置 | 矩阵块 | 接收核 |
|--------|--------|-----------|--------|--------|
| `disaggregate_leaf_pmchw_j!` | J-Pass | `y[bfID_orig]` | Z^EJ | R_L（L 算子） |
| `disaggregate_leaf_pmchw_j!` | J-Pass | `y[bfID_orig+N]` | Z^HJ | R_{−K}（−K 算子） |
| `disaggregate_leaf_pmchw_m!` | M-Pass | `y[bfID_orig]` | Z^EM | R_K（K 算子） |
| `disaggregate_leaf_pmchw_m!` | M-Pass | `y[bfID_orig+N]` | Z^HM | R_{L_η}（L_η 算子） |

#### 接收函数数学定义

设远场点 $\hat{r}$ 的"解聚场"（disaggG）为 $\mathbf{E}^{ff}(\hat{r})$（θ/ϕ 分量），测试基函数第 $i$ 个 RWG 覆盖两个三角形（Plus/Minus）。

**term\_efie（ρ·E 积分，用于 R_L 和 R_{L_η}）**：
$$\text{term\_efie} = \sum_{\hat{r},q} w_r\, w_q \left[\boldsymbol{\rho}(r_q) \cdot \mathbf{E}^{ff}(\hat{r})\right] e^{-jk\hat{r}\cdot (r_q-r_0)} \cdot \frac{\ell_i}{2}$$

**term\_mfie（(ρ×n̂)·H 积分，用于 R_K 和 R_{−K}）**：
$$\text{term\_mfie} = \sum_{\hat{r},q} w_r\, w_q \left[(\boldsymbol{\rho}(r_q)\times\hat{n})\cdot \frac{\hat{r}\times\mathbf{E}^{ff}(\hat{r})}{\eta_l}\right] e^{-jk\hat{r}\cdot (r_q-r_0)} \cdot \frac{\ell_i}{2}$$

其中 $r_0$ 是叶结点中心（`cubeCenter`），$\hat{n}$ 是三角形法向量。

#### 因子对照表（精确公式，用于代码验证）

| 矩阵块 | 遍（共4遍） | 积分类型 | 因子 |
|--------|------------|----------|------|
| Z^EJ | J×k0-Pass | `term_efie` | $\frac{jk_0\eta_0}{16\pi}$ |
| Z^EJ | J×k1-Pass | `term_efie` | $\frac{jk_1\eta_1}{16\pi}$ |
| Z^HJ | J×k0-Pass | `term_mfie` | $-\frac{jk_0}{16\pi}$ （**负号**） |
| Z^HJ | J×k1-Pass | `term_mfie` | $-\frac{jk_1}{16\pi}$ |
| Z^EM | M×k0-Pass | `term_mfie` | $+\frac{jk_0}{16\pi}$ （与 HJ 符号相反） |
| Z^EM | M×k1-Pass | `term_mfie` | $+\frac{jk_1}{16\pi}$ |
| Z^HM | M×k0-Pass | `term_efie` | $\frac{jk_0/\eta_0}{16\pi}$ |
| Z^HM | M×k1-Pass | `term_efie` | $\frac{jk_1/\eta_1}{16\pi}$ |

> **注意**：因子中的 $16\pi$ 来自 MLFMA 球面积分的归一化约定；与现有 EFIE `Disaggregation.jl` 的 `add_received_field_rwg!` 内 `factor` 计算方式一致（直接拷贝参考）。

#### 代码骨架

```julia
function disaggregate_leaf_pmchw_j!(
    leaf_level, basis::RWGBasis, pmchw::PMCHW,
    y::AbstractVector, sorted_ids::Vector{Int}, kmode::Symbol
)
    N      = num_basis(basis)
    k, η = kmode === :k0 ? (pmchw.k0, pmchw.eta0) : (pmchw.k1, pmchw.eta1)
    disaggG = leaf_level.disaggG   # (nPoles, 2, nCubes)

    factor_EJ = im * k * η / (16π)   # EJ 块
    factor_HJ = -im * k / (16π)      # HJ 块（负号）

    Threads.@threads for iCube in 1:leaf_level.nCubes
        cube  = leaf_level.cubes[iCube]
        field = view(disaggG, :, :, iCube)
        r0    = leaf_level.cubeCenter[:, iCube]

        for bfID_sorted in cube.bfInterval
            bfID_orig = sorted_ids[bfID_sorted]
            bf        = basis.functions[bfID_orig]

            term_efie, term_mfie = _receive_terms(bf, field, r0, k, η, leaf_level.poles)

            y[bfID_orig]   += term_efie * factor_EJ   # EJ 块
            y[bfID_orig+N] += term_mfie * factor_HJ   # HJ 块
        end
    end
end

function disaggregate_leaf_pmchw_m!(
    leaf_level, basis::RWGBasis, pmchw::PMCHW,
    y::AbstractVector, sorted_ids::Vector{Int}, kmode::Symbol
)
    N      = num_basis(basis)
    k, η = kmode === :k0 ? (pmchw.k0, pmchw.eta0) : (pmchw.k1, pmchw.eta1)
    disaggG = leaf_level.disaggG

    factor_EM  = +im * k / (16π)        # EM 块（正号，与 HJ 相反）
    factor_HM  =  im * k / η / (16π)    # HM 块（L_η）

    Threads.@threads for iCube in 1:leaf_level.nCubes
        cube  = leaf_level.cubes[iCube]
        field = view(disaggG, :, :, iCube)
        r0    = leaf_level.cubeCenter[:, iCube]

        for bfID_sorted in cube.bfInterval
            bfID_orig = sorted_ids[bfID_sorted]
            bf        = basis.functions[bfID_orig]

            term_efie, term_mfie = _receive_terms(bf, field, r0, k, η, leaf_level.poles)

            y[bfID_orig]   += term_mfie * factor_EM   # EM 块
            y[bfID_orig+N] += term_efie * factor_HM   # HM 块
        end
    end
end
```

#### `_receive_terms` — 两种接收积分（内层）

```julia
"""
对基函数 bf（RWG），在 disaggG 的球面极点"远场"中计算：
  term_efie = Σ_r Σ_q  w_r * w_q * (ρ(r_q) · E^ff(r̂_r)) * exp(-jk r̂_r·(r_q-r0)) * ℓ/2
  term_mfie = Σ_r Σ_q  w_r * w_q * ((ρ(r_q)×n̂)·H^ff(r̂_r)) * exp(-jk r̂_r·(r_q-r0)) * ℓ/2

其中：
  E^ff(r̂_r) = (disaggG[r,1,iCube], disaggG[r,2,iCube])  = (Eθ, Eϕ)
  H^ff = k̂×E/η，故 Hθ = -Eϕ/η，Hϕ = Eθ/η
  ρ(r_q)  — RWG 基函数面内矢量（从自由顶点指向积分点，分 Plus/Minus）
  n̂       — 三角形法向量（来自顶点叉积，需归一化）
  r0      — 叶结点中心（cubeCenter[:,iCube]）

实现参考 Disaggregation.jl 的 add_received_field_rwg! 内层循环；
将其 EFIE 和 MFIE 两个分支合并为同一函数返回两个值即可。
"""
function _receive_terms(bf::RWGFunction, field, r0, k, η, poles)
    term_efie = zero(ComplexF64)
    term_mfie = zero(ComplexF64)
    # 实现：遍历 Plus/Minus 三角形 → 高斯积分点 → 球面极点加权求和
    # （参考 Disaggregation.jl 中 add_received_field_rwg! 的内层循环）
    return term_efie, term_mfie
end
```

### 5.7 assemble_near_field_pmchw — 2N×2N 近场稀疏矩阵

**职责**：遍历八叉树所有近邻三角形对，对每对 (test_bf_i, src_bf_j)（均为 RWG，编号 1:N），计算 4 个矩阵元素并写入 2N×2N 稀疏矩阵。

```
对每个近邻对 (i, j)（i, j ∈ 1:N）：

  Z[i,     j    ] += efie_interaction_k0(tri_test, tri_src)
                   + efie_interaction_k1(tri_test, tri_src)           # EJ 块

  Z[i,     j+N  ] += k_interaction_k0(tri_test, tri_src)
                   + k_interaction_k1(tri_test, tri_src)              # EM 块

  Z[i+N,   j    ]  = -Z[i, j+N]                                     # HJ = -EM（精确关系）

  Z[i+N,   j+N  ] += efie_eta_interaction_k0(tri_test, tri_src)
                   + efie_eta_interaction_k1(tri_test, tri_src)       # HM 块
```

其中：
- `efie_interaction(tri_test, tri_src; k, η)` — 标准 EFIE L 算子（PMCHW.jl 中 EJ 块已有）
- `k_interaction(tri_test, tri_src; k)` — K 算子（PMCHW.jl 中 EM 块已有）
- `efie_eta_interaction(tri_test, tri_src; k, η)` — L 算子但因子 $\eta \to 1/\eta$（HM 块）

**实现策略**：
1. 从 `src/Operators/PMCHW.jl` 的 `assemble_impedance_matrix` 中提取 EJ/EM/HM 三个块的单元交互函数
2. 近邻对从 `octree` 的近场对列表获取（与现有 `assemble_near_field` 流程相同）
3. 近场矩阵验证：`Z_near[1:N, N+1:2N]` 应与 `assemble_impedance_matrix(pmchw, basis)` 的 EM 块一致（仅近邻对，不包含远场对）

---

## 6. 子任务分解与优先级

| 步骤 | 内容 | 前置 | 优先级 |
|------|------|------|-------|
| **15.1** | TDD: `test_pmchw_excitation.jl` (DeltaGap + input_impedance_pmchw) | — | P0 🔴 |
| **15.2** | 实现 `excitation_vector(PMCHW, DeltaGapSource, RWGBasis)` | 15.1 RED | P0 🔴 |
| **15.3** | 实现 `input_impedance(op::PMCHW, source, I_2N, basis)` | 15.1 RED | P0 🔴 |
| **15.4** | TDD: `test_scfie_delta_gap.jl` | — | P0 🔴 |
| **15.5** | 实现 `excitation_vector(SCFIE, DeltaGapSource, rwg, swg)` | 15.4 RED | P0 🔴 |
| **15.6** | 基准脚本 `benchmark/accuracy/run_B1_B5_antenna.jl` | 15.3, 15.5 | P1 🟠 |
| **15.7** | TDD: `test_pmchw_mlfma_operator.jl` | — | P1 🟠 |
| **15.8** | 实现 `assemble_near_field_pmchw`（2N×2N，4 块：EJ/EM/HJ/HM） | 15.7 RED | P1 🟠 |
| **15.9** | 实现 `aggregate_leaf_pmchw!`（单函数，x_range 参数区分 J/M） | 15.7 RED | P1 🟠 |
| **15.10** | 实现 `disaggregate_leaf_pmchw_j!` 和 `_m!`（四块接收核函数） | 15.9 | P1 🟠 |
| **15.11** | 组装 `PMCHWMLFMAOperator` struct（两棵 N 点八叉树 octree0/octree1）+ 构造函数 + `mul!`（4 遍远场） | 15.8, 15.10 | P1 🟠 |
| **15.12** | 更新 `generate_report.jl` 加入 B1–B5 | 15.6 | P2 🟡 |
| **15.13** | 检视迭代 × 2 轮 | 所有 | P2 🟡 |

### 双八叉树 + 四遍远场设计亮点（与旧方案对比）

| 旧方案（已废弃） | 新方案（Gibson Algorithm 14） |
|----------------|------------------------------|
| `MagneticRWGBasis` 标签类型（需新文件） | 无新类型，`PMCHWMLFMAOperator.jl` 自包含 |
| 2N 点八叉树（J/M 坐标重复） | 两个 N 点八叉树（octree0/k0, octree1/k1）；J/M 共享几何 |
| 混入同一 `aggS`，解聚时无法分离 | J×k0/J×k1/M×k0/M×k1 四遍独立 `aggS`；每遍前强制清零 |
| `bfID ≤ N` 判断（排序后失效） | `kmode` Symbol 参数明确指定 k；`x_range` 明确指定 J/M 系数范围 |
| 解聚靠 sorted bfID 判断行类型（脆弱） | sorted_ids 仅 N 长，原始 bfID 直接写 `y[bfID]` 和 `y[bfID+N]` |
| 需扩展 Aggregation.jl / Disaggregation.jl | 不修改现有文件，PMCHWMLFMAOperator.jl 自包含 |

---

## 7. DoD（完成定义）

| 检查项 | 通过标准 |
|--------|---------|
| 新增激励 API 单元测试 | `test_pmchw_excitation.jl` 全通过（含 zero-feed 错误捕获） |
| SCFIE delta-gap 测试 | `test_scfie_delta_gap.jl` 全通过 |
| B1 PMCHW Direct | Re(Z_in) > 0，εᵣ→1 时 Z_in 趋近 EFIE 结果（相对误差 <10%） |
| B2 PMCHW MLFMA | GMRES 收敛，Z_in_B2 与 Z_in_B1 误差 <5%（使用 PMCHWMLFMAOperator） |
| B3 VS-EFIE Direct | Z_in Re 误差 (vs EFIE-only) <10% |
| B4 VS-CFIE Direct | Z_in 与 B3 差异 <5Ω |
| B5 VS-CFIE MLFMA | Z_in 与 B4 差异 <5% |
| MagneticRWGBasis 设计 | ~~已废弃~~ — 不需要该类型 |
| `disaggregate_leaf_pmchw_j!/m!` | 4 个块的因子（EJ/EM/HJ/HM）与 Direct solve 的矩阵元素对齐（数量级 + 符号） |
| `assemble_near_field_pmchw` | 2N×2N 近场矩阵；`Z_near[1:N,N+1:2N]` 与 `assemble_impedance_matrix(pmchw,basis)[1:N,N+1:2N]`（近邻部分）一致 |
| 检视迭代 | ≥ 2 轮连续 clean |

---

## 8. 技术风险与缓解

| 风险 | 概率 | 缓解措施 |
|------|------|---------|
| PMCHW delta-gap 符号/因子错误 | 中 | TDD: 先验证 ε→1 极限等于 EFIE，再跑实际介质 |
| `disaggregate_leaf_pmchw_j/m!` 的 EJ/EM/HJ/HM 因子搞错符号或 η | 高 | 单独验证每块：用 Z_near（直接法等效）与 Direct 矩阵元素逐项对比 |
| 四遍远场 aggS 清零遗漏（前一遍残留污染后一遍） | 中 | 在每遍聚合前显式 `fill!(lv.aggS, 0)`；单元测试：零输入 → 零输出 |
| octree0 和 octree1 的叶结构不一致（Lebedev 阶数不同）导致近邻对不匹配 | 低 | `assemble_near_field_pmchw` 只用 octree0 的近邻对；两树几何相同，near-pair 一致 |
| 介质天线无解析参考值 | 高 | 使用自洽检查（Re>0、freq sweep 趋势、ε→1 极限）代替绝对精度 |
| k1 复数时 Lebedev 展开阶数不足 | 中 | octree1 用 `λ1 = 2π/|k1|` 自动选极点阶数；先用无损介质验证 |
| `assemble_near_field_pmchw` K 块因子 | 中 | 与 `PMCHW.jl` 里的 K 块近场计算逐项对比校验 |

---

## 9. 进度跟踪入口

- ROADMAP: [REFACTORING_ROADMAP.md](REFACTORING_ROADMAP.md) `## Phase 15` 节
- PROGRESS: [REFACTORING_PROGRESS.md](REFACTORING_PROGRESS.md) `Phase 15` 节
- 基准脚本: `benchmark/accuracy/run_B1_B5_antenna.jl`
- 测试文件: `test/test_pmchw_excitation.jl`, `test/test_scfie_delta_gap.jl`, `test/test_pmchw_mlfma_operator.jl`
- 实现文件: `src/FastAlgorithms/MLFMA/PMCHWMLFMAOperator.jl`（自包含：struct/构造函数/mul!/聚合/解聚/近场装配）

---

## 附录: VS-EFIE 与 VS-CFIE 命名约定

本文档遵循以下命名：

| 名称 | EMSuite 代码 | α 值 | 说明 |
|------|------------|------|------|
| **VS-EFIE** | `SCFIE(freq, perms; alpha=0.0)` | 0.0 | 表面 EFIE + 体积 EFIE |
| **VS-CFIE** | `SCFIE(freq, perms; alpha=0.5)` | 0.5 | 表面 CFIE (EFIE+MFIE) + 体积 EFIE |

> α=0.5 为 CFIE 的标准混合参数（50% EFIE + 50% MFIE），通常比纯 EFIE 条件数更好。

---

*本计划文档依照 `copilot-instructions.md` 原则 9 规范撰写。*
