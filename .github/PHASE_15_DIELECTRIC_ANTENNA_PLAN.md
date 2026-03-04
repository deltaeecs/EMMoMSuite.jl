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

### 5.0 设计原则（修订）

> **指导思想**：参照 VS-EFIE（SCFIE）的混合基函数 MLFMA 模式。VS-EFIE 中 `bases = [RWGBasis, SWGBasis]`，两类基函数按空间位置在同一八叉树中交错（不是按类型连续分区）。PMCHW 也应如此——J 和 M 电磁流共享同一网格，八叉树自然将 (J_i, M_i) 打包进同一叶结点。
>
> **不做的事**：不创建两个独立的 `MLFMAOperator` 包装器，不在 `mul!` 外层手动切割 `y_E = y[1:N]` / `y_H = y[N+1:2N]` 索引块。
>
> **做的事**：引入轻量 `MagneticRWGBasis` 标签类型；在底层 `aggregate_leaf!`/`disaggregate_leaf!`/`assemble_near_field` 里按 operator isa PMCHW 分派正确核函数；`PMCHWMLFMAOperator.mul!` 直接操控八叉树内部，两趟聚合-平移-解聚覆盖全部 4 个矩阵块。

### 5.1 文件结构

```
src/BasisFunctions/
└── MagneticRWG.jl          ← 新建：MagneticRWGBasis 包装类型

src/FastAlgorithms/MLFMA/
├── PMCHWMLFMAOperator.jl   ← 新建：struct + mul! + 自定义 disaggregate_leaf_pmchw!
├── Aggregation.jl          ← 扩展：支持 MagneticRWGBasis（调用同一 add_radiation_pattern_rwg!）
├── Disaggregation.jl       ← 扩展：新增 disaggregate_leaf_pmchw! 处理 4 块测试
└── MLFMAOperator.jl        ← 扩展：assemble_near_field 识别 PMCHW 的 4 种交叉块
```

### 5.2 MagneticRWGBasis — M 电流标签类型

**位置**: `src/BasisFunctions/MagneticRWG.jl`

```julia
"""
    MagneticRWGBasis{IT,FT} <: AbstractBasisFunction

轻量包装，将 RWGBasis 标记为"磁流（M）"类型。
空间几何与原始 RWGBasis 完全相同；唯一作用：在 MLFMA
聚合/解聚/近场装配时触发 PMCHW 专用的 K 算子核函数分支。
"""
struct MagneticRWGBasis{IT,FT} <: AbstractBasisFunction
    basis::RWGBasis{IT,FT}
end

# 委托接口：从底层 RWGBasis 获取所有几何属性
CoreModule.num_basis(b::MagneticRWGBasis) = num_basis(b.basis)
Base.getproperty(b::MagneticRWGBasis, s::Symbol) =
    s === :basis ? getfield(b, :basis) : getproperty(b.basis, s)
```

> **关键效果**：`bases = [rwg_basis, MagneticRWGBasis(rwg_basis)]` 时，`basis_offsets = [N, 2N]`，两类基函数中心坐标完全重合；八叉树对 2N 个点排序后，叶结点中 J_i 和 M_i 自然相邻（"混着的"，正是用户要求的效果）。

### 5.3 PMCHWMLFMAOperator 结构体

```julia
struct PMCHWMLFMAOperator{FT,CT} <: AbstractIntegralOperator
    pmchw     :: PMCHW{FT,CT}
    basis     :: RWGBasis              # N 个 DOF 的 RWGBasis（J 和 M 共用）
    Z_near    :: SparseMatrixCSC{CT}   # 2N×2N 近场（含 4 个交叉块）
    octree    :: OctreeInfo            # 从 2N 中心点建立（J 和 M 坐标重叠）
    sorted_ids     :: Vector{Int}      # 2N 长排列
    inv_sorted_ids :: Vector{Int}
    freq      :: FT
end

Base.size(op::PMCHWMLFMAOperator)    = size(op.Z_near)
Base.eltype(::PMCHWMLFMAOperator{FT,CT}) where {FT,CT} = CT
```

### 5.4 mul! — 两趟聚合共享八叉树

```julia
function LinearAlgebra.mul!(y::AbstractVector, A::PMCHWMLFMAOperator, x::AbstractVector)
    N = num_basis(A.basis)

    # ① 近场（2N×2N 稀疏矩阵，包含全部 4 块近邻交互）
    mul!(y, A.Z_near, x)

    # ② 远场 — 趟 1：J 源（x[1:N]）
    #    聚合：standard RWGBasis radiation pattern
    x_J = @view x[1:N]
    x_pad_J = zeros(eltype(x), 2N); x_pad_J[1:N] .= x_J   # 只含 J 系数
    aggregate!(A.octree, AbstractBasisFunction[A.basis, MagneticRWGBasis(A.basis)],
               [N, 2N], A.pmchw, x_pad_J, A.sorted_ids)
    #    注：aggregate_leaf! 对 J 和 M 用同一 add_radiation_pattern_rwg!
    #        （M 辐射花样积分形式 ∫ρ e^{jkr̂·r} 与 J 相同；区别在解聚）

    _pmchw_translate_disaggregate!(y, A, :J)

    # ③ 远场 — 趟 2：M 源（x[N+1:2N]）
    x_M = @view x[N+1:2N]
    x_pad_M = zeros(eltype(x), 2N); x_pad_M[N+1:2N] .= x_M
    aggregate!(A.octree, AbstractBasisFunction[A.basis, MagneticRWGBasis(A.basis)],
               [N, 2N], A.pmchw, x_pad_M, A.sorted_ids)

    _pmchw_translate_disaggregate!(y, A, :M)

    return y
end

function _pmchw_translate_disaggregate!(y, A::PMCHWMLFMAOperator, src_type::Symbol)
    # 平移（水平传递）
    for levelID = 2:A.octree.nLevels
        translate!(A.octree.levels[levelID])
    end
    # 向下传递
    for levelID = 2:(A.octree.nLevels-1)
        disaggregate_downward!(A.octree.levels[levelID], A.octree.levels[levelID+1])
    end
    # 叶结点解聚 —— PMCHW 专用：根据 src_type + test_type 选核函数
    N = num_basis(A.basis)
    y_far = zeros(eltype(y), 2N)
    disaggregate_leaf_pmchw!(
        A.octree.levels[A.octree.nLevels],
        A.basis, A.pmchw, y_far, A.sorted_ids, src_type,
    )
    y .+= y_far
end
```

### 5.5 disaggregate_leaf_pmchw! — 4 块测试公式

PMCHW 的远场测试逻辑（在 `Disaggregation.jl` 新增，或写入 `PMCHWMLFMAOperator.jl`）：

| src_type | test_type | 核函数 | 因子 |
|----------|-----------|--------|------|
| J | E 方程（bfID ≤ N） | L（EFIE 测试 ρ·E） | η₀/k₀ + η₁/k₁ 相应因子 |
| J | H 方程（bfID > N） | −K（MFIE 测试 (ρ×n̂)·H，负号） | −1/k 相关 |
| M | E 方程（bfID ≤ N） | K（MFIE 测试 (ρ×n̂)·H） | +1/k 相关 |
| M | H 方程（bfID > N） | Lη（EFIE/η 测试 ρ·E） | 1/η 调整 |

实现参考 `add_received_field_rwg!` 中已有的 EFIE/MFIE 分支：

```julia
function disaggregate_leaf_pmchw!(
    leaf_level, basis::RWGBasis, pmchw::PMCHW, y_far, sorted_ids, src_type::Symbol
)
    N = num_basis(basis)
    k0, eta0, k1, eta1 = pmchw.k0, pmchw.eta0, pmchw.k1, pmchw.eta1
    elem_info = get_triangles_info(basis.mesh, basis)
    gq        = GaussQuadratureInfo(:Triangle, 3, Float64)
    poles     = leaf_level.poles
    nPoles    = length(poles.r̂sθsϕs)
    disaggG   = leaf_level.disaggG

    Threads.@threads for iCube = 1:leaf_level.nCubes
        cube = leaf_level.cubes[iCube]
        field = view(disaggG, :, :, iCube)
        for bfID_sorted in cube.bfInterval
            bfID   = sorted_ids[bfID_sorted]
            is_E_row = bfID ≤ N                     # E 方程（前 N 行）
            local_id = is_E_row ? bfID : bfID - N   # 本地基函数编号
            bf       = basis.functions[local_id]

            # 计算 term_efie（ρ·E 积分）和 term_mfie（(ρ×n̂)·H 积分）
            # … (积分循环，复用 add_received_field_rwg! 的内层逻辑)

            val = zero(ComplexF64)
            if src_type === :J && is_E_row
                # EJ 块：L 算子，factor = jk0η0/(16π) + jk1η1/(16π)
                val = term_efie_k0 * factor_EJ_k0 + term_efie_k1 * factor_EJ_k1
            elseif src_type === :J && !is_E_row
                # HJ 块：-K 算子，factor = -(jk0/(16π) + jk1/(16π))
                val = -(term_mfie_k0 * factor_K_k0 + term_mfie_k1 * factor_K_k1)
            elseif src_type === :M && is_E_row
                # EM 块：+K 算子
                val = term_mfie_k0 * factor_K_k0 + term_mfie_k1 * factor_K_k1
            else   # :M, H 方程
                # HM 块：Lη 算子，factor = jk0/(η0·16π) + jk1/(η1·16π)
                val = term_efie_k0 * factor_HM_k0 + term_efie_k1 * factor_HM_k1
            end

            y_far[bfID] += val
        end
    end
end
```

### 5.6 assemble_near_field 近场扩展（PMCHW 4 块）

在 `MLFMAOperator.jl` 的 `assemble_near_field` 内，**新增 PMCHW 分支**：

```julia
# 在现有 if operator isa SCFIE 分支后添加：
elseif operator isa PMCHW
    # RWG×RWG → EJ 块（L 核）
    if !isempty(my_tris) && !isempty(neigh_tris)
        # (与现有 efie_interaction! 类似，但 operator 用 pmchw 的 L 子算子)
    end
    # RWG×MagneticRWG → EM 块（K 核）
    # MagneticRWG×RWG → HJ 块（-K 核）
    # MagneticRWG×MagneticRWG → HM 块（Lη 核）
end
```

> **注**：可直接调用 `pmchw_block_interaction!(Z_local, pmchw, tri_test, tri_src, test_type, src_type)` 辅助函数，封装 EJ/EM/HJ/HM 四种核函数选择逻辑，代码清晰，符合 Legacy PMCHW 装配逻辑（`assemble_impedance_matrix` 中已有对应子块计算）。

### 5.7 构造函数

```julia
function PMCHWMLFMAOperator(pmchw::PMCHW, basis::RWGBasis, leaf_size::Float64)
    # 1. 构建八叉树：从 2N 个中心点（N J + N M，位置相同）
    m_basis = MagneticRWGBasis(basis)
    bases   = AbstractBasisFunction[basis, m_basis]
    offsets = cumsum([num_basis(b) for b in bases])   # [N, 2N]
    centers_J = reduce(hcat, [bf.center for bf in basis.functions])
    centers_M = reduce(hcat, [bf.center for bf in m_basis.functions])  # 同 J
    all_centers = hcat(centers_J, centers_M)
    λ = 299792458.0 / pmchw.freq
    octree, sorted_ids = build_octree(all_centers, leaf_size; λ = λ)
    N = num_basis(basis)
    inv_sorted_ids = zeros(Int, 2N)
    for i = 1:2N; inv_sorted_ids[sorted_ids[i]] = i; end

    # 2. 装配 Z_near（2N×2N，含全部 4 块近邻交互）
    Z_near = assemble_near_field(pmchw, bases, offsets, octree, sorted_ids, inv_sorted_ids)

    FT = typeof(pmchw.freq)
    CT = Complex{FT}
    return PMCHWMLFMAOperator{FT,CT}(pmchw, basis, Z_near, octree, sorted_ids, inv_sorted_ids, pmchw.freq)
end
```

**关键优势**：
- 八叉树建在 2N 点上，J_i 与 M_i 坐标相同 → 叶结点内 J/M 自然交错（"混着的"）
- 两趟 `aggregate!` 复用完全相同的八叉树结构，无重复构建
- 没有 `y_E = y[1:N]` 等外部索引切割；J/M 的区别由 `bfID ≤ N` 在解聚内部处理
- 与 VS-EFIE 设计完全对称（VS-EFIE 也是 `bases = [RWG, SWG]`，这里是 `[RWG, MagneticRWG]`）

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
| **15.8** | 新建 `MagneticRWGBasis` 标签类型（`src/BasisFunctions/MagneticRWG.jl`） | 15.7 RED | P1 🟠 |
| **15.9** | 扩展 `aggregate_leaf!`：支持 `MagneticRWGBasis`（复用 RWG 辐射花样积分） | 15.8 | P1 🟠 |
| **15.10** | 新建 `PMCHWMLFMAOperator` + `disaggregate_leaf_pmchw!`（4 块测试公式） | 15.9 | P1 🟠 |
| **15.11** | 扩展 `assemble_near_field`：PMCHW 的 4 种交叉块近场（K 核 + Lη 核） | 15.10 | P1 🟠 |
| **15.12** | 更新 `generate_report.jl` 加入 B1–B5 | 15.6 | P2 🟡 |
| **15.13** | 检视迭代 × 2 轮 | 所有 | P2 🟡 |

### 受益于 MagneticRWGBasis 的设计亮点

与旧方案（"阶段1仅近场 → 阶段2加远场"）相比，新方案：

| 旧方案 | 新方案 |
|--------|--------|
| `PMCHWMLFMAOperator` 包含 2 个独立 `MLFMAOperator` 对象 | 只含 1 个 `OctreeInfo` + `Z_near` |
| `mul!` 外层手动 `y_E = y[1:N]` / `y_H = y[N+1:2N]` | J/M 区分在 `disaggregate_leaf_pmchw!` 内部以 `bfID ≤ N` 判断 |
| 两套八叉树（k0 和 k1 分别构建） | 一套八叉树，两趟聚合-平移-解聚，复用树结构 |
| K 块远场另建 CFIE 算子 | 直接在 `disaggregate_leaf_pmchw!` 里按 `term_mfie` 公式计算，零额外算子 |
| 不可扩展为混合 J/M 空间排列 | 八叉树对 2N 点排序，叶内 (J_i, M_i) 天然交错，与 VS-EFIE 对称 |

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
| MagneticRWGBasis 设计 | 八叉树叶结点中 (J_i, M_i) 按空间位置交错；与 VS-EFIE 的 (RWG,SWG) 模式对称 |
| disaggregate_leaf_pmchw! | 4 个块的因子（EJ/EM/HJ/HM）与 Direct solve 的矩阵元素对齐（数量级 + 符号） |
| `assemble_near_field` PMCHW 分支 | K 块近场矩阵元素与 `assemble_impedance_matrix(pmchw,basis)[1:N,N+1:2N]` 一致 |
| 检视迭代 | ≥ 2 轮连续 clean |

---

## 8. 技术风险与缓解

| 风险 | 概率 | 缓解措施 |
|------|------|---------|
| PMCHW delta-gap 符号/因子错误 | 中 | TDD: 先验证 ε→1 极限等于 EFIE，再跑实际介质 |
| `disaggregate_leaf_pmchw!` 的 EJ/EM/HJ/HM 因子搞错符号或 η | 高 | 单独验证每块：用 Z_near（直接法等效）与 Direct 矩阵元素逐项对比 |
| `MagneticRWGBasis` 在 `aggregate_leaf!` 中被遗漏（error 分支） | 低 | 扩展 `aggregate_leaf!` 前先写单元测试，确认 `basis isa MagneticRWGBasis` 被识别 |
| 两趟聚合重置 `aggS` 不彻底（第2趟叠加第1趟残留） | 中 | 在 `aggregate!` 入口处显式 `fill!(level.aggS, 0)` |
| 介质天线无解析参考值 | 高 | 使用自洽检查（Re>0、freq sweep 趋势、ε→1 极限）代替绝对精度 |
| k1 复数时 Lebedev 展开阶数不足 | 中 | 以 max(|k0|, |k1|) 选展开阶数；先用无损介质验证 |
| `assemble_near_field` PMCHW K 块因子 | 中 | 与 `PMCHW.jl` 里的 `pmchw_em_interaction!` 函数返回值对比校验 |

---

## 9. 进度跟踪入口

- ROADMAP: [REFACTORING_ROADMAP.md](REFACTORING_ROADMAP.md) `## Phase 15` 节
- PROGRESS: [REFACTORING_PROGRESS.md](REFACTORING_PROGRESS.md) `Phase 15` 节
- 基准脚本: `benchmark/accuracy/run_B1_B5_antenna.jl`
- 测试文件: `test/test_pmchw_excitation.jl`, `test/test_scfie_delta_gap.jl`, `test/test_pmchw_mlfma_operator.jl`
- 新类型: `src/BasisFunctions/MagneticRWG.jl`
- 实现文件: `src/FastAlgorithms/MLFMA/PMCHWMLFMAOperator.jl`
- 修改文件: `src/FastAlgorithms/MLFMA/Aggregation.jl`（MagneticRWGBasis 分支）, `src/FastAlgorithms/MLFMA/Disaggregation.jl`（disaggregate_leaf_pmchw!）, `src/FastAlgorithms/MLFMA/MLFMAOperator.jl`（assemble_near_field PMCHW 分支）

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
