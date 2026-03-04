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

## 5. PMCHWMLFMAOperator 实现计划

### 5.1 文件结构

```
src/FastAlgorithms/MLFMA/
├── PMCHWMLFMAOperator.jl   ← 新建
│   ├── struct PMCHWMLFMAOperator{FT,CT}
│   ├── PMCHWMLFMAOperator(pmchw, basis, leaf_size)  # 构造函数
│   ├── Base.size, Base.eltype
│   ├── mul!(y, A::PMCHWMLFMAOperator, x)            # 核心计算
│   └── get_leaf_intervals(op::PMCHWMLFMAOperator)   # 用于 Block Jacobi
└── MLFMA.jl  ← 新增 include 和 export
```

### 5.2 结构体定义

```julia
struct PMCHWMLFMAOperator{FT,CT} <: AbstractIntegralOperator
    pmchw     :: PMCHW{FT,CT}          # 底层 PMCHW 算子（含 k0, k1, η0, η1）
    basis     :: RWGBasis              # 共用 RWGBasis（J 和 M 都在同一网格上）
    Z_near    :: SparseMatrixCSC{CT}   # 2N×2N 近场稀疏矩阵（从 PMCHW 直接装配）
    # 两个 EFIE-型 MLFMAOperator（分别处理 k0 和 k1 的 L 贡献）
    mlfma_k0  :: MLFMAOperator{FT,CT}  # EFIE(k0, η0)
    mlfma_k1  :: MLFMAOperator{FT,CT}  # EFIE(k1, η1)
    # 可选: K-型算子（阶段 2 添加）
    # mlfma_k0_K :: MLFMAOperator{FT,CT}
    # mlfma_k1_K :: MLFMAOperator{FT,CT}
    sorted_ids     :: Vector{Int}
    inv_sorted_ids :: Vector{Int}
    freq      :: FT                    # 供 AbstractIntegralOperator 接口
end
```

### 5.3 mul! 的阶段实现

**阶段 1（仅近场，可立即验证接口）**:
```julia
function LinearAlgebra.mul!(y, A::PMCHWMLFMAOperator, x)
    mul!(y, A.Z_near, x)   # 纯近场：等价于直接法但用稀疏矩阵
end
```

**阶段 2（完整 MLFMA）**:
```julia
function LinearAlgebra.mul!(y, A::PMCHWMLFMAOperator, x)
    N = size(A.mlfma_k0, 1)   # RWGBasis 的 N
    x_J = x[1:N];   x_M = x[N+1:2N]
    y_E = view(y, 1:N);   y_H = view(y, N+1:2N)

    # 近场贡献
    mul!(y, A.Z_near, x)   # sparse mat-vec (2N×2N)

    # 远场贡献（从 Z_far_col[J] 和 Z_far_col[M] 分别计算）
    # EJ block: L(k0)*J + L(k1)*J
    y_EJ_k0 = zeros(CT, N); mul!(y_EJ_k0, A.mlfma_k0, x_J)
    y_EJ_k1 = zeros(CT, N); mul!(y_EJ_k1, A.mlfma_k1, x_J)
    y_E .+= y_EJ_k0 .+ y_EJ_k1
    # ... (K 块和 Lη 块需要额外的 K-type MLFMA 实例)
end
```

> 阶段 1 先完成，验证 B1/B2 接口正确性后，再在阶段 2 扩展完整 far-field。

### 5.4 近场矩阵构建

参考 `assemble_impedance_matrix(pmchw, basis)` 已返回完整 2N×2N 矩阵。  
PMCHWMLFMAOperator 构建时直接截取 Z_near_sparse 部分（近邻块）：

```julia
function PMCHWMLFMAOperator(pmchw::PMCHW, basis::RWGBasis, leaf_size::Float64)
    # 1. 构建单一 RWGBasis 的八叉树
    # 2. 装配 Z_near (2N×2N 近邻稀疏矩阵，从完整矩阵筛选近邻对)
    # 3. 构建 mlfma_k0 = MLFMAOperator(EFIE(freq, η0_equiv), basis, leaf_size)
    #    其中 η0_equiv = η0  → 统一远场聚合
    # 4. 构建 mlfma_k1 = MLFMAOperator(EFIE(freq_equiv, η1), basis, leaf_size)
    #    注意: k1 ≠ k0，需要分别构建八叉树（不同的 Lebedev 展开阶数）
end
```

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
| **15.8** | `PMCHWMLFMAOperator` 阶段 1（仅近场） | 15.7 RED | P1 🟠 |
| **15.9** | `PMCHWMLFMAOperator` 阶段 2（添加 Z_far，K 块） | 15.8 | P2 🟡 |
| **15.10** | 更新 `generate_report.jl` 加入 B1–B5 | 15.6 | P2 🟡 |
| **15.11** | 检视迭代 × 2 轮 | 所有 | P2 🟡 |

---

## 7. DoD（完成定义）

| 检查项 | 通过标准 |
|--------|---------|
| 新增激励 API 单元测试 | `test_pmchw_excitation.jl` 全通过（含 zero-feed 错误捕获） |
| SCFIE delta-gap 测试 | `test_scfie_delta_gap.jl` 全通过 |
| B1 PMCHW Direct | Re(Z_in) > 0，εᵣ→1 时 Z_in 趋近 EFIE 结果（相对误差 <10%） |
| B2 PMCHW MLFMA (阶段1) | GMRES 收敛，Z_in_B2 与 Z_in_B1 误差 <5% |
| B3 VS-EFIE Direct | Z_in Re误差 (vs EFIE-only) <10% |
| B4 VS-CFIE Direct | Z_in 与 B3 差异 <5Ω |
| B5 VS-CFIE MLFMA | Z_in 与 B4 差异 <5% |
| PMCHWMLFMAOperator 阶段 2 | B2 MLFMA far-field residual <1e-3，Z_in 与 B1 误差 <5% |
| 检视迭代 | ≥ 2 轮连续 clean |

---

## 8. 技术风险与缓解

| 风险 | 概率 | 缓解措施 |
|------|------|---------|
| PMCHW delta-gap 符号/因子错误 | 中 | TDD: 先验证 ε→1 极限等于 EFIE，再跑实际介质 |
| PMCHWMLFMAOperator K-块 far-field 复杂 | 高 | 阶段 1 先跳过 K 块远场（仅近场），验证收敛后再加 |
| 介质天线无解析参考值 | 高 | 使用自洽检查（Re>0、freq sweep 趋势、ε→1 极限）代替绝对精度 |
| sphere_600MHz.nas 网格不适合天线馈电（弯曲面难找馈边） | 中 | 程序生成介质球（generate_cylinder_mesh 改为球面网格），或用平面偶极子形状 |
| k1 复数时 MLFMA Lebedev 展开需要更高阶 | 中 | 以 max(|k0|, |k1|) 为准选展开阶数；有损介质先用无损版本验证 |

---

## 9. 进度跟踪入口

- ROADMAP: [REFACTORING_ROADMAP.md](REFACTORING_ROADMAP.md) `## Phase 15` 节
- PROGRESS: [REFACTORING_PROGRESS.md](REFACTORING_PROGRESS.md) `Phase 15` 节
- 基准脚本: `benchmark/accuracy/run_B1_B5_antenna.jl`
- 测试文件: `test/test_pmchw_excitation.jl`, `test/test_scfie_delta_gap.jl`, `test/test_pmchw_mlfma_operator.jl`
- 实现文件: `src/FastAlgorithms/MLFMA/PMCHWMLFMAOperator.jl`

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
