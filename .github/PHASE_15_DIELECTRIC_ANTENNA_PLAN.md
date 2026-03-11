# Phase 15: 介质与金属-介质混合天线精度测试 + PMCHW Transmission Block/Operator 架构

> 创建日期: 2026-03-07  
> 状态: 计划中  
> 前置: Phase 14 完成（A1–A4, F1–F9, P1–P3 基准脚本已创建）

> 2026-03-07 下一子流补记：B1–B5 与 PMCHW shell / Gate S 已完成后，Phase 15 的下一正式开发焦点已切换为 alternative formulation 基线，优先项为 N-Muller dense baseline。对应执行文档：`.github/plans/phase_15_nmuller_dense_baseline.md`。
> 2026-03-07 子流续记：`NMuller` 已补齐 `DeltaGapSource` 激励与 `input_impedance` 最小链路；`PMCHWMLFMAOperator` 已新增 `PMCHWMLFMAErrorBudget`，把 `near_range / L_min / leaf_size_eff` 从隐式启发式提升为显式接口。
> 2026-03-07 预算子流补记：已新增 `benchmark/compare_pmchw_mlfma_budget.jl`，用于在固定 small / medium 夹具上量化比较不同 `PMCHWMLFMAErrorBudget` 配置对 `nnz_near`、matvec 误差、强形式 `Z_in` 误差与耗时的影响。
> 2026-03-07 预算子流续记：默认预算集合现已冻结为 `default / loose_near / fixed_leaf_0p04_nr9` 三组代表性配置；并已新增 `benchmark/compare_pmchw_mlfma_budget_krylov.jl`，把这三组预算直接接入 medium 夹具的 `short/long` Krylov 对照。
> 2026-03-08 主线修正：PMCHW K-type receive 修复后，旧的 `7`–`16Ω` long-Krylov budget 大分叉已被推翻。现行 `BG2`、Arnoldi 子空间诊断与 GMRES 轨迹诊断一致表明：`MLFMA strong` 与 `dense strong` 在 medium 夹具上已基本贴合；剩余问题应继续收敛到 **更深 Krylov 轨迹 / formulation 条件数放大**，而不是继续停留在 budget 主导的 backend fidelity 归因。

## 当前总计划

### A. 当前主线

1. **PMCHW 仍是主交付 formulation**：继续围绕 PMCHW shell + strong-form + MLFMA backend 做 medium 级精度与收敛治理。
2. **N-Muller 是对照组，不是当前阻塞主线**：它只用于回答“剩余问题更像 formulation 还是 backend”这个归因问题。
3. **当前剩余问题已经从 receive bug / budget 大分叉 收缩到更深 Krylov / conditioning 机制**。

### B. 当前已完成的主干

1. B1–B5 天线基准已打通。
2. PMCHW block/operator shell 与 Gate S 四路对照已完成。
3. PMCHW M-pass receive 链错误已修复，BF1 / BG2 / budget benchmark 已同步更新。
4. Medium 级 Arnoldi 子空间与 GMRES 轨迹诊断已落地。

### C. N-Muller 的明确角色

1. **已完成**：dense baseline、small/medium dense 对照、conditioning / GMRES 行为对照。
2. **未作为当前主线**：N-Muller 的 DeltaGap / input_impedance 端口语义仍未校准，因此不进入正式天线阻抗验收门。
3. **后续使用方式**：只在需要区分 PMCHW formulation 问题与 PMCHW MLFMA backend 问题时作为 dense 对照基线使用。

### D. 下一阶段执行顺序

1. 先继续推进 PMCHW 的更深 Krylov 诊断，而不是扩展 N-Muller 天线语义。
2. 已完成：在同一 medium 夹具上，已分别用 random RHS 与 `PlaneWave` 物理激励把 PMCHW / N-Muller 的 dense GMRES 轨迹做成正式对照；但该结论现已被 restart 扫描进一步细化。
3. 已完成：PMCHW 自身的 medium dense `weak/strong` plane-wave 轨迹也已固化；当前结果表明 strong-form 只有边际改善，不能单独解除主问题。
4. 已完成：plane-wave dense restart 扫描已证明默认 `restart=20` 会显著夸大 PMCHW 的坏轨迹；把 `restart` 提升到 `250` 后，PMCHW 可把相对残差压到 `~1e-5`，但相对 LU 误差仍明显落后于 `NMuller`。
5. 下一步若要继续收缩剩余边界，应进入更深 Arnoldi 子空间或显式 full-restart / full-GMRES 轨迹诊断，而不是扩展 N-Muller 天线语义。

---

## 0. 开发原则遵循声明

> 2026-03-07 实施状态补记：PMCHW/SCFIE 的 DeltaGap 激励与 PMCHW 输入阻抗 API 已在源码中落地，并已通过 `test/test_pmchw_excitation.jl` 与 `test/test_scfie_delta_gap.jl` 验证且纳入默认 `runtests.jl`。Dense shell 最小实现已落地，新增 `PMCHWBlockOperator`、`DensePMCHWBackend` 与 `test/test_pmchw_operator_shell.jl`；默认 `strong_form` 现已接入真实 RWG surface Gram 的 2N block pairing。另已新增 `test/test_pmchw_gate_s_dense.jl`，完成 Gate S 的 dense 半边回归；`PMCHWMLFMAOperator` 也已通过 `MatrixFreePMCHWBackend` 接入 shell，并由 `test/test_pmchw_operator_shell_mlfma.jl` 在小夹具上验证 weak/strong 链路可执行。最新已新增 `test/test_pmchw_gate_s_mlfma_medium.jl`，在 `N=540` 夹具上正式完成 `dense weak / dense strong / MLFMA weak / MLFMA strong` 四路 Gate S 对照并通过；该回归运行约 9 分钟，暂不并入默认 `runtests.jl`。Phase 15 天线基准脚本 `benchmark/accuracy/run_B1_B5_antenna.jl` 现已完整跑通 `B1`–`B5`，对应 `generate_report.jl` 汇总也已对齐。另：`NMuller` 已补齐 `DeltaGapSource` + `input_impedance` 最小天线接口链，并由 `test/test_nmuller_excitation.jl` 锁定；PMCHW MLFMA backend 现已新增 `PMCHWMLFMAErrorBudget`，可显式控制 `near_range`、`L_min` 与 `leaf_size_eff`，旧构造入口保持兼容。

本 Phase 适用以下核心原则（见 `copilot-instructions.md`）：

| 原则 | 适用内容 |
|------|---------|
| **原则 1 (TDD)** | 新增激励 API、PMCHW block/operator shell、backend 适配层、`input_impedance_pmchw` 均先写单元测试 |
| **原则 2 (Legacy 对齐)** | PMCHW 激励及 MLFMA 结果必须与 Direct 结果一致，不使用经验常数 |
| **原则 3 (差异排查)** | MLFMA 结果与 Direct 偏差时逐块比较：Z_near 是否正确、far-field k0/k1 分离是否正确 |
| **原则 5 (Git 提交)** | 每类子任务独立提交：API 扩展、block/operator 外壳、backend 适配、单元测试、基准脚本 |
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
| PMCHW block/operator shell | `PMCHW.jl` `assemble_impedance_matrix` | 四块 `EJ/EM/HJ/HM` 组合必须与 Direct PMCHW 一致 |
| PMCHW MLFMA backend | `PMCHW.jl` `assemble_impedance_matrix` + 现有 MLFMA shared core | 只负责 nonlocal evaluator，不再承载顶层 transmission 语义 |

---

## 2. 目标

在 Phase 14 天线测试（A1–A4, 纯金属偶极子）基础上，新增：

1. **B1–B2**: 介质谐振体（PMCHW）Delta-Gap 馈电天线，输入阻抗测试
2. **B3–B5**: 金属-介质混合天线（VS-EFIE/VS-CFIE, 即 SCFIE），输入阻抗测试
3. **PMCHW block/operator shell**: 先实现 Bempp 风格 2×2 块 transmission 外壳，支持 weak/strong form 与 backend 切换
4. **PMCHW MLFMA backend**: 将现有 PMCHWMLFMA 路径收编为 backend，服务于 B2 的 fast 求解路径

> 2026-03-07 子流补记：`15.14` 与 `15.15` 已在小球 dense 基线层面完成；当前已有 `test/test_nmuller.jl`、`test/test_nmuller_comparison.jl`、`test/test_nmuller_excitation.jl` 与 `benchmark/compare_pmchw_nmuller_sphere.jl` 作为正式回归/对照入口。与此同时，PMCHW fast backend 已新增 `PMCHWMLFMAErrorBudget` 作为显式误差预算接口。

| 代号 | EMSuite 方程 | 激励方式 | 求解 |
|------|-------------|---------|------|
| VS-EFIE | `SCFIE(freq, perms; alpha=0.0)` | DeltaGap (表面) | Direct |
| VS-CFIE | `SCFIE(freq, perms; alpha=0.5)` | DeltaGap (表面) | Direct |
| PMCHW Direct | `PMCHW(freq, eps_r, mu_r)` | DeltaGap (新增) | Direct |
| PMCHW Dense Shell | PMCHW block/operator shell + Dense backend | DeltaGap | Direct / GMRES |
| PMCHW MLFMA | PMCHW block/operator shell + MLFMA backend | DeltaGap | GMRES |

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

### 3.3 PMCHW transmission block/operator shell（新主线）

**目标位置**: transmission operator 层（顶层） + backend 层（Dense / MLFMA / 未来 H-matrix）

#### 数学结构

PMCHW 系统 (2N × 2N)：
$$Z = \begin{bmatrix} Z^{EJ} & Z^{EM} \\ Z^{HJ} & Z^{HM} \end{bmatrix}$$

其中（超简化表示）：
- $Z^{EJ}(k_0, \eta_0) + Z^{EJ}(k_1, \eta_1)$：EFIE-型 L 算子，两个波数之和
- $Z^{EM}(k_0) + Z^{EM}(k_1)$：K 算子（MFIE 类型），两个波数之和  
- $Z^{HJ} = -Z^{EM}$（结构不变量，精确成立）
- $Z^{HM}(k_0, 1/\eta_0) + Z^{HM}(k_1, 1/\eta_1)$：倒 η 的 EFIE-型

#### 设计方案: Bempp 风格 block/operator 外壳 + backend

```
PMCHWBlockOperator
├── EJ block operator
├── EM block operator
├── HJ block operator
├── HM block operator
├── weak_form()
├── strong_form()
└── backend binding
        ├── DenseBackend
        ├── MLFMABackend
        └── HMatrixBackend (deferred)
```

**backend 责任划分**：

```
顶层 shell:
- 保留四块语义
- 定义 weak_form / strong_form
- 负责 block algebra 与 mass-matrix-aware solve

DenseBackend:
- 为四块提供 dense block matvec / matrix

MLFMABackend:
- 为四块提供 nonlocal evaluator
- 近场 singular / near part 仍由 shell 显式组合

HMatrixBackend:
- 后续阶段再接入
```

**实现策略**：
- 阶段 1: 先完成 PMCHW block/operator shell + Dense backend
- 阶段 2: 把现有 `PMCHWMLFMAOperator` 收编成 MLFMABackend，并保留兼容 facade
- 阶段 3: 在统一 shell 下完成 weak/strong 对照与 Gate S
- 阶段 4: 未来再接入 HMatrixBackend

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

## 5. PMCHW MLFMA backend 收编方案

> **第四次修订**  
> Gibson Algorithm 14 仍然是本地 MLFMA backend 的理论依据；但它不再定义顶层 transmission 架构。  
> 顶层架构已切换为 **Bempp 风格 block/operator shell**，而“双八叉树 + 四遍远场”只作为 MLFMABackend 的内部实现保留。

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

### 5.1 文件结构（修订后主线）

```
src/
├── IntegralEquations/
│   └── PMCHW*.jl                  ← shell 层: 四块 block/operator 语义 + weak/strong form
├── Solvers/
│   └── ...                        ← shell 层强形式与 mass-matrix-aware solve 的调用点
└── FastAlgorithms/
    ├── MLFMA/
    │   └── PMCHWMLFMAOperator.jl  ← MLFMA backend / 兼容 facade
    └── HMatrix/                   ← 未来阶段再引入
```

**冻结规则**：
- 顶层 transmission 语义必须留在 shell 层，不再由单个 `PMCHWMLFMAOperator.jl` 独占。
- `PMCHWMLFMAOperator.jl` 允许继续存在，但角色改为 backend 内核或兼容 facade。
- H-matrix 目录/后端在本 Phase 不落地，只保留接口预留。

### 5.2 shell / backend 分工

| 层级 | 责任 | 当前阶段 |
|------|------|---------|
| PMCHW block/operator shell | 定义 `EJ/EM/HJ/HM` 四块、`weak_form()`、`strong_form()`、block algebra | **先实现** |
| Dense backend | 作为 shell 的第一条可验证后端，建立 direct/GMRES/strong-form 基线 | **先实现** |
| MLFMA backend | 收编现有双八叉树 + 四遍远场实现，提供 fast matvec | **第二阶段** |
| H-matrix backend | 参考 FreeFEM/Htool，提供压缩矩阵 backend | **后续阶段** |

### 5.3 MLFMA backend 保留内容

> 双波数设计仍然保留，但只服务于 backend，不再决定顶层接口。

- PMCHW 每块都仍是 `k0 + k1` 的叠加，因此 MLFMA backend 内部继续保留：
  - `octree0` for `k0`
  - `octree1` for `k1`
  - 四遍远场 `J×k0`, `J×k1`, `M×k0`, `M×k1`
- `assemble_near_field_pmchw` 仍然有效，但其产物应由 shell 统一接入四块 block 语义，而不是直接把 backend 暴露成顶层 transmission 对象。
- `aggS` 清零、leaf/translation 预算、near/far 划分等规则继续作为 backend 内部正确性合同存在。

### 5.4 当前 Phase 的接口目标

```julia
abstract type AbstractPMCHWBackend end

struct PMCHWBlockOperator{B<:AbstractPMCHWBackend}
    pmchw
    basis
    backend::B
end

weak_form(op::PMCHWBlockOperator)
strong_form(op::PMCHWBlockOperator)
```

**当前约束**：
- Dense backend 必须先在 shell 下跑通 B1/B2/Gate S。
- 现有 `PMCHWMLFMAOperator(...)` 若保留，只能作为构造 `PMCHWBlockOperator(..., MLFMABackend(...))` 的兼容入口。
- 禁止继续把新的 transmission 语义、strong-form 逻辑或 block algebra 直接塞回 `PMCHWMLFMAOperator.jl`。
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
| **15.7** | TDD: `test_pmchw_operator_shell.jl`（四块 block/operator + weak/strong form） | — | P1 🟠 |
| **15.8** | 实现 PMCHW block/operator shell + Dense backend | 15.7 RED | P1 🟠 |
| **15.9** | 在 shell 下完成 Gate S: `dense weak/strong` 对照 | 15.8 | P1 🟠 |
| **15.10** | TDD: `test_pmchw_mlfma_operator.jl`（backend 收编与 shell 绑定） | 15.8, 15.9 | P1 🟠 |
| **15.11** | 收编现有 `PMCHWMLFMAOperator` 为 MLFMA backend（双八叉树 + 四遍远场保留为内部实现） | 15.10 | P1 🟠 |
| **15.12** | 更新 `generate_report.jl` 加入 B1–B5 | 15.6 | P2 🟡 |
| **15.13** | 检视迭代 × 2 轮 | 所有 | P2 🟡 |
| **15.14** | TDD + 实现 `NMuller` dense baseline（alternative formulation） | 15.13 前可并行启动 | P1 🟠 |
| **15.15** | 同一 dielectric sphere 上比较 PMCHW / N-Muller conditioning 与 GMRES 行为 | 15.14 | P1 🟠 |

### 新架构亮点（与旧方案对比）

| 旧方案（已废弃） | 新方案（Bempp 风格 shell + backend） |
|----------------|------------------------------|
| transmission 顶层由单体 `PMCHWMLFMAOperator` 承载 | 顶层固定为 `EJ/EM/HJ/HM` block/operator shell |
| strong-form 无明确归属 | `strong_form` 固定在 shell 层 |
| Dense/MLFMA 共享接口边界不清晰 | Dense 先行，MLFMA 作为 backend 收编 |
| H-matrix 可能与当前主线混杂推进 | H-matrix 正式延后，待 shell 稳定后再引入 |
| Gibson 双八叉树决定整个架构 | Gibson 双八叉树仅决定 MLFMA backend 内部实现 |

---

## 7. DoD（完成定义）

| 检查项 | 通过标准 |
|--------|---------|
| 新增激励 API 单元测试 | `test_pmchw_excitation.jl` 全通过（含 zero-feed 错误捕获） |
| SCFIE delta-gap 测试 | `test_scfie_delta_gap.jl` 全通过 |
| B1 PMCHW Direct | Re(Z_in) > 0，εᵣ→1 时 Z_in 趋近 EFIE 结果（相对误差 <10%） |
| B2 PMCHW MLFMA | GMRES 收敛，Z_in_B2 与 Z_in_B1 误差 <5%（通过 shell 调用 MLFMA backend） |
| B3 VS-EFIE Direct | Z_in Re 误差 (vs EFIE-only) <10% |
| B4 VS-CFIE Direct | Z_in 与 B3 差异 <5Ω |
| B5 VS-CFIE MLFMA | Z_in 与 B4 差异 <5% |
| MagneticRWGBasis 设计 | ~~已废弃~~ — 不需要该类型 |
| PMCHW block/operator shell | 四块 `EJ/EM/HJ/HM` 语义、`weak_form`、`strong_form` 在 Dense backend 下全部可验证 |
| MLFMA backend 收编 | backend 经 shell 调用后，保持现有 B2 / Gate C 验收能力 |
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
- 实现主线: PMCHW shell 层 + Dense backend + MLFMA backend 绑定
- 兼容入口: `src/FastAlgorithms/MLFMA/PMCHWMLFMAOperator.jl`（backend/facade）

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
