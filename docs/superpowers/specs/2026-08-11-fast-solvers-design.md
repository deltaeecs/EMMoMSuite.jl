# EMMoMSuite 快速求解扩展设计文档（ACA / MLACA / MLFMA 加速）

**日期：** 2026-08-11
**状态：** 待用户确认（brainstorming 设计门）
**目标：** 在现有 MLFMA 之上，新增 ACA 与 MLACA 低秩压缩算子并与 GMRES 集成，以稠密参照做精度/内存/耗时/压缩率门控；随后按实测补齐 MLFMA 加速项（预条件接线、参数配置化、低频 ACA 互补）。

---

## 1. 背景与现状证据

### 1.1 仓库现状（2026-08-11 取证）

- `EMMoMSuite` v0.1.2（Julia，分支 `codex/rel-v0.1.2`，工作树干净）。
- 已有完整 MLFMA：`src/FastAlgorithms/MLFMA/`（Octree、Lebedev 插值、聚合/转移/解聚、`MLFMAOperator`、`PMCHWMLFMAOperator`、MPI 版）。
- `MLFMAOperator <: AbstractIntegralOperator`，实现 `mul!`，可直接进 GMRES；近场稀疏装配 `assemble_near_field` 可复用。
- `src/Core/Interfaces.jl` 只有 `AbstractFastAlgorithm` 占位；ACA/MLACA 无实现。
- `docs/src/theory/fast_algorithms.md` §3 已写明"低频 MLFMA 与 ACA"路线（ACA 为核无关代数压缩，适合复杂介质核与低频）。
- 预条件器已存在：`DiagonalPreconditioner`、`ILUPreconditioner`、`SPAIPreconditioner`、`BlockJacobiPreconditioner`（`src/Solvers/Preconditioners.jl`），M4 重点是"验证 + 接线"而非从零实现。
- 硬编码参数：`src/FastAlgorithms/MLFMA/Precomputations.jl:73` `nInterp = 6`（截断精度等亦为常量），需配置化。
- M0 基线（本机实测，N=792 EFIE 球体 @300MHz）：稠密装配 0.826 s；MLFMA 建立 5.58 s；近场与稠密逐元素最大差 1.01e-15；ILU(0.01) 预条件 GMRES 2.67 s；相对残差 2.18e-6；解相对误差 7.83e-6。

### 1.2 CEM 知识库依据

- Gibson《Method of Moments in Electromagnetics, 3rd》Ch9：部分主元 ACA（Algorithm 6）、Frobenius 范数递推估计（式 9.4 及 Algorithm 6 step 12-13）、早期终止（step 7）、QR/SVD 再压缩（τ_SVD ≈ 10·τ_ACA，额外压缩 20-30%）、聚类（K-means/octree）、单层分组近/远分离。
- Gibson Ch10：MLACA 递归子块压缩（L 层将分组再分为 2^L 子组，逐级对 U 块再压缩）、对角块最细层稠密、U/V 型用于直接 LU；参考压缩数据（PEC 球/再入体，L=0..5）。
- Ergül & Gürel《MLFMA》Ch2/Ch4：MLFMA 积分方程背景与并行化。

## 2. 设计决策

### 2.1 总体路线：方案 A（仓库内自研，先迭代后直接）

1. **ACA 核心模块**：与算子/基函数解耦的通用低秩压缩（`aca(getrow, getcol, m, n; tol, maxrank)`），返回 `LowRankBlock(U, V)`（Z ≈ U·Vᵀ，转置约定无共轭，已数值验证）；可选 QR/SVD 再压缩。
2. **块求值器**：`BlockEvaluator(op, basis)` 预计算三角形/四面体信息与基函数映射，提供 `eval_block(rows, cols)`，供 ACA 按行/列采样，避免装配整块。
3. **ACAOperator**：复用 MLFMA 八叉树聚类与 `assemble_near_field`；叶层非邻盒子对按 ACA 压缩为低秩块；实现 `mul!` + GMRES。复数对称矩阵用转置语义（`V*(Uᵀ*x)`）应用 (j,i) 块，避免重复压缩。
4. **MLACAOperator**（M3）：在八叉树多层结构上做 H-矩阵风格递归块压缩（可容许远块 ACA、近块下钻、叶层稠密），实现 `mul!` + GMRES。
5. **MLFMA 加速**（M4）：预条件器（ILU/SPAI/BlockJacobi）与 MLFMA/ACA 算子接线并基准验证；截断/插值参数配置化；低频由 ACA 互补。
6. **基准/报告/文档**（M5）：稠密 vs MLFMA vs ACA vs MLACA 的 MatVec 误差、求解残差、内存、耗时、压缩率；`benchmark/` 脚本 + CSV 报告；docs/theory 与 API 刷新。

### 2.2 关键接口

```julia
struct LowRankBlock{T}
    U::Matrix{T}   # m×k
    V::Matrix{T}   # n×k，Z ≈ U * transpose(V)
end

aca(getrow, getcol, m, n; tol=1e-4, maxrank=min(m,n), recompress=true) -> LowRankBlock
recompress!(B::LowRankBlock; tol=10*tol_aca) -> LowRankBlock

struct ACAOperator{FT,CT} <: AbstractIntegralOperator
    octree, bases, basis_offsets, Z_near, blocks, operator, sorted_ids, inv_sorted_ids, params
end
LinearAlgebra.mul!(y, A::ACAOperator, x)   # 近场稀疏 + 远场低秩块
```

### 2.3 精度与门控

- ACA 容差默认 `tol=1e-4`（Gibson 常见取值），QR/SVD 再压缩容差 `10·tol`。
- 门控：`‖Z_block − U·V′‖_F / ‖Z_block‖_F ≤ 5·tol`（测试）；整体 MatVec 相对误差 ≤ 1e-3；GMRES 相对残差 ≤ 1e-5 或按问题规模放宽；内存/耗时/压缩率分开报告。
- 近场一律复用现有精确 `assemble_near_field`，不引入 fallback 或启发式掩膜。

### 2.4 不做（本期范围外）

- Gibson U/V 型 MLACA 直接块 LU（多 RHS 直接求解）→ 后续计划。
- K-means 聚类（复用八叉树即可满足聚类目标）。
- 外部 H-matrix 库引入。
