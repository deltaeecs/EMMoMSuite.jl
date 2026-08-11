# EMMoMSuite 快速求解扩展实现计划（ACA → MLACA → MLFMA 加速）

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在现有 MLFMA 之上新增 ACA 与 MLACA 低秩压缩算子并与 GMRES 集成，用稠密参照做精度/内存/耗时/压缩率门控，随后补齐 MLFMA 加速项与基准报告。

**Architecture:** 复用现有 MLFMA 八叉树与近场稀疏装配；新增 `src/FastAlgorithms/ACA/` 模块（ACA 核心 → 块求值器 → `ACAOperator` → 后续 `MLACAOperator`）。所有算子实现 `AbstractIntegralOperator` + `mul!`，与 IterativeSolvers/GMRES 直接兼容。

**Tech Stack:** Julia 1.10+（本机 1.12.3）、IterativeSolvers、LinearAlgebra、SparseArrays、IncompleteLU、现有 FastAlgorithms/MLFMA、IntegralEquations。

---

## Task 1: M0 基线复跑与证据记录（已完成）

**Files:**
- Run: `benchmark/benchmark_mlfma_vs_direct.jl`
- Record: `docs/superpowers/specs/2026-08-11-fast-solvers-design.md` §1.1

- [x] **Step 1: 运行基线**

Run: `julia --project=. benchmark\benchmark_mlfma_vs_direct.jl`

Expected/实际输出（2026-08-11 本机）：
```
Unknowns: 792
Direct Assembly: 0.8262 s
MLFMA Setup: 5.5817 s
Max Difference (near vs direct): 1.01e-15
Relative Residual: 2.1831e-06
Relative Error (Norm): 7.83e-06
```

结论：近场精确、MLFMA 解与稠密一致（7.8e-6），作为后续 ACA/MLACA 的对比基线。

---

## Task 2: ACA 核心模块（TDD）

**Files:**
- Create: `src/FastAlgorithms/ACA/ACA.jl`
- Create: `test/test_aca.jl`
- Modify: `src/FastAlgorithms/FastAlgorithms.jl`（include + using + export）
- Modify: `src/EMMoMSuite.jl`（FastAlgorithms 导出区）
- Modify: `test/runtests.jl`（include test_aca.jl）

- [ ] **Step 1: 写失败测试（部分主元 ACA 重构低秩块）**

```julia
# test/test_aca.jl
using Test
using EMMoMSuite.FastAlgorithms.ACA
using LinearAlgebra
using Random

@testset "ACA core" begin
    Random.seed!(42)

    @testset "exact low-rank reconstruction" begin
        m, n, k0 = 60, 50, 4
        U0 = randn(ComplexF64, m, k0) .+ im .* randn(ComplexF64, m, k0) ./ 10
        V0 = randn(ComplexF64, n, k0) .+ im .* randn(ComplexF64, n, k0) ./ 10
        Z = U0 * V0'
        B = aca(Z; tol=1e-12)
        err = norm(Z - B.U * transpose(B.V)) / norm(Z)
        @test size(B.U, 2) <= k0
        @test err < 1e-10
    end

    @testset "Green's function well-separated block" begin
        m = n = 100
        rng = Random.MersenneTwister(7)
        src = 3 .* rand(rng, 3, n)
        dst = 3 .* rand(rng, 3, m) .+ 8.0
        k = 1.0
        Z = zeros(ComplexF64, m, n)
        for i in 1:m, j in 1:n
            R = norm(dst[:, i] - src[:, j])
            Z[i, j] = exp(-im * k * R) / (4π * R)
        end
        B = aca(Z; tol=1e-4)
        err = norm(Z - B.U * transpose(B.V)) / norm(Z)
        @test err < 5e-4
        @test size(B.U, 2) < 40
    end

    @testset "symmetric block transpose application" begin
        m = n = 40
        U0 = randn(ComplexF64, m, 3)
        V0 = randn(ComplexF64, n, 3)
        Z = U0 * V0'
        B = aca(Z; tol=1e-12)
        x = randn(ComplexF64, m)
        y1 = B.U * (transpose(B.V) * x)
        # 转置块应用：transpose(Z) * x == V * (transpose(U) * x)
        y2 = B.V * (transpose(B.U) * x)
        @test norm(y1 - Z * x) / norm(Z * x) < 1e-10
        @test norm(y2 - transpose(Z) * x) / norm(Z * x) < 1e-10
    end

    @testset "zero block early termination" begin
        Z = zeros(ComplexF64, 30, 30)
        B = aca(Z; tol=1e-4)
        @test size(B.U, 2) == 0
    end
end
```

- [ ] **Step 2: 运行测试确认失败**

Run: `julia --project=. -e 'using Test; include("test/test_aca.jl")'`

Expected: `ERROR: UndefVarError: aca not defined`（模块不存在）。

- [ ] **Step 3: 实现 ACA 核心**

```julia
# src/FastAlgorithms/ACA/ACA.jl
module ACA

using LinearAlgebra

export LowRankBlock, aca, recompress!, compression_stats

"""
    LowRankBlock{T}

低秩块：`Z ≈ U * transpose(V)`（转置约定，无共轭），其中 `U` 为 m×k、`V` 为 n×k。
"""
struct LowRankBlock{T}
    U::Matrix{T}
    V::Matrix{T}
end

Base.size(B::LowRankBlock) = (size(B.U, 1), size(B.V, 1))
Base.size(B::LowRankBlock, i::Int) = size(B)[i]

"""
    aca(Z::AbstractMatrix; tol=1e-4, maxrank=min(size(Z)...), recompress=true)

对稠密块 Z 做部分主元自适应交叉近似（Gibson Ch9 Algorithm 6，复数采用转置约定）。
"""
function aca(Z::AbstractMatrix{T}; tol::Real=1e-4, maxrank::Int=min(size(Z)...),
             recompress::Bool=true) where {T}
    m, n = size(Z)
    getrow(i) = Z[i, :]
    getcol(j) = Z[:, j]
    return aca(T, getrow, getcol, m, n; tol=tol, maxrank=maxrank, recompress=recompress)
end

"""
    aca(::Type{T}, getrow, getcol, m, n; tol=1e-4, maxrank=min(m,n), recompress=true)

自适应交叉近似：通过按行/列采样构造 `Z ≈ U*transpose(V)`，收敛判据为
`‖u_k‖‖v_k‖ ≤ tol * ‖Z̃‖_F`，其中 `‖Z̃‖_F²` 用递推估计（Algorithm 6 step 12）。
行残差 `R(I_k,:) = Z(I_k,:) − Σ_l (u_l)_{I_k} v_l`（无共轭）。
"""
function aca(::Type{T}, getrow::Function, getcol::Function, m::Int, n::Int;
             tol::Real=1e-4, maxrank::Int=min(m, n), recompress::Bool=true) where {T}
    # 1. 寻找首个非零行（Gibson 9.2.1.1）
    Iprev = 0
    for i in 1:m
        r = getrow(i)
        if norm(r) > 0
            Iprev = i
            break
        end
    end
    Iprev == 0 && return LowRankBlock(zeros(T, m, 0), zeros(T, n, 0))

    U = Matrix{T}(undef, m, 0)
    V = Matrix{T}(undef, n, 0)
    rows_used = falses(m)
    cols_used = falses(n)
    normZ2 = 0.0
    RT = real(T)
    eps_t = eps(RT)
    k = 0

    while k < maxrank
        # 6. 更新第 Iprev 行残差
        Rrow = k == 0 ? getrow(Iprev) : getrow(Iprev) .- (V * U[Iprev, :])

        # 7. 早期终止：剩余列残差近似为 0
        Jk = 0
        best = 0.0
        for j in 1:n
            cols_used[j] && continue
            a = abs(Rrow[j])
            if a > best
                best = a
                Jk = j
            end
        end
        best < eps_t * max(norm(Rrow), one(RT)) && break

        # 8-9. 主元列 + v_k
        pivot = Rrow[Jk]
        v = Rrow ./ pivot

        # 10-11. 列残差 + u_k
        Rcol = k == 0 ? getcol(Jk) : getcol(Jk) .- (U * V[Jk, :])
        u = Rcol

        # 12. Frobenius 范数递推估计
        s = 0.0
        for l in 1:k
            s += 2 * real(dot(U[:, l], u) * dot(V[:, l], v))
        end
        normZ2 += s + real(dot(u, u) * dot(v, v))

        U = hcat(U, u)
        V = hcat(V, v)
        rows_used[Iprev] = true
        cols_used[Jk] = true
        k += 1

        # 13. 收敛判据
        if norm(u) * norm(v) <= tol * sqrt(max(normZ2, 0.0))
            break
        end

        # 14. 下一主元行
        Iprev = 0
        best = 0.0
        for i in 1:m
            rows_used[i] && continue
            a = abs(Rcol[i])
            if a > best
                best = a
                Iprev = i
            end
        end
        Iprev == 0 && break
    end

    B = LowRankBlock(U, V)
    recompress && (B = recompress!(B; tol=10 * tol))
    return B
end

"""
    recompress!(B::LowRankBlock; tol=1e-3)

QR/SVD 再压缩（Gibson 9.2.2，转置约定）：`Z = Q_u R_u R_vᵀ Q_vᵀ`，对 `R_u*R_vᵀ` 做 SVD，
按 `σ ≥ tol*σ_max` 截断，通常可再获得 20-30% 额外压缩。
"""
function recompress!(B::LowRankBlock{T}; tol::Real=1e-3) where {T}
    U, V = B.U, B.V
    k = size(U, 2)
    k == 0 && return B
    Fu = qr(U)
    Fv = qr(V)
    Ru = Matrix(Fu.R)
    Rv = Matrix(Fv.R)
    S = Ru * transpose(Rv)
    F = svd(S)
    σmax = F.S[1]
    r = count(s -> s >= tol * σmax, F.S)
    r < 1 && (r = 1)
    U_new = Matrix(Fu.Q) * (F.U[:, 1:r] * Diagonal(F.S[1:r]))
    V_new = Matrix(Fv.Q) * conj(F.V[:, 1:r])
    return LowRankBlock(U_new, V_new)
end

"""
    compression_stats(B, m, n) -> (rank, ratio)

压缩统计：`ratio = 1 - k*(m+n)/(m*n)`。
"""
function compression_stats(B::LowRankBlock, m::Int, n::Int)
    k = size(B.U, 2)
    return k, 1 - k * (m + n) / (m * n)
end

end # module ACA
```

- [ ] **Step 4: 运行测试确认通过**

Run: `julia --project=. -e 'using Test; include("test/test_aca.jl")'`

Expected: `Test Summary: | Pass ...` 全绿（4 个 testset）。

- [ ] **Step 5: 注册模块并导出**

```julia
# src/FastAlgorithms/FastAlgorithms.jl
include("ACA/ACA.jl")
using .ACA
export ACA, LowRankBlock, aca, recompress!
```

并在 `src/EMMoMSuite.jl` 的 FastAlgorithms 导出区补充 `LowRankBlock, aca, recompress!`。

- [ ] **Step 6: 回归**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`（或先跑 `test/runtests.jl` 的快照）。

Expected: 既有测试全绿，新增 test_aca.jl 通过。

- [ ] **Step 7: 提交**

```bash
git add src/FastAlgorithms/ACA/ACA.jl src/FastAlgorithms/FastAlgorithms.jl src/EMMoMSuite.jl test/test_aca.jl test/runtests.jl
git commit -m "feat(ACA): partial-pivoting ACA with Frobenius estimate and QR/SVD recompression"
```

---

## Task 3: 块求值器（TDD）

**Files:**
- Create: `src/FastAlgorithms/ACA/BlockEvaluator.jl`
- Create: `test/test_block_evaluator.jl`
- Modify: `src/FastAlgorithms/FastAlgorithms.jl`、`test/runtests.jl`

- [ ] **Step 1: 写失败测试**

```julia
# test/test_block_evaluator.jl
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.ACA
using LinearAlgebra

@testset "BlockEvaluator" begin
    mesh = generate_sphere_mesh(0.5, 6, 10)
    basis = RWGBasis(mesh)
    efie = EFIE(300e6)
    Z = assemble_impedance_matrix(efie, basis)
    ev = BlockEvaluator(efie, basis)
    rows = [1, 5, 7, 42, 100]
    cols = [2, 3, 8, 99]
    B = eval_block(ev, rows, cols)
    @test size(B) == (length(rows), length(cols))
    for (ii, i) in enumerate(rows), (jj, j) in enumerate(cols)
        @test isapprox(B[ii, jj], Z[i, j]; atol=1e-10)
    end
end
```

- [ ] **Step 2: 运行确认失败**

Expected: `ERROR: UndefVarError: BlockEvaluator not defined`。

- [ ] **Step 3: 实现**

```julia
# src/FastAlgorithms/ACA/BlockEvaluator.jl
module BlockEvaluatorModule

using ...IntegralEquations
using ...IntegralEquations.Impedance: get_triangles_info
using ...IntegralEquations.EFIEModule: efie_interaction!, EFIE
using ...IntegralEquations.MFIEModule: mfie_interaction!, MFIE
using ...IntegralEquations.CFIEModule: CFIE
using ...BasisFunctions: RWGBasis

export BlockEvaluator, eval_block

"""
    BlockEvaluator(op, basis)

预计算三角形信息与基函数→三角形映射，供 `eval_block` 按行/列集合求值阻抗块。
当前支持 RWG 单基函数表面算子（EFIE/MFIE/CFIE）。
"""
struct BlockEvaluator{OP,BF}
    op::OP
    basis::BF
    tri_info::Vector
    tri_to_rwg::Vector{Vector{Tuple{Int,Int,Float64}}}
    basis_tris::Vector{Vector{Int}}
    basis_local::Vector{Vector{Tuple{Int,Int,Float64}}}  # (tri_id, local_edge, sign)
end

function BlockEvaluator(op::AbstractIntegralOperator, basis::RWGBasis)
    nt = num_elements(basis.mesh)
    tri_to_rwg = [Vector{Tuple{Int,Int,Float64}}() for _ in 1:nt]
    basis_tris = [Int[] for _ in 1:num_basis(basis)]
    basis_local = [Tuple{Int,Int,Float64}[] for _ in 1:num_basis(basis)]
    for (i, f) in enumerate(basis.functions)
        for k in 1:2
            t = f.support[k]
            if t > 0
                push!(tri_to_rwg[t], (f.local_edge_idx[k], i, f.signs[k]))
                push!(basis_tris[i], t)
                push!(basis_local[i], (t, f.local_edge_idx[k], f.signs[k]))
            end
        end
    end
    tri_info = get_triangles_info(basis.mesh, basis)
    return BlockEvaluator(op, basis, tri_info, tri_to_rwg, basis_tris, basis_local)
end

"""
    eval_block(ev::BlockEvaluator, rows::Vector{Int}, cols::Vector{Int}) -> Matrix

求值 `Z[rows, cols]`（全局基函数索引）。按三角形对去重计算后分发，支持
EFIE/MFIE/CFIE 组合（含 EFIE 近场规范序转置修复）。
"""
function eval_block(ev::BlockEvaluator, rows::Vector{Int}, cols::Vector{Int})
    op = ev.op
    m, n = length(rows), length(cols)
    Z = zeros(ComplexF64, m, n)
    rowmap = Dict(r => i for (i, r) in enumerate(rows))
    colmap = Dict(c => j for (j, c) in enumerate(cols))

    row_tris = unique(vcat([ev.basis_tris[i] for i in rows]...))
    col_tris = unique(vcat([ev.basis_tris[j] for j in cols]...))

    Zloc = zeros(ComplexF64, 3, 3)
    needs_mfie = op isa CFIE
    Zmfie = needs_mfie ? zeros(ComplexF64, 3, 3) : Zloc

    for t1 in row_tris
        for t2 in col_tris
            fill!(Zloc, 0)
            if op isa CFIE
                efie_interaction!(Zloc, op.efie, ev.tri_info[t1], ev.tri_info[t2])
                fill!(Zmfie, 0)
                mfie_interaction!(Zmfie, op.mfie, ev.tri_info[t1], ev.tri_info[t2])
                @inbounds for q in eachindex(Zloc)
                    Zloc[q] = op.alpha * Zloc[q] + (1 - op.alpha) * Zmfie[q]
                end
            elseif op isa MFIE
                mfie_interaction!(Zloc, op, ev.tri_info[t1], ev.tri_info[t2])
            else
                efie_interaction!(Zloc, op, ev.tri_info[t1], ev.tri_info[t2])
            end
            # 分发到 (rows, cols) 块
            for (lt1, g1, s1) in ev.tri_to_rwg[t1]
                haskey(rowmap, g1) || continue
                for (lt2, g2, s2) in ev.tri_to_rwg[t2]
                    haskey(colmap, g2) || continue
                    val = Zloc[lt1, lt2] * s1 * s2
                    Z[rowmap[g1], colmap[g2]] += val
                end
            end
        end
    end
    return Z
end

end # module BlockEvaluatorModule
```

- [ ] **Step 4: 运行确认通过**（`julia --project=. -e 'using Test; include("test/test_block_evaluator.jl")'`）

- [ ] **Step 5: 注册模块并导出**

```julia
# src/FastAlgorithms/FastAlgorithms.jl（追加）
include("ACA/BlockEvaluator.jl")
using .BlockEvaluatorModule
export ACA, LowRankBlock, aca, recompress!, BlockEvaluator, eval_block
```

并在 `test/runtests.jl` 中加入 `include("test_aca.jl")`、`include("test_block_evaluator.jl")`。

- [ ] **Step 6: 回归**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`

Expected: 既有测试全绿，新增 block evaluator 测试通过。

- [ ] **Step 7: 提交**

```bash
git add src/FastAlgorithms/ACA/BlockEvaluator.jl test/test_block_evaluator.jl src/FastAlgorithms/FastAlgorithms.jl test/runtests.jl
git commit -m "feat(ACA): block evaluator for EFIE/MFIE/CFIE RWG blocks"
```

---

## Task 4: ACAOperator（TDD）

**Files:**
- Create: `src/FastAlgorithms/ACA/ACAOperator.jl`
- Create: `test/test_aca_operator.jl`
- Modify: `src/FastAlgorithms/FastAlgorithms.jl`、`src/EMMoMSuite.jl`、`test/runtests.jl`

- [ ] **Step 1: 写失败测试**

```julia
# test/test_aca_operator.jl
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.ACA
using EMMoMSuite.FastAlgorithms.ACAOperatorModule: ACAOperator
using LinearAlgebra
using IterativeSolvers

@testset "ACAOperator" begin
    mesh = generate_sphere_mesh(0.5, 5, 8)
    basis = RWGBasis(mesh)
    efie = EFIE(300e6)
    Z = assemble_impedance_matrix(efie, basis)
    N = size(Z, 1)

    leaf = 0.5 * (299792458.0 / 300e6)
    op = ACAOperator(efie, basis, leaf; tol=1e-4)

    @test size(op) == (N, N)
    @test !isempty(op.blocks)

    x = randn(ComplexF64, N)
    y = op * x
    err = norm(y - Z * x) / norm(Z * x)
    @test err < 1e-2   # 远场 ACA 容差 1e-4 下的整体 MatVec 误差（门控按实测收紧/放宽）

    src = PlaneWave(300e6, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(efie, src, basis)
    I_aca = gmres(op, V; abstol=1e-6, reltol=1e-8, maxiter=200, restart=50)
    I_dir = Z \ V
    @test norm(I_aca - I_dir) / norm(I_dir) < 1e-2
end
```

- [ ] **Step 2: 运行确认失败**（`UndefVarError: ACAOperator not defined`）

- [ ] **Step 3: 实现 ACAOperator**

```julia
# src/FastAlgorithms/ACA/ACAOperator.jl
module ACAOperatorModule

using LinearAlgebra
using SparseArrays
using ...CoreModule
using ...Geometry
using ...BasisFunctions
using ...IntegralEquations
using ..Octree
using ..OctreeBuilder
using ..MLFMAOperatorModule: assemble_near_field
using ..ACA: LowRankBlock, aca
using ..BlockEvaluatorModule: BlockEvaluator, eval_block

export ACAOperator, ACAParams, get_leaf_intervals

struct ACAParams
    tol::Float64
    maxrank::Int
    recompress::Bool
    symmetric::Bool
end
ACAParams(; tol=1e-4, maxrank=512, recompress=true, symmetric=true) =
    ACAParams(tol, maxrank, recompress, symmetric)

struct ACABlock{CT}
    rows::Vector{Int}   # 全局行索引
    cols::Vector{Int}   # 全局列索引
    U::Matrix{CT}
    V::Matrix{CT}
end

struct ACAOperator{FT,CT} <: AbstractIntegralOperator
    octree::OctreeInfo
    bases::Vector{AbstractBasisFunction}
    basis_offsets::Vector{Int}
    Z_near::SparseMatrixCSC{CT,Int}
    blocks::Vector{ACABlock{CT}}
    operator::AbstractIntegralOperator
    sorted_ids::Vector{Int}
    inv_sorted_ids::Vector{Int}
    params::ACAParams
end

Base.eltype(::ACAOperator{FT,CT}) where {FT,CT} = CT
Base.size(A::ACAOperator) = size(A.Z_near)
Base.size(A::ACAOperator, i::Int) = size(A.Z_near, i)

function ACAOperator(operator::AbstractIntegralOperator, basis::AbstractBasisFunction,
                     leafCubeEdgel::Float64; tol::Real=1e-4, maxrank::Int=512,
                     recompress::Bool=true, symmetric::Bool=true,
                     near_range::Int=4, interp_method::Val=Val(:Lagrange2Step))
    return ACAOperator(operator, [basis], leafCubeEdgel; tol=tol, maxrank=maxrank,
                       recompress=recompress, symmetric=symmetric,
                       near_range=near_range, interp_method=interp_method)
end

function ACAOperator(operator::AbstractIntegralOperator,
                     bases::Vector{<:AbstractBasisFunction}, leafCubeEdgel::Float64;
                     tol::Real=1e-4, maxrank::Int=512, recompress::Bool=true,
                     symmetric::Bool=true, near_range::Int=4,
                     interp_method::Val=Val(:Lagrange2Step))
    params = ACAParams(; tol=Float64(tol), maxrank=maxrank, recompress=recompress,
                       symmetric=symmetric)
    bf_centers_list = [reduce(hcat, [bf.center for bf in b.functions]) for b in bases]
    bf_centers = reduce(hcat, bf_centers_list)
    λ = Constants.c0 / operator.freq
    octree, sorted_ids = build_octree(bf_centers, leafCubeEdgel; λ=λ,
                                      interp_method=interp_method, near_range=near_range)
    N = length(sorted_ids)
    inv_sorted_ids = zeros(Int, N)
    for i in 1:N
        inv_sorted_ids[sorted_ids[i]] = i
    end
    basis_offsets = cumsum([num_basis(b) for b in bases])
    abstract_bases = Vector{AbstractBasisFunction}(bases)
    Z_near = assemble_near_field(operator, abstract_bases, basis_offsets,
                                 octree, sorted_ids, inv_sorted_ids)

    # 叶层非邻盒子对 → ACA 压缩
    leaf_level = octree.levels[octree.nLevels]
    cubes = leaf_level.cubes
    n_cubes = length(cubes)
    ev = BlockEvaluator(operator, bases[1])   # 当前支持单 RWG 基函数
    blocks = Vector{ACABlock{ComplexF64}}()
    for i in 1:n_cubes, j in (i+1):n_cubes
        isempty(cubes[i].bfInterval) && continue
        isempty(cubes[j].bfInterval) && continue
        j in cubes[i].neighbors && continue
        rows = sorted_ids[cubes[i].bfInterval]
        cols = sorted_ids[cubes[j].bfInterval]
        (length(rows) <= 1 || length(cols) <= 1) && continue
        B = aca((r) -> vec(eval_block(ev, [rows[r]], cols)),
                (c) -> vec(eval_block(ev, rows, [cols[c]])),
                length(rows), length(cols);
                tol=params.tol, maxrank=params.maxrank, recompress=params.recompress)
        size(B.U, 2) > 0 && push!(blocks, ACABlock(rows, cols, B.U, B.V))
    end

    FT = eltype(bases[1].mesh.node)
    CT = ComplexF64
    return ACAOperator{FT,CT}(octree, abstract_bases, basis_offsets, Z_near,
                              blocks, operator, sorted_ids, inv_sorted_ids, params)
end

function get_leaf_intervals(A::ACAOperator)
    leaf_level = A.octree.levels[A.octree.nLevels]
    return [cube.bfInterval for cube in leaf_level.cubes]
end

function Base.:*(A::ACAOperator, x::AbstractVector)
    y = similar(x)
    mul!(y, A, x)
    return y
end

function LinearAlgebra.mul!(y::AbstractVector, A::ACAOperator, x::AbstractVector)
    mul!(y, A.Z_near, x)
    for blk in A.blocks
        xc = view(x, blk.cols)
        t = blk.U * (transpose(blk.V) * xc)
        @inbounds for (idx, r) in enumerate(blk.rows)
            y[r] += t[idx]
        end
        if A.params.symmetric
            xi = view(x, blk.rows)
            t2 = blk.V * (transpose(blk.U) * xi)
            @inbounds for (idx, c) in enumerate(blk.cols)
                y[c] += t2[idx]
            end
        end
    end
    return y
end

end # module ACAOperatorModule
```

注意：`build_octree` 的近邻范围与 `assemble_near_field` 必须一致（都用 `near_range`），否则近/远场会重叠或漏算——测试中会显式核对 `Z_near + Σ 远场块 ≈ Z_dense`。

- [ ] **Step 4: 注册模块并导出**

```julia
# src/FastAlgorithms/FastAlgorithms.jl（追加）
include("ACA/ACAOperator.jl")
using .ACAOperatorModule
export ACA, LowRankBlock, aca, recompress!, BlockEvaluator, eval_block, ACAOperator, ACAParams
```

并在 `src/EMMoMSuite.jl` 的 FastAlgorithms 导出区补充 `ACAOperator`；在 `test/runtests.jl` 加入 `include("test_aca_operator.jl")`。

- [ ] **Step 5: 运行确认通过**

Expected: `test_aca_operator.jl` 全绿；若整体 MatVec 误差 > 1e-2，检查近/远场划分一致性或调低 `tol`（门控以稠密参照为准，不放松参照）。

- [ ] **Step 6: 回归 + 提交**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
git add src/FastAlgorithms/ACA/ACAOperator.jl test/test_aca_operator.jl src/FastAlgorithms/FastAlgorithms.jl src/EMMoMSuite.jl test/runtests.jl
git commit -m "feat(ACA): ACAOperator with octree clustering, near-field reuse, and GMRES integration"
```

---

## Task 5: ACA 基准与报告（M0 对比）

**Files:**
- Create: `benchmark/benchmark_aca_vs_direct.jl`（复用 `benchmark_mlfma_vs_direct.jl` 的球体夹具）
- Create: `benchmark/results/aca_benchmark_<date>.csv`

- [ ] **Step 1: 写基准脚本**（参数：freq、radius、leaf、tol=1e-4；输出 N、装配/建立/求解时间、MatVec 误差、相对残差、解误差、压缩率、内存估计）
- [ ] **Step 2: 运行并记录 CSV**
- [ ] **Step 3: 提交** `git add benchmark/benchmark_aca_vs_direct.jl benchmark/results/... && git commit -m "bench(ACA): ACA vs direct baseline report"`

---

## Task 6+: 后续里程碑（另行计划文档细化）

> 状态：M3/M4/M5 与后续方向（非对称 MLACA、PMCHW、直接块 LU、N=1.1 万
> 本地大规模基准）均已完成，见设计文档 §2.4 与 docs/theory/fast_algorithms.md §7-9。

### M3 MLACAOperator ✅
- 八叉树多层簇结构上的 H-矩阵风格递归块压缩：可容许远块 ACA、近块下钻、叶层稠密（`Z_near`）。
- MatVec 递归应用；对照稠密 / ACAOperator / MLFMA；N=11352 压缩率 71.5%（EFIE）。

### M4 MLFMA 加速 ✅
- 预条件接线：`ILUPreconditioner(op)` / `SPAIPreconditioner(op)` / `BlockJacobiPreconditioner(op)`
  在 MLFMA/ACA/MLACA 算子上的迭代数基准（Identity/ILU/SPAI/BlockJacobi 扫描）。
- 参数配置化：`nInterp`、`precision_digits` 构造参数（默认保持现状），MPI 构造器同步支持。
- 低频：ACAOperator 覆盖低频场景（30 MHz、0.03λ 叶层），与 MLFMA 互补。

### M5 综合基准与文档 ✅
- `benchmark/run_full_fast_solvers_benchmark.jl`、`benchmark/run_large_fast_solvers_benchmark.jl`
  （N=792..11352，多配方，finite/NaN 检查，本地手动运行、不进 CI）。
- 报告：精度/内存/耗时/压缩率分列 CSV + docs 表格。
- `docs/src/theory/fast_algorithms.md`、API 文档更新。

### 后续方向（已完成）
- 非对称 MLACA（CFIE）✅；PMCHW 2N 多基函数 ✅；直接块 LU（多 RHS）✅。
