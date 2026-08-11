module ACAOperatorModule

using LinearAlgebra
using SparseArrays
using ...CoreModule
using ...CoreModule: Constants
using ...Geometry
using ...BasisFunctions
using ...IntegralEquations
using ...IntegralEquations.PMCHWModule: PMCHW
using ..MLFMA: OctreeInfo, build_octree
import ..MLFMA: get_leaf_intervals
using ..MLFMA.MLFMAOperatorModule: assemble_near_field
using ..MLFMA.PMCHWMLFMAOperatorModule: assemble_near_field_pmchw
using ..ACA: LowRankBlock, aca
using ..BlockEvaluatorModule: BlockEvaluator, PMCHWBlockEvaluator, eval_block

export ACAOperator, ACAParams

"""
    ACAParams

ACA 算子参数：
- `tol`：ACA 压缩容差（默认 `1e-4`）。
- `maxrank`：单块最大秩（默认 `512`）。
- `recompress`：是否执行 QR/SVD 再压缩（默认 `true`，`τ_SVD = 10*tol`）。
- `symmetric`：矩阵是否复数对称（默认 `true`）。对称时只压缩 (i,j) 块并用
  转置语义应用 (j,i) 块；非对称时两个方向分别压缩。
"""
struct ACAParams
    tol::Float64
    maxrank::Int
    recompress::Bool
    symmetric::Bool
end

ACAParams(; tol::Real = 1e-4, maxrank::Int = 512, recompress::Bool = true, symmetric::Bool = true) =
    ACAParams(Float64(tol), maxrank, recompress, symmetric)

"""
    ACABlock{CT}

叶层非邻盒子对的低秩块：`Z[rows, cols] ≈ U * transpose(V)`。
`rows`/`cols` 为全局基函数索引。
"""
struct ACABlock{CT}
    rows::Vector{Int}
    cols::Vector{Int}
    U::Matrix{CT}
    V::Matrix{CT}
end

"""
    ACAOperator{FT,CT} <: AbstractIntegralOperator

自适应交叉近似（ACA）算子：复用 MLFMA 八叉树聚类与近场稀疏装配，
叶层非邻盒子对（远场）按 ACA 压缩为低秩块，实现 `mul!` 后可直接与
IterativeSolvers（GMRES 等）配合。

# 构造
    ACAOperator(operator, basis, leafCubeEdgel; tol=1e-4, maxrank=512,
                recompress=true, symmetric=true, near_range=1,
                interp_method=Val(:Lagrange2Step))

# 参数
- `operator`：EFIE/MFIE/CFIE（当前支持 RWG 单基函数表面算子）。
- `basis`：RWG 基函数。
- `leafCubeEdgel`：叶层盒子边长（米）。
- `near_range`：近场 Chebyshev 邻域半径（默认 1；小算例下远大于 1 会使全部
  块成为近场、无远场压缩）。

# 参考
- Gibson, The Method of Moments in Electromagnetics, 3rd, Ch9。
"""
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

function ACAOperator(
    operator::AbstractIntegralOperator,
    basis::AbstractBasisFunction,
    leafCubeEdgel::Float64;
    tol::Real = 1e-4,
    maxrank::Int = 512,
    recompress::Bool = true,
    symmetric::Bool = true,
    near_range::Int = 1,
    interp_method::Val = Val(:Lagrange2Step),
    nInterp::Int = 6,
    precision_digits::Real = 9.0,
)
    return ACAOperator(
        operator,
        [basis],
        leafCubeEdgel;
        tol = tol,
        maxrank = maxrank,
        recompress = recompress,
        symmetric = symmetric,
        near_range = near_range,
        interp_method = interp_method,
        nInterp = nInterp,
        precision_digits = precision_digits,
    )
end

function ACAOperator(
    operator::AbstractIntegralOperator,
    bases::Vector{<:AbstractBasisFunction},
    leafCubeEdgel::Float64;
    tol::Real = 1e-4,
    maxrank::Int = 512,
    recompress::Bool = true,
    symmetric::Bool = true,
    near_range::Int = 1,
    interp_method::Val = Val(:Lagrange2Step),
    nInterp::Int = 6,
    precision_digits::Real = 9.0,
)
    params = ACAParams(; tol = tol, maxrank = maxrank, recompress = recompress, symmetric = symmetric)

    # 1. 八叉树聚类（复用 MLFMA build_octree）
    bf_centers_list = [reduce(hcat, [bf.center for bf in b.functions]) for b in bases]
    bf_centers = reduce(hcat, bf_centers_list)
    λ = Constants.c0 / operator.freq
    octree, sorted_ids = build_octree(
        bf_centers,
        leafCubeEdgel;
        λ = λ,
        interp_method = interp_method,
        near_range = near_range,
        nInterp = nInterp,
        precision_digits = precision_digits,
    )
    N = length(sorted_ids)
    inv_sorted_ids = zeros(Int, N)
    for i in 1:N
        inv_sorted_ids[sorted_ids[i]] = i
    end
    basis_offsets = cumsum([num_basis(b) for b in bases])
    abstract_bases = Vector{AbstractBasisFunction}(bases)

    # 2. 近场稀疏装配（复用 MLFMA assemble_near_field）
    Z_near = assemble_near_field(
        operator,
        abstract_bases,
        basis_offsets,
        octree,
        sorted_ids,
        inv_sorted_ids,
    )

    # 3. 远场：叶层非邻盒子对 → ACA 压缩
    leaf_level = octree.levels[octree.nLevels]
    cubes = leaf_level.cubes
    n_cubes = length(cubes)
    ev = BlockEvaluator(operator, bases[1]) # 当前支持单 RWG 基函数
    blocks = Vector{ACABlock{ComplexF64}}()

    for i in 1:n_cubes
        isempty(cubes[i].bfInterval) && continue
        for j in (i + 1):n_cubes
            isempty(cubes[j].bfInterval) && continue
            j in cubes[i].neighbors && continue

            rows_i = sorted_ids[cubes[i].bfInterval]
            cols_j = sorted_ids[cubes[j].bfInterval]

            # 块 (i, j)
            B = aca(
                ComplexF64,
                (r) -> vec(eval_block(ev, [rows_i[r]], collect(cols_j))),
                (c) -> vec(eval_block(ev, collect(rows_i), [cols_j[c]])),
                length(rows_i),
                length(cols_j);
                tol = params.tol,
                maxrank = params.maxrank,
                recompress = params.recompress,
            )
            size(B.U, 2) > 0 &&
                push!(blocks, ACABlock(collect(rows_i), collect(cols_j), B.U, B.V))

            # 非对称：块 (j, i) 单独压缩
            if !params.symmetric
                B2 = aca(
                    ComplexF64,
                    (r) -> vec(eval_block(ev, [cols_j[r]], collect(rows_i))),
                    (c) -> vec(eval_block(ev, collect(cols_j), [rows_i[c]])),
                    length(cols_j),
                    length(rows_i);
                    tol = params.tol,
                    maxrank = params.maxrank,
                    recompress = params.recompress,
                )
                size(B2.U, 2) > 0 &&
                    push!(blocks, ACABlock(collect(cols_j), collect(rows_i), B2.U, B2.V))
            end
        end
    end

    FT = eltype(bases[1].mesh.node)
    CT = eltype(Z_near)
    return ACAOperator{FT,CT}(
        octree,
        abstract_bases,
        basis_offsets,
        Z_near,
        blocks,
        operator,
        sorted_ids,
        inv_sorted_ids,
        params,
    )
end

"""
    ACAOperator(pmchw::PMCHW, basis::RWGBasis, leafCubeEdgel; tol=1e-4, ...)

PMCHW 系统的 ACA 算子（2N×2N）：复用 MLFMA 的 `assemble_near_field_pmchw`
近场装配与 `PMCHWBlockEvaluator` 远场块求值。PMCHW 非对称（Z^HJ = -Z^EM），
对每个叶层非邻盒子对的两个方向分别压缩。
"""
function ACAOperator(
    pmchw::PMCHW,
    basis::RWGBasis,
    leafCubeEdgel::Float64;
    tol::Real = 1e-4,
    maxrank::Int = 512,
    recompress::Bool = true,
    symmetric::Bool = false,
    near_range::Int = 1,
    interp_method::Val = Val(:Lagrange2Step),
    nInterp::Int = 6,
    precision_digits::Real = 9.0,
)
    symmetric && error("PMCHW 系统非对称（Z^HJ = -Z^EM），ACAOperator 必须使用 symmetric=false")
    params = ACAParams(; tol = tol, maxrank = maxrank, recompress = recompress, symmetric = false)
    N = num_basis(basis)

    bf_centers = reduce(hcat, [bf.center for bf in basis.functions])
    λ = Constants.c0 / pmchw.freq
    octree, sorted_ids = build_octree(
        bf_centers,
        leafCubeEdgel;
        λ = λ,
        interp_method = interp_method,
        near_range = near_range,
        nInterp = nInterp,
        precision_digits = precision_digits,
    )
    inv_sorted_ids = zeros(Int, N)
    for i in 1:N
        inv_sorted_ids[sorted_ids[i]] = i
    end
    basis_offsets = cumsum([num_basis(basis)])
    abstract_bases = Vector{AbstractBasisFunction}([basis])

    Z_near = assemble_near_field_pmchw(pmchw, basis, octree, sorted_ids, inv_sorted_ids)
    ev = PMCHWBlockEvaluator(pmchw, basis)

    leaf_level = octree.levels[octree.nLevels]
    cubes = leaf_level.cubes
    n_cubes = length(cubes)
    blocks = Vector{ACABlock{ComplexF64}}()

    for i in 1:n_cubes
        isempty(cubes[i].bfInterval) && continue
        for j in (i + 1):n_cubes
            isempty(cubes[j].bfInterval) && continue
            j in cubes[i].neighbors && continue

            rowsJ = sorted_ids[cubes[i].bfInterval]
            colsJ = sorted_ids[cubes[j].bfInterval]
            rows = vcat(rowsJ, N .+ rowsJ)
            cols = vcat(colsJ, N .+ colsJ)

            B = aca(
                ComplexF64,
                (r) -> vec(eval_block(ev, [rows[r]], cols)),
                (c) -> vec(eval_block(ev, rows, [cols[c]])),
                length(rows),
                length(cols);
                tol = params.tol,
                maxrank = params.maxrank,
                recompress = params.recompress,
            )
            size(B.U, 2) > 0 && push!(blocks, ACABlock(rows, cols, B.U, B.V))

            # 转置方向 (j, i)
            rowsT = vcat(colsJ, N .+ colsJ)
            colsT = vcat(rowsJ, N .+ rowsJ)
            B2 = aca(
                ComplexF64,
                (r) -> vec(eval_block(ev, [rowsT[r]], colsT)),
                (c) -> vec(eval_block(ev, rowsT, [colsT[c]])),
                length(rowsT),
                length(colsT);
                tol = params.tol,
                maxrank = params.maxrank,
                recompress = params.recompress,
            )
            size(B2.U, 2) > 0 && push!(blocks, ACABlock(rowsT, colsT, B2.U, B2.V))
        end
    end

    FT = eltype(basis.mesh.node)
    CT = eltype(Z_near)
    return ACAOperator{FT,CT}(
        octree,
        abstract_bases,
        basis_offsets,
        Z_near,
        blocks,
        pmchw,
        sorted_ids,
        inv_sorted_ids,
        params,
    )
end

function get_leaf_intervals(A::ACAOperator)
    leaf_level = A.octree.levels[A.octree.nLevels]
    return [cube.bfInterval for cube in leaf_level.cubes]
end

import ....Solvers: ILUPreconditioner, SPAIPreconditioner, BlockJacobiPreconditioner
ILUPreconditioner(op::ACAOperator; τ::Real = 0.01) = ILUPreconditioner(op.Z_near; τ = τ)
SPAIPreconditioner(op::ACAOperator) = SPAIPreconditioner(op.Z_near)
BlockJacobiPreconditioner(op::ACAOperator) = BlockJacobiPreconditioner(op.Z_near, _leaf_block_indices(op))
BlockJacobiPreconditioner(op::ACAOperator, ::Any) = BlockJacobiPreconditioner(op)

function _leaf_block_indices(op::ACAOperator)
    leaf_level = op.octree.levels[op.octree.nLevels]
    blocks = Vector{Vector{Int}}()
    for cube in leaf_level.cubes
        isempty(cube.bfInterval) && continue
        push!(blocks, collect(op.sorted_ids[cube.bfInterval]))
    end
    return blocks
end

function Base.:*(A::ACAOperator, x::AbstractVector)
    y = similar(x)
    mul!(y, A, x)
    return y
end

"""
    mul!(y, A::ACAOperator, x)

计算 `y = A*x = Z_near*x + Σ 远场低秩块`。
对称矩阵的 (j,i) 块用转置语义应用：`Z[cols, rows] ≈ V * (Uᵀ * x[rows])`
（转置约定无共轭），避免重复压缩。
"""
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
