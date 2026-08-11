module MLACAOperatorModule

using LinearAlgebra
using SparseArrays
using ...CoreModule
using ...CoreModule: Constants
using ...Geometry
using ...BasisFunctions
using ...IntegralEquations
using ..MLFMA: OctreeInfo, build_octree
using ..MLFMA.MLFMAOperatorModule: assemble_near_field
using ..ACA: LowRankBlock, aca
using ..BlockEvaluatorModule: BlockEvaluator, eval_block
using ..ACAOperatorModule: ACAParams

export MLACAOperator

"""
    MLACABlock{CT}

多层低秩块：`Z[rows, cols] ≈ U * transpose(V)`，`rows`/`cols` 为全局基函数索引。
同一块的行/列可跨越多个叶层盒子（多层压缩证据）。
"""
struct MLACABlock{CT}
    rows::Vector{Int}
    cols::Vector{Int}
    U::Matrix{CT}
    V::Matrix{CT}
end

"""
    MLACAOperator{FT,CT} <: AbstractIntegralOperator

多层自适应交叉近似（MLACA，H-矩阵风格）：在八叉树多层结构上递归块压缩。
- 可容许对（非邻盒子）→ ACA 压缩为低秩块（可跨多个叶层盒子）。
- 近邻对 → 下钻到子层；叶层近邻对由 `Z_near`（复用 MLFMA 近场装配）覆盖。
- 对角块 → 递归到子层，叶层自/邻对在 `Z_near` 中。

当前仅支持复数对称算子（如 EFIE，`symmetric=true`）；非对称问题请使用
`ACAOperator(symmetric=false)`。

# 参考
- Gibson, The Method of Moments in Electromagnetics, 3rd, Ch10（MLACA 递归压缩）。
"""
struct MLACAOperator{FT,CT} <: AbstractIntegralOperator
    octree::OctreeInfo
    bases::Vector{AbstractBasisFunction}
    basis_offsets::Vector{Int}
    Z_near::SparseMatrixCSC{CT,Int}
    blocks::Vector{MLACABlock{CT}}
    operator::AbstractIntegralOperator
    sorted_ids::Vector{Int}
    inv_sorted_ids::Vector{Int}
    params::ACAParams
end

Base.eltype(::MLACAOperator{FT,CT}) where {FT,CT} = CT
Base.size(A::MLACAOperator) = size(A.Z_near)
Base.size(A::MLACAOperator, i::Int) = size(A.Z_near, i)

function MLACAOperator(
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
    symmetric || error("MLACAOperator 当前仅支持对称算子（symmetric=true）")
    params = ACAParams(; tol = tol, maxrank = maxrank, recompress = recompress, symmetric = true)

    # 1. 八叉树聚类（复用 MLFMA build_octree）
    bf_centers = reduce(hcat, [bf.center for bf in basis.functions])
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
    basis_offsets = cumsum([num_basis(basis)])
    abstract_bases = Vector{AbstractBasisFunction}([basis])

    # 2. 近场稀疏装配（复用 MLFMA assemble_near_field）
    Z_near = assemble_near_field(
        operator,
        abstract_bases,
        basis_offsets,
        octree,
        sorted_ids,
        inv_sorted_ids,
    )

    # 3. 多层递归块压缩
    ev = BlockEvaluator(operator, basis)
    blocks = Vector{MLACABlock{ComplexF64}}()
    levels = octree.levels
    nLevels = octree.nLevels
    top = levels[1]
    for i in 1:length(top.cubes), j in i:length(top.cubes)
        _compress_pair!(
            blocks,
            ev,
            levels,
            nLevels,
            1,
            i,
            j,
            sorted_ids,
            params,
        )
    end

    FT = eltype(basis.mesh.node)
    CT = eltype(Z_near)
    return MLACAOperator{FT,CT}(
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

function _collect_subtree!(
    ids::Vector{Int},
    levels,
    nLevels::Int,
    levelID::Int,
    cube_idx::Int,
    sorted_ids::Vector{Int},
)
    level = levels[levelID]
    cube = level.cubes[cube_idx]
    if levelID == nLevels
        append!(ids, sorted_ids[cube.bfInterval])
    else
        for k in cube.kidsInterval
            _collect_subtree!(ids, levels, nLevels, levelID + 1, k, sorted_ids)
        end
    end
    return ids
end

subtree_ids(levels, nLevels::Int, levelID::Int, cube_idx::Int, sorted_ids::Vector{Int}) =
    _collect_subtree!(Int[], levels, nLevels, levelID, cube_idx, sorted_ids)

"""
    _compress_pair!(blocks, ev, levels, nLevels, levelID, iA, iB, sorted_ids, params)

递归压缩盒子对 (iA, iB)：
- `iA == iB`：对角块下钻到子层（叶层由 `Z_near` 覆盖）。
- 非邻（可容许）→ ACA 压缩整个子树块。
- 近邻且非叶 → 下钻到全部子盒子对；近邻且叶层 → `Z_near` 覆盖。
对称矩阵只处理 `iA ≤ iB`，下三角由 `mul!` 转置语义应用。
"""
function _compress_pair!(
    blocks::Vector{MLACABlock{CT}},
    ev::BlockEvaluator,
    levels,
    nLevels::Int,
    levelID::Int,
    iA::Int,
    iB::Int,
    sorted_ids::Vector{Int},
    params::ACAParams,
) where {CT}
    level = levels[levelID]
    cubeA = level.cubes[iA]
    cubeB = level.cubes[iB]

    if iA == iB
        levelID == nLevels && return
        kids = cubeA.kidsInterval
        nk = length(kids)
        for a in 1:nk, b in a:nk
            _compress_pair!(blocks, ev, levels, nLevels, levelID + 1, kids[a], kids[b], sorted_ids, params)
        end
        return
    end

    if !(iB in cubeA.neighbors)
        # 可容许 → 压缩整个子树块
        rows = subtree_ids(levels, nLevels, levelID, iA, sorted_ids)
        cols = subtree_ids(levels, nLevels, levelID, iB, sorted_ids)
        (isempty(rows) || isempty(cols)) && return
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
        size(B.U, 2) > 0 && push!(blocks, MLACABlock(rows, cols, B.U, B.V))
        return
    end

    # 近邻：下钻（叶层邻对已在 Z_near）
    levelID == nLevels && return
    for a in cubeA.kidsInterval, b in cubeB.kidsInterval
        _compress_pair!(blocks, ev, levels, nLevels, levelID + 1, a, b, sorted_ids, params)
    end
    return
end

function Base.:*(A::MLACAOperator, x::AbstractVector)
    y = similar(x)
    mul!(y, A, x)
    return y
end

import ....Solvers: ILUPreconditioner, SPAIPreconditioner, BlockJacobiPreconditioner
ILUPreconditioner(op::MLACAOperator; τ::Real = 0.01) = ILUPreconditioner(op.Z_near; τ = τ)
SPAIPreconditioner(op::MLACAOperator) = SPAIPreconditioner(op.Z_near)
BlockJacobiPreconditioner(op::MLACAOperator) = BlockJacobiPreconditioner(op.Z_near, _leaf_block_indices(op))
BlockJacobiPreconditioner(op::MLACAOperator, ::Any) = BlockJacobiPreconditioner(op)

function _leaf_block_indices(op::MLACAOperator)
    leaf_level = op.octree.levels[op.octree.nLevels]
    blocks = Vector{Vector{Int}}()
    for cube in leaf_level.cubes
        isempty(cube.bfInterval) && continue
        push!(blocks, collect(op.sorted_ids[cube.bfInterval]))
    end
    return blocks
end

"""
    mul!(y, A::MLACAOperator, x)

`y = A*x = Z_near*x + Σ 多层低秩块`。对称矩阵的下三角块用转置语义
（`V * (Uᵀ * x)`，无共轭）应用。
"""
function LinearAlgebra.mul!(y::AbstractVector, A::MLACAOperator, x::AbstractVector)
    mul!(y, A.Z_near, x)
    for blk in A.blocks
        xc = view(x, blk.cols)
        t = blk.U * (transpose(blk.V) * xc)
        @inbounds for (idx, r) in enumerate(blk.rows)
            y[r] += t[idx]
        end
        xi = view(x, blk.rows)
        t2 = blk.V * (transpose(blk.U) * xi)
        @inbounds for (idx, c) in enumerate(blk.cols)
            y[c] += t2[idx]
        end
    end
    return y
end

end # module MLACAOperatorModule
