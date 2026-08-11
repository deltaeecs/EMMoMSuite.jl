module HMatrixModule

using LinearAlgebra
using ...CoreModule: num_basis
using ...IntegralEquations.PMCHWModule: PMCHW
using ..BlockLUModule: extract_block

export HMatrixNode, hmatrix_from_mlaca, materialize, hmatrix_size, hmatrix_count

"""
    HMatrixNode{T}

H 矩阵树节点：
- `:dense`：稠密块（`dense`）；
- `:lowrank`：低秩块（`Z ≈ U * transpose(V)`，转置约定）；
- `:split`：按行/列子块递归划分（`children` 行主序，`nrowblocks × ncolblocks`）。

`rows`/`cols` 为全局索引；H-LU 阶段 `factor` 存放分解因子。
"""
mutable struct HMatrixNode{T}
    rows::Vector{Int}
    cols::Vector{Int}
    kind::Symbol
    dense::Matrix{T}
    U::Matrix{T}
    V::Matrix{T}
    children::Vector{HMatrixNode{T}}
    nrowblocks::Int
    ncolblocks::Int
    factor::Any
end

HMatrixNode(rows, cols, kind, dense, U, V, children, nrb, ncb) =
    HMatrixNode(rows, cols, kind, dense, U, V, children, nrb, ncb, nothing)

hmatrix_size(H::HMatrixNode) = (length(H.rows), length(H.cols))
hmatrix_count(H::HMatrixNode) = H.kind == :split ? sum(hmatrix_count, H.children) : 1

function _subtree_ids(levels, nLevels::Int, levelID::Int, cube_idx::Int, sorted_ids::Vector{Int})
    level = levels[levelID]
    cube = level.cubes[cube_idx]
    ids = Int[]
    if levelID == nLevels
        append!(ids, sorted_ids[cube.bfInterval])
    else
        for k in cube.kidsInterval
            append!(ids, _subtree_ids(levels, nLevels, levelID + 1, k, sorted_ids))
        end
    end
    return ids
end

function _build_hnode(op, levels, nLevels::Int, levelID::Int, iA::Int, iB::Int,
                      sorted_ids::Vector{Int}, pmchw::Bool, Npass::Int,
                      lookup)
    level = levels[levelID]
    cubeA = level.cubes[iA]
    cubeB = level.cubes[iB]
    idsA = _subtree_ids(levels, nLevels, levelID, iA, sorted_ids)
    idsB = _subtree_ids(levels, nLevels, levelID, iB, sorted_ids)
    rows = pmchw ? vcat(idsA, Npass .+ idsA) : idsA
    cols = pmchw ? vcat(idsB, Npass .+ idsB) : idsB
    CT = ComplexF64

    if iA == iB
        if levelID == nLevels
            d = extract_block(op, rows, cols)
            return HMatrixNode(rows, cols, :dense, d, Matrix{CT}(undef, 0, 0), Matrix{CT}(undef, 0, 0), HMatrixNode{CT}[], 0, 0)
        end
        kids = cubeA.kidsInterval
        nk = length(kids)
        children = HMatrixNode{CT}[
            _build_hnode(op, levels, nLevels, levelID + 1, kids[a], kids[b], sorted_ids, pmchw, Npass, lookup)
            for a in 1:nk for b in 1:nk
        ]
        return HMatrixNode(rows, cols, :split, Matrix{CT}(undef, 0, 0), Matrix{CT}(undef, 0, 0), Matrix{CT}(undef, 0, 0), children, nk, nk)
    end

    if !(iB in cubeA.neighbors)
        key = (Tuple(rows), Tuple(cols))
        if haskey(lookup, key)
            U, V = lookup[key]
            return HMatrixNode(rows, cols, :lowrank, Matrix{CT}(undef, 0, 0), U, V, HMatrixNode{CT}[], 0, 0)
        end
    end

    if levelID == nLevels
        d = extract_block(op, rows, cols)
        return HMatrixNode(rows, cols, :dense, d, Matrix{CT}(undef, 0, 0), Matrix{CT}(undef, 0, 0), HMatrixNode{CT}[], 0, 0)
    end

    kidsA = cubeA.kidsInterval
    kidsB = cubeB.kidsInterval
    children = HMatrixNode{CT}[
        _build_hnode(op, levels, nLevels, levelID + 1, a, b, sorted_ids, pmchw, Npass, lookup)
        for a in kidsA for b in kidsB
    ]
    return HMatrixNode(rows, cols, :split, Matrix{CT}(undef, 0, 0), Matrix{CT}(undef, 0, 0), Matrix{CT}(undef, 0, 0), children, length(kidsA), length(kidsB))
end

"""
    hmatrix_from_mlaca(op) -> HMatrixNode

将 MLACA/ACA 算子的分层低秩结构重建为显式 H 矩阵树。可容许对 → `:lowrank`
（对称算子反方向用转置因子 `V*(Uᵀ)`），叶层近对 → `:dense`，近对非叶 → `:split`。
PMCHW（2N 系统）自动展开 J/M 双通道。
"""
function hmatrix_from_mlaca(op)
    levels = op.octree.levels
    nLevels = op.octree.nLevels
    sorted_ids = op.sorted_ids
    pmchw = op.operator isa PMCHW
    Npass = pmchw ? num_basis(op.bases[1]) : 0

    # 低秩块查找：键 = (rows, cols) 元组 → (U, V)
    lookup = Dict{Tuple{Tuple{Vararg{Int}},Tuple{Vararg{Int}}},Tuple{Matrix{ComplexF64},Matrix{ComplexF64}}}()
    for blk in op.blocks
        lookup[(Tuple(blk.rows), Tuple(blk.cols))] = (blk.U, blk.V)
        if op.params.symmetric
            lookup[(Tuple(blk.cols), Tuple(blk.rows))] = (blk.V, blk.U)
        end
    end

    top = levels[1]
    n_top = length(top.cubes)
    CT = ComplexF64
    if n_top == 1 && nLevels == 1
        ids = sorted_ids[top.cubes[1].bfInterval]
        rows = pmchw ? vcat(ids, Npass .+ ids) : ids
        d = extract_block(op, rows, rows)
        return HMatrixNode(rows, rows, :dense, d, Matrix{CT}(undef, 0, 0), Matrix{CT}(undef, 0, 0), HMatrixNode{CT}[], 0, 0)
    end

    children = HMatrixNode{CT}[
        _build_hnode(op, levels, nLevels, 1, i, j, sorted_ids, pmchw, Npass, lookup)
        for i in 1:n_top for j in 1:n_top
    ]
    root_rows = vcat([_subtree_ids(levels, nLevels, 1, i, sorted_ids) for i in 1:n_top]...)
    if pmchw
        root_rows = vcat(root_rows, Npass .+ root_rows)
    end
    return HMatrixNode(root_rows, root_rows, :split, Matrix{CT}(undef, 0, 0), Matrix{CT}(undef, 0, 0), Matrix{CT}(undef, 0, 0), children, n_top, n_top)
end

"""
    materialize(H::HMatrixNode) -> Matrix

稠密化 H 矩阵树（测试与门控用）。
"""
function materialize(H::HMatrixNode{T}) where {T}
    m, n = length(H.rows), length(H.cols)
    if H.kind == :dense
        return H.dense
    elseif H.kind == :lowrank
        return H.U * transpose(H.V)
    end
    Z = zeros(T, m, n)
    rowpos = Dict(g => i for (i, g) in enumerate(H.rows))
    colpos = Dict(g => j for (j, g) in enumerate(H.cols))
    for ch in H.children
        r = [rowpos[g] for g in ch.rows]
        c = [colpos[g] for g in ch.cols]
        Z[r, c] = materialize(ch)
    end
    return Z
end

end # module HMatrixModule
