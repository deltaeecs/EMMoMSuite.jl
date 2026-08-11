module HMatrixModule

using LinearAlgebra
using ...CoreModule: num_basis
using ...IntegralEquations.PMCHWModule: PMCHW
using ..ACA: LowRankBlock, aca
using ..BlockLUModule: extract_block

export HMatrixNode, hmatrix_from_mlaca, materialize, hmatrix_size, hmatrix_count,
    h_lu!, h_lu_solve

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

# =============================================================================
# H-LU 分解与多 RHS 求解（标准分层块 LU；默认不截断保证精确分解）
# =============================================================================

function _child(A::HMatrixNode, i::Int, j::Int)
    return A.children[(i - 1) * A.ncolblocks + j]
end

_mul_dense(X::HMatrixNode, Y::HMatrixNode) = materialize(X) * materialize(Y)

"""
    h_lu!(A::HMatrixNode; tol=1e-4, recompress=false)

原地 H-LU 分解（标准分层块 LU）：对角递归分解；离对角
`L_sb = (Z_sb − Σ L_sp U_pb) U_bb⁻¹`、`U_bs = L_bb⁻¹ (Z_bs − Σ L_bp U_ps)`。
默认 `recompress=false` 保证精确分解；`recompress=true` 时用 ACA 截断离对角
因子块（校验误差 ≤ 10·tol，否则回退稠密）。
"""
function h_lu!(A::HMatrixNode{T}; tol::Real = 1e-4, recompress::Bool = false) where {T}
    if A.kind == :dense
        A.factor = lu(A.dense, NoPivot())
        return A
    end
    A.kind == :lowrank && return A

    nb = A.nrowblocks
    for b in 1:nb
        Ab = _child(A, b, b)
        # 对角更新：A_bb = Z_bb - Σ_{p<b} L_bp U_pb（先更新再分解）
        Dd = materialize(Ab)
        for p in 1:(b - 1)
            Dd .-= _mul_dense(_child(A, b, p), _child(A, p, b))
        end
        _put_back!(Ab, Dd; tol = tol, recompress = recompress)
        h_lu!(Ab; tol = tol, recompress = recompress)

        for s in (b + 1):nb
            # L_sb = (Z_sb - Σ_{p<b} L_sp U_pb) * U_bb⁻¹
            R = _child(A, s, b)
            D = materialize(R)
            for p in 1:(b - 1)
                D .-= _mul_dense(_child(A, s, p), _child(A, p, b))
            end
            _h_rdiv_U!(D, Ab)
            _put_back!(R, D; tol = tol, recompress = recompress)

            # U_bs = L_bb⁻¹ * (Z_bs - Σ_{p<b} L_bp U_ps)
            C = _child(A, b, s)
            E = materialize(C)
            for p in 1:(b - 1)
                E .-= _mul_dense(_child(A, b, p), _child(A, p, s))
            end
            _h_ldiv_L!(E, Ab)
            _put_back!(C, E; tol = tol, recompress = recompress)
        end
    end
    return A
end

# 右除 U：X = X * U_bb⁻¹（U 为上三角块）
function _h_rdiv_U!(X::Matrix, D::HMatrixNode)
    if D.kind == :dense
        F = D.factor
        X .= X / UpperTriangular(Matrix(F.U))
        return X
    end
    nb = D.nrowblocks
    colpos = Dict(g => i for (i, g) in enumerate(D.rows))
    for b in 1:nb
        cb = [colpos[g] for g in _child(D, b, b).rows]
        rhs = X[:, cb]
        for t in 1:(b - 1)
            Utb = _child(D, t, b)
            ct = [colpos[g] for g in _child(D, t, t).rows]
            rhs .-= X[:, ct] * materialize(Utb)
        end
        X[:, cb] = _h_rdiv_U!(rhs, _child(D, b, b))
    end
    return X
end

# 左除 L：X = L_bb⁻¹ * X（L 为下三角块）
function _h_ldiv_L!(X::Matrix, D::HMatrixNode)
    if D.kind == :dense
        F = D.factor
        X .= LowerTriangular(Matrix(F.L)) \ X
        return X
    end
    nb = D.nrowblocks
    rowpos = Dict(g => i for (i, g) in enumerate(D.rows))
    for b in 1:nb
        rb = [rowpos[g] for g in _child(D, b, b).rows]
        rhs = X[rb, :]
        for p in 1:(b - 1)
            Lbp = _child(D, b, p)
            rp = [rowpos[g] for g in _child(D, p, p).rows]
            rhs .-= materialize(Lbp) * X[rp, :]
        end
        X[rb, :] = _h_ldiv_L!(rhs, _child(D, b, b))
    end
    return X
end

# 左除 U：X = U_bb⁻¹ * X（U 为上三角块）
function _h_ldiv_U!(X::Matrix, D::HMatrixNode)
    if D.kind == :dense
        F = D.factor
        X .= UpperTriangular(Matrix(F.U)) \ X
        return X
    end
    nb = D.nrowblocks
    rowpos = Dict(g => i for (i, g) in enumerate(D.rows))
    for b in nb:-1:1
        rb = [rowpos[g] for g in _child(D, b, b).rows]
        rhs = X[rb, :]
        for p in (b + 1):nb
            Ubp = _child(D, b, p)
            rp = [rowpos[g] for g in _child(D, p, p).rows]
            rhs .-= materialize(Ubp) * X[rp, :]
        end
        X[rb, :] = _h_ldiv_U!(rhs, _child(D, b, b))
    end
    return X
end

"""
    _put_back!(node, D; tol, recompress)

将稠密块按节点的树结构写回：`:dense` 存稠密；`:lowrank` 用 ACA 再压缩
（误差校验后保留低秩，否则回退稠密）；`:split` 递归写回子块。
"""
function _put_back!(node::HMatrixNode, D::Matrix; tol::Real = 1e-4, recompress::Bool = false)
    if node.kind == :dense
        node.dense = D
        node.factor = nothing
    elseif node.kind == :lowrank
        if recompress
            B = aca(D; tol = tol, maxrank = min(size(D)...), recompress = true)
            if size(B.U, 2) * (length(node.rows) + length(node.cols)) < length(D)
                err = norm(D - B.U * transpose(B.V)) / norm(D)
                if err <= 10 * tol
                    node.U = B.U
                    node.V = B.V
                    return node
                end
            end
        end
        node.kind = :dense
        node.dense = D
    else
        rowpos = Dict(g => i for (i, g) in enumerate(node.rows))
        colpos = Dict(g => j for (j, g) in enumerate(node.cols))
        for ch in node.children
            r = [rowpos[g] for g in ch.rows]
            c = [colpos[g] for g in ch.cols]
            _put_back!(ch, D[r, c]; tol = tol, recompress = recompress)
        end
    end
    return node
end

"""
    h_lu_solve(H::HMatrixNode, B::AbstractMatrix) -> X

用已分解的 H 树做多 RHS 直接求解（前代 + 回代）。输入/输出为全局索引顺序。
"""
function h_lu_solve(H::HMatrixNode, B::AbstractMatrix)
    n = length(H.rows)
    nrhs = size(B, 2)
    CT = eltype(B)
    if H.kind == :dense
        return Matrix(H.factor \ B[H.rows, :])
    end

    Bp = B[H.rows, :]   # 置换到树序
    X = zeros(CT, n, nrhs)
    Y = zeros(CT, n, nrhs)
    nb = H.nrowblocks
    rowpos = Dict(g => i for (i, g) in enumerate(H.rows))

    # 前代：Y_b = L_bb⁻¹ (b_b - Σ_{p<b} L_bp Y_p)
    for b in 1:nb
        rb = [rowpos[g] for g in _child(H, b, b).rows]
        rhs = Bp[rb, :]
        for p in 1:(b - 1)
            Lbp = _child(H, b, p)
            rp = [rowpos[g] for g in _child(H, p, p).rows]
            rhs .-= materialize(Lbp) * Y[rp, :]
        end
        Y[rb, :] = _h_ldiv_L!(rhs, _child(H, b, b))
    end
    # 回代：X_b = U_bb⁻¹ (Y_b - Σ_{p>b} U_bp X_p)
    for b in nb:-1:1
        rb = [rowpos[g] for g in _child(H, b, b).rows]
        rhs = Y[rb, :]
        for p in (b + 1):nb
            Ubp = _child(H, b, p)
            rp = [rowpos[g] for g in _child(H, p, p).rows]
            rhs .-= materialize(Ubp) * X[rp, :]
        end
        X[rb, :] = _h_ldiv_U!(rhs, _child(H, b, b))
    end

    Xg = zeros(CT, size(B, 1), nrhs)
    Xg[H.rows, :] = X
    return Xg
end

h_lu_solve(H::HMatrixNode, b::AbstractVector) = vec(h_lu_solve(H, reshape(b, :, 1)))

Base.:\(H::HMatrixNode, b::AbstractVector) = h_lu_solve(H, b)
Base.:\(H::HMatrixNode, B::AbstractMatrix) = h_lu_solve(H, B)

end # module HMatrixModule
