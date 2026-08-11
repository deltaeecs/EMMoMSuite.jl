module BlockLUModule

using LinearAlgebra
using SparseArrays
using ..ACA: LowRankBlock, aca

export BlockLUFactorization, block_lu, block_lu_solve

"""
    BlockLUFactorization{CT}

ACA/MLACA 压缩矩阵的直接块 LU 分解（Gibson《Method of Moments》Ch9 Algorithm 7）。
按八叉树叶层盒子分块：
- 对角块 `diag[b]`：`lu(A_bb)`，其中 `A_bb = Z_bb - Σ_{p<b} L_bp U_pb`；
- 下三角 `L[s][b]`（s > b）：`L_sb = (Z_sb - Σ_{p<b} L_sp U_pb) * U_bb⁻¹`；
- 上三角 `U[b][s]`（s > b）：`U_bs = L_bb⁻¹ * (Z_bs - Σ_{p<b} L_bp U_ps)`。

分解一次后支持多 RHS 直接求解（前代 + 回代）。分解得到的离对角块可
用 ACA 再压缩（`recompress=true`）控制存储。
"""
struct BlockLUFactorization{CT}
    N::Int
    blocks::Vector{Vector{Int}}       # 每个对角块的全局索引
    diag::Vector{Any}                 # 对角块 lu(A_bb, NoPivot()) 分解
    L::Vector{Vector{Any}}            # L[s][b] for s > b
    U::Vector{Vector{Any}}            # U[b][s] for s > b
end

Base.size(F::BlockLUFactorization) = (F.N, F.N)

"""
    block_lu(op; tol=1e-4, recompress=true) -> BlockLUFactorization

对 ACA/MLACA 算子做叶层块 LU 分解。`Z_bb` 与离对角块由近场稀疏矩阵与
低秩块重建（`extract_block`）。
"""
function block_lu(op; tol::Real = 1e-4, recompress::Bool = true)
    leaf_level = op.octree.levels[op.octree.nLevels]
    cubes = leaf_level.cubes
    idxs = [collect(op.sorted_ids[c.bfInterval]) for c in cubes]
    filter!(!isempty, idxs)
    M = length(idxs)
    CT = eltype(op)

    diag = Vector{Any}(undef, M)
    L = [Vector{Any}(undef, M) for _ in 1:M]
    U = [Vector{Any}(undef, M) for _ in 1:M]

    for b in 1:M
        # A_bb = Z_bb - Σ_{p<b} L_bp U_pb
        A = extract_block(op, idxs[b], idxs[b])
        for p in 1:(b-1)
            A .-= mul_dense(L[b][p], U[p][b])
        end
        # 块 LU 要求 A_bb = L_bb * U_bb 精确成立（无行交换），否则块更新公式失效
        diag[b] = lu(A, NoPivot())

        for s in (b+1):M
            # L_sb = (Z_sb - Σ_{p<b} L_sp U_pb) * U_bb⁻¹
            R = extract_block(op, idxs[s], idxs[b])
            for p in 1:(b-1)
                R .-= mul_dense(L[s][p], U[p][b])
            end
            Lblk = R / UpperTriangular(Matrix(diag[b].U))
            L[s][b] = recompress ? _maybe_aca(Lblk; tol = tol) : Lblk

            # U_bs = L_bb⁻¹ * (Z_bs - Σ_{p<b} L_bp U_ps)
            C = extract_block(op, idxs[b], idxs[s])
            for p in 1:(b-1)
                C .-= mul_dense(L[b][p], U[p][s])
            end
            Ublk = LowerTriangular(Matrix(diag[b].L)) \ C
            U[b][s] = recompress ? _maybe_aca(Ublk; tol = tol) : Ublk
        end
    end

    return BlockLUFactorization{CT}(size(op, 1), idxs, diag, L, U)
end

"""
    extract_block(op, Ia::Vector{Int}, Ib::Vector{Int}) -> Matrix{ComplexF64}

从近场稀疏矩阵与低秩块重建稠密块 `Z[Ia, Ib]`。对称算子的反向块由转置
贡献补齐（`V * (Uᵀ)`）；非对称算子两个方向均已显式存储。
"""
function extract_block(op, Ia::Vector{Int}, Ib::Vector{Int})
    CT = eltype(op)
    Z = zeros(CT, length(Ia), length(Ib))
    Z .+= op.Z_near[Ia, Ib]
    rowmap = Dict(g => i for (i, g) in enumerate(Ia))
    colmap = Dict(g => j for (j, g) in enumerate(Ib))
    for blk in op.blocks
        # 正向贡献：blk.rows → Ia，blk.cols → Ib
        rp = Int[]; tp = Int[]
        for (i, r) in enumerate(blk.rows)
            haskey(rowmap, r) && (push!(rp, i); push!(tp, rowmap[r]))
        end
        cp = Int[]; sp = Int[]
        for (j, c) in enumerate(blk.cols)
            haskey(colmap, c) && (push!(cp, j); push!(sp, colmap[c]))
        end
        if !isempty(rp) && !isempty(cp)
            Z[tp, sp] .+= blk.U[rp, :] * transpose(blk.V[cp, :])
        end
        # 反向贡献（仅对称算子：下三角由转置补齐；非对称算子已显式存储两个方向）
        if op.params.symmetric
            rp2 = Int[]; sp2 = Int[]
            for (i, r) in enumerate(blk.rows)
                haskey(colmap, r) && (push!(rp2, i); push!(sp2, colmap[r]))
            end
            cp2 = Int[]; tp2 = Int[]
            for (j, c) in enumerate(blk.cols)
                haskey(rowmap, c) && (push!(cp2, j); push!(tp2, rowmap[c]))
            end
            if !isempty(rp2) && !isempty(cp2)
                Z[tp2, sp2] .+= blk.V[cp2, :] * transpose(blk.U[rp2, :])
            end
        end
    end
    return Z
end

# 块乘法（稠密化；矩阵尺寸小，正确性优先）
mul_dense(X::AbstractMatrix, Y::AbstractMatrix) = X * Y
mul_dense(X::LowRankBlock, Y::AbstractMatrix) = X.U * (transpose(X.V) * Y)
mul_dense(X::AbstractMatrix, Y::LowRankBlock) = (X * Y.U) * transpose(Y.V)
mul_dense(X::LowRankBlock, Y::LowRankBlock) =
    (X.U * (transpose(X.V) * Y.U)) * transpose(Y.V)

function _maybe_aca(A::Matrix{CT}; tol::Real = 1e-4) where {CT}
    m, n = size(A)
    k_full = min(m, n)
    B = aca(A; tol = tol, maxrank = k_full, recompress = true)
    if size(B.U, 2) * (m + n) < m * n
        return B
    end
    return A
end

"""
    block_lu_solve(F::BlockLUFactorization, B::AbstractMatrix) -> X

多 RHS 直接求解 `F * X = B`（前代 + 回代）。
"""
function block_lu_solve(F::BlockLUFactorization, B::AbstractMatrix)
    M = length(F.blocks)
    X = zeros(eltype(B), size(B))
    Y = zeros(eltype(B), size(B))

    # 前代：y_b = L_bb⁻¹ (b_b - Σ_{s<b} L_bs y_s)
    for b in 1:M
        rhs = B[F.blocks[b], :]
        for s in 1:(b-1)
            rhs .-= apply_block(F.L[b][s], Y[F.blocks[s], :])
        end
        Y[F.blocks[b], :] = LowerTriangular(Matrix(F.diag[b].L)) \ rhs
    end

    # 回代：x_b = U_bb⁻¹ (y_b - Σ_{s>b} U_bs x_s)
    for b in M:-1:1
        rhs = Y[F.blocks[b], :]
        for s in (b+1):M
            rhs .-= apply_block(F.U[b][s], X[F.blocks[s], :])
        end
        X[F.blocks[b], :] = UpperTriangular(Matrix(F.diag[b].U)) \ rhs
    end
    return X
end

block_lu_solve(F::BlockLUFactorization, b::AbstractVector) = vec(block_lu_solve(F, reshape(b, :, 1)))

apply_block(X::AbstractMatrix, v::AbstractVecOrMat) = X * v
apply_block(X::LowRankBlock, v::AbstractVecOrMat) = X.U * (transpose(X.V) * v)

Base.:\(F::BlockLUFactorization, b::AbstractVector) = block_lu_solve(F, b)
Base.:\(F::BlockLUFactorization, B::AbstractMatrix) = block_lu_solve(F, B)

end # module BlockLUModule
