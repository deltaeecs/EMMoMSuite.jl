# Preconditioners.jl — MPI 分布式预条件（Phase 15 目标：预条件支持 MPI 并行）
#
# 设计：
#   - `DistributedBlockJacobiPreconditioner`：块按叶 cube 归属分到各秩（rank 拥有
#     `(i_cube-1)%P` 号 cube 的 J/M 行块），构造只提取本秩块的行；施加时各秩
#     对本秩块做 LU 求解 → 1 次 Allreduce 汇聚完整 y。
#   - `DistributedDiagonalPreconditioner`：对角逆向量每秩复制（O(N)），施加无通信。
#   - `apply_mpi_preconditioner!(y, P, x)`：统一入口；`nothing` 与串行预条件
#     （BlockJacobi/Diagonal/...）作为复制式回退（每秩完整施加，结果一致）。
#
# 内存：预条件块不再每秩全量复制（旧 BlockJacobiPreconditioner 每秩保存全部块的
# LU 分解），分布式版本每秩只存 P 分之一。

using LinearAlgebra
using SparseArrays
import MPI

import ..Solvers:
    BlockJacobiPreconditioner,
    DiagonalPreconditioner,
    IdentityPreconditioner,
    SPAIPreconditioner,
    ILUPreconditioner

import ..FastAlgorithms.MLFMA.MLFMAOperatorModule: MLFMAOperatorMPI
import ..FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperatorMPI

"""
    DistributedBlockJacobiPreconditioner{CT,FT}

MPI 分布式块 Jacobi 预条件：块（叶 cube 的基函数行集）按 cube 归属分到各秩，
每秩只保存并求解自己拥有的块。`apply_mpi_preconditioner!` 施加后 Allreduce。
"""
struct DistributedBlockJacobiPreconditioner{CT}
    blocks::Vector{LU{CT,Matrix{CT},Vector{Int}}}  # 每块 LU 分解（本秩拥有，类型稳定）
    block_rows::Vector{Vector{Int}}    # 每块的全局行号（本秩拥有）
    comm
end

"""
    DistributedDiagonalPreconditioner{T}

MPI 分布式对角（Jacobi）预条件：对角逆全量复制（O(N)/秩），施加为逐元素乘法，无通信。
"""
struct DistributedDiagonalPreconditioner{T}
    diag_inv::Vector{T}
end

function DistributedBlockJacobiPreconditioner(
    Z_near_local::SparseMatrixCSC{CT,Int},
    block_rows::Vector{Vector{Int}},
    comm,
) where {CT}
    blocks = Vector{LU{CT,Matrix{CT},Vector{Int}}}(undef, length(block_rows))
    for (ib, idx) in enumerate(block_rows)
        B = Matrix{CT}(Z_near_local[idx, idx])
        blocks[ib] = lu(B)
    end
    return DistributedBlockJacobiPreconditioner{CT}(blocks, block_rows, comm)
end

"""
    DistributedBlockJacobiPreconditioner(op::MLFMAOperatorMPI)

从 MPI MLFMA 算子构造分布式块 Jacobi：块 = 本秩拥有的叶 cube 基函数行
（与 `Z_near_local` 的 cube 分区一致，行提取完整）。
"""
function DistributedBlockJacobiPreconditioner(op::MLFMAOperatorMPI)
    comm = op.comm
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    leaf = op.octree.levels[op.octree.nLevels]
    block_rows = Vector{Vector{Int}}()
    for (ic, cube) in enumerate(leaf.cubes)
        isempty(cube.bfInterval) && continue
        (ic - 1) % P == rank || continue
        push!(block_rows, collect(op.sorted_ids[cube.bfInterval]))
    end
    return DistributedBlockJacobiPreconditioner(op.Z_near_local, block_rows, comm)
end

"""
    DistributedBlockJacobiPreconditioner(op::PMCHWMLFMAOperatorMPI)

PMCHW 2N×2N 系统：块 = 本秩拥有的叶 cube 的 J 行 ∪ M 行（i 与 i+N）。
"""
function DistributedBlockJacobiPreconditioner(op::PMCHWMLFMAOperatorMPI)
    comm = op.comm
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    N = size(op, 1) ÷ 2
    leaf = op.octree0.levels[op.octree0.nLevels]
    block_rows = Vector{Vector{Int}}()
    for (ic, cube) in enumerate(leaf.cubes)
        isempty(cube.bfInterval) && continue
        (ic - 1) % P == rank || continue
        ids = collect(op.sorted_ids0[cube.bfInterval])
        n_leaf = length(ids)
        rows = Vector{Int}(undef, 2 * n_leaf)
        copyto!(rows, 1, ids, 1, n_leaf)
        for i = 1:n_leaf
            rows[n_leaf + i] = ids[i] + N
        end
        push!(block_rows, rows)
    end
    return DistributedBlockJacobiPreconditioner(op.Z_near_local, block_rows, comm)
end

"""
    apply_mpi_preconditioner!(y, P, x) → y

MPI 预条件统一施加入口（左预条件 M⁻¹）：
- `nothing`：y = x
- `DistributedBlockJacobiPreconditioner`：本秩块 LU 求解 + Allreduce
- `DistributedDiagonalPreconditioner`：逐元素除法（无通信）
- 串行预条件：每秩完整施加（复制式回退，结果一致但内存每秩全量）
"""
function apply_mpi_preconditioner!(y::AbstractVector{CT}, P::Nothing, x::AbstractVector{CT}) where {CT}
    copyto!(y, x)
    return y
end

function apply_mpi_preconditioner!(
    y::AbstractVector{CT},
    P::DistributedBlockJacobiPreconditioner,
    x::AbstractVector{CT},
) where {CT}
    fill!(y, zero(CT))
    @inbounds for ib in eachindex(P.blocks)
        idx = P.block_rows[ib]
        xb = view(x, idx)
        yb = view(y, idx)
        ldiv!(yb, P.blocks[ib], xb)
    end
    MPI.Allreduce!(y, +, P.comm)
    return y
end

function apply_mpi_preconditioner!(
    y::AbstractVector{CT},
    P::DistributedDiagonalPreconditioner,
    x::AbstractVector{CT},
) where {CT}
    @inbounds for i in eachindex(y)
        y[i] = P.diag_inv[i] * x[i]
    end
    return y
end

function apply_mpi_preconditioner!(
    y::AbstractVector{CT},
    P::Union{BlockJacobiPreconditioner,DiagonalPreconditioner,IdentityPreconditioner,SPAIPreconditioner,ILUPreconditioner},
    x::AbstractVector{CT},
) where {CT}
    ldiv!(y, P, x)
    return y
end

"""
    DistributedDiagonalPreconditioner(A, comm)

从全量（或本地行）对角提取构造。`A` 需支持 `A[i,i]`。
"""
function DistributedDiagonalPreconditioner(A::AbstractMatrix{T}, comm) where {T}
    n = size(A, 1)
    diag_inv = Vector{T}(undef, n)
    for i = 1:n
        d = A[i, i]
        diag_inv[i] = abs(d) < 1e-15 ? one(T) : inv(d)
    end
    return DistributedDiagonalPreconditioner(diag_inv)
end

function DistributedDiagonalPreconditioner(op::MLFMAOperatorMPI)
    N = size(op, 1)
    diag_inv = zeros(ComplexF64, N)
    # 从本地近场行提取对角（各秩只填自己的行，Allreduce 汇聚）
    Z = op.Z_near_local
    I, J, V = findnz(Z)
    for k in eachindex(I)
        if I[k] == J[k]
            diag_inv[I[k]] = V[k]
        end
    end
    MPI.Allreduce!(diag_inv, +, op.comm)
    for i = 1:N
        diag_inv[i] = abs(diag_inv[i]) < 1e-15 ? 1.0 : inv(diag_inv[i])
    end
    return DistributedDiagonalPreconditioner(diag_inv)
end

function DistributedDiagonalPreconditioner(op::PMCHWMLFMAOperatorMPI)
    N = size(op, 1)
    diag_inv = zeros(ComplexF64, N)
    I, J, V = findnz(op.Z_near_local)
    for k in eachindex(I)
        if I[k] == J[k]
            diag_inv[I[k]] = V[k]
        end
    end
    MPI.Allreduce!(diag_inv, +, op.comm)
    for i = 1:N
        diag_inv[i] = abs(diag_inv[i]) < 1e-15 ? 1.0 : inv(diag_inv[i])
    end
    return DistributedDiagonalPreconditioner(diag_inv)
end
