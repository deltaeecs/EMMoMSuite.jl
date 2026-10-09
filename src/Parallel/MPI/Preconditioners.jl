# Preconditioners.jl — MPI 分布式预条件（Phase 15 目标：预条件支持 MPI 并行）
#
# 设计：
#   - `DistributedBlockJacobiPreconditioner`：块按叶 cube 归属分到各秩（rank 拥有
#     `(i_cube-1)%P` 号 cube 的 J/M 行块），构造只提取本秩块的行；施加时各秩
#     对本秩块做 LU 求解 → 1 次 Allreduce 汇聚完整 y。
#   - `DistributedDiagonalPreconditioner`：对角逆向量每秩复制（O(N)），施加无通信。
#   - `DistributedSPAIPreconditioner`：逐八叉树块 SAI 行块按 cube 归属分到各秩，
#     构造只向各块所有者请求所需行（不做全量 Z_near 汇总）；施加 SpMV + Allreduce。
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
- `DistributedSPAIPreconditioner`：本秩 SAI 行块 SpMV + Allreduce
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
    DistributedSPAIPreconditioner{CT}

MPI 分布式逐八叉树块 SAI 预条件（左预条件）：块按叶 cube 归属分到各秩
（rank 拥有 `(i_cube-1)%P` 号 cube 的 M 行块，与 `DistributedBlockJacobiPreconditioner`
同一约定），每秩只构造并保存自己拥有 cube 的 SAI 行块。

构造不交换整个 `Z_near`：每秩只向各块所有者请求本秩 cube 扩展邻域
（neibfs ∪ neisNeibfs）涉及的**那几行** `Z_near_local`（Alltoall 请求 +
Alltoallv 回传行非零元），避免逐行/全量 Allgatherv 交换。
施加：`y_local = M_local * x`（x 每秩复制）+ 1 次 Allreduce 汇聚完整 y。
"""
struct DistributedSPAIPreconditioner{CT}
    M_local::SparseMatrixCSC{CT,Int}   # 本秩拥有的行块（仅 owned cube 的 M 行非零）
    comm
end

"""
    DistributedSPAIPreconditioner(op::MLFMAOperatorMPI)

从 MPI MLFMA 算子构造分布式逐块 SAI：按块所有者分布构造（仅交换所需行），
施加走 `apply_mpi_preconditioner!`（SpMV + Allreduce）。
"""
function DistributedSPAIPreconditioner(op::MLFMAOperatorMPI)
    comm = op.comm
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    Z_local = op.Z_near_local
    CT = eltype(Z_local)
    sorted_ids = op.sorted_ids
    N = size(Z_local, 1)
    leaf = op.octree.levels[op.octree.nLevels].cubes

    owned = Int[]   # 本秩拥有的非空 cube 编号
    for (ic, c) in enumerate(leaf)
        isempty(c.bfInterval) && continue
        (ic - 1) % P == rank && push!(owned, ic)
    end

    # 1. 排序位置 → 拥有者秩（一次遍历叶层建立查找表）
    pos_owner = fill(Int32(-1), N)
    for (ic, c) in enumerate(leaf)
        isempty(c.bfInterval) && continue
        o = Int32((ic - 1) % P)
        for s in c.bfInterval
            pos_owner[s] = o
        end
    end

    # 2. 本秩所需的全局行 = 所有 owned cube 的 neibfs ∪ neisNeibfs
    needed_pos = Int[]                       # 排序位置（用于查拥有者）
    cube_sets = Dict{Int,Tuple{Vector{Int},Vector{Int},Vector{Int}}}()  # ic => (cbfs, neibfs, neisNeibfs) 排序位置
    for ic in owned
        cube = leaf[ic]
        cbfs_s = collect(cube.bfInterval)
        neibfs_s = sort!(unique!(vcat(cbfs_s,
            [collect(leaf[i].bfInterval) for i in cube.neighbors]...)))
        # 列集 = 本块 ∪ N(N(c))（本块显式包含，保证 cbfsInCnnei 总能找到）
        neis_s = sort!(unique!(vcat(
            collect(cube.bfInterval),
            [collect(leaf[j].bfInterval) for i in cube.neighbors
             for j in leaf[i].neighbors]...)))
        cube_sets[ic] = (cbfs_s, neibfs_s, neis_s)
        append!(needed_pos, neibfs_s)
        append!(needed_pos, neis_s)
    end
    sort!(needed_pos)
    unique!(needed_pos)

    # 3. 按拥有者分组请求全局行号；Alltoall 请求 + Alltoallv 回传行非零元
    req_rows = [Int[] for _ = 1:P]           # 向 rank q 请求的全局行号
    for s in needed_pos
        push!(req_rows[pos_owner[s] + 1], sorted_ids[s])
    end
    foreach(req_rows) do r
        sort!(r)
        unique!(r)
    end

    scounts = Int32[length(req_rows[q]) for q = 1:P]
    rcounts = MPI.Alltoallv!(MPI.VBuffer(scounts, fill(Int32(1), P)), MPI.VBuffer(Vector{Int32}(undef, P), fill(Int32(1), P)), comm)
    send_rows = reduce(vcat, req_rows; init = Int[])
    recv_rows = MPI.Alltoallv!(MPI.VBuffer(send_rows, scounts),
                               MPI.VBuffer(Vector{Int}(undef, sum(rcounts)), rcounts), comm)

    # 打包应答：每行 [全局行号, nnz] + 列号(Int32) + 值(ComplexF64)
    hdr_out = [Int32[] for _ = 1:P]
    col_out = [Int32[] for _ = 1:P]
    val_out = [ComplexF64[] for _ = 1:P]
    off = 0
    for q = 1:P
        for k = 1:rcounts[q]
            g = recv_rows[off + k]
            col, val = findnz(Z_local[g, :])
            push!(hdr_out[q], Int32(g), Int32(length(col)))
            append!(col_out[q], col)
            append!(val_out[q], val)
        end
        off += rcounts[q]
    end

    hcounts = Int32[length(hdr_out[q]) for q = 1:P]
    ccounts = Int32[length(col_out[q]) for q = 1:P]
    vcounts = Int32[length(val_out[q]) for q = 1:P]
    rh = MPI.Alltoallv!(MPI.VBuffer(hcounts, fill(Int32(1), P)), MPI.VBuffer(Vector{Int32}(undef, P), fill(Int32(1), P)), comm)
    rc = MPI.Alltoallv!(MPI.VBuffer(ccounts, fill(Int32(1), P)), MPI.VBuffer(Vector{Int32}(undef, P), fill(Int32(1), P)), comm)
    rv = MPI.Alltoallv!(MPI.VBuffer(vcounts, fill(Int32(1), P)), MPI.VBuffer(Vector{Int32}(undef, P), fill(Int32(1), P)), comm)
    hdr_in = MPI.Alltoallv!(MPI.VBuffer(reduce(vcat, hdr_out; init = Int32[]), hcounts),
                            MPI.VBuffer(Vector{Int32}(undef, sum(rh)), rh), comm)
    col_in = MPI.Alltoallv!(MPI.VBuffer(reduce(vcat, col_out; init = Int32[]), ccounts),
                            MPI.VBuffer(Vector{Int32}(undef, sum(rc)), rc), comm)
    val_in = MPI.Alltoallv!(MPI.VBuffer(reduce(vcat, val_out; init = ComplexF64[]), vcounts),
                            MPI.VBuffer(Vector{ComplexF64}(undef, sum(rv)), rv), comm)

    # 4. 用收到的行非零元拼出 Z_sub（仅含所需行的稀疏子块，全局编号）
    I_sub = Int[]
    J_sub = Int[]
    V_sub = ComplexF64[]
    hoff = 0
    coff = 0
    voff = 0
    for q = 1:P
        for _ = 1:rh[q]÷2
            g = hdr_in[hoff+1]
            nnz = hdr_in[hoff+2]
            for t = 1:nnz
                push!(I_sub, g)
                push!(J_sub, col_in[coff+t])
                push!(V_sub, val_in[voff+t])
            end
            hoff += 2
            coff += nnz
            voff += nnz
        end
    end
    Z_sub = sparse(I_sub, J_sub, V_sub, N, N)

    # 5. 逐 owned cube 做块最小二乘（与串行 SPAIPreconditioner 同款算法）
    IM = Int[]
    JM = Int[]
    VM = ComplexF64[]
    for ic in owned
        cbfs_s, neibfs_s, neis_s = cube_sets[ic]
        cbfs = sorted_ids[cbfs_s]
        neibfs = sorted_ids[neibfs_s]
        neis = sorted_ids[neis_s]
        Znn = Matrix{CT}(Z_sub[neibfs, neis])
        ZnnHZnn = Znn * Znn'
        # 列下标：cbfs 的排序位置在 neis_s（排序位置向量）中的局部列号
        PH = lu!(ZnnHZnn) \ view(Znn, :, [searchsortedfirst(neis_s, b) for b in cbfs_s])
        PHt = PH'
        for (a, r) in enumerate(cbfs)
            for (bidx, c) in enumerate(neibfs)
                v = PHt[a, bidx]
                if abs(v) > 1e-12
                    push!(IM, r)
                    push!(JM, c)
                    push!(VM, v)
                end
            end
        end
    end

    M_local = sparse(IM, JM, VM, N, N)
    return DistributedSPAIPreconditioner{CT}(M_local, comm)
end

function apply_mpi_preconditioner!(
    y::AbstractVector{CT},
    P::DistributedSPAIPreconditioner,
    x::AbstractVector{CT},
) where {CT}
    mul!(y, P.M_local, x)
    MPI.Allreduce!(y, +, P.comm)
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
