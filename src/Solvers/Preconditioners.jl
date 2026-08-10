module Preconditioners

using LinearAlgebra
using SparseArrays
using IncompleteLU

export AbstractPreconditioner,
    DiagonalPreconditioner,
    IdentityPreconditioner,
    ILUPreconditioner,
    SPAIPreconditioner,
    BlockJacobiPreconditioner

"""
    AbstractPreconditioner

Abstract base type for preconditioners in iterative solvers.

Preconditioners accelerate iterative solution of linear systems ``\\mathbf{Ax} = \\mathbf{b}``
by approximating ``\\mathbf{A}^{-1}`` with a cheaper-to-apply operator ``\\mathbf{M}^{-1}``.

The preconditioned system is:
```math
\\mathbf{M}^{-1}\\mathbf{Ax} = \\mathbf{M}^{-1}\\mathbf{b}
```

# Interface Requirements

Preconditioners must implement:
- `ldiv!(y, P, x)`: Compute `y = P \\ x` (in-place)
- `ldiv!(P, x)`: Compute `x = P \\ x` (in-place)
- `\\(P, x)`: Compute `P \\ x` (allocating)

# See Also

- [`DiagonalPreconditioner`](@ref): Jacobi preconditioner
- [`ILUPreconditioner`](@ref): Incomplete LU factorization
- [`BlockJacobiPreconditioner`](@ref): Block-diagonal approximation
"""
abstract type AbstractPreconditioner end

# --- Identity Preconditioner ---
struct IdentityPreconditioner <: AbstractPreconditioner end
LinearAlgebra.ldiv!(y, P::IdentityPreconditioner, x) = copyto!(y, x)
LinearAlgebra.ldiv!(P::IdentityPreconditioner, x) = x
Base.:\(P::IdentityPreconditioner, x) = x

# --- Diagonal (Jacobi) Preconditioner ---
"""
    DiagonalPreconditioner{T} <: AbstractPreconditioner

Diagonal (Jacobi) preconditioner: ``M = \\text{diag}(A)``.

The simplest preconditioner, using only the diagonal of the matrix:
```math
\\mathbf{M}^{-1} = \\text{diag}(A_{11}^{-1}, A_{22}^{-1}, \\ldots, A_{nn}^{-1})
```

# Construction

    DiagonalPreconditioner(A::AbstractMatrix)

# Characteristics

- **Memory**: O(N) — stores only diagonal
- **Setup time**: O(N) — extract diagonal
- **Application time**: O(N) — element-wise division
- **Effectiveness**: Poor for ill-conditioned systems, good for diagonally dominant

# Example

```julia
# GMRES with diagonal preconditioning
Z = assemble_impedance_matrix(efie)
P = DiagonalPreconditioner(Z)
V = compute_excitation_vector(efie, plane_wave)

I = gmres(Z, V, Pl=P, abstol=1e-6)
```

# When to Use

- ✅ Quick setup, low memory
- ✅ Diagonally dominant systems
- ❌ Poor for EFIE/MFIE (off-diagonal dominant)
- ❌ Better options available (ILU, Block-Jacobi)

# See Also

- [`ILUPreconditioner`](@ref): Better for sparse systems
- [`BlockJacobiPreconditioner`](@ref): Better for MoM matrices
"""
struct DiagonalPreconditioner{T} <: AbstractPreconditioner
    diag_inv::Vector{T}
end

function DiagonalPreconditioner(A::AbstractMatrix{T}) where {T}
    n = size(A, 1)
    diag_inv = Vector{T}(undef, n)
    for i = 1:n
        d = A[i, i]
        if abs(d) < 1e-15
            # Handle zero diagonal? 
            # For EFIE, diagonal should not be zero.
            diag_inv[i] = 1.0
        else
            diag_inv[i] = 1.0 / d
        end
    end
    return DiagonalPreconditioner(diag_inv)
end

function LinearAlgebra.ldiv!(y::AbstractVector, P::DiagonalPreconditioner, x::AbstractVector)
    @inbounds for i in eachindex(y)
        y[i] = P.diag_inv[i] * x[i]
    end
    return y
end

function LinearAlgebra.ldiv!(P::DiagonalPreconditioner, x::AbstractVector)
    @inbounds for i in eachindex(x)
        x[i] *= P.diag_inv[i]
    end
    return x
end

Base.:\(P::DiagonalPreconditioner, x) = P.diag_inv .* x

# --- ILU Preconditioner ---
"""
    ILUPreconditioner{T, F} <: AbstractPreconditioner

Incomplete LU factorization preconditioner.

Approximates ``\\mathbf{A} \\approx \\mathbf{L}\\mathbf{U}`` with sparse lower/upper 
triangular factors, using drop tolerance to control fill-in.

# Construction

    ILUPreconditioner(A::SparseMatrixCSC; τ=0.01)

# Arguments

- `A::SparseMatrixCSC`: Sparse matrix to precondition
- `τ::Float64=0.01`: Drop tolerance (smaller → more accurate, more fill)
  - τ = 0: No dropping (exact LU if no pivoting)
  - τ = 0.01: Typical (1% threshold)
  - τ = 0.1: Aggressive dropping (less accurate, faster)

# Characteristics

- **Memory**: O(N × fill_factor) — depends on τ
- **Setup time**: O(N² × fill_factor) — factorization
- **Application time**: O(N × fill_factor) — triangular solves
- **Effectiveness**: Excellent for sparse systems

# Example

```julia
# GMRES with ILU(0.01) preconditioning
Z = assemble_impedance_matrix(efie)
P = ILUPreconditioner(Z, τ=0.01)
V = compute_excitation_vector(efie, plane_wave)

I = gmres(Z, V, Pl=P, abstol=1e-6, restart=50)
```

# When to Use

- ✅ Sparse matrices (EFIE/MFIE direct assembly)
- ✅ Medium-sized problems (N < 50k)
- ❌ Dense matrices (use Block-Jacobi)
- ❌ Very large problems (setup cost high)

# Notes

- Uses IncompleteLU.jl package
- Drop tolerance trades accuracy for speed/memory
- For MLFMA operators, use Block-Jacobi instead

# See Also

- [`DiagonalPreconditioner`](@ref): Simpler, faster, less effective
- [`BlockJacobiPreconditioner`](@ref): For MLFMA operators
"""
struct ILUPreconditioner{T,F} <: AbstractPreconditioner
    ilu_factor::F
end

function ILUPreconditioner(A::SparseMatrixCSC{T}; τ = 0.01) where {T}
    ilu_factor = ilu(A, τ = τ)
    return ILUPreconditioner{T,typeof(ilu_factor)}(ilu_factor)
end

function LinearAlgebra.ldiv!(y::AbstractVector, P::ILUPreconditioner, x::AbstractVector)
    ldiv!(y, P.ilu_factor, x)
    return y
end

function LinearAlgebra.ldiv!(P::ILUPreconditioner, x::AbstractVector)
    ldiv!(P.ilu_factor, x)
    return x
end

Base.:\(P::ILUPreconditioner, x) = P.ilu_factor \ x

# --- SPAI Preconditioner ---
"""
    SPAIPreconditioner(A::SparseMatrixCSC)

稀疏近似逆预条件（论文式 (2-63)~(2-68)）：构造稀疏矩阵 `M ≈ A^{-1}`，
逐列求解最小二乘问题

```math
\\min \\|\\bm{A} \\bm{M} - \\overline{\\bm{I}}\\|_F^2, \\qquad
\\min_{\\bm{m}_k} \\|\\bm{A}\\, \\bm{m}_k - \\bm{e}_k\\|_2
```

实现按列独立构造：对第 `k` 列取 `A[:, k]` 的非零模式 `J`，收集相关行
`I = {i : A[i, J] ≠ 0}`，组装稠密最小二乘系统 `A_hat * x = e_hat`
并调用 `qr!` 求解（经典 QR 路线，论文式 (2-66)）。每列问题相互独立，
天然可并行（`Threads.@threads` 按行循环，每线程独立工作区）。

论文推荐改为 LU 分解路线（论文式 (2-67)~(2-68)）：

```math
\\bm{Z}_{near}^{(C_{in}, C_{inn})}\\bm{Z}_{near}^{(C_{in}, C_{inn})H} = \\bm{L}\\bm{U}, \\qquad
(\\bm{P}_l^i)^H = (\\bm{L}\\bm{U})^{-1} \\bm{Z}_{near}^{(C_{in}, C_i)}
```

该路线分解矩阵规模更小且利用单位矩阵的稀疏性，可作为后续优化方向。

# Arguments
- `A`: 近场稀疏矩阵（一般取 `Z_near`）。

# Returns
- `SPAIPreconditioner`：应用时执行稀疏矩阵向量乘 `y = M * x`。
"""
struct SPAIPreconditioner{T} <: AbstractPreconditioner
    M::SparseMatrixCSC{T,Int}
    function SPAIPreconditioner{T}(M::SparseMatrixCSC{T,Int}) where {T}
        new{T}(M)
    end
end

mutable struct SPAIWorkspace{T}
    I_vec::Vector{Int}
    row_map::Dict{Int,Int}
    A_hat::Matrix{T}
    e_hat::Vector{T}
end

function SPAIPreconditioner(A::SparseMatrixCSC{T,Int}) where {T}
    n = size(A, 1)

    # Per-thread storage
    # Use maxthreadid() to ensure we cover all possible thread IDs
    max_tid = Threads.maxthreadid()
    I_M_thread = [Int[] for _ = 1:max_tid]
    J_M_thread = [Int[] for _ = 1:max_tid]
    V_M_thread = [T[] for _ = 1:max_tid]

    # Per-thread workspace
    workspaces = [
        SPAIWorkspace(
            Int[],
            Dict{Int,Int}(),
            Matrix{T}(undef, 64, 64), # Initial capacity
            Vector{T}(undef, 64),
        ) for _ = 1:max_tid
    ]

    Threads.@threads for k = 1:n
        tid = Threads.threadid()
        ws = workspaces[tid]

        # 1. Define sparsity pattern J for column k of M
        # We use the sparsity pattern of A[:, k]
        r_start = A.colptr[k]
        r_end = A.colptr[k+1] - 1
        J = view(A.rowval, r_start:r_end)

        if isempty(J)
            continue
        end

        # 2. Identify rows I where A[:, J] is non-zero
        empty!(ws.I_vec)
        for j in J
            rj_start = A.colptr[j]
            rj_end = A.colptr[j+1] - 1
            for idx = rj_start:rj_end
                push!(ws.I_vec, A.rowval[idx])
            end
        end
        sort!(ws.I_vec)
        unique!(ws.I_vec)
        I = ws.I_vec

        # 3. Build dense LS system A_hat * x = e_hat
        m_hat = length(I)
        n_hat = length(J)

        # Resize buffers if needed
        if size(ws.A_hat, 1) < m_hat || size(ws.A_hat, 2) < n_hat
            new_m = max(m_hat, size(ws.A_hat, 1) * 2)
            new_n = max(n_hat, size(ws.A_hat, 2) * 2)
            ws.A_hat = Matrix{T}(undef, new_m, new_n)
        end
        if length(ws.e_hat) < m_hat
            resize!(ws.e_hat, max(m_hat, length(ws.e_hat) * 2))
        end

        # Map global row indices to local 1..m_hat
        empty!(ws.row_map)
        for (idx, r) in enumerate(I)
            ws.row_map[r] = idx
        end

        # Fill A_hat (only the top-left m_hat x n_hat part)
        # Zero out the used part
        for c = 1:n_hat
            for r = 1:m_hat
                ws.A_hat[r, c] = zero(T)
            end
        end

        for (j_local, j_global) in enumerate(J)
            rj_start = A.colptr[j_global]
            rj_end = A.colptr[j_global+1] - 1
            for idx = rj_start:rj_end
                r_global = A.rowval[idx]
                if haskey(ws.row_map, r_global)
                    r_local = ws.row_map[r_global]
                    ws.A_hat[r_local, j_local] = A.nzval[idx]
                end
            end
        end

        # Fill e_hat
        for r = 1:m_hat
            ws.e_hat[r] = zero(T)
        end
        if haskey(ws.row_map, k)
            ws.e_hat[ws.row_map[k]] = one(T)
        end

        # 4. Solve LS
        A_sub = view(ws.A_hat, 1:m_hat, 1:n_hat)
        e_sub = view(ws.e_hat, 1:m_hat)

        try
            # Use qr! which is in-place on A_sub
            F = qr!(A_sub)
            x = F \ e_sub

            # 5. Store result
            for (j_local, val) in enumerate(x)
                if abs(val) > 1e-12
                    push!(I_M_thread[tid], J[j_local]) # Row index in M
                    push!(J_M_thread[tid], k)          # Col index in M
                    push!(V_M_thread[tid], val)
                end
            end
        catch
            # Ignore singular cases
        end
    end

    # Concatenate
    I_M = reduce(vcat, I_M_thread)
    J_M = reduce(vcat, J_M_thread)
    V_M = reduce(vcat, V_M_thread)

    M = sparse(I_M, J_M, V_M, n, n)
    return SPAIPreconditioner{T}(M)
end

function LinearAlgebra.ldiv!(y::AbstractVector, P::SPAIPreconditioner, x::AbstractVector)
    mul!(y, P.M, x)
    return y
end

function LinearAlgebra.ldiv!(P::SPAIPreconditioner, x::AbstractVector)
    x_new = P.M * x
    copyto!(x, x_new)
    return x
end

Base.:\(P::SPAIPreconditioner, x) = P.M * x

# --- Block Jacobi Preconditioner ---
"""
    BlockJacobiPreconditioner

Block Jacobi preconditioner for MLFMA. Extracts diagonal blocks of Z_near
corresponding to each octree cube and factorizes them independently.

Construction is embarrassingly parallel (independent LU per block), much
cheaper than sparse LU of the full Z_near.

# Constructor
    BlockJacobiPreconditioner(Z_near, block_indices)

- `Z_near::SparseMatrixCSC`: The near-field sparse impedance matrix.
- `block_indices::Vector{<:AbstractVector{Int}}`: One index set per cube, giving
  the matrix rows and columns belonging to that block. Indices need not be contiguous.
"""
struct BlockJacobiPreconditioner{T} <: AbstractPreconditioner
    blocks::Vector{LU{T,Matrix{T},Vector{Int}}}  # LU factorizations per block
    indices::Vector{Vector{Int}}                    # basis function indices
end

function BlockJacobiPreconditioner(
    Z_near::SparseMatrixCSC{T,Int},
    block_indices::AbstractVector{<:AbstractVector{Int}},
) where {T}
    # Filter out empty blocks and normalize to dense index vectors.
    non_empty = Vector{Vector{Int}}()
    for idxs in block_indices
        isempty(idxs) && continue
        push!(non_empty, collect(idxs))
    end
    n_blocks = length(non_empty)

    # Pre-allocate block storage
    blocks = Vector{LU{T,Matrix{T},Vector{Int}}}(undef, n_blocks)

    # Build per-block LU in parallel
    Threads.@threads for ib = 1:n_blocks
        idxs = non_empty[ib]
        n = length(idxs)
        # Extract dense diagonal block from sparse Z_near
        B = Matrix{T}(undef, n, n)
        @inbounds for (jj, j_global) in enumerate(idxs)
            for (ii, i_global) in enumerate(idxs)
                B[ii, jj] = Z_near[i_global, j_global]
            end
        end
        blocks[ib] = lu(B)
    end

    return BlockJacobiPreconditioner{T}(blocks, non_empty)
end

function LinearAlgebra.ldiv!(y::AbstractVector, P::BlockJacobiPreconditioner, x::AbstractVector)
    # Apply each block's LU independently (can be parallelized)
    Threads.@threads for ib in eachindex(P.blocks)
        idxs = P.indices[ib]
        x_block = view(x, idxs)
        y_block = view(y, idxs)
        ldiv!(y_block, P.blocks[ib], x_block)
    end
    return y
end

function LinearAlgebra.ldiv!(P::BlockJacobiPreconditioner, x::AbstractVector)
    Threads.@threads for ib in eachindex(P.blocks)
        idxs = P.indices[ib]
        x_block = x[idxs]  # copy needed since ldiv! overwrites
        ldiv!(view(x, idxs), P.blocks[ib], x_block)
    end
    return x
end

Base.:\(P::BlockJacobiPreconditioner, x) = begin
    y = similar(x)
    ldiv!(y, P, x)
    return y
end

end
