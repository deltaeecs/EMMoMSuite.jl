module Preconditioners

using LinearAlgebra
using SparseArrays
using IncompleteLU

export AbstractPreconditioner, DiagonalPreconditioner, IdentityPreconditioner, ILUPreconditioner, SPAIPreconditioner

abstract type AbstractPreconditioner end

# --- Identity Preconditioner ---
struct IdentityPreconditioner <: AbstractPreconditioner end
LinearAlgebra.ldiv!(y, P::IdentityPreconditioner, x) = copyto!(y, x)
LinearAlgebra.ldiv!(P::IdentityPreconditioner, x) = x
Base.:\(P::IdentityPreconditioner, x) = x

# --- Diagonal (Jacobi) Preconditioner ---
struct DiagonalPreconditioner{T} <: AbstractPreconditioner
    diag_inv::Vector{T}
end

function DiagonalPreconditioner(A::AbstractMatrix{T}) where T
    n = size(A, 1)
    diag_inv = Vector{T}(undef, n)
    for i in 1:n
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
struct ILUPreconditioner{T, F} <: AbstractPreconditioner
    ilu_factor::F
end

function ILUPreconditioner(A::SparseMatrixCSC{T}; τ=0.01) where T
    ilu_factor = ilu(A, τ=τ)
    return ILUPreconditioner{T, typeof(ilu_factor)}(ilu_factor)
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
struct SPAIPreconditioner{T} <: AbstractPreconditioner
    M::SparseMatrixCSC{T, Int}
    function SPAIPreconditioner{T}(M::SparseMatrixCSC{T, Int}) where T
        new{T}(M)
    end
end

mutable struct SPAIWorkspace{T}
    I_vec::Vector{Int}
    row_map::Dict{Int,Int}
    A_hat::Matrix{T}
    e_hat::Vector{T}
end

function SPAIPreconditioner(A::SparseMatrixCSC{T, Int}) where T
    n = size(A, 1)
    
    # Per-thread storage
    # Use maxthreadid() to ensure we cover all possible thread IDs
    max_tid = Threads.maxthreadid()
    I_M_thread = [Int[] for _ in 1:max_tid]
    J_M_thread = [Int[] for _ in 1:max_tid]
    V_M_thread = [T[] for _ in 1:max_tid]
    
    # Per-thread workspace
    workspaces = [
        SPAIWorkspace(
            Int[], 
            Dict{Int,Int}(),
            Matrix{T}(undef, 64, 64), # Initial capacity
            Vector{T}(undef, 64)
        ) 
        for _ in 1:max_tid
    ]
    
    Threads.@threads for k in 1:n
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
            for idx in rj_start:rj_end
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
        for c in 1:n_hat
            for r in 1:m_hat
                ws.A_hat[r, c] = zero(T)
            end
        end
        
        for (j_local, j_global) in enumerate(J)
            rj_start = A.colptr[j_global]
            rj_end = A.colptr[j_global+1] - 1
            for idx in rj_start:rj_end
                r_global = A.rowval[idx]
                if haskey(ws.row_map, r_global)
                    r_local = ws.row_map[r_global]
                    ws.A_hat[r_local, j_local] = A.nzval[idx]
                end
            end
        end
        
        # Fill e_hat
        for r in 1:m_hat
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

end
