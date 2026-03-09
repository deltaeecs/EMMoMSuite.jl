module Impedance

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using StaticArrays
using LinearAlgebra
using Base.Threads

export assemble_generic, assemble_generic!

"""
    assemble_generic(operator, basis, interaction_func; symmetric=false)

Generic assembly function for impedance matrix.
`interaction_func` should have signature:
`interaction_func(Z_local, operator, tri_test::TriangleInfo, tri_source::TriangleInfo)`

## Threading Strategy: Per-row SpinLock Array

**Problem**: A single global SpinLock serializes all Z[m,n] writes across threads,
causing ~15-25% overhead at N≈14559 (100% contention probability).

**Solution**: Replace the single global lock with N per-row SpinLocks.
- Contention probability: ~nthreads/N ≈ 0.03% — effectively zero.
- Each lock is held for only 3 array writes (~30ns critical section).
- Memory overhead: N × sizeof(SpinLock) ≈ negligible.
"""
function assemble_generic!(
    Z::AbstractMatrix{CT},
    operator,
    basis::RWGBasis{IT,FT},
    interaction_func;
    symmetric::Bool = false,
    accumulate::Bool = false,
) where {IT,FT,CT}
    N = num_basis(basis)
    size(Z) == (N, N) || throw(DimensionMismatch("目标矩阵尺寸必须为 $(N)×$(N)，当前为 $(size(Z))"))
    accumulate || fill!(Z, zero(CT))

    mesh = basis.mesh
    nt = num_elements(mesh)

    # Precompute TriangleInfo
    tris_info = Vector{TriangleInfo{IT,FT}}(undef, nt)
    Threads.@threads for t = 1:nt
        tris_info[t] = get_triangle_info(mesh, basis, t)
    end

    n_threads = Threads.nthreads()

    # Per-row SpinLocks: N locks instead of 1 global lock.
    # Contention probability ≈ nthreads / N ≈ 0.03% — effectively zero.
    row_locks = [SpinLock() for _ = 1:N]

    # Cyclic scheduling for balanced load in symmetric (triangular) loops
    Threads.@threads for tid = 1:n_threads
        # Thread-local 3×3 interaction buffer
        Z_local = zeros(CT, 3, 3)

        for t_test = tid:n_threads:nt
            tri_test = tris_info[t_test]

            start_source = symmetric ? t_test : 1

            for t_source = start_source:nt
                tri_source = tris_info[t_source]

                fill!(Z_local, zero(CT))
                interaction_func(Z_local, operator, tri_test, tri_source)

                # --- Forward writes: Z[test_bf, source_bf] ---
                # Lock per test-BF row (3 locks per triangle pair)
                @inbounds for i = 1:3
                    row_idx = tri_test.inBfsID[i]
                    row_idx == 0 && continue
                    sign_test = tri_test.bfsSign[i]

                    lock(row_locks[row_idx])
                    for j = 1:3
                        col_idx = tri_source.inBfsID[j]
                        if col_idx != 0
                            Z[row_idx, col_idx] += Z_local[i, j] * sign_test * tri_source.bfsSign[j]
                        end
                    end
                    unlock(row_locks[row_idx])
                end

                # --- Symmetric transpose writes: Z[source_bf, test_bf] ---
                if symmetric && t_test != t_source
                    @inbounds for j = 1:3
                        col_idx = tri_source.inBfsID[j]
                        col_idx == 0 && continue
                        sign_src = tri_source.bfsSign[j]

                        lock(row_locks[col_idx])
                        for i = 1:3
                            row_idx = tri_test.inBfsID[i]
                            if row_idx != 0
                                Z[col_idx, row_idx] +=
                                    Z_local[i, j] * tri_test.bfsSign[i] * sign_src
                            end
                        end
                        unlock(row_locks[col_idx])
                    end
                end
            end
        end
    end

    return Z
end

function assemble_generic(
    operator,
    basis::RWGBasis{IT,FT},
    interaction_func;
    symmetric::Bool = false,
) where {IT,FT}
    N = num_basis(basis)
    CT = Complex{FT}
    Z = zeros(CT, N, N)
    return assemble_generic!(Z, operator, basis, interaction_func; symmetric)
end

end
