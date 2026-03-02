module Assembly

using MPI
using Base.Threads: SpinLock
using ...CoreModule
using ...Geometry
using ...BasisFunctions
using ...IntegralEquations
using ..MPIArrays

using ...IntegralEquations.Impedance: get_triangle_info
using ...IntegralEquations.EFIEModule: efie_interaction!
using ...IntegralEquations.MFIEModule: mfie_interaction!

export assemble_impedance_matrix_parallel

# ─── Type-stable interaction kernel ──────────────────────────────────────────
@inline function _fill_local!(Z_loc::Matrix, op::EFIE, tri_t, tri_s, ::Matrix, ::Matrix)
    efie_interaction!(Z_loc, op, tri_t, tri_s)
end
@inline function _fill_local!(Z_loc::Matrix, op::MFIE, tri_t, tri_s, ::Matrix, ::Matrix)
    mfie_interaction!(Z_loc, op, tri_t, tri_s)
end
@inline function _fill_local!(Z_loc::Matrix, op::CFIE, tri_t, tri_s, buf_e::Matrix, buf_m::Matrix)
    fill!(buf_e, zero(eltype(buf_e))); fill!(buf_m, zero(eltype(buf_m)))
    efie_interaction!(buf_e, op.efie, tri_t, tri_s)
    mfie_interaction!(buf_m, op.mfie, tri_t, tri_s)
    α = op.alpha
    @. Z_loc = α * buf_e + (1 - α) * buf_m
end

# ─── Main function ────────────────────────────────────────────────────────────
function assemble_impedance_matrix_parallel(operator, basis::RWGBasis{IT,FT}) where {IT,FT}
    comm    = MPI.COMM_WORLD
    rank    = MPI.Comm_rank(comm)
    n_procs = MPI.Comm_size(comm)

    N  = num_basis(basis)
    CT = Complex{FT}

    # Distributed matrix (column partition)
    Z          = mpiarray(CT, N, N; comm = comm)
    local_cols = Z.indices[2]

    rank == 0 && println(
        "RWG MPI Assembly (column partition): N=$N, procs=$n_procs, threads=$(Threads.nthreads())",
    )
    MPI.Barrier(comm)
    println("  Rank $rank owns columns $(first(local_cols)):$(last(local_cols))")
    MPI.Barrier(comm)
    flush(stdout)

    mesh      = basis.mesh
    nt        = num_elements(mesh)
    tris_info = [get_triangle_info(mesh, basis, t) for t = 1:nt]

    # Collect source triangles whose DOFs fall in local_cols
    # (no Set{Int} — direct UnitRange membership check)
    src_tris = Int[]
    for it = 1:nt
        tri = tris_info[it]
        for j = 1:3
            m = tri.inBfsID[j]
            if m != 0 && m in local_cols
                push!(src_tris, it)
                break
            end
        end
    end

    row_locks = [SpinLock() for _ = 1:N]
    n_threads = Threads.nthreads()
    next_idx  = Threads.Atomic{Int}(1)

    Threads.@threads for _ = 1:n_threads
        Z_loc  = zeros(CT, 3, 3)
        buf_e  = zeros(CT, 3, 3)
        buf_m  = zeros(CT, 3, 3)

        while true
            tidx = Threads.atomic_add!(next_idx, 1)
            tidx > length(src_tris) && break

            t_src    = src_tris[tidx]
            tri_src  = tris_info[t_src]

            for t_tst = 1:nt
                tri_tst = tris_info[t_tst]
                fill!(Z_loc, zero(CT))
                _fill_local!(Z_loc, operator, tri_tst, tri_src, buf_e, buf_m)

                @inbounds for j = 1:3
                    col = tri_src.inBfsID[j]
                    col == 0 && continue
                    col in local_cols || continue
                    bf_s   = basis.functions[col]
                    sign_s = (bf_s.support[1] == t_src) ? bf_s.signs[1] : bf_s.signs[2]
                    lock(row_locks[col])
                    @inbounds for i = 1:3
                        row = tri_tst.inBfsID[i]
                        row == 0 && continue
                        bf_t   = basis.functions[row]
                        sign_t = (bf_t.support[1] == t_tst) ? bf_t.signs[1] : bf_t.signs[2]
                        Z[row, col] += Z_loc[i, j] * sign_t * sign_s
                    end
                    unlock(row_locks[col])
                end
            end
        end
    end

    sync!(Z)

    rank == 0 && println("RWG MPI Assembly Completed.")
    rank == 0 && flush(stdout)

    return Z
end

end
