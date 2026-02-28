module Assembly

using MPI
using ...CoreModule
using ...Geometry
using ...BasisFunctions
using ...IntegralEquations
using ..MPIArrays

using ...IntegralEquations.Impedance: get_triangle_info
using ...IntegralEquations.EFIEModule: efie_interaction!
using ...IntegralEquations.MFIEModule: mfie_interaction!

export assemble_impedance_matrix_parallel

function assemble_impedance_matrix_parallel(operator, basis::RWGBasis{IT,FT}) where {IT,FT}
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    n_procs = MPI.Comm_size(comm)

    N = num_basis(basis)
    CT = Complex{FT}

    # Create distributed matrix (partitioned by columns)
    Z = mpiarray(CT, N, N; comm = comm)

    # Get local column range
    # Z.indices is (1:N, col_start:col_end)
    local_cols = Z.indices[2]

    # Identify source triangles relevant to local columns
    # We use a Set to avoid duplicates
    src_tris = Set{Int}()
    for col in local_cols
        bf = basis.functions[col]
        if bf.support[1] != 0
            push!(src_tris, bf.support[1])
        end
        if bf.support[2] != 0
            push!(src_tris, bf.support[2])
        end
    end

    # Convert to sorted vector for deterministic iteration
    src_tris_vec = sort(collect(src_tris))

    mesh = basis.mesh
    nt = num_elements(mesh)

    # Precompute TriangleInfo
    tris_info = [get_triangle_info(mesh, basis, t) for t = 1:nt]

    Z_local = zeros(CT, 3, 3)
    Z_efie_local = zeros(CT, 3, 3)  # Workspace for CFIE
    Z_mfie_local = zeros(CT, 3, 3)  # Workspace for CFIE

    # Loop over RELEVANT source triangles
    for t_source in src_tris_vec
        tri_source = tris_info[t_source]

        # Loop over ALL test triangles (since we own full rows)
        for t_test = 1:nt
            tri_test = tris_info[t_test]

            # Compute interaction
            fill!(Z_local, zero(CT))

            if isa(operator, CFIE)
                # CFIE = alpha * EFIE + (1-alpha) * MFIE
                fill!(Z_efie_local, zero(CT))
                fill!(Z_mfie_local, zero(CT))
                efie_interaction!(Z_efie_local, operator.efie, tri_test, tri_source)
                mfie_interaction!(Z_mfie_local, operator.mfie, tri_test, tri_source)
                α = operator.alpha
                @. Z_local = α * Z_efie_local + (1 - α) * Z_mfie_local
            elseif isa(operator, EFIE)
                efie_interaction!(Z_local, operator, tri_test, tri_source)
            elseif isa(operator, MFIE)
                mfie_interaction!(Z_local, operator, tri_test, tri_source)
            else
                error("Unknown operator type for parallel assembly")
            end

            # Distribute to global matrix
            for j = 1:3 # Source (Columns)
                col_idx = tri_source.inBfsID[j]

                # Check if this column is local
                if col_idx in local_cols
                    if col_idx != 0
                        bf_src = basis.functions[col_idx]
                        sign_src =
                            (bf_src.support[1] == t_source) ? bf_src.signs[1] : bf_src.signs[2]

                        for i = 1:3 # Test (Rows)
                            row_idx = tri_test.inBfsID[i]
                            if row_idx != 0
                                bf = basis.functions[row_idx]
                                sign_test = (bf.support[1] == t_test) ? bf.signs[1] : bf.signs[2]

                                # Add to Z
                                # Z is MPIArray. Z[row, col] works if col is local.
                                Z[row_idx, col_idx] += Z_local[i, j] * sign_test * sign_src
                            end
                        end
                    end
                end
            end
        end
    end

    sync!(Z) # Ensure ghost data is synced (if any)
    return Z
end

end
