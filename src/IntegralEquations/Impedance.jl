module Impedance

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using StaticArrays
using LinearAlgebra
using Base.Threads

export assemble_generic

"""
    assemble_generic(operator, basis, interaction_func; symmetric=false)

Generic assembly function for impedance matrix.
`interaction_func` should have signature:
`interaction_func(Z_local, operator, tri_test::TriangleInfo, tri_source::TriangleInfo)`
"""
function assemble_generic(operator, basis::RWGBasis{IT, FT}, interaction_func; symmetric::Bool=false) where {IT, FT}
    N = num_basis(basis)
    CT = Complex{FT}
    Z = zeros(CT, N, N)
    
    mesh = basis.mesh
    nt = num_elements(mesh)
    
    # Precompute TriangleInfo
    tris_info = Vector{TriangleInfo{IT, FT}}(undef, nt)
    Threads.@threads for t in 1:nt
        tris_info[t] = get_triangle_info(mesh, basis, t)
    end
    
    lockZ = SpinLock()
    n_threads = Threads.nthreads()
    use_threads = n_threads > 1

    # Use cyclic scheduling to balance load for symmetric assembly (triangular loop)
    # Thread k processes indices k, k+N, k+2N...
    Threads.@threads for tid in 1:n_threads
        for t_test in tid:n_threads:nt
            tri_test = tris_info[t_test]
            
            # Thread-local storage
            Z_local = zeros(CT, 3, 3)
            
            # If symmetric, only loop t_source >= t_test
            start_source = symmetric ? t_test : 1
            
            for t_source in start_source:nt
                tri_source = tris_info[t_source]
                
                fill!(Z_local, zero(CT))
                interaction_func(Z_local, operator, tri_test, tri_source)
                
                # Distribute to global matrix
                if use_threads lock(lockZ) end
                try
                    for i in 1:3
                        row_idx = tri_test.inBfsID[i]
                        if row_idx != 0
                            sign_test = tri_test.bfsSign[i]
                            
                            for j in 1:3
                                col_idx = tri_source.inBfsID[j]
                                if col_idx != 0
                                    sign_src = tri_source.bfsSign[j]
                                    
                                    val = Z_local[i, j] * sign_test * sign_src
                                    Z[row_idx, col_idx] += val
                                    
                                    if symmetric && t_test != t_source
                                        Z[col_idx, row_idx] += val
                                    end
                                end
                            end
                        end
                    end
                finally
                    if use_threads unlock(lockZ) end
                end
            end
        end
    end
    
    return Z
end

end
