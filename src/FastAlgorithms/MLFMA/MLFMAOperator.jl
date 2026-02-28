module MLFMAOperatorModule

using LinearAlgebra
using StaticArrays
using SparseArrays
using ....CoreModule
using ....Geometry
using ....BasisFunctions
using ....IntegralEquations
using ..Octree
using ..OctreeBuilder
using ..Aggregation
using ..Translation
using ..Disaggregation
using ..Level

using ....IntegralEquations.Impedance: get_triangle_info, get_triangles_info
using ....IntegralEquations.EFIEModule: efie_interaction!, EFIE
using ....IntegralEquations.MFIEModule: mfie_interaction!, MFIE
using ....IntegralEquations.CFIEModule: CFIE
using ....IntegralEquations.VEFIEModule: vefie_element_interaction, vefie_element_interaction_kernel, precompute_vefie_basis, vefie_mass_matrix_cached, TetBasisCache, get_tetrahedra_info, VEFIE
using ....IntegralEquations.SCFIEModule: SCFIE, scfie_coupling_interaction, assemble_fss_boundary_correction_sparse

export MLFMAOperator, mul!

"""
    MLFMAOperator

A linear operator representing the impedance matrix Z, computed using MLFMA.
Z = Z_near + Z_far
Z_far is computed via MLFMA (Aggregation -> Translation -> Disaggregation).
Z_near is computed directly and stored as a sparse matrix.
"""
struct MLFMAOperator{FT, CT} <: AbstractIntegralOperator
    octree::OctreeInfo
    bases::Vector{AbstractBasisFunction}
    basis_offsets::Vector{Int} # Cumulative counts: [N1, N1+N2, ...]
    Z_near::SparseMatrixCSC{CT, Int}
    operator::AbstractIntegralOperator # Underlying EFIE/MFIE/VEFIE/SCFIE
    sorted_ids::Vector{Int} # Permutation to map original basis to sorted basis
    inv_sorted_ids::Vector{Int} # Inverse permutation
end

function MLFMAOperator(operator::AbstractIntegralOperator, basis::AbstractBasisFunction, leafCubeEdgel::Float64)
    return MLFMAOperator(operator, [basis], leafCubeEdgel)
end

function MLFMAOperator(operator::AbstractIntegralOperator, bases::Vector{<:AbstractBasisFunction}, leafCubeEdgel::Float64)
    # 1. Build Octree
    # Concatenate centers from all bases
    bf_centers_list = [reduce(hcat, [bf.center for bf in b.functions]) for b in bases]
    bf_centers = reduce(hcat, bf_centers_list)
    
    lambda = 299792458.0 / operator.freq
    octree, sorted_ids = build_octree(bf_centers, leafCubeEdgel; λ=lambda)
    
    # Inverse permutation
    N = length(sorted_ids)
    inv_sorted_ids = zeros(Int, N)
    for i in 1:N
        inv_sorted_ids[sorted_ids[i]] = i
    end
    
    # Basis offsets
    basis_offsets = cumsum([num_basis(b) for b in bases])
    
    # 2. Compute Near Field Matrix
    println("Assembling Near Field Matrix...")
    # Cast bases to Vector{AbstractBasisFunction}
    abstract_bases = Vector{AbstractBasisFunction}(bases)
    Z_near = assemble_near_field(operator, abstract_bases, basis_offsets, octree, sorted_ids, inv_sorted_ids)
    
    # Add Fss boundary correction for SCFIE operators
    if operator isa SCFIE && length(bases) >= 2
        surf_basis = nothing
        vol_basis = nothing
        for b in bases
            if b isa RWGBasis
                surf_basis = b
            elseif b isa SWGBasis
                vol_basis = b
            end
        end
        if surf_basis !== nothing && vol_basis !== nothing
            Z_fss = assemble_fss_boundary_correction_sparse(operator, surf_basis, vol_basis)
            Z_near = Z_near + Z_fss
        end
    end
    
    FT = eltype(bases[1].mesh.node)
    CT = eltype(Z_near)
    
    return MLFMAOperator{FT, CT}(octree, abstract_bases, basis_offsets, Z_near, operator, sorted_ids, inv_sorted_ids)
end

Base.eltype(op::MLFMAOperator{FT, CT}) where {FT, CT} = CT
Base.size(op::MLFMAOperator) = size(op.Z_near)
Base.size(op::MLFMAOperator, i::Int) = size(op.Z_near, i)

function Base.:*(A::MLFMAOperator, x::AbstractVector)
    y = similar(x)
    mul!(y, A, x)
    return y
end

"""
    mul!(y::AbstractVector, A::MLFMAOperator, x::AbstractVector)

Compute y = A * x using MLFMA.
"""
function LinearAlgebra.mul!(y::AbstractVector, A::MLFMAOperator, x::AbstractVector)
    # 1. Near Field: y = Z_near * x
    mul!(y, A.Z_near, x)
    
    # 2. Far Field (MLFMA)
    
    # 2.1 Aggregation (Upward Pass)
    # Computes radiation patterns at leaf level and aggregates up to level 2
    aggregate!(A.octree, A.bases, A.basis_offsets, A.operator, x, A.sorted_ids)
    
    # 2.2 Translation (Horizontal Pass)
    # Computes received fields at each level from far neighbors
    for levelID in 2:A.octree.nLevels
        level = A.octree.levels[levelID]
        translate!(level)
    end
    
    # 2.3 Disaggregation (Downward Pass)
    # Propagates received fields from parent to child
    for levelID in 2:(A.octree.nLevels - 1)
        parentLevel = A.octree.levels[levelID]
        childLevel = A.octree.levels[levelID + 1]
        disaggregate_downward!(parentLevel, childLevel)
    end
    
    # 2.4 Final Collection (Leaf Level)
    # Projects received fields onto test functions
    y_far = zeros(ComplexF64, length(x))
    leafLevel = A.octree.levels[A.octree.nLevels]
    disaggregate_leaf!(leafLevel, A.bases, A.basis_offsets, A.operator, y_far, A.sorted_ids)
    
    # Apply EFIE factor to Far Field part if needed
    # Note: The factor is multiplied by 4 because the aggregation/disaggregation
    # each use l/2 (giving l²/4), but efie.factor = jkη/(16π) already includes
    # a 1/4 from the RWG normalization in the Direct solver code path.
    # Without the ×4 correction, the far-field would be 4× too small.
    # Legacy code avoids this by using -jk/(16π²) in translation + jkη in disagg.
    if hasfield(typeof(A.operator), :factor)
        y_far .*= (4 * A.operator.factor)
    end
    
    # Add to y
    y .+= y_far
    
    return y
end

function get_basis_index(global_idx::Int, offsets::Vector{Int})
    idx = searchsortedfirst(offsets, global_idx)
    local_idx = idx == 1 ? global_idx : global_idx - offsets[idx-1]
    return idx, local_idx
end

function assemble_near_field(operator, bases::Vector{AbstractBasisFunction}, offsets::Vector{Int}, octree, sorted_ids, inv_sorted_ids)
    N = offsets[end]
    leaf_level = octree.levels[octree.nLevels]
    n_cubes = length(leaf_level.cubes)
    
    # 1. Build Element -> Basis Maps
    max_tri_id = 0
    max_tet_id = 0
    
    for b in bases
        if b isa RWGBasis
            max_tri_id = max(max_tri_id, num_elements(b.mesh))
        elseif b isa SWGBasis
            max_tet_id = max(max_tet_id, num_elements(b.mesh))
        end
    end
    
    tri_to_rwg = [Vector{Tuple{Int, Int, Float64}}() for _ in 1:max_tri_id]
    tet_to_swg = [Vector{Tuple{Int, Int, Float64}}() for _ in 1:max_tet_id]
    
    all_tris = Vector{TriangleInfo{Int, Float64}}()
    all_tets = Vector{TetrahedraInfo{Int, Float64, ComplexF64}}()
    
    for (b_idx, b) in enumerate(bases)
        offset = b_idx == 1 ? 0 : offsets[b_idx-1]
        
        if b isa RWGBasis
            for (i, f) in enumerate(b.functions)
                global_id = offset + i
                for k in 1:2
                    tri_id = f.support[k]
                    if tri_id > 0
                        push!(tri_to_rwg[tri_id], (f.local_edge_idx[k], global_id, f.signs[k]))
                    end
                end
            end
            append!(all_tris, get_triangles_info(b.mesh, b))
        elseif b isa SWGBasis
            for (i, f) in enumerate(b.functions)
                global_id = offset + i
                for k in 1:2
                    tet_id = f.support[k]
                    if tet_id > 0
                        push!(tet_to_swg[tet_id], (f.local_face_idx[k], global_id, f.signs[k]))
                    end
                end
            end
            if hasfield(typeof(operator), :permittivities)
                append!(all_tets, get_tetrahedra_info(b.mesh, b, operator.permittivities))
            elseif operator isa VEFIE
                append!(all_tets, get_tetrahedra_info(b.mesh, b, operator.permittivities))
            end
        end
    end
    
    # 2. Precompute Cube -> Unique Elements
    cube_tris_vec = [Int[] for _ in 1:n_cubes]
    cube_tets_vec = [Int[] for _ in 1:n_cubes]
    
    for i_cube in 1:n_cubes
        cube = leaf_level.cubes[i_cube]
        if isempty(cube.bfInterval) continue end
        
        tris_set = Set{Int}()
        tets_set = Set{Int}()
        
        for sorted_idx in cube.bfInterval
            global_id = sorted_ids[sorted_idx]
            b_idx, local_idx = get_basis_index(global_id, offsets)
            b = bases[b_idx]
            f = b.functions[local_idx]
            
            if b isa RWGBasis
                for k in 1:2
                    if f.support[k] > 0 push!(tris_set, f.support[k]) end
                end
            elseif b isa SWGBasis
                for k in 1:2
                    if f.support[k] > 0 push!(tets_set, f.support[k]) end
                end
            end
        end
        cube_tris_vec[i_cube] = collect(tris_set)
        cube_tets_vec[i_cube] = collect(tets_set)
    end
    
    # 3. Pre-create operators
    efie_op = nothing
    mfie_op = nothing
    vefie_op = nothing
    vefie_caches = nothing
    
    if operator isa SCFIE
        efie_op = EFIE(operator.freq)
        mfie_op = MFIE(operator.freq)
        vefie_op = VEFIE(operator.freq, operator.permittivities)
    elseif operator isa VEFIE
        vefie_op = operator
    end
    
    # Pre-compute VEFIE caches for all tets (needed for vefie_element_interaction_kernel)
    if vefie_op !== nothing && !isempty(all_tets)
        vefie_caches = precompute_vefie_basis(vefie_op, all_tets)
    end
    
    # 4. Assembly Loop
    n_threads = Threads.nthreads()
    max_tid = Threads.maxthreadid()
    
    # Pre-estimate nnz: each triangle pair contributes up to 9 entries.
    # Total near-field pairs ≈ sum(|cube_tris| × |neigh_tris| × |neighbors|)
    estimated_nnz_per_thread = 0
    non_empty_cubes = 0
    for i_c in 1:n_cubes
        cube_c = leaf_level.cubes[i_c]
        if isempty(cube_c.bfInterval) continue end
        non_empty_cubes += 1
        n_my = length(cube_tris_vec[i_c]) + length(cube_tets_vec[i_c])
        for ni in cube_c.neighbors
            n_neigh = length(cube_tris_vec[ni]) + length(cube_tets_vec[ni])
            estimated_nnz_per_thread += n_my * n_neigh * 9
        end
    end
    estimated_nnz_per_thread = max(1024, div(estimated_nnz_per_thread, max(n_threads, 1)))
    
    Is = [Vector{Int}(undef, estimated_nnz_per_thread) for _ in 1:max_tid]
    Js = [Vector{Int}(undef, estimated_nnz_per_thread) for _ in 1:max_tid]
    Vs = [Vector{ComplexF64}(undef, estimated_nnz_per_thread) for _ in 1:max_tid]
    counts = zeros(Int, max_tid)  # actual count per thread
    
    println("Assembling Near Field Matrix (Optimized) with $n_threads threads...")
    println("  Non-empty Cubes: $non_empty_cubes, est nnz/thread: $estimated_nnz_per_thread")
    
    # Progress counter
    counter = Threads.Atomic{Int}(0)
    total_cubes = n_cubes
    
    # Determine if CFIE/SCFIE needs Z_mfie buffer
    needs_mfie = operator isa CFIE || operator isa SCFIE
    
    Threads.@threads for i_cube in 1:n_cubes
        tid = Threads.threadid()
        
        # Progress update (every 200 cubes to reduce overhead)
        c = Threads.atomic_add!(counter, 1)
        if c % 200 == 0 || c == total_cubes
            print("\rProgress: $c / $total_cubes cubes")
        end
        
        cube = leaf_level.cubes[i_cube]
        if isempty(cube.bfInterval) continue end
        
        my_tris = cube_tris_vec[i_cube]
        my_tets = cube_tets_vec[i_cube]
        
        # Thread-local 3×3 interaction buffers (reused per pair)
        Z_local = @MMatrix zeros(ComplexF64, 3, 3)
        Z_mfie_buf = needs_mfie ? @MMatrix(zeros(ComplexF64, 3, 3)) : Z_local
        
        neighbors = cube.neighbors
        for neighbor_idx in neighbors
            neighbor_cube = leaf_level.cubes[neighbor_idx]
            if isempty(neighbor_cube.bfInterval) continue end
            
            neigh_tris = cube_tris_vec[neighbor_idx]
            neigh_tets = cube_tets_vec[neighbor_idx]
            
            # 1. Surface-Surface (Tri-Tri)
            if !isempty(my_tris) && !isempty(neigh_tris)
                for t_test in my_tris
                    tri_test = all_tris[t_test]
                    for t_src in neigh_tris
                        tri_src = all_tris[t_src]
                        
                        fill!(Z_local, zero(ComplexF64))
                        if operator isa SCFIE
                             efie_interaction!(Z_local, efie_op, tri_test, tri_src)
                             fill!(Z_mfie_buf, zero(ComplexF64))
                             mfie_interaction!(Z_mfie_buf, mfie_op, tri_test, tri_src)
                             alpha = operator.alpha
                             # Note: mfie_interaction! already includes eta internally
                             # Match CFIE formula: alpha * Z_efie + (1-alpha) * Z_mfie
                             @inbounds for idx in eachindex(Z_local)
                                 Z_local[idx] = alpha * Z_local[idx] + (1.0 - alpha) * Z_mfie_buf[idx]
                             end
                        elseif operator isa CFIE
                             efie_interaction!(Z_local, operator.efie, tri_test, tri_src)
                             fill!(Z_mfie_buf, zero(ComplexF64))
                             mfie_interaction!(Z_mfie_buf, operator.mfie, tri_test, tri_src)
                             @inbounds for idx in eachindex(Z_local)
                                 Z_local[idx] = operator.alpha * Z_local[idx] + (1.0 - operator.alpha) * Z_mfie_buf[idx]
                             end
                        elseif operator isa MFIE
                             mfie_interaction!(Z_local, operator, tri_test, tri_src)
                        else
                             efie_interaction!(Z_local, operator, tri_test, tri_src)
                        end
                        
                        distribute_term!(Is, Js, Vs, counts, tid,
                                         Z_local, 
                                         tri_to_rwg[t_test], tri_to_rwg[t_src], 
                                         cube.bfInterval, neighbor_cube.bfInterval, 
                                         inv_sorted_ids)
                    end
                end
            end
            
            # 2. Volume-Volume (Tet-Tet)
            # Note: vefie_element_interaction_kernel already includes bfsSign in its output,
            #       so we use distribute_term_nosign! to avoid double-counting.
            if !isempty(my_tets) && !isempty(neigh_tets)
                for t_test in my_tets
                    tet_test = all_tets[t_test]
                    cache_test = vefie_caches[t_test]
                    for t_src in neigh_tets
                        tet_src = all_tets[t_src]
                        cache_src = vefie_caches[t_src]
                        
                        Z_ts = vefie_element_interaction_kernel(vefie_op, tet_test, tet_src, cache_test, cache_src)
                        
                        # Add mass matrix for self-interaction (same tet)
                        if t_test == t_src
                            M_t = vefie_mass_matrix_cached(vefie_op, tet_test, cache_test)
                            Z_ts = Z_ts .+ M_t
                        end
                        
                        distribute_term_nosign!(Is, Js, Vs, counts, tid,
                                         Z_ts, 
                                         tet_to_swg[t_test], tet_to_swg[t_src], 
                                         cube.bfInterval, neighbor_cube.bfInterval, 
                                         inv_sorted_ids)
                    end
                end
            end
            
            # 3. Surface-Volume (Tri-Tet)
            # Note: scfie_coupling_interaction already includes bfsSign for both
            #       tri (via tri.bfsSign) and tet (via tet.bfsSign), so use nosign.
            if operator isa SCFIE
                if !isempty(my_tris) && !isempty(neigh_tets)
                    for t_test in my_tris
                        tri_test = all_tris[t_test]
                        for t_src in neigh_tets
                            tet_src = all_tets[t_src]
                            
                            Z_sv, _ = scfie_coupling_interaction(operator, tri_test, tet_src)
                            distribute_term_nosign!(Is, Js, Vs, counts, tid,
                                             Z_sv, 
                                             tri_to_rwg[t_test], tet_to_swg[t_src], 
                                             cube.bfInterval, neighbor_cube.bfInterval, 
                                             inv_sorted_ids)
                        end
                    end
                end
                
                # 4. Volume-Surface (Tet-Tri)
                if !isempty(my_tets) && !isempty(neigh_tris)
                    for t_test in my_tets
                        tet_test = all_tets[t_test]
                        for t_src in neigh_tris
                            tri_src = all_tris[t_src]
                            
                            _, Z_vs = scfie_coupling_interaction(operator, tri_src, tet_test)
                            
                            distribute_term_nosign!(Is, Js, Vs, counts, tid,
                                             Z_vs, 
                                             tet_to_swg[t_test], tri_to_rwg[t_src], 
                                             cube.bfInterval, neighbor_cube.bfInterval, 
                                             inv_sorted_ids)
                        end
                    end
                end
            end
        end
    end
    
    # Merge results — trim to actual counts
    for tid in 1:max_tid
        ct = counts[tid]
        resize!(Is[tid], ct)
        resize!(Js[tid], ct)
        resize!(Vs[tid], ct)
    end
    I_total = reduce(vcat, Is)
    J_total = reduce(vcat, Js)
    V_total = reduce(vcat, Vs)
    
    return sparse(I_total, J_total, V_total, N, N)
end

# Specialized interaction functions to ensure type stability

"""
    _ensure_capacity!(Is, Js, Vs, counts, tid, needed)

Ensure COO arrays have at least `counts[tid] + needed` capacity.
Uses amortized doubling to minimize reallocations.
"""
@inline function _ensure_capacity!(Is::Vector{Vector{Int}}, Js::Vector{Vector{Int}}, 
                                    Vs::Vector{Vector{ComplexF64}}, counts::Vector{Int}, 
                                    tid::Int, needed::Int)
    required = counts[tid] + needed
    cap = length(Is[tid])
    if required > cap
        new_cap = max(required, cap * 2)
        resize!(Is[tid], new_cap)
        resize!(Js[tid], new_cap)
        resize!(Vs[tid], new_cap)
    end
    nothing
end

function distribute_term!(Is, Js, Vs, counts, tid,
                          Z_local, 
                          test_bases, src_bases, 
                          test_interval, src_interval, 
                          inv_sorted_ids)
    
    # test_bases: Vector of (local_idx, global_basis_id, sign)
    # src_bases: Vector of (local_idx, global_basis_id, sign)
    
    # Pre-check capacity (up to 9 entries per pair)
    _ensure_capacity!(Is, Js, Vs, counts, tid, length(test_bases) * length(src_bases))
    
    @inbounds for (loc_test, glob_test, sign_test) in test_bases
        sorted_idx_test = inv_sorted_ids[glob_test]
        if !(sorted_idx_test in test_interval)
            continue
        end
        
        for (loc_src, glob_src, sign_src) in src_bases
            sorted_idx_src = inv_sorted_ids[glob_src]
            
            if !(sorted_idx_src in src_interval)
                continue
            end
            
            val = Z_local[loc_test, loc_src] * sign_test * sign_src
            
            if abs(val) > 1e-16
                ct = counts[tid] + 1
                counts[tid] = ct
                Is[tid][ct] = glob_test
                Js[tid][ct] = glob_src
                Vs[tid][ct] = val
            end
        end
    end
end

"""
    distribute_term_nosign!(Is, Js, Vs, counts, tid, Z_local, test_bases, src_bases, ...)

Same as `distribute_term!` but does NOT multiply by basis function signs.
Used when the element-level interaction function already includes bfsSign
(e.g., VEFIE, SCFIE coupling).
"""
function distribute_term_nosign!(Is, Js, Vs, counts, tid,
                          Z_local, 
                          test_bases, src_bases, 
                          test_interval, src_interval, 
                          inv_sorted_ids)
    
    _ensure_capacity!(Is, Js, Vs, counts, tid, length(test_bases) * length(src_bases))
    
    @inbounds for (loc_test, glob_test, _) in test_bases
        sorted_idx_test = inv_sorted_ids[glob_test]
        if !(sorted_idx_test in test_interval)
            continue
        end
        
        for (loc_src, glob_src, _) in src_bases
            sorted_idx_src = inv_sorted_ids[glob_src]
            
            if !(sorted_idx_src in src_interval)
                continue
            end
            
            val = Z_local[loc_test, loc_src]
            
            if abs(val) > 1e-16
                ct = counts[tid] + 1
                counts[tid] = ct
                Is[tid][ct] = glob_test
                Js[tid][ct] = glob_src
                Vs[tid][ct] = val
            end
        end
    end
end

end
