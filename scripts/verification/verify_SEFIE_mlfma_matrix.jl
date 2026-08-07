using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.MLFMA
using LinearAlgebra
using Printf
using Random
using SparseArrays

function verify_mlfma_matrix()
    println("==================================================")
    println("   Verification: SEFIE MLFMA vs Direct Matrix     ")
    println("==================================================")

    # 1. Parameters
    freq = 300e6
    lambda = 299792458.0 / freq
    radius = 2.0 # Larger radius
    
    # 2. Mesh
    n_theta = 12
    n_phi = 24
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    println("Mesh: $(num_vertices(mesh)) vertices, $(num_elements(mesh)) elements")

    # 3. Basis
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")

    # 4. Operator
    efie = EFIE(freq)

    # 5. Direct Assembly
    println("\n[Direct] Assembling Impedance Matrix...")
    t_direct = @elapsed begin
        Z_direct = assemble_impedance_matrix(efie, basis)
    end
    println("        Done in $(t_direct) s")

    # 6. MLFMA Setup
    println("\n[MLFMA] Setting up Operator...")
    # Use small box size to force multiple levels
    leafCubeEdgel = 0.25 * lambda
    mlfma_op = MLFMAOperator(efie, basis, leafCubeEdgel)
    println("        Levels: $(mlfma_op.octree.nLevels)")
    println("        L (Leaf): $(mlfma_op.octree.levels[mlfma_op.octree.nLevels].L)")

    # 7. Comparison
    println("\n[Comparison] Computing Matrix-Vector Product...")
    Random.seed!(1234)
    x = randn(ComplexF64, N)
    # x = zeros(ComplexF64, N)
    # x[1] = 1.0
    # println("Using x = e_1 (First Basis Function)")
    
    # Direct MVP
    y_direct = Z_direct * x
    
    # MLFMA MVP
    y_mlfma = mlfma_op * x
    
    # Analyze components
    y_near = mlfma_op.Z_near * x
    y_far = y_mlfma - y_near
    
    println(@sprintf("        Near Field Norm: %.4e", norm(y_near)))
    println(@sprintf("        Far Field Norm:  %.4e", norm(y_far)))
    
    # Check Z_near consistency
    println("\n[Consistency Check] Z_near vs Z_direct (subset)...")
    rows, cols, vals = findnz(mlfma_op.Z_near)
    n_check = min(10, length(vals))
    
    for i in 1:n_check
        r = rows[i]
        c = cols[i]
        v_near = vals[i]
        v_direct = Z_direct[r, c]
        
        println(@sprintf("    (%d, %d): Near=%.2e, Direct=%.2e, Ratio=%.4f", 
            r, c, abs(v_near), abs(v_direct), abs(v_near)/abs(v_direct)))
    end
    
    # Error
    diff = y_direct - y_mlfma
    rel_err = norm(diff) / norm(y_direct)
    
    println(@sprintf("        Direct Norm: %.4e", norm(y_direct)))
    println(@sprintf("        MLFMA Norm:  %.4e", norm(y_mlfma)))
    println(@sprintf("        Diff Norm:   %.4e", norm(diff)))
    println(@sprintf("        Rel Error:   %.4e", rel_err))
    
    # Find max error index
    diff_mag = abs.(diff)
    max_err, idx = findmax(diff_mag)
    println("\n[Error Analysis]")
    println("Max Error at Index $idx: $max_err")
    println("Z_direct[$idx, 1]: ", y_direct[idx])
    println("Z_mlfma[$idx, 1]:  ", y_mlfma[idx])
    
    # Check distance
    c1 = basis.functions[1].center
    ci = basis.functions[idx].center
    dist = norm(c1 - ci)
    println("Distance: $dist")
    
    # Check if Near
    is_near = mlfma_op.Z_near[idx, 1] != 0.0
    println("Is Near: $is_near")
    if is_near
        println("Z_near[$idx, 1]: ", mlfma_op.Z_near[idx, 1])
    end
    
    # Check Cubes
    s1 = mlfma_op.inv_sorted_ids[1]
    s_idx = mlfma_op.inv_sorted_ids[idx]
    
    leaf_level = mlfma_op.octree.levels[mlfma_op.octree.nLevels]
    cube1 = nothing
    cube_idx = nothing
    
    for c in leaf_level.cubes
        if s1 in c.bfInterval
            cube1 = c
        end
        if s_idx in c.bfInterval
            cube_idx = c
        end
    end
    
    if cube1 !== nothing && cube_idx !== nothing
        println("Basis 1 in Cube: ", cube1.ID3D)
        println("Basis $idx in Cube: ", cube_idx.ID3D)
        
        # Check if neighbors
        is_neighbor = false
        # Check if cube_idx is in cube1.neighbors
        # We need the index of cube_idx in the cubes list
        # But cube.neighbors stores indices into the cubes list.
        # We need to find the index of cube_idx.
        
        # Let's find the index of the cubes in the list
        c1_list_idx = findfirst(c -> c === cube1, leaf_level.cubes)
        ci_list_idx = findfirst(c -> c === cube_idx, leaf_level.cubes)
        
        println("Cube 1 List Index: $c1_list_idx")
        println("Cube $idx List Index: $ci_list_idx")
        
        if ci_list_idx in cube1.neighbors
            println("Cubes ARE neighbors in Octree.")
        else
            println("Cubes are NOT neighbors in Octree.")
            println("Cube 1 Neighbors: ", cube1.neighbors)
        end
    end

    if rel_err < 0.05
        println("SUCCESS: MLFMA matches Direct Matrix.")
    else
        println("FAILURE: MLFMA mismatch.")
        
        # Analyze scaling factor
        ratio = norm(y_mlfma) / norm(y_direct)
        println(@sprintf("        Ratio (MLFMA/Direct): %.4f", ratio))
        println(@sprintf("        20*log10(Ratio): %.4f dB", 20*log10(ratio)))
    end
end

verify_mlfma_matrix()
