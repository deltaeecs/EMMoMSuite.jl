using EMSuite
using LinearAlgebra
using StaticArrays
using Test
using Printf

function debug_scaling()
    println("==================================================")
    println("   Debug: MLFMA Scaling Factor                    ")
    println("==================================================")

    freq = 300e6
    k = 2π * freq / 299792458.0
    eta = 376.73
    lambda = 299792458.0 / freq
    
    # 1. Define two pairs of triangles far apart (z-axis)
    # Group 1 (Source) at z=0
    # Square composed of 2 triangles
    # v1=(0,0), v2=(0.1,0), v3=(0,0.1), v4=(0.1,0.1)
    g1_verts = [0.0 0.1 0.0 0.1; 0.0 0.0 0.1 0.1; 0.0 0.0 0.0 0.0]
    g1_elems = [1 2; 2 4; 3 3] # Tri 1: 1-2-3, Tri 2: 2-4-3
    
    # Group 2 (Test) just outside neighbor list
    # leaf_size = 0.25.
    # Neighbors are within +/- 1 index.
    # We want index +2.
    # dist = 2.0 * leaf_size = 0.5.
    
    leaf_size = 0.1
    dist = 6.25 * leaf_size # 0.625m.
    dist_vec = [dist, 0.0, 0.0]
    g2_verts = g1_verts .+ dist_vec
    g2_elems = g1_elems .+ 4
    
    vertices = hcat(g1_verts, g2_verts)
    elements = hcat(g1_elems, g2_elems)
    
    num_tris = size(elements, 2)
    mesh = TriangleMesh(num_tris, vertices, elements)
    basis = RWGBasis(mesh)
    
    # Find the internal edge in Group 1 and Group 2
    # Edge 2-3 is common in Group 1.
    # Edge 6-7 is common in Group 2.
    
    # We can just search for basis functions
    bfs = basis.functions
    println("Number of BFs: $(length(bfs))")
    println("Basis Map Type: $(typeof(basis.basis_map))")
    println("Basis Map Size: $(size(basis.basis_map))")
    
    if length(bfs) < 2
        error("Not enough basis functions!")
    end
    
    bf_src_idx = 1
    bf_test_idx = 2
    
    # Ensure they are far apart
    c1 = bfs[bf_src_idx].center
    c2 = bfs[bf_test_idx].center
    println("Distance: $(norm(c1 - c2))")
    
    # 2. Compute Direct Z
    efie = EFIE(freq)
    Z_direct = assemble_impedance_matrix(efie, basis)
    val_direct = Z_direct[bf_test_idx, bf_src_idx]
    println("Direct Z: $val_direct")
    
    # 3. Compute MLFMA Z
    # We need to manually run the MLFMA pipeline for these two basis functions.
    # Or use MLFMAOperator with a small leaf size so they are in different cubes.
    
    leaf_size = 0.25 # 0.25 meter. Distance is 10m.
    mlfma_op = MLFMAOperator(efie, basis, leaf_size)
    
    # Create input vector x (1 at src, 0 elsewhere)
    x = zeros(ComplexF64, num_basis(basis))
    x[bf_src_idx] = 1.0
    
    # Compute y = A * x
    y = mlfma_op * x
    val_mlfma = y[bf_test_idx]
    println("MLFMA Z: $val_mlfma")
    
    # 4. Compare
    ratio = abs(val_mlfma) / abs(val_direct)
    println("Ratio (MLFMA / Direct): $ratio")
    println("Expected: 1.0")
    
    # Check if it's a factor of 4 or 16 or pi
    println("Ratio / 4: $(ratio / 4)")
    println("Ratio * 4: $(ratio * 4)")
    println("Ratio / 16: $(ratio / 16)")
    println("Ratio * 16: $(ratio * 16)")
    println("Ratio / pi: $(ratio / pi)")
    println("Ratio * pi: $(ratio * pi)")
end

debug_scaling()
