using EMSuite
using LinearAlgebra
using StaticArrays
using Test
using Printf
using SparseArrays

function debug_mlfma()
    println("==================================================")
    println("   Debug: MLFMA Near Field vs Direct Solver       ")
    println("==================================================")

    # 1. Parameters
    freq = 300e6
    lambda = 299792458.0 / freq
    radius = 0.5
    
    # 2. Mesh & Basis
    mesh = generate_sphere_mesh(radius, 12, 24)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")

    # 3. Direct Solver
    efie = EFIE(freq)
    println("Assembling Direct Matrix...")
    Z_direct = assemble_impedance_matrix(efie, basis)
    
    # 4. MLFMA (Near Field Only)
    # Set leaf size to be very large so everything is in one cube (or few neighbors)
    large_leaf_size = 10.0 * radius 
    println("Building MLFMA with leaf_size = $large_leaf_size (Should be all near field)")
    
    mlfma_op = MLFMAOperator(efie, basis, large_leaf_size)
    
    # Check Octree stats
    n_levels = mlfma_op.octree.nLevels
    n_cubes = length(mlfma_op.octree.levels[n_levels].cubes)
    println("Levels: $n_levels")
    println("Leaf Cubes: $n_cubes")
    
    # Extract Z_near
    Z_near = mlfma_op.Z_near
    
    # Compare Z_near with Z_direct
    # Note: Z_near is sparse, Z_direct is dense
    diff = norm(Z_near - Z_direct) / norm(Z_direct)
    println("Relative Error (Z_near vs Z_direct): $diff")
    
    if diff < 1e-10
        println("SUCCESS: Near Field assembly matches Direct Solver.")
    else
        println("FAILURE: Near Field assembly mismatch.")
        
        # Analyze mismatch
        rows, cols, vals = findnz(Z_near)
        err_max = 0.0
        for i in 1:N, j in 1:N
            val_near = Z_near[i,j]
            val_direct = Z_direct[i,j]
            if abs(val_near - val_direct) > 1e-6
                # println("Mismatch at ($i, $j): Near=$val_near, Direct=$val_direct")
                err_max = max(err_max, abs(val_near - val_direct))
            end
        end
        println("Max element error: $err_max")
    end
end

debug_mlfma()
