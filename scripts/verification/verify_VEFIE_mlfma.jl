using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Solvers
using EMSuite.CoreModule.Sources
using EMSuite.FastAlgorithms.MLFMA
using LinearAlgebra
using Printf
using Random

function verify_vefie_mlfma()
    println("==================================================")
    println("   Verification: VEFIE MLFMA                      ")
    println("==================================================")

    # 1. Parameters
    # Increase frequency to ensure Far Field interactions exist
    # Object size ~0.4m. At 300MHz, lambda=1m -> 0.4 lambda (Small)
    # At 1.5GHz, lambda=0.2m -> 2.0 lambda (Moderate)
    freq = 1.5e9 
    c0 = 299792458.0
    lambda = c0 / freq
    
    # 2. Mesh
    mesh_file = joinpath(@__DIR__, "../../../MoM_AllinOne/meshfiles/Tetra.nas")
    if !isfile(mesh_file)
        println("Error: Mesh file not found: $mesh_file")
        return
    end
    
    println("Loading mesh from $mesh_file...")
    # Scale mesh by 1e-3 (mm to m)
    mesh = read_nas_mesh(mesh_file, scale=1e-3)
    
    # Center mesh
    println("Centering mesh...")
    min_c = minimum(mesh.node, dims=2)
    max_c = maximum(mesh.node, dims=2)
    center = (min_c + max_c) / 2
    mesh.node .-= center
    
    println("Mesh: $(num_vertices(mesh)) vertices, $(num_elements(mesh)) tetrahedra")
    
    # 3. Basis
    println("Setting up SWG Basis...")
    basis = SWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")

    # 4. Permittivity
    eps_r = 2.0 + 0.0im
    permittivities = fill(eps_r, num_elements(mesh))

    # 5. Operator
    println("Setting up VEFIE...")
    # Use new constructor with permittivities
    vefie = VEFIE(freq, permittivities)
    
    # 6. Direct Assembly (Reference)
    println("Assembling Direct Matrix (Reference)...")
    Z_direct = assemble_impedance_matrix(vefie, basis)
    
    # 7. MLFMA Setup
    println("Setting up MLFMA...")
    leaf_box_size = lambda / 4 # Standard MLFMA box size
    mlfma_op = MLFMAOperator(vefie, basis, leaf_box_size)
    
    # 8. Verification
    println("Verifying Matrix-Vector Product...")
    Random.seed!(1234)
    x = rand(ComplexF64, N)
    
    # Direct M-V
    y_direct = Z_direct * x
    
    # MLFMA M-V
    y_mlfma = mlfma_op * x
    
    # Compare
    diff = norm(y_mlfma - y_direct)
    rel_error = diff / norm(y_direct)
    
    println("    Direct Norm: $(norm(y_direct))")
    println("    MLFMA Norm:  $(norm(y_mlfma))")
    println("    Diff Norm:   $diff")
    println("    Rel Error:   $rel_error")
    
    ratio = norm(y_mlfma) / norm(y_direct)
    println("    Ratio:       $ratio")
    
    if rel_error < 0.05
        println("SUCCESS: MLFMA matches Direct Solver.")
    else
        println("FAILURE: MLFMA mismatch.")
    end
end

verify_vefie_mlfma()
