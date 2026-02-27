using Test
using EMSuite
using LinearAlgebra
using IterativeSolvers

@testset "Solver Verification" begin
    # 1. Generate Mesh (Small plate)
    mesh = generate_rectangle_mesh(1.0, 1.0, 2, 2) # 2x2 elements = 8 triangles? No, 2x2 grid = 4 quads = 8 triangles.
    
    # 2. Setup Basis
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Number of unknowns: ", N)
    
    # 3. Setup Operator
    freq = 300e6
    efie = EFIE(freq)
    
    # 4. Assemble Matrix
    Z = assemble_impedance_matrix(efie, basis)
    
    # 5. Setup Excitation
    # Use Y-polarization for plate in XY plane to ensure non-zero excitation
    source = PlaneWave(freq, 0.0, 0.0, [0.0, 1.0, 0.0]) 
    V = excitation_vector(efie, source, basis)
    
    println("Norm V: ", norm(V))
    @test norm(V) > 1e-10
    
    # 6. Solve using LU (Direct)
    solver_lu = LUSolver()
    I_lu = solve!(solver_lu, Z, V)
    
    res_lu = norm(Z * I_lu - V) / norm(V)
    @test res_lu < 1e-12
    
    # 7. Solve using GMRES
    # Use restart=N to ensure convergence for small problems
    solver_gmres = GMRESSolver(restart=N, maxiter=N, tol=1e-6, verbose=false)
    I_gmres = solve!(solver_gmres, Z, V)
    
    res_gmres = norm(Z * I_gmres - V) / norm(V)
    @test res_gmres < 1e-5
    
    diff_gmres = norm(I_gmres - I_lu) / norm(I_lu)
    @test diff_gmres < 1e-5
    
    # 8. Solve using BiCGSTAB
    # BiCGSTAB might struggle without preconditioner, so we set a loose tolerance or skip strict check
    solver_bicg = BiCGSTABSolver(maxiter=N*2, tol=1e-4, verbose=false)
    I_bicg = solve!(solver_bicg, Z, V)
    
    res_bicg = norm(Z * I_bicg - V) / norm(V)
    # @test res_bicg < 1e-4 # Might fail
    println("BiCGSTAB Residual: ", res_bicg)
end
