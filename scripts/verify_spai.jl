using EMSuite.CoreModule
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.FastAlgorithms.MLFMA
using EMSuite.Solvers
using LinearAlgebra
using SparseArrays

function verify_spai()
    println("=== Verifying SPAI Preconditioner ===")
    
    # 1. Setup
    radius = 1.0
    freq = 150e6 # ka approx pi
    k = 2π * freq / 299792458.0
    println("Frequency: $freq Hz")
    println("Wavenumber k: $k")
    println("ka: $(k*radius)")
    
    # Mesh
    n_theta = 16
    n_phi = 32
    println("Generating mesh ($n_theta x $n_phi)...")
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    println("Number of elements: ", num_elements(mesh))
    
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Number of unknowns: ", N)
    
    operator = EFIE(freq)
    
    # 2. MLFMA Operator
    println("Setting up MLFMA...")
    # Box size 1.0 m
    mlfma = MLFMAOperator(operator, basis, 1.0)
    
    # 3. RHS (Plane Wave)
    println("Computing RHS...")
    pw = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    b = excitation_vector(operator, pw, basis)
    
    # 4. Compare Preconditioners
    println("\n--- Computing Preconditioners ---")
    
    println("1. Diagonal...")
    t_diag = @elapsed P_diag = DiagonalPreconditioner(mlfma.Z_near)
    println("   Time: $t_diag s")
    
    println("2. ILU (tau=0.01)...")
    t_ilu = @elapsed P_ilu = ILUPreconditioner(mlfma.Z_near, τ=0.01)
    println("   Time: $t_ilu s")
    
    println("3. SPAI (Static Pattern)...")
    t_spai = @elapsed P_spai = SPAIPreconditioner(mlfma.Z_near)
    println("   Time: $t_spai s")
    
    # 5. Solve with GMRES
    solver = GMRESSolver(restart=50, maxiter=100, tol=1e-4, verbose=false)
    
    println("\n--- Solving with GMRES ---")
    
    println("1. Diagonal Preconditioner:")
    t_solve_diag = @elapsed x_diag = solve!(solver, mlfma, b; Pl=P_diag)
    res_diag = norm(b - mlfma * x_diag) / norm(b)
    println("   Time: $t_solve_diag s, Residual: $res_diag")
    
    println("2. ILU Preconditioner:")
    t_solve_ilu = @elapsed x_ilu = solve!(solver, mlfma, b; Pl=P_ilu)
    res_ilu = norm(b - mlfma * x_ilu) / norm(b)
    println("   Time: $t_solve_ilu s, Residual: $res_ilu")
    
    println("3. SPAI Preconditioner:")
    t_solve_spai = @elapsed x_spai = solve!(solver, mlfma, b; Pl=P_spai)
    res_spai = norm(b - mlfma * x_spai) / norm(b)
    println("   Time: $t_solve_spai s, Residual: $res_spai")
    
    # 6. Solve with BiCGSTAB
    solver_bicg = BiCGSTABSolver(maxiter=100, tol=1e-4, verbose=false)
    
    println("\n--- Solving with BiCGSTAB ---")
    
    println("1. ILU Preconditioner:")
    t_solve_ilu_bicg = @elapsed x_ilu_bicg = solve!(solver_bicg, mlfma, b; Pl=P_ilu)
    res_ilu_bicg = norm(b - mlfma * x_ilu_bicg) / norm(b)
    println("   Time: $t_solve_ilu_bicg s, Residual: $res_ilu_bicg")
    
    println("2. SPAI Preconditioner:")
    t_solve_spai_bicg = @elapsed x_spai_bicg = solve!(solver_bicg, mlfma, b; Pl=P_spai)
    res_spai_bicg = norm(b - mlfma * x_spai_bicg) / norm(b)
    println("   Time: $t_solve_spai_bicg s, Residual: $res_spai_bicg")
    
end

verify_spai()
