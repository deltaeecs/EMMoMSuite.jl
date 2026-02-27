using EMSuite.CoreModule
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.FastAlgorithms.MLFMA
using EMSuite.Solvers
using LinearAlgebra
using SparseArrays

function verify_solvers()
    println("=== Verifying Solvers ===")
    
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
    # Increase box size to include more near terms in preconditioner
    # Box size 1.0 m (Radius)
    mlfma = MLFMAOperator(operator, basis, 1.0)
    
    # 3. RHS (Plane Wave)
    println("Computing RHS...")
    # Incident from +z (theta=0), polarized in x
    pw = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    b = excitation_vector(operator, pw, basis)
    
    # 4. Solvers
    println("\n--- Testing GMRES (ILU Preconditioner) ---")
    println("Computing ILU Preconditioner...")
    P_ilu = ILUPreconditioner(mlfma.Z_near, τ=0.01)
    
    solver_gmres = GMRESSolver(restart=50, maxiter=100, tol=1e-4, verbose=true)
    
    println("Solving with GMRES...")
    x_gmres = solve!(solver_gmres, mlfma, b; Pl=P_ilu)
    
    res_gmres = norm(b - mlfma * x_gmres) / norm(b)
    println("GMRES Relative Residual: $res_gmres")
    
    # 5. BiCGSTAB
    println("\n--- Testing BiCGSTAB (ILU Preconditioner) ---")
    
    solver_bicg = BiCGSTABSolver(maxiter=100, tol=1e-4, verbose=true)
    
    println("Solving with BiCGSTAB...")
    x_bicg = solve!(solver_bicg, mlfma, b; Pl=P_ilu)
    
    res_bicg = norm(b - mlfma * x_bicg) / norm(b)
    println("BiCGSTAB Relative Residual: $res_bicg")
    
    # Compare solutions
    diff_sol = norm(x_gmres - x_bicg) / norm(x_gmres)
    println("Difference between GMRES and BiCGSTAB solutions: $diff_sol")
    
end

verify_solvers()
