using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.Solvers
using EMMoMSuite.CoreModule.Sources
using LinearAlgebra
using Printf

# 1. Setup
radius = 0.5
n_theta = 20
n_phi = 40
println("Generating sphere mesh...")
mesh = generate_sphere_mesh(radius, n_theta, n_phi)
basis = RWGBasis(mesh)
n_unknowns = num_basis(basis)
println("Unknowns: $n_unknowns")

freq = 300e6
efie = EFIE(freq)
println("Assembling Z...")
Z = assemble_impedance_matrix(efie, basis)

source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
V = excitation_vector(efie, source, basis)

# 2. Solve without Preconditioner
println("\n--- GMRES (No Preconditioner) ---")
solver_gmres = GMRESSolver(tol=1e-4, maxiter=1000, verbose=true)
t1 = @elapsed begin
    I1 = solve!(solver_gmres, Z, V)
end
println("Time: $t1 s")

# 3. Solve with Diagonal Preconditioner
println("\n--- GMRES (Diagonal Preconditioner) ---")
P = DiagonalPreconditioner(Z)
t2 = @elapsed begin
    I2 = solve!(solver_gmres, Z, V; Pl=P)
end
println("Time: $t2 s")

# 4. Solve with BiCGSTAB (No Preconditioner)
println("\n--- BiCGSTAB (No Preconditioner) ---")
solver_bicg = BiCGSTABSolver(tol=1e-4, maxiter=1000, verbose=true)
t3 = @elapsed begin
    I3 = solve!(solver_bicg, Z, V)
end
println("Time: $t3 s")

# 5. Solve with BiCGSTAB (Diagonal Preconditioner)
println("\n--- BiCGSTAB (Diagonal Preconditioner) ---")
t4 = @elapsed begin
    I4 = solve!(solver_bicg, Z, V; Pl=P)
end
println("Time: $t4 s")

# Check agreement
diff = norm(I2 - I1) / norm(I1)
println("\nDifference (GMRES vs GMRES+Precond): $diff")
