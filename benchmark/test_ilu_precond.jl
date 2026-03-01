"""Quick test: ILU vs Block Jacobi preconditioner for MLFMA."""

using EMSuite
using LinearAlgebra
using SparseArrays

# Build Plate EFIE MLFMA (small test case)
mesh_file = joinpath(@__DIR__, "plate_benchmark.nas")
freq = 3e8
mesh = read_nas_mesh(mesh_file, scale=1.0)
set_frequency!(freq)
basis = RWGBasis(mesh)
efie = EFIE(freq)
source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
V = excitation_vector(efie, source, basis)
lambda = EMSuite.Constants.c0 / freq
leaf_size = 0.35 * lambda
Z_mlfma = MLFMAOperator(efie, basis, leaf_size)
V_sorted = V[Z_mlfma.sorted_ids]

println("\n--- Sparse LU Preconditioner ---")
t_lu = @elapsed P_lu = lu(Z_mlfma.Z_near)
println("  Build: $(round(t_lu, digits=3))s")

struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
t_solve_lu = @elapsed I_lu = solve!(solver, Z_mlfma, V_sorted; Pl=LUPreconditioner(P_lu))
println("  Solve: $(round(t_solve_lu, digits=3))s")
println("  Total: $(round(t_lu + t_solve_lu, digits=3))s")

println("\n--- ILU Preconditioner (τ=0.01) ---")
t_ilu = @elapsed P_ilu = ILUPreconditioner(Z_mlfma.Z_near; τ=0.01)
println("  Build: $(round(t_ilu, digits=3))s")

solver2 = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
t_solve_ilu = @elapsed I_ilu = solve!(solver2, Z_mlfma, V_sorted; Pl=P_ilu)
println("  Solve: $(round(t_solve_ilu, digits=3))s")
println("  Total: $(round(t_ilu + t_solve_ilu, digits=3))s")

println("\n--- Block Jacobi Preconditioner ---")
t_bj = @elapsed begin
    intervals = get_leaf_intervals(Z_mlfma)
    P_bj = BlockJacobiPreconditioner(Z_mlfma.Z_near, intervals)
end
println("  Build: $(round(t_bj, digits=3))s")

solver3 = GMRESSolver(restart=50, maxiter=200, tol=1e-3, verbose=true)
t_solve_bj = @elapsed I_bj = solve!(solver3, Z_mlfma, V_sorted; Pl=P_bj)
println("  Solve: $(round(t_solve_bj, digits=3))s")
println("  Total: $(round(t_bj + t_solve_bj, digits=3))s")
