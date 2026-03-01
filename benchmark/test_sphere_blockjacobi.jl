"""Test Block Jacobi on Sphere CFIE MLFMA (Case 4) — well-conditioned problem."""

using EMSuite
using LinearAlgebra
using SparseArrays

const MESH_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles")

mesh_file = joinpath(MESH_DIR, "sphere_600MHz.nas")
freq = 6e8

println("Loading mesh...")
mesh = read_nas_mesh(mesh_file, scale=1.0)
set_frequency!(freq)
basis = RWGBasis(mesh)
println("N = $(num_basis(basis))")

cfie = CFIE(freq, 0.5)
source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
V = excitation_vector(cfie, source, basis)

lambda = EMSuite.Constants.c0 / freq
leaf_size = 0.35 * lambda

println("Building MLFMA...")
Z_mlfma = MLFMAOperator(cfie, basis, leaf_size)
V_sorted = V[Z_mlfma.sorted_ids]

# --- Sparse LU baseline ---
println("\n--- Sparse LU Preconditioner ---")

struct LUPreconditioner2
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner2, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner2, x) = (x .= P.F \ x)

t_lu = @elapsed F_lu = lu(Z_mlfma.Z_near)
P_lu = LUPreconditioner2(F_lu)
println("  Build: $(round(t_lu, digits=3))s")

solver_lu = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
t_solve_lu = @elapsed I_lu = solve!(solver_lu, Z_mlfma, V_sorted; Pl=P_lu)
println("  Solve: $(round(t_solve_lu, digits=3))s")
println("  Total: $(round(t_lu + t_solve_lu, digits=3))s")

# --- Block Jacobi ---
println("\n--- Block Jacobi Preconditioner ---")
t_bj = @elapsed begin
    intervals = get_leaf_intervals(Z_mlfma)
    P_bj = BlockJacobiPreconditioner(Z_mlfma.Z_near, intervals)
end
println("  Build: $(round(t_bj, digits=3))s")

solver_bj = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
t_solve_bj = @elapsed I_bj = solve!(solver_bj, Z_mlfma, V_sorted; Pl=P_bj)
println("  Solve: $(round(t_solve_bj, digits=3))s")
println("  Total: $(round(t_bj + t_solve_bj, digits=3))s")
