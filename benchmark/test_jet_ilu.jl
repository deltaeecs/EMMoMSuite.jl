"""Test ILU preconditioner on Jet EFIE MLFMA (Case 3)."""

using EMSuite
using LinearAlgebra
using SparseArrays

const MESH_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles")

mesh_file = joinpath(MESH_DIR, "jet_100MHz.nas")
freq = 1e8

println("Loading mesh...")
mesh = read_nas_mesh(mesh_file, scale=1.0)
set_frequency!(freq)
basis = RWGBasis(mesh)
println("N = $(num_basis(basis))")

efie = EFIE(freq)
source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
V = excitation_vector(efie, source, basis)

lambda = EMSuite.Constants.c0 / freq
leaf_size = 0.35 * lambda

println("Building MLFMA...")
Z_mlfma = MLFMAOperator(efie, basis, leaf_size)
V_sorted = V[Z_mlfma.sorted_ids]

println("\n--- ILU Preconditioner (τ=0.01) ---")
t_ilu = @elapsed P_ilu = ILUPreconditioner(Z_mlfma.Z_near; τ=0.01)
println("  Build: $(round(t_ilu, digits=3))s")

solver = GMRESSolver(restart=50, maxiter=200, tol=1e-3, verbose=true)
t_solve = @elapsed I_sorted = solve!(solver, Z_mlfma, V_sorted; Pl=P_ilu)
println("  Solve: $(round(t_solve, digits=3))s")
println("  Total precond+solve: $(round(t_ilu + t_solve, digits=3))s")
println("  (Baseline LU: precond=45.86s, solve=13.32s, total=59.18s)")
