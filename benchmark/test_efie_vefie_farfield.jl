using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.FastAlgorithms.MLFMA
using LinearAlgebra, SparseArrays

println("=== Isolated EFIE/VEFIE Far-Field Test ===")

mesh_file = joinpath(@__DIR__, "../../MoM_Basics/meshfiles/TriTetra.nas")
surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale=0.001)
basis_surf = RWGBasis(surf_mesh)
basis_vol = SWGBasis(vol_mesh)
n_surf = num_basis(basis_surf)
n_vol = num_basis(basis_vol)

freq = 4e9  # 4 GHz: λ=0.075m, mesh=0.1m≈1.33λ, edge≈0.13λ (better sampled)
eps_r = 2.0 * (1 - 0.001im)
perms = fill(eps_r, num_elements(vol_mesh))
lambda = 299792458.0 / freq
leaf_size = 0.25 * lambda

# Test 1: Pure EFIE (surface only)
println("\n--- EFIE (surface only) ---")
efie = EFIE(freq)
Z_efie_direct = assemble_impedance_matrix(efie, basis_surf)
mlfma_efie = MLFMAOperator(efie, [basis_surf], leaf_size)
near_pct = round(100 * nnz(mlfma_efie.Z_near) / n_surf^2, digits=1)
println("Near-field: $near_pct%")

x_s = randn(ComplexF64, n_surf)
y_d = Z_efie_direct * x_s
y_m = mlfma_efie * x_s
println("EFIE rel err: $(norm(y_d - y_m) / norm(y_d))")

# Test 2: Pure VEFIE (volume only)
println("\n--- VEFIE (volume only) ---")
vefie = VEFIE(freq, perms)
Z_vefie_direct = assemble_impedance_matrix(vefie, basis_vol)
mlfma_vefie = MLFMAOperator(vefie, [basis_vol], leaf_size)
near_pct2 = round(100 * nnz(mlfma_vefie.Z_near) / n_vol^2, digits=1)
println("Near-field: $near_pct2%")

x_v = randn(ComplexF64, n_vol)
y_d2 = Z_vefie_direct * x_v
y_m2 = mlfma_vefie * x_v
println("VEFIE rel err: $(norm(y_d2 - y_m2) / norm(y_d2))")
