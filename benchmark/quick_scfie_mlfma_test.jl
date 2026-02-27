using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.FastAlgorithms.MLFMA
using EMSuite.Solvers
using LinearAlgebra
using SparseArrays

println("=== Quick SCFIE MLFMA Test (TriTetra mesh) ===")

# Use simple TriTetra mesh (coordinates in mm, need scale=0.001)
mesh_file = joinpath(@__DIR__, "../../MoM_Basics/meshfiles/TriTetra.nas")
surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale=0.001)
println("Surface: $(num_elements(surf_mesh)) tri, Volume: $(num_elements(vol_mesh)) tet")
println("Mesh extent: x=$(extrema(vol_mesh.node[1,:])), y=$(extrema(vol_mesh.node[2,:])), z=$(extrema(vol_mesh.node[3,:]))")

basis_surf = RWGBasis(surf_mesh)
basis_vol = SWGBasis(vol_mesh)
n_surf = num_basis(basis_surf)
n_vol = num_basis(basis_vol)
n_total = n_surf + n_vol
println("Unknowns: surf=$n_surf, vol=$n_vol, total=$n_total")

# Match Legacy: freq=2GHz, eps_r=2(1-0.001j)
freq = 2e9
eps_r = 2.0 * (1 - 0.001im)
perms = fill(eps_r, num_elements(vol_mesh))

# Direct solver
println("\n--- Direct Solver ---")
scfie = SCFIE(freq, perms; alpha=0.5)
Z_direct = assemble_impedance_matrix(scfie, basis_surf, basis_vol)
println("Z_direct norm: $(norm(Z_direct))")

# Check block norms
Z_SS = Z_direct[1:n_surf, 1:n_surf]
Z_SV = Z_direct[1:n_surf, n_surf+1:end]
Z_VS = Z_direct[n_surf+1:end, 1:n_surf]
Z_VV = Z_direct[n_surf+1:end, n_surf+1:end]
println("Block norms: SS=$(norm(Z_SS)), SV=$(norm(Z_SV)), VS=$(norm(Z_VS)), VV=$(norm(Z_VV))")

# MLFMA with 0.25λ leaf size
println("\n--- MLFMA ---")
lambda = 299792458.0 / freq
leaf_size = 0.25 * lambda
println("lambda=$lambda m, leaf_size=$leaf_size m")

mlfma_op = MLFMAOperator(scfie, [basis_surf, basis_vol], leaf_size)
println("Octree Levels: $(mlfma_op.octree.nLevels)")
near_nnz = nnz(mlfma_op.Z_near)
println("Z_near nnz: $near_nnz / $(n_total^2) ($(round(100*near_nnz/n_total^2, digits=1))%)")

# MatVec comparison
x = randn(ComplexF64, n_total)
y_direct = Z_direct * x
y_mlfma = mlfma_op * x

rel_err = norm(y_direct - y_mlfma) / norm(y_direct)
println("\nRelative Error (MLFMA vs Direct): $rel_err")

y_d_s = y_direct[1:n_surf]; y_m_s = y_mlfma[1:n_surf]
y_d_v = y_direct[n_surf+1:end]; y_m_v = y_mlfma[n_surf+1:end]
norm(y_d_s) > 0 && println("  Surface block rel err: $(norm(y_d_s - y_m_s) / norm(y_d_s))")
norm(y_d_v) > 0 && println("  Volume block rel err: $(norm(y_d_v - y_m_v) / norm(y_d_v))")

near_nnz == n_total^2 && println("\nNote: All near-field (mesh < 1λ). Tests assembly only.")
println(rel_err < 0.05 ? "VERIFICATION PASSED (near-field)" : "VERIFICATION FAILED")

# === Far-Field Test (higher freq) ===
# Use higher frequency to make mesh electrically large (>1λ) → ensures far-field interactions
println("\n\n=== Far-Field Test (4 GHz, mesh≈1.33λ) ===")
freq2 = 4e9  # λ=0.075m, mesh=0.1m≈1.33λ, edge≈0.13λ
eps_r2 = 2.0 * (1 - 0.001im)
perms2 = fill(eps_r2, num_elements(vol_mesh))

scfie2 = SCFIE(freq2, perms2; alpha=0.5)
Z_direct2 = assemble_impedance_matrix(scfie2, basis_surf, basis_vol)
println("Z_direct2 norm: $(norm(Z_direct2))")

lambda2 = 299792458.0 / freq2
leaf_size2 = 0.25 * lambda2
println("lambda=$lambda2 m, leaf_size=$leaf_size2 m")

mlfma_op2 = MLFMAOperator(scfie2, [basis_surf, basis_vol], leaf_size2)
near_nnz2 = nnz(mlfma_op2.Z_near)
println("Octree: $(mlfma_op2.octree.nLevels) levels, Z_near: $near_nnz2 / $(n_total^2) ($(round(100*near_nnz2/n_total^2, digits=1))%)")

x2 = randn(ComplexF64, n_total)
y_d2 = Z_direct2 * x2
y_m2 = mlfma_op2 * x2

rel_err2 = norm(y_d2 - y_m2) / norm(y_d2)
println("\nRelative Error (far-field): $rel_err2")

y_d2_s = y_d2[1:n_surf]; y_m2_s = y_m2[1:n_surf]
y_d2_v = y_d2[n_surf+1:end]; y_m2_v = y_m2[n_surf+1:end]
norm(y_d2_s) > 0 && println("  Surface block rel err: $(norm(y_d2_s - y_m2_s) / norm(y_d2_s))")
norm(y_d2_v) > 0 && println("  Volume block rel err: $(norm(y_d2_v - y_m2_v) / norm(y_d2_v))")

println(rel_err2 < 0.1 ? "FAR-FIELD VERIFICATION PASSED" : "FAR-FIELD VERIFICATION FAILED")
