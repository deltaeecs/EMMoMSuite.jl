# Direct comparison: Legacy vs EMSuite full Z matrix on plate mesh
using MoM_Basics, MoM_Kernels
using EMSuite
using LinearAlgebra, Statistics
using Printf

mesh_file = joinpath(@__DIR__, "plate_benchmark.nas")
println("Mesh file: $mesh_file")
println("Exists: $(isfile(mesh_file))")

freq = 3e8

# ========== Legacy Z matrix ==========
println("\n===== Legacy Assembly =====")
MoM_Basics.setPrecision!(Float64)
MoM_Basics.SimulationParams.SHOWIMAGE = false
MoM_Kernels.inputParameters(frequency=freq, ieT=:EFIE)

meshData, _ = MoM_Basics.getMeshData(mesh_file; meshUnit=:m)
ngeo_l, nbf_l, geosInfo_l, _ = MoM_Basics.getBFsFromMeshData(meshData; sbfT=:RWG)
println("Legacy: $ngeo_l triangles, $nbf_l unknowns")

Z_leg = MoM_Kernels.impedancemat4EFIE4PEC(geosInfo_l, nbf_l, MoM_Basics.RWG)
println("Z_leg size: $(size(Z_leg)),  Z_leg[1,1] = $(Z_leg[1,1])")

# ========== EMSuite Z matrix ==========
println("\n===== EMSuite Assembly =====")
mesh_ems = read_nas_mesh(mesh_file, scale=1.0)
set_frequency!(freq)
basis_ems = RWGBasis(mesh_ems)
N_ems = num_basis(basis_ems)
println("EMSuite: $(EMSuite.CoreModule.num_elements(mesh_ems)) triangles, $N_ems unknowns")

efie = EFIE(freq)
Z_ems = assemble_impedance_matrix(efie, basis_ems)
println("Z_ems size: $(size(Z_ems)),  Z_ems[1,1] = $(Z_ems[1,1])")

# ========== Check basis function ordering ==========
# The ordering may differ. Let's compare magnitudes of diagonals.
println("\n===== Diagonal Comparison =====")
N = min(size(Z_ems, 1), size(Z_leg, 1))
println("Compared dimension: $N")

if N > 0
    for k in 1:min(20, N)
        re = real(Z_ems[k,k])
        ie = imag(Z_ems[k,k])
        rl = real(Z_leg[k,k])
        il = imag(Z_leg[k,k])
        ratio = Z_ems[k,k] / Z_leg[k,k]
        @printf("Z[%3d,%3d] | EMSuite=%+10.4f%+10.4fi | Legacy=%+10.4f%+10.4fi | ratio=%.4f%+.4fi (|r|=%.6f)\n",
            k, k, re, ie, rl, il, real(ratio), imag(ratio), abs(ratio))
    end

    # Diagonal stats
    diag_e = diag(Z_ems[1:N, 1:N])
    diag_l = diag(Z_leg[1:N, 1:N])
    diag_ratio = abs.(diag_e) ./ abs.(diag_l)
    
    println("\nDiagonal |ratio| stats:")
    @printf("  Min:  %.6f\n", minimum(diag_ratio))
    @printf("  Max:  %.6f\n", maximum(diag_ratio))
    @printf("  Mean: %.6f\n", mean(diag_ratio))
    @printf("  Std:  %.6f\n", std(diag_ratio))

    # Frobenius norm comparison
    fn_e = norm(Z_ems)
    fn_l = norm(Z_leg)
    @printf("\nFrobenius norm: EMSuite=%.4f, Legacy=%.4f, ratio=%.6f\n", fn_e, fn_l, fn_e/fn_l)

    # Element-wise max diff
    diff_mat = abs.(Z_ems[1:N,1:N] .- Z_leg[1:N,1:N])
    @printf("Max |diff|: %.6e\n", maximum(diff_mat))
    @printf("Mean |diff|: %.6e\n", mean(diff_mat))
    
    # Relative error
    scale = max.(abs.(Z_leg[1:N,1:N]), 1e-30)
    rel_err = diff_mat ./ scale
    @printf("Max rel error: %.6e\n", maximum(rel_err))
    @printf("Mean rel error: %.6e\n", mean(rel_err))
    
    # Sample off-diagonal
    println("\nOff-diagonal comparison (first 10 elements with |Z| > 0.01):")
    count = 0
    for i in 1:N, j in 1:N
        if i != j && abs(Z_leg[i,j]) > 0.01 && count < 10
            ratio = Z_ems[i,j] / Z_leg[i,j]
            @printf("  Z[%d,%d]: |ratio|=%.6f  phase_diff=%.2f°\n",
                i, j, abs(ratio), angle(ratio)*180/π)
            count += 1
        end
    end
end
