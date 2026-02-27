# Direct comparison: Legacy vs EMSuite full Z matrix on small mesh
using MoM_Basics
using MoM_Kernels
using EMSuite
using LinearAlgebra
using Printf

# Use the plate benchmark mesh (small)
mesh_file = joinpath(@__DIR__, "plate_benchmark.nas")
println("Mesh file: $mesh_file")
println("Exists: $(isfile(mesh_file))")

# ========== EMSuite Z matrix ==========
println("\n===== EMSuite Assembly =====")
mesh_ems = read_nas_mesh(mesh_file, scale=1.0)
freq = 3e8
set_frequency!(freq)
basis_ems = RWGBasis(mesh_ems)
N = num_basis(basis_ems)
println("Unknowns: $N")

efie = EFIE(freq)
Z_ems = assemble_impedance_matrix(efie, basis_ems)

println("Z_ems[1,1] = $(Z_ems[1,1])")
println("Z_ems size: $(size(Z_ems))")

# ========== Legacy Z matrix ==========
println("\n===== Legacy Assembly =====")
# Initialize Legacy parameters
# Need to load mesh with Legacy API
using MoM_Basics: SimulationParams, inputBasicParameters, getMeshData
using MoM_Basics: setCal4Param!

# Set Legacy frequencies
inputBasicParameters(freq)

# Load mesh with Legacy system  
meshData, _ = getMeshData(mesh_file; rscale=1.0)
# Setup triangles
triinfos = meshData.trianglesInfo
ntri = length(triinfos)
println("Legacy triangles: $ntri")

# Get number of basis functions
nbf = meshData.nbf
println("Legacy unknowns: $nbf")

# Call Legacy impedance matrix assembly
Zmat_leg = zeros(ComplexF32, nbf, nbf)
MoM_Kernels.impedancemat4EFIE4PEC!(Zmat_leg, triinfos, MoM_Basics.BasisFunctionType.RWG)

println("Z_leg[1,1] = $(Zmat_leg[1,1])")

# ========== Compare ==========
println("\n===== Comparison =====")
Z_leg64 = Complex{Float64}.(Zmat_leg)  # Convert to Float64

# Diagonal comparison
println("\nDiagonal comparison (first 10 elements):")
for k in 1:min(10, N)
    ratio = Z_ems[k,k] / Z_leg64[k,k]
    @printf("Z[%d,%d]: EMSuite=%.4f%+.4fi  Legacy=%.4f%+.4fi  ratio=%.6f\n",
        k, k, real(Z_ems[k,k]), imag(Z_ems[k,k]),
        real(Z_leg64[k,k]), imag(Z_leg64[k,k]), abs(ratio))
end

# Full matrix comparison
diffs = abs.(Z_ems .- Z_leg64)
ratios = abs.(Z_ems) ./ max.(abs.(Z_leg64), 1e-30)
println("\nMatrix-wide statistics:")
println("Max |diff|: $(maximum(diffs))")
println("Mean |diff|: $(mean(diffs))")
println("Frobenius norm EMSuite: $(norm(Z_ems))")
println("Frobenius norm Legacy:  $(norm(Z_leg64))")
println("Frobenius norm ratio:   $(norm(Z_ems)/norm(Z_leg64))")

# Check diagonal ratio
diag_ratios = abs.(diag(Z_ems)) ./ abs.(diag(Z_leg64))
println("\nDiagonal ratio stats:")
println("  Min: $(minimum(diag_ratios))")
println("  Max: $(maximum(diag_ratios))")
println("  Mean: $(mean(diag_ratios))")
println("  Std: $(std(diag_ratios))")

# Check off-diagonal ratio (sample)
println("\nOff-diagonal sample (first 5×5 block):")
for i in 1:min(5,N), j in 1:min(5,N)
    if i != j && abs(Z_leg64[i,j]) > 1e-10
        ratio = Z_ems[i,j] / Z_leg64[i,j]
        @printf("  Z[%d,%d]: ratio=%.4f+%.4fi (|ratio|=%.6f)\n",
            i, j, real(ratio), imag(ratio), abs(ratio))
    end
end
