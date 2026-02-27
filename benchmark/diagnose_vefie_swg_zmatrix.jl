# Diagnose VEFIE SWG Z matrix differences between EMSuite and Legacy
# Compare element by element to find the factor difference
using Pkg

# First run Legacy to save Z matrix
println("=== Step 1: Legacy Z matrix ===")
Pkg.activate(joinpath(@__DIR__, "../../LegacyBenchmark"))

using MoM_Basics
using MoM_Kernels
using LinearAlgebra
using Printf
using JLD2

setPrecision!(Float64)
SimulationParams.SHOWIMAGE = false

filename = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles/Tetra.nas")
frequency = 300e6
ieT = :EFIE
sbfT = :nothing
vbfT = :SWG

inputParameters(; frequency=frequency, ieT=ieT)
updateVSBFTParams!(; sbfT=sbfT, vbfT=vbfT)

meshData, εᵣs = getMeshData(filename; meshUnit=:mm)
ngeo, nbf, geosInfo, bfsInfo = getBFsFromMeshData(meshData; sbfT=sbfT, vbfT=vbfT)
setGeosPermittivity!(geosInfo, 2.0 + 0.0im)

println("Legacy: ngeo=$ngeo, nbf=$nbf")
Z_legacy = getImpedanceMatrix(geosInfo, nbf)
println("Legacy Z[1,1] = $(Z_legacy[1,1])")
println("Legacy ||Z|| = $(norm(Z_legacy))")

# Excitation
source_legacy = PlaneWave(Float64(π), 0.0, 0.0, 1.0)
V_legacy = getExcitationVector(geosInfo, nbf, source_legacy)
println("Legacy ||V|| = $(norm(V_legacy))")

# Save for comparison
outdir = joinpath(@__DIR__, "../test_results")
mkpath(outdir)
@save joinpath(outdir, "legacy_vefie_swg_Zmatrix.jld2") Z_legacy V_legacy

println("\nLegacy data saved.")

# Key matrix statistics
println("\n=== Legacy Z matrix statistics ===")
println("  diag mean:  $(mean(abs.(diag(Z_legacy))))")
println("  offdiag mean: $(mean(abs.(Z_legacy - Diagonal(diag(Z_legacy)))))") 
println("  symmetric: $(isapprox(Z_legacy, transpose(Z_legacy), rtol=1e-6))")
println("  Frobenius: $(norm(Z_legacy))")

# Print first 5x5 block
println("\n  Z_legacy[1:5, 1:5] =")
for i in 1:min(5, nbf)
    for j in 1:min(5, nbf)
        @printf("    (%+.6e, %+.6e)  ", real(Z_legacy[i,j]), imag(Z_legacy[i,j]))
    end
    println()
end

using Statistics: mean
