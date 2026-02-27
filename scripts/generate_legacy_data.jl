using Pkg
# Activate MoM_Kernels environment
Pkg.activate(joinpath(@__DIR__, "../../MoM_Kernels"))
using MoM_Kernels
using MoM_Basics
using Serialization
using LinearAlgebra

# Parameters
# Correct path to mesh file relative to EMSuite/scripts
filename = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles/sphere_600MHz.nas")
frequency = 300e6
meshUnit = :m

println("Loading mesh from $filename...")

# --- EFIE ---
println("Generating Legacy EFIE Data...")
inputParameters(;frequency = frequency, ieT = :EFIE)
updateVSBFTParams!(;sbfT = :RWG, vbfT = :nothing)
meshData, εᵣs = getMeshData(filename; meshUnit=meshUnit)
ngeo, nbf, geosInfo, bfsInfo = getBFsFromMeshData(meshData; sbfT = :RWG, vbfT = :nothing)

# Force precision to Float64 for comparison
setPrecision!(Float64)

Z_EFIE = getImpedanceMatrix(geosInfo, nbf)
# PlaneWave(theta, phi, pol_theta, pol_phi)
source = PlaneWave(π/2, 0.0, 0.0, 1.0) 
V_EFIE = getExcitationVector(geosInfo, nbf, source)

# --- MFIE ---
println("Generating Legacy MFIE Data...")
inputParameters(;frequency = frequency, ieT = :MFIE)
# Re-calculate Z for MFIE
Z_MFIE = getImpedanceMatrix(geosInfo, nbf)
V_MFIE = getExcitationVector(geosInfo, nbf, source)

# Save using Serialization
output_file = joinpath(@__DIR__, "legacy_data.jls")
println("Saving matrices to $output_file...")
serialize(output_file, Dict(
    "Z_EFIE" => Z_EFIE,
    "Z_MFIE" => Z_MFIE,
    "V_EFIE" => V_EFIE,
    "V_MFIE" => V_MFIE
))
println("Done.")
