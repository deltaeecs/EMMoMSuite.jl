using MoM_AllinOne
using LinearAlgebra
using StaticArrays
using Printf

# Setup
setPrecision!(Float64)
frequency = 300e6
ieT  = :EFIE
sbfT = :RWG
vbfT = :nothing

inputParameters(;frequency = frequency, ieT = ieT)
updateVSBFTParams!(;sbfT = sbfT, vbfT = vbfT)

# Mesh
filename = joinpath(@__DIR__, "..", "temp_sphere.nas")
meshData, εᵣs   =  getMeshData(filename; meshUnit=:m)
ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData; sbfT = sbfT, vbfT = vbfT)

# Get first triangle
tri = geosInfo[1]
println("Triangle 1:")
println("  Edges: $(tri.edgel)")
println("  Area: $(tri.area)")

# Compute Singular Values
C4divk2 = Params.C4divk²
println("  C4divk2: $C4divk2")

a, b, c = abs.(tri.edgel)
sF1 = MoM_Kernels.singularF1(a, b, c)
F1 = C4divk2 * sF1
println("  sF1: $sF1")
println("  F1: $F1")

sF21 = MoM_Kernels.singularF21(a, b, c, tri.area^2)
println("  sF21: $sF21")

diff = sF21 - F1
println("  Diff (sF21 - F1): $diff")

# Check F22
# For m=1, n=2. k=3.
# a, b, c = edgel[3], edgel[1], edgel[2]
a2, b2, c2 = abs(tri.edgel[3]), abs(tri.edgel[1]), abs(tri.edgel[2])
sF22 = MoM_Kernels.singularF22(a2, b2, c2, tri.area^2)
println("  sF22 (1,2): $sF22")
println("  Diff (sF22 - F1): $(sF22 - F1)")
