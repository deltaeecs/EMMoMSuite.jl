using MoM_AllinOne
using LinearAlgebra
using JLD2
using StaticArrays

# 1. Setup
setPrecision!(Float64)
SimulationParams.SHOWIMAGE = false

# Mesh
# Ensure we use the same mesh as the EMSuite tests
filename = joinpath(@__DIR__, "..", "temp_sphere.nas")
if !isfile(filename)
    error("Mesh file not found: $filename. Please run setup.jl or gen_mesh_file.jl first.")
end
meshUnit = :m

# Physics
frequency = 300e6
ieT  = :EFIE # Start with EFIE for baseline
sbfT = :RWG
vbfT = :nothing

# 2. Initialize Legacy System
inputParameters(;frequency = frequency, ieT = ieT)
updateVSBFTParams!(;sbfT = sbfT, vbfT = vbfT)

# 3. Load Mesh & Generate Basis
println("Loading Mesh and Generating Basis Functions...")
meshData, εᵣs   =  getMeshData(filename; meshUnit=meshUnit)
ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData; sbfT = sbfT, vbfT = vbfT)

println("Basis Functions Generated: $nbf")

# 4. Compute Z Matrix
println("Computing Impedance Matrix...")
Z = getImpedanceMatrix(geosInfo, nbf)

# 5. Extract Data for Verification

# 5.1 Basis Function Data
# We want to export a simplified struct or dictionary for each BF
# Note: bfsInfo is a Vector{RWG}
bf_data = []
for (i, bf) in enumerate(bfsInfo)
    # bf is of type RWG
    # We convert StaticArrays to standard Arrays for easier JLD2 compatibility/inspection
    push!(bf_data, Dict(
        "id" => bf.bfID,
        "isbd" => bf.isbd,
        "edge_length" => bf.edgel,
        "center" => Vector(bf.center),
        "triangles" => Vector(bf.inGeo),      # Triangle IDs
        "local_edge_ids" => Vector(bf.inGeoID) # Local Edge IDs (1,2,3)
    ))
end

# 5.2 Triangle Data
tri_data = []
for (i, tri) in enumerate(geosInfo)
    # tri is of type TriangleInfo
    push!(tri_data, Dict(
        "id" => tri.triID,
        "area" => tri.area,
        "center" => Vector(tri.center),
        "vertices" => Matrix(tri.vertices),
        "edge_lengths" => Vector(tri.edgel), # Vector of 3 edge lengths
        "bf_ids" => Vector(tri.inBfsID)      # Vector of 3 BF IDs
    ))
end

# 6. Save to JLD2
output_file = joinpath(@__DIR__, "..", "legacy_golden_data.jld2")
println("Saving Golden Data to $output_file...")
@save output_file Z bf_data tri_data frequency

println("Done.")
