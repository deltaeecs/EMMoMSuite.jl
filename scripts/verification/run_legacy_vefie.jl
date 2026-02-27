using Pkg
# Activate MoM_AllinOne environment
project_path = joinpath(@__DIR__, "../../../MoM_AllinOne")
Pkg.activate(project_path)

# Add local packages to LOAD_PATH
push!(LOAD_PATH, joinpath(@__DIR__, "../../../MoM_Basics"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../MoM_Kernels"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../MoM_MPI"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../MPIArray4MoMs"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../../MoM_Visualizing"))
push!(LOAD_PATH, project_path)

using MoM_AllinOne
using DataFrames, CSV, Printf
using LinearAlgebra

# Parameters
filename = joinpath(project_path, "meshfiles/Tetra.nas")
meshUnit = :mm 
frequency = 300e6
ieT = :EFIE # Based on examples_VEFIE_direct.jl, it uses :EFIE but with vbfT=:SWG
sbfT = :nothing
vbfT = :SWG
solverT = :direct

# Source: Incident from +z (theta=pi), x-pol (alpha=0)
source = PlaneWave(Float64(π), 0.0, 0.0, 1.0)

# Observation: E-plane (phi=0)
θs_obs = collect(0:1.0:180.0) .* (π/180.0)
ϕs_obs = [0.0]

println("Running Legacy VEFIE...")
println("Mesh: $filename")
println("Freq: $frequency")

# 1. Update Parameters
inputParameters(;frequency = frequency, ieT = ieT)
updateVSBFTParams!(;sbfT = sbfT, vbfT = vbfT)

# 2. Read Mesh
println("Reading Mesh...")
meshData, εᵣs = getMeshData(filename; meshUnit=meshUnit)

# 3. Basis Functions
println("Generating Basis...")
ngeo, nbf, geosInfo, bfsInfo = getBFsFromMeshData(meshData; sbfT = sbfT, vbfT = vbfT)

# 4. Set Permittivity
println("Setting Permittivity...")
setGeosPermittivity!(geosInfo, 2.0 + 0.0im)

# 5. Impedance Matrix
println("Assembling Matrix...")
Zmat = getImpedanceMatrix(geosInfo, nbf)

# 6. Excitation
println("Excitation...")
V = getExcitationVector(geosInfo, size(Zmat, 1), source)

println("Legacy Z norm: ", norm(Zmat))
println("Legacy V norm: ", norm(V))

# 7. Solve
println("Solving...")
ICoeff, ch = solve(Zmat, V; solverT = solverT)

# 8. RCS
println("Calculating RCS...")
RCSθsϕs, RCSθsϕsdB, RCS, RCSdB = radarCrossSection(θs_obs, ϕs_obs, ICoeff, geosInfo)

# 9. Save Results
output_file = joinpath(@__DIR__, "legacy_vefie_rcs.txt")
open(output_file, "w") do io
    println(io, "Theta(deg) RCS(dBsm)")
    for i in 1:length(θs_obs)
        deg = θs_obs[i] * 180/π
        val = RCSdB[i, 1] # phi index 1
        @printf(io, "%10.4f %15.4f\n", deg, val)
    end
end
println("Saved Legacy RCS to $output_file")
