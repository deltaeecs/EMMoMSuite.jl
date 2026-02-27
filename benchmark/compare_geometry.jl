# Compare triangle geometry between Legacy and EMSuite
using MoM_Basics, MoM_Kernels
using EMSuite
using EMSuite.BasisFunctions: get_triangle_info
using LinearAlgebra, Printf

mesh_file = joinpath(@__DIR__, "plate_benchmark.nas")
freq = 3e8

# ========== Legacy ==========
MoM_Basics.setPrecision!(Float64)
MoM_Basics.SimulationParams.SHOWIMAGE = false
MoM_Kernels.inputParameters(frequency=freq, ieT=:EFIE)
meshData, _ = MoM_Basics.getMeshData(mesh_file; meshUnit=:m)
ngeo_l, nbf_l, geosInfo_l, _ = MoM_Basics.getBFsFromMeshData(meshData; sbfT=:RWG)

# ========== EMSuite ==========
mesh_ems = read_nas_mesh(mesh_file, scale=1.0)
set_frequency!(freq)
basis_ems = RWGBasis(mesh_ems)

println("Legacy triangles: $ngeo_l, EMSuite triangles: $(EMSuite.CoreModule.num_elements(mesh_ems))")
println("Legacy unknowns: $nbf_l, EMSuite unknowns: $(num_basis(basis_ems))")

# Compare first few triangles
println("\n===== Triangle 1 Comparison =====")
tri_l = geosInfo_l[1]
tri_e = get_triangle_info(mesh_ems, basis_ems, 1)

println("\nLegacy Triangle 1:")
println("  Vertices: $(tri_l.vertices)")
println("  Area: $(tri_l.area)")
println("  edgel: $(tri_l.edgel)")
println("  inBfsID: $(tri_l.inBfsID)")

println("\nEMSuite Triangle 1:")
println("  Vertices: $(tri_e.vertices)")
println("  Area: $(tri_e.area)")
println("  edgel: $(tri_e.edgel)")
println("  bfsSign: $(tri_e.bfsSign)")
println("  inBfsID: $(tri_e.inBfsID)")

# Compare Legacy edgel (signed) vs EMSuite edgel * bfsSign
println("\nLegacy signed edgel: $(tri_l.edgel)")
println("EMSuite edgel*sign:  $(tri_e.edgel .* tri_e.bfsSign)")

# Compare MORE triangles
println("\n===== Triangle area comparison (first 10) =====")
for t in 1:min(10, ngeo_l)
    area_l = geosInfo_l[t].area
    tri_e_t = get_triangle_info(mesh_ems, basis_ems, t)
    area_e = tri_e_t.area
    @printf("Tri %3d: Legacy area=%.8f, EMSuite area=%.8f, ratio=%.6f\n",
        t, area_l, area_e, area_e/area_l)
end

# Compare basis function mapping
println("\n===== Basis function mapping comparison =====")
println("Legacy BF1: inBfsID on each tri")
for t in 1:ngeo_l
    for k in 1:3
        if geosInfo_l[t].inBfsID[k] == 1
            println("  Found on Legacy tri $t, local edge $k, edgel=$(geosInfo_l[t].edgel[k])")
        end
    end
end

println("EMSuite BF1: inBfsID on each tri")
ntri_e = EMSuite.CoreModule.num_elements(mesh_ems)
for t in 1:ntri_e
    tri = get_triangle_info(mesh_ems, basis_ems, t)
    for k in 1:3
        if tri.inBfsID[k] == 1
            println("  Found on EMSuite tri $t, local edge $k, edgel=$(tri.edgel[k]), sign=$(tri.bfsSign[k])")
        end
    end
end

# Check Gauss quadrature comparison
println("\n===== Gauss Quadrature =====")
efie = EFIE(freq)
println("EMSuite far GQ points: $(length(efie.gq_far.weight))")
println("EMSuite far GQ weights sum: $(sum(efie.gq_far.weight))")
println("EMSuite near GQ points: $(length(efie.gq_near.weight))")
println("EMSuite near GQ weights sum: $(sum(efie.gq_near.weight))")

# Legacy GQ info
println("\nLegacy GQPNTri (far): $(MoM_Basics.GQPNTri)")
println("Legacy GQPNTriSglr (near): $(MoM_Basics.GQPNTriSglr)")
println("Legacy far weights sum: $(sum(MoM_Basics.TriGQInfo.weight))")
println("Legacy near weights sum: $(sum(MoM_Basics.TriGQInfoSglr.weight))")
