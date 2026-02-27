# Compare faceSingularityIgIvecg between EMSuite and Legacy
using MoM_Basics, MoM_Kernels
using EMSuite
using EMSuite.BasisFunctions: get_triangle_info
using EMSuite.IntegralEquations.EFIEModule: EFIE
using EMSuite.IntegralEquations.EFIEModule.Impedance.Geometry: get_global_quad_points
using EMSuite.IntegralEquations.EFIEModule.Singularities: faceSingularityIgIvecg as ems_faceSingularityIgIvecg
using LinearAlgebra, Printf

mesh_file = joinpath(@__DIR__, "plate_benchmark.nas")
freq = 3e8

# Setup Legacy
MoM_Basics.setPrecision!(Float64)
MoM_Basics.SimulationParams.SHOWIMAGE = false
MoM_Kernels.inputParameters(frequency=freq, ieT=:EFIE)
meshData, _ = MoM_Basics.getMeshData(mesh_file; meshUnit=:m)
ngeo_l, nbf_l, geosInfo_l, _ = MoM_Basics.getBFsFromMeshData(meshData; sbfT=:RWG)

# Setup EMSuite
mesh_ems = read_nas_mesh(mesh_file, scale=1.0)
set_frequency!(freq)
basis_ems = RWGBasis(mesh_ems)
efie = EFIE(freq)

# ===== Compare on first observation point =====
tri_e1 = get_triangle_info(mesh_ems, basis_ems, 1)
tri_e2 = get_triangle_info(mesh_ems, basis_ems, 2)

# Use first near GQ point of tri 1 as observation
gq_near = efie.gq_near
r_test = get_global_quad_points(tri_e1, gq_near)
robs = r_test[1]
println("Observation point: $robs")

# ===== EMSuite call =====
println("\n===== EMSuite faceSingularityIgIvecg =====")
println("  Source tri vertices: $(tri_e2.vertices)")
println("  Source tri edgel: $(tri_e2.edgel)")
println("  Source tri edgev̂: $(tri_e2.edgev̂)")
println("  Source tri edgen̂: $(tri_e2.edgen̂)")
println("  Source tri area: $(tri_e2.area)")
println("  Source tri facen̂: $(tri_e2.facen̂)")

Ig_ems, IvecSg_ems = ems_faceSingularityIgIvecg(
    robs, tri_e2.vertices, tri_e2.edgel, tri_e2.edgev̂,
    tri_e2.edgen̂, tri_e2.area, tri_e2.facen̂, efie.SSCg)

println("\nEMSuite Ig = $Ig_ems")
println("EMSuite IvecSg = $IvecSg_ems")

# ===== Legacy call =====
println("\n===== Legacy faceSingularityIgIvecg =====")
tris_l = geosInfo_l[2]
println("  Legacy tri vertices: $(tris_l.vertices)")
println("  Legacy tri edgel: $(tris_l.edgel)")
println("  Legacy tri area: $(tris_l.area)")

# Legacy faceSingularityIgIvecg(robs, tri, abs_area, facen)
Ig_leg, Ivecg_leg = MoM_Basics.faceSingularityIgIvecg(
    robs, tris_l, abs(tris_l.area), tris_l.facen̂)

println("\nLegacy Ig = $Ig_leg")
println("Legacy Ivecg = $Ivecg_leg")

# ===== Compare =====
println("\n===== Comparison =====")
println("Ig ratio: $(Ig_ems / Ig_leg)")
println("|Ig ratio|: $(abs(Ig_ems / Ig_leg))")

for i in 1:3
    ratio = IvecSg_ems[i] / Ivecg_leg[i]
    @printf("IvecSg[%d] ratio: %.6f + %.6fi (|r|=%.6f)\n", i, real(ratio), imag(ratio), abs(ratio))
end

# ===== Also check with Legacy abs edgel =====
println("\n===== Legacy with abs(edgel) =====")
# Manually call EMSuite's function with Legacy's data (unsigned edgel)
abs_edgel = abs.(tris_l.edgel)
println("Legacy abs edgel: $abs_edgel")
println("EMSuite edgel:    $(tri_e2.edgel)")
println("Same? $(abs_edgel ≈ tri_e2.edgel)")
