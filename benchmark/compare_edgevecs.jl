# Compare edge vectors and normals between Legacy and EMSuite
using MoM_Basics, MoM_Kernels
using EMSuite
using EMSuite.BasisFunctions: get_triangle_info
using EMSuite.IntegralEquations.EFIEModule: EFIE
using EMSuite.IntegralEquations.EFIEModule.Impedance.Geometry: get_global_quad_points
using EMSuite.IntegralEquations.EFIEModule.Singularities: faceSingularityIgIvecg as ems_fSIgIvecg, compute_SSCg
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

# Triangle 2 (source for near-interaction)
tri_l = geosInfo_l[2]
tri_e = get_triangle_info(mesh_ems, basis_ems, 2)

println("===== Triangle 2 edgev̂ =====")
for i in 1:3
    ev_l = tri_l.edgev̂[:, i]
    ev_e = tri_e.edgev̂[:, i]
    @printf("Edge %d:\n  Legacy: [%.6f, %.6f, %.6f]\n  EMSuite:[%.6f, %.6f, %.6f]\n  Same? %s\n",
        i, ev_l[1], ev_l[2], ev_l[3], ev_e[1], ev_e[2], ev_e[3], ev_l ≈ ev_e)
end

println("\n===== Triangle 2 edgen̂ =====")
for i in 1:3
    en_l = tri_l.edgen̂[:, i]
    en_e = tri_e.edgen̂[:, i]
    @printf("Edge %d:\n  Legacy: [%.6f, %.6f, %.6f]\n  EMSuite:[%.6f, %.6f, %.6f]\n  Same? %s\n",
        i, en_l[1], en_l[2], en_l[3], en_e[1], en_e[2], en_e[3], en_l ≈ en_e)
end

println("\n===== Triangle 2 facen̂ =====")
fn_l = tri_l.facen̂
fn_e = tri_e.facen̂
@printf("Legacy: [%.6f, %.6f, %.6f]\nEMSuite:[%.6f, %.6f, %.6f]\nSame? %s\n",
    fn_l[1], fn_l[2], fn_l[3], fn_e[1], fn_e[2], fn_e[3], fn_l ≈ fn_e)

# Compare SSCg coefficients
println("\n===== SSCg coefficients =====")
SSCg_ems = efie.SSCg
# Legacy SSCg is an OffsetArray with index 0:SglrOrder-1
# We'll extract it manually
k = MoM_Basics.Params.K_0
jk = im * k
for n in 0:5
    ems_val = SSCg_ems[n+1]  # EMSuite 1-indexed
    # Legacy: coeffgreen(n) = (-jk)^n / n!
    # Legacy SSCg[n] = coeffgreen(n) (OffsetArray with offset -1)
    leg_val = (-jk)^n / factorial(n)
    @printf("n=%d: EMSuite SSCg[%d] = %.6f%+.6fi, Legacy = %.6f%+.6fi, Same? %s\n",
        n, n+1, real(ems_val), imag(ems_val), real(leg_val), imag(leg_val), ems_val ≈ leg_val)
end

# Now call BOTH faceSingularityIgIvecg with the SAME inputs (EMSuite's data)
println("\n===== Call both singularity functions with SAME inputs =====")
tri_e1 = get_triangle_info(mesh_ems, basis_ems, 1)
gq_near = efie.gq_near
r_test = get_global_quad_points(tri_e1, gq_near)
robs = r_test[1]

# EMSuite call
Ig_ems, Ivec_ems = ems_fSIgIvecg(
    robs, tri_e.vertices, tri_e.edgel, tri_e.edgev̂,
    tri_e.edgen̂, tri_e.area, tri_e.facen̂, efie.SSCg)
println("EMSuite: Ig = $Ig_ems")
println("EMSuite: Ivec = $Ivec_ems")

# Legacy call with Legacy's triangle info (has signed edgel)
Ig_leg, Ivec_leg = MoM_Kernels.faceSingularityIgIvecg(
    robs, tri_l, abs(tri_l.area), tri_l.facen̂)
println("\nLegacy:  Ig = $Ig_leg")
println("Legacy:  Ivec = $Ivec_leg")

println("\nIg ratio: $(Ig_ems / Ig_leg)")
println("|Ig ratio|: $(abs(Ig_ems / Ig_leg))")

for i in 1:3
    if abs(Ivec_leg[i]) > 1e-20
        ratio = Ivec_ems[i] / Ivec_leg[i]
        @printf("Ivec[%d] ratio: %.6f + %.6fi\n", i, real(ratio), imag(ratio))
    end
end

# Also test: use Legacy's actual edgev and edgen in EMSuite's function
println("\n===== EMSuite function with Legacy's edgev̂/edgen̂ =====")
Ig_ems2, Ivec_ems2 = ems_fSIgIvecg(
    robs, tri_l.vertices, abs.(tri_l.edgel), tri_l.edgev̂,
    tri_l.edgen̂, abs(tri_l.area), tri_l.facen̂, efie.SSCg)
println("EMSuite(legacy data): Ig = $Ig_ems2")
println("EMSuite(legacy data): Ivec = $Ivec_ems2")
println("Ig ratio vs Legacy: $(Ig_ems2 / Ig_leg)")
