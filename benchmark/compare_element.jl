# Direct element-by-element comparison using BOTH Legacy and EMSuite APIs
using MoM_Basics, MoM_Kernels
using EMSuite
using EMSuite.BasisFunctions: get_triangle_info
using EMSuite.IntegralEquations.EFIEModule: EFIE, efie_interaction!, calc_interaction!, calc_self_interaction!, calc_near_interaction!, is_adjacent
using EMSuite.IntegralEquations.EFIEModule.Impedance.Geometry: get_global_quad_points, GaussQuadratureInfo
using LinearAlgebra, Printf

mesh_file = joinpath(@__DIR__, "plate_benchmark.nas")
freq = 3e8

# ========== Setup Legacy ==========
MoM_Basics.setPrecision!(Float64)
MoM_Basics.SimulationParams.SHOWIMAGE = false
MoM_Kernels.inputParameters(frequency=freq, ieT=:EFIE)
meshData, _ = MoM_Basics.getMeshData(mesh_file; meshUnit=:m)
ngeo_l, nbf_l, geosInfo_l, _ = MoM_Basics.getBFsFromMeshData(meshData; sbfT=:RWG)

# ========== Setup EMSuite ==========
mesh_ems = read_nas_mesh(mesh_file, scale=1.0)
set_frequency!(freq)
basis_ems = RWGBasis(mesh_ems)
efie = EFIE(freq)

println("EMSuite factor: $(efie.factor)")
println("Legacy JKηdiv16π: $(MoM_Basics.Params.JKηdiv16π)")
println("Factor ratio: $(efie.factor / MoM_Basics.Params.JKηdiv16π)")

# ========== Test 1: Self interaction (tri 1) ==========
println("\n===== Self Interaction (Triangle 1) =====")
tri_e1 = get_triangle_info(mesh_ems, basis_ems, 1)

Z_self_ems = zeros(ComplexF64, 3, 3)
calc_self_interaction!(Z_self_ems, efie, tri_e1)
Z_self_ems .*= efie.factor

Z_self_leg = MoM_Kernels.EFIEOnTris(geosInfo_l[1])

println("EMSuite Z_self[1,1] = $(Z_self_ems[1,1])")
println("Legacy  Z_self[1,1] = $(Z_self_leg[1,1])")
println("Ratio: $(Z_self_ems[1,1] / Z_self_leg[1,1])")
println("|ratio| = $(abs(Z_self_ems[1,1] / Z_self_leg[1,1]))")

# Show all elements
for mi in 1:3, ni in 1:3
    if abs(Z_self_leg[mi,ni]) > 1e-20
        ratio = Z_self_ems[mi,ni] / Z_self_leg[mi,ni]
        @printf("  [%d,%d] ratio = %.6f + %.6fi (|r|=%.6f)\n", mi, ni, real(ratio), imag(ratio), abs(ratio))
    end
end

# ========== Test 2: Far interaction (tri 1 vs distant tri) ==========
# Find a distant triangle
distant_tri = 100
println("\n===== Far Interaction (Tri 1 vs Tri $distant_tri) =====")
tri_e2 = get_triangle_info(mesh_ems, basis_ems, distant_tri)
println("Adjacent? $(is_adjacent(tri_e1, tri_e2))")

Z_far_ems = zeros(ComplexF64, 3, 3)
r_test = get_global_quad_points(tri_e1, efie.gq_far)
r_src = get_global_quad_points(tri_e2, efie.gq_far)
calc_interaction!(Z_far_ems, efie, tri_e1, tri_e2, r_test, r_src)
Z_far_ems .*= efie.factor

Z_far_leg = MoM_Kernels.EFIEOnTris(geosInfo_l[1], geosInfo_l[distant_tri])

println("EMSuite Z_far[1,1] = $(Z_far_ems[1,1])")
println("Legacy  Z_far[1,1] = $(Z_far_leg[1,1])")
if abs(Z_far_leg[1,1]) > 1e-20
    println("Ratio: $(Z_far_ems[1,1] / Z_far_leg[1,1])")
    println("|ratio| = $(abs(Z_far_ems[1,1] / Z_far_leg[1,1]))")
end

for mi in 1:3, ni in 1:3
    if abs(Z_far_leg[mi,ni]) > 1e-20
        ratio = Z_far_ems[mi,ni] / Z_far_leg[mi,ni]
        @printf("  [%d,%d] ratio = %.6f + %.6fi (|r|=%.6f)\n", mi, ni, real(ratio), imag(ratio), abs(ratio))
    end
end

# ========== Test 3: Near interaction (tri 1 vs tri 2, if adjacent) ==========
tri_e2_near = get_triangle_info(mesh_ems, basis_ems, 2)
println("\n===== Near Interaction (Tri 1 vs Tri 2) =====")
println("Adjacent? $(is_adjacent(tri_e1, tri_e2_near))")
println("Same tri? $(tri_e1.triID == tri_e2_near.triID)")

# Check Legacy distance
Rts = norm(geosInfo_l[1].center - geosInfo_l[2].center)
println("Legacy center distance: $Rts, Rsglr: $(MoM_Basics.Params.Rsglr)")
println("Is near in Legacy? $(Rts < MoM_Basics.Params.Rsglr)")

Z_near_ems = zeros(ComplexF64, 3, 3)
calc_near_interaction!(Z_near_ems, efie, tri_e1, tri_e2_near)
Z_near_ems .*= efie.factor

Z_near_leg = MoM_Kernels.EFIEOnNearTris(geosInfo_l[1], geosInfo_l[2])

println("EMSuite Z_near[1,1] = $(Z_near_ems[1,1])")
println("Legacy  Z_near[1,1] = $(Z_near_leg[1,1])")

for mi in 1:3, ni in 1:3
    if abs(Z_near_leg[mi,ni]) > 1e-20
        ratio = Z_near_ems[mi,ni] / Z_near_leg[mi,ni]
        @printf("  [%d,%d] ratio = %.6f + %.6fi (|r|=%.6f)\n", mi, ni, real(ratio), imag(ratio), abs(ratio))
    end
end

# ========== Test 4: Check assembly contribution to Z[1,1] ==========
println("\n===== Assembly Contribution to Z[1,1] =====")
# BF1 exists on tri 1 (local edge 1, sign +1) and tri 2 (local edge 2, sign -1)
# Z[1,1] = Self(tri1)[1,1]*sign1*sign1 + Self(tri2)[2,2]*sign2*sign2 + Near(tri1,tri2)[1,2]*sign1*sign2*2

# Self from tri 1
Z_s1 = zeros(ComplexF64, 3, 3)
calc_self_interaction!(Z_s1, efie, tri_e1)
Z_s1 .*= efie.factor
contrib_self1 = Z_s1[1,1] * 1 * 1  # sign^2 = 1
println("Self(tri1) contribution: $contrib_self1 (|Z|=$(abs(contrib_self1)))")

# Self from tri 2
Z_s2 = zeros(ComplexF64, 3, 3)
calc_self_interaction!(Z_s2, efie, tri_e2_near)
Z_s2 .*= efie.factor
# BF1 is local edge 2 on tri 2 
contrib_self2 = Z_s2[2,2] * (-1) * (-1)  # sign^2 = 1
println("Self(tri2) contribution: $contrib_self2 (|Z|=$(abs(contrib_self2)))")

# Near from (tri1, tri2) - tri1 has BF1 at local 1, tri2 has BF1 at local 2
Z_n12 = zeros(ComplexF64, 3, 3)
calc_near_interaction!(Z_n12, efie, tri_e1, tri_e2_near)
Z_n12 .*= efie.factor
# Contribution: Z_n12[1,2] * sign_test(+1) * sign_src(-1) * 2 (symmetric fill)
contrib_near = Z_n12[1,2] * (+1) * (-1) * 2  # doubled by symmetric fill
println("Near(tri1,tri2) contribution: $contrib_near (|Z|=$(abs(contrib_near)))")

total_ems = contrib_self1 + contrib_self2 + contrib_near
println("\nTotal EMSuite Z[1,1] from these: $total_ems (|Z|=$(abs(total_ems)))")
println("Actual EMSuite Z[1,1]: $(Z_self_ems)")  # This won't be right, let me compute properly

# Now compare with Legacy
Z_leg_full = MoM_Kernels.impedancemat4EFIE4PEC(geosInfo_l, nbf_l, MoM_Basics.RWG)
println("Actual Legacy Z[1,1]: $(Z_leg_full[1,1])")

# But also, there might be more adjacent triangles contributing to Z[1,1]
# Find ALL triangles that contain BF1 edges
println("\n--- All near-interaction pairs that contribute to Z[1,1] ---")
ntri = EMSuite.CoreModule.num_elements(mesh_ems)
for t1 in 1:ntri
    tri1 = get_triangle_info(mesh_ems, basis_ems, t1)
    # Check if BF1 is on this triangle
    local_idx1 = 0
    sign1 = 0
    for k in 1:3
        if tri1.inBfsID[k] == 1
            local_idx1 = k
            sign1 = tri1.bfsSign[k]
            break
        end
    end
    local_idx1 == 0 && continue
    
    for t2 in t1:ntri
        tri2 = get_triangle_info(mesh_ems, basis_ems, t2)
        local_idx2 = 0
        sign2 = 0
        for k in 1:3
            if tri2.inBfsID[k] == 1
                local_idx2 = k
                sign2 = tri2.bfsSign[k]
                break
            end
        end
        local_idx2 == 0 && continue
        
        if t1 == t2
            # Self
            Z_tmp = zeros(ComplexF64, 3, 3)
            calc_self_interaction!(Z_tmp, efie, tri1)
            Z_tmp .*= efie.factor
            c = Z_tmp[local_idx1, local_idx2] * sign1 * sign2
            println("  Self tri $t1: Z_local[$local_idx1,$local_idx2]*$sign1*$sign2 = $c")
        elseif is_adjacent(tri1, tri2)
            Z_tmp = zeros(ComplexF64, 3, 3)
            calc_near_interaction!(Z_tmp, efie, tri1, tri2)
            Z_tmp .*= efie.factor
            c = Z_tmp[local_idx1, local_idx2] * sign1 * sign2
            c2 = c * 2  # symmetric fill doubles diagonal
            println("  Near tri $t1-$t2: Z_local[$local_idx1,$local_idx2]*$sign1*$sign2 = $c (×2 = $c2)")
        else
            Z_tmp = zeros(ComplexF64, 3, 3)
            r_t = get_global_quad_points(tri1, efie.gq_far)
            r_s = get_global_quad_points(tri2, efie.gq_far)
            calc_interaction!(Z_tmp, efie, tri1, tri2, r_t, r_s)
            Z_tmp .*= efie.factor
            c = Z_tmp[local_idx1, local_idx2] * sign1 * sign2
            c2 = c * 2
            println("  Far tri $t1-$t2: Z_local[$local_idx1,$local_idx2]*$sign1*$sign2 = $c (×2 = $c2)")
        end
    end
end
