# Comprehensive Z element comparison: fix near interactino with /A_src, and verify far field
using EMSuite
using LinearAlgebra
using Printf

using EMSuite.IntegralEquations.EFIEModule: EFIE, calc_self_interaction!, calc_near_interaction!, calc_interaction!, is_adjacent
using EMSuite.IntegralEquations.EFIEModule.Impedance.Geometry: TriangleInfo, GaussQuadratureInfo, get_global_quad_points
using EMSuite.IntegralEquations.EFIEModule.Kernels: green_function_free_space
using EMSuite.IntegralEquations.EFIEModule.Singularities: faceSingularityIgIvecg, singularF1, singularF21, singularF22, greenfunc_star, compute_SSCg
using EMSuite.CoreModule: num_elements

const MOM_ALLINONE_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne")

mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
mesh = read_nas_mesh(mesh_file, scale=1.0)
freq = 1e8
set_frequency!(freq)
basis = RWGBasis(mesh)
efie = EFIE(freq)

println("k = $(efie.k), eta = $(efie.eta), factor = $(efie.factor)")

bf1 = basis.functions[1]
t_plus = bf1.support[1]
t_minus = bf1.support[2]
tri_plus = TriangleInfo(mesh, t_plus)
tri_minus = TriangleInfo(mesh, t_minus)
println("BF1: T+=$t_plus (area=$(tri_plus.area)), T-=$t_minus (area=$(tri_minus.area))")

# ==========================================
# Test A: Self interaction - EMSuite vs Legacy formula
# ==========================================
println("\n" * "="^60)
println("TEST A: Self interaction on T+ (id=$t_plus)")
println("="^60)

# EMSuite self interaction
Z_self_ems = zeros(ComplexF64, 3, 3)
calc_self_interaction!(Z_self_ems, efie, tri_plus)
Z_self_ems .*= efie.factor

# Legacy self interaction (manual implementation)
Z_self_leg = zeros(ComplexF64, 3, 3)
C4divk2 = efie.C4divk2
gq_near = efie.gq_near
r_pts = get_global_quad_points(tri_plus, gq_near)
w_pts = gq_near.weight
k = efie.k
area2 = tri_plus.area^2

# singularF1 value
sF1 = singularF1(tri_plus.edgel[1], tri_plus.edgel[2], tri_plus.edgel[3])
F1_val = C4divk2 * sF1

for ni in 1:3
    for mi in 1:3
        Ztemp = 0.0im
        # Smooth part: (e^{-jkR} - 1)/R
        for gj in 1:length(w_pts), gi in 1:length(w_pts)
            rgi = r_pts[gi]
            rgj = r_pts[gj]
            R = norm(rgi - rgj)
            G_star = greenfunc_star(R, k)
            rho_m = rgi - tri_plus.vertices[:, mi]
            rho_n = rgj - tri_plus.vertices[:, ni]
            Ztemp += (dot(rho_m, rho_n) - C4divk2) * G_star * w_pts[gi] * w_pts[gj]
        end
        # Singular part
        if mi == ni
            a = tri_plus.edgel[mi]
            b = tri_plus.edgel[mod1(mi+1,3)]
            c = tri_plus.edgel[mod1(mi+2,3)]
            Ztemp += singularF21(a, b, c, area2) - F1_val
        else
            k_idx = 6 - mi - ni
            a = tri_plus.edgel[k_idx]
            b = tri_plus.edgel[mod1(k_idx+1,3)]
            c = tri_plus.edgel[mod1(k_idx+2,3)]
            Ztemp += singularF22(a, b, c, area2) - F1_val
        end
        Z_self_leg[mi, ni] = Ztemp * tri_plus.edgel[mi] * tri_plus.edgel[ni] * efie.factor
    end
end

println("Self Z[1,1]: EMSuite=$(Z_self_ems[1,1]), Legacy=$(Z_self_leg[1,1])")
if abs(Z_self_leg[1,1]) > 1e-15
    r = Z_self_ems[1,1] / Z_self_leg[1,1]
    println("Self ratio: $(abs(r))")
end

# ==========================================
# Test B: Far interaction - EMSuite vs Legacy formula
# ==========================================
println("\n" * "="^60)
println("TEST B: Far-field interaction comparison")
println("="^60)

ntri = num_elements(mesh)
local t_far_id = 0
for t in 1:ntri
    tri_t = TriangleInfo(mesh, t)
    R_dist = norm(tri_plus.center - tri_t.center)
    if R_dist > 2.0 && R_dist < 3.0 && !is_adjacent(tri_plus, tri_t) && t != t_plus
        global t_far_id = t
        break
    end
end

if t_far_id > 0
    tri_far = TriangleInfo(mesh, t_far_id)
    println("Far triangle $t_far_id, dist=$(norm(tri_plus.center - tri_far.center))")
    
    gq_far = efie.gq_far
    r_test = get_global_quad_points(tri_plus, gq_far)
    r_src = get_global_quad_points(tri_far, gq_far)
    
    # EMSuite far
    Z_far_ems = zeros(ComplexF64, 3, 3)
    calc_interaction!(Z_far_ems, efie, tri_plus, tri_far, r_test, r_src)
    Z_far_ems .*= efie.factor
    
    # Legacy far (manual)
    Z_far_leg = zeros(ComplexF64, 3, 3)
    w_far = gq_far.weight
    for gj in 1:length(w_far), gi in 1:length(w_far)
        G = green_function_free_space(r_test[gi], r_src[gj], efie.k)
        gw = G * w_far[gi] * w_far[gj]
        for ni in 1:3, mi in 1:3
            rho_m = r_test[gi] - tri_plus.vertices[:, mi]
            rho_n = r_src[gj] - tri_far.vertices[:, ni]
            Z_far_leg[mi, ni] += (dot(rho_m, rho_n) - C4divk2) * gw
        end
    end
    for ni in 1:3, mi in 1:3
        Z_far_leg[mi, ni] *= tri_plus.edgel[mi] * tri_far.edgel[ni] * efie.factor
    end
    
    println("Far Z[1,1]: EMSuite=$(Z_far_ems[1,1]), Legacy=$(Z_far_leg[1,1])")
    if abs(Z_far_leg[1,1]) > 1e-15
        r = Z_far_ems[1,1] / Z_far_leg[1,1]
        println("Far ratio: $(abs(r))")
    end
else
    println("No suitable far triangle found!")
end

# ==========================================
# Test C: Excitation vector comparison
# ==========================================
println("\n" * "="^60)
println("TEST C: Excitation vector comparison")
println("="^60)

source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
V_ems = excitation_vector(efie, source, basis)

# Legacy excitation formula (manual):
# V[n] = Σ_k sign_k * l_n/2 * Σ_i w_i * ρ(r_i) · E_inc(r_i)
# For EFIE: V = �?f_n · E_inc dS

# Compute V[1] manually
V1_manual = 0.0im
gq3 = GaussQuadratureInfo(:Triangle, 3, Float64)
for sc in 1:2
    t_i = bf1.support[sc]
    tri_i = TriangleInfo(mesh, t_i)
    sign_val = bf1.signs[sc]
    local_edge = bf1.local_edge_idx[sc]
    area_i = tri_i.area
    
    r_quad = get_global_quad_points(tri_i, gq3)
    w_quad = gq3.weight
    
    for q in 1:length(w_quad)
        rqi = r_quad[q]
        v_free = tri_i.vertices[:, local_edge]
        rho = rqi - v_free
        f_val = sign_val * (bf1.edge_length / (2 * area_i)) * rho
        E_inc = EMSuite.CoreModule.incident_field(source, rqi)
        V1_manual += area_i * w_quad[q] * dot(f_val, E_inc)
    end
end

println("V[1] EMSuite: $(V_ems[1])")
println("V[1] Manual:  $(V1_manual)")
println("V ratio: $(abs(V_ems[1] / V1_manual))")

println("\nDone.")
