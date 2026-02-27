# Deep diagnosis: compare near-interaction Z contributions with Legacy formula
using EMSuite
using LinearAlgebra
using Printf

# Access internal modules
using EMSuite.IntegralEquations.EFIEModule: EFIE, efie_interaction!, calc_self_interaction!, calc_near_interaction!, calc_interaction!, is_adjacent
using EMSuite.IntegralEquations.EFIEModule.Impedance.Geometry: TriangleInfo, GaussQuadratureInfo, get_global_quad_points
using EMSuite.IntegralEquations.EFIEModule.Kernels: green_function_free_space
using EMSuite.IntegralEquations.EFIEModule.Singularities: faceSingularityIgIvecg
using EMSuite.CoreModule: num_elements

const MOM_ALLINONE_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne")

# Load mesh
mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
mesh = read_nas_mesh(mesh_file, scale=1.0)
freq = 1e8
set_frequency!(freq)
basis = RWGBasis(mesh)

efie = EFIE(freq)
println("k = $(efie.k), eta = $(efie.eta)")
println("factor = $(efie.factor)")
println("C4divk2 = $(efie.C4divk2)")

# Find basis function 1's support triangles
bf1 = basis.functions[1]
t_plus = bf1.support[1]
t_minus = bf1.support[2]
println("\nBasis function 1: T+=$t_plus, T-=$t_minus")

tri_plus = TriangleInfo(mesh, t_plus)
tri_minus = TriangleInfo(mesh, t_minus)

println("T+ area: $(tri_plus.area), T- area: $(tri_minus.area)")

# ==========================================
# Test 1: Self interaction on T+
# ==========================================
println("\n" * "=" ^ 50)
println("Test 1: Self interaction Z_local on T+ (id=$t_plus)")
println("=" ^ 50)

Z_self = zeros(ComplexF64, 3, 3)
calc_self_interaction!(Z_self, efie, tri_plus)
# Note: efie.factor not yet applied
Z_self_pre = copy(Z_self)  # before factor
Z_self .*= efie.factor
println("Z_self (before factor):")
for m in 1:3, n in 1:3
    @printf("  [%d,%d] = %+.6e %+.6ej\n", m, n, real(Z_self_pre[m,n]), imag(Z_self_pre[m,n]))
end
println("Z_self (after factor):")
for m in 1:3, n in 1:3
    @printf("  [%d,%d] = %+.6e %+.6ej\n", m, n, real(Z_self[m,n]), imag(Z_self[m,n]))
end

# ==========================================
# Test 2: Near interaction T+ → T-
# ==========================================
println("\n" * "=" ^ 50)
println("Test 2: Near interaction T+($t_plus) → T-($t_minus)")
println("=" ^ 50)

Z_near = zeros(ComplexF64, 3, 3)
calc_near_interaction!(Z_near, efie, tri_plus, tri_minus)
Z_near_pre = copy(Z_near)
Z_near .*= efie.factor
println("Z_near (before factor):")
for m in 1:3, n in 1:3
    @printf("  [%d,%d] = %+.6e %+.6ej\n", m, n, real(Z_near_pre[m,n]), imag(Z_near_pre[m,n]))
end
println("Z_near (after factor):")
for m in 1:3, n in 1:3
    @printf("  [%d,%d] = %+.6e %+.6ej\n", m, n, real(Z_near[m,n]), imag(Z_near[m,n]))
end

# ==========================================
# Test 3: Recompute near interaction using Legacy formula
# (manually, not calling calc_near_interaction!)
# ==========================================
println("\n" * "=" ^ 50)
println("Test 3: Near interaction using LEGACY formula (manual)")
println("=" ^ 50)

# Legacy formula:
# Ztemp += ((ρm·ρn - C4/k²)*Ig - ρm·Ivec) * weight[gi] / A_src
# Ztemp *= lm * ln * JKηdiv16π

gq_near = efie.gq_near
r_test = get_global_quad_points(tri_plus, gq_near)
w_test = gq_near.weight
C4divk2 = efie.C4divk2

Z_legacy = zeros(ComplexF64, 3, 3)

# Print first Ig value to compare
r0_test = r_test[1]
Ig0, IvecSg0 = faceSingularityIgIvecg(
    r0_test, tri_minus.vertices,
    tri_minus.edgel, tri_minus.edgev̂,
    tri_minus.edgen̂, tri_minus.area,
    tri_minus.facen̂, efie.SSCg)
println("\nAt test point 1:")
println("  Ig = $Ig0")
println("  IvecSg = $IvecSg0")
println("  |Ig| = $(abs(Ig0))")
println("  |IvecSg| = $(norm(IvecSg0))")
println("  A_src = $(tri_minus.area)")
println("  Ig / A_src = $(Ig0 / tri_minus.area)")

for gi in 1:length(w_test)
    ri = r_test[gi]
    
    Ig, IvecSg = faceSingularityIgIvecg(
        ri, tri_minus.vertices,
        tri_minus.edgel, tri_minus.edgev̂,
        tri_minus.edgen̂, tri_minus.area,
        tri_minus.facen̂, efie.SSCg)
    
    for ni in 1:3
        rho_n = ri - tri_minus.vertices[:, ni]  # ρn at test point
        for mi in 1:3
            rho_m = ri - tri_plus.vertices[:, mi]   # ρm at test point
            
            # Legacy formula (exactly):
            Ztemp = ((dot(rho_m, rho_n) - C4divk2) * Ig - dot(rho_m, IvecSg)) * w_test[gi] / tri_minus.area
            
            Z_legacy[mi, ni] += Ztemp
        end
    end
end

# Apply edge lengths and factor
for ni in 1:3
    ln = tri_minus.edgel[ni]
    for mi in 1:3
        lm = tri_plus.edgel[mi]
        Z_legacy[mi, ni] *= lm * ln * efie.factor
    end
end

println("\nZ_legacy (Legacy formula):")
for m in 1:3, n in 1:3
    @printf("  [%d,%d] = %+.6e %+.6ej\n", m, n, real(Z_legacy[m,n]), imag(Z_legacy[m,n]))
end

# Ratios
println("\nRatio Z_near/Z_legacy (EMSuite/Legacy):")
for m in 1:3, n in 1:3
    if abs(Z_legacy[m,n]) > 1e-15
        ratio = Z_near[m,n] / Z_legacy[m,n]
        @printf("  [%d,%d] ratio = %+.6f %+.6fj  (|ratio|=%.6f)\n", m, n, real(ratio), imag(ratio), abs(ratio))
    end
end

# ==========================================
# Test 4: Find a far-field pair and compare
# ==========================================
println("\n" * "=" ^ 50)
println("Test 4: Far-field interaction comparison")
println("=" ^ 50)

# Find a far triangle (not adjacent to T+)
t_far = 0
for t in 1:num_elements(mesh)
    if t != t_plus
        tri_far = TriangleInfo(mesh, t)
        R = norm(tri_plus.center - tri_far.center)
        if R > 2.0 && R < 3.0  # Find one at moderate distance
            if !is_adjacent(tri_plus, tri_far)
                t_far = t
                break
            end
        end
    end
end

if t_far > 0
    tri_far = TriangleInfo(mesh, t_far)
    R_dist = norm(tri_plus.center - tri_far.center)
    println("Far triangle: $t_far, distance=$R_dist")
    
    Z_far_emsuite = zeros(ComplexF64, 3, 3)
    r_src = get_global_quad_points(tri_far, efie.gq_far)
    r_test_far = get_global_quad_points(tri_plus, efie.gq_far)
    calc_interaction!(Z_far_emsuite, efie, tri_plus, tri_far, r_test_far, r_src)
    Z_far_emsuite .*= efie.factor
    
    # Now compute using direct Legacy-style formula
    gq_far = efie.gq_far
    w_far = gq_far.weight
    Z_far_legacy = zeros(ComplexF64, 3, 3)
    
    for gj in 1:length(w_far)
        rgj = r_src[gj]
        for gi in 1:length(w_far)
            rgi = r_test_far[gi]
            G = green_function_free_space(rgi, rgj, efie.k)
            gw = G * w_far[gi] * w_far[gj]
            for ni in 1:3
                rho_n = rgj - tri_far.vertices[:, ni]
                for mi in 1:3
                    rho_m = rgi - tri_plus.vertices[:, mi]
                    Z_far_legacy[mi, ni] += (dot(rho_m, rho_n) - C4divk2) * gw
                end
            end
        end
    end
    
    for ni in 1:3
        ln = tri_far.edgel[ni]
        for mi in 1:3
            lm = tri_plus.edgel[mi]
            Z_far_legacy[mi, ni] *= lm * ln * efie.factor
        end
    end
    
    println("Ratio Z_far_emsuite/Z_far_legacy:")
    for m in 1:3, n in 1:3
        if abs(Z_far_legacy[m,n]) > 1e-15
            ratio = Z_far_emsuite[m,n] / Z_far_legacy[m,n]
            @printf("  [%d,%d] ratio = %+.6f %+.6fj  (|ratio|=%.6f)\n", m, n, real(ratio), imag(ratio), abs(ratio))
        end
    end
end

println("\nDone.")
