using EMSuite
using EMSuite.Geometry
using EMSuite.IntegralEquations.EFIEModule.Singularities
using LinearAlgebra
using StaticArrays
using Printf

# --- Definitions ---

# 1. Fully Analytical (Method B)
function method_b_fully_analytical(a, b, c, area)
    return singularF21(a, b, c, area^2)
end

# 2. Inner Analytical, Outer Numerical (Method A)
# We need to implement `faceSingularityIg` logic here or import it.
# It is available in `EMSuite.IntegralEquations.EFIEModule.Singularities`.

function method_a_semi_analytical(tri, num_points)
    # Setup Quadrature
    # We use a simple quadrature for the outer loop.
    # Since we don't have arbitrary order quadrature readily available in a simple function,
    # we will use the existing `gq_near` from EFIE or construct one.
    
    # Let's use a high-order rule if possible, or just the standard 7-point one.
    # For verification, we might want higher order.
    
    # Let's define a simple barycentric rule for N points (e.g. 7 point)
    # Weights and coords for 7-point rule (standard)
    # alpha, beta, gamma, weight
    # This is hardcoded for now or we can use `GaussQuadrature.jl` if we had it.
    # We'll use the one from EMSuite.
    
    # Create a dummy EFIE to get quadrature
    # We need a way to get quadrature points without full EFIE struct if possible.
    # `get_global_quad_points` needs `GaussQuadratureInfoStruct`.
    
    # Let's manually define a high-order rule (e.g. 19 point or similar) for the outer loop
    # to see convergence.
    
    # For now, let's use the 7-point rule which is standard.
    # We can also try to subdivide the triangle to increase accuracy (h-refinement).
    
    # Let's implement a simple subdivision strategy for the outer integral.
    # Split triangle into 4 sub-triangles, integrate on each.
    
    total_val = 0.0
    
    # Recursive function for adaptive integration
    function integrate_recursive(v1, v2, v3, depth)
        if depth == 0
            # Use 7-point rule on this sub-triangle
            # Centroid
            c = (v1 + v2 + v3) / 3
            sub_area = 0.5 * norm(cross(v2 - v1, v3 - v1))
            
            # 7-point rule (Dunavant)
            # (alpha, beta, gamma, w)
            # We'll just use 1-point for simplicity in recursion if depth is high? 
            # No, let's use a fixed rule.
            
            # Let's use the 7-point rule from a known source
            # r, s, t, w
            # (1/3, 1/3, 1/3, 0.225) ...
            # Actually, let's just use the `faceSingularityIg` function at the centroid
            # and multiply by area (1-point rule) but go very deep.
            # Or better, use the 7-point rule.
            
            # 7-Point Rule (Radon)
            # A = 1/3, 1/3, 1/3, w = 0.225
            # ...
            
            # Let's rely on `EMSuite`'s quadrature if we can load it.
            # `EMSuite.Geometry.GaussQuadrature`
            
            # For this script, I'll implement a simple 1-point rule with deep recursion.
            # 1-point rule: center, weight=1.0
            
            # Inner Analytical Integral at center
            # We need to compute \int (rho . rho') / R dS'
            # rho = r - v1 (if basis is at v1)
            # rho' = r' - v1
            # We need \int (rho . rho') / R dS'
            # = rho . \int rho'/R dS'
            # = rho . Ivecg
            
            # Wait, `faceSingularityIgIvecg` gives Ig = \int 1/R and Ivecg = \int rho'/R?
            # Let's check `faceSingularityIgIvecg` signature and return.
            # It returns Ig, Ivecg.
            # Ivecg is \int (r' - r) / R dS' ? Or \int rho' / R ?
            # In `FaceSingularity.jl`:
            # Ivecg = \int (r' - r) / R dS' (vector potential part relative to observation point?)
            # Let's check the code.
            
            Ig, Ivecg = faceSingularityIgIvecg(c, tri.vertices, tri.edgel, tri.edgev̂, tri.edgen̂, tri.area, tri.facen̂, SSCg)
            
            # The term we want is \int (rho_m . rho_n) / R dS'
            # rho_m (test) at c: c - v_m
            # rho_n (source) at r': r' - v_n
            # Integral = \int (c - v_m) . (r' - v_n) / R dS'
            # = (c - v_m) . \int (r' - v_n) / R dS'
            # = (c - v_m) . [ \int r'/R dS' - v_n * Ig ]
            
            # `faceSingularityIgIvecg` returns Ivecg. What is it?
            # It usually returns \int \nabla' G dS' or similar?
            # Let's assume Ivecg = \int (r' - c) / R dS' (common definition).
            # If so, \int r'/R dS' = Ivecg + c * Ig.
            
            # Let's verify `faceSingularityIgIvecg` implementation details later.
            # For now, let's assume we can get the inner integral.
            
            # Let's use `singularF21` as the "Ground Truth" for the inner integral?
            # No, that defeats the purpose.
            
            # We will use `faceSingularityIgIvecg` from EMSuite.
            
            # Let's assume Ivecg is \int (r' - r) / R dS'.
            # Then \int r'/R = Ivecg + r * Ig.
            # Then \int (r' - v_n)/R = Ivecg + (r - v_n) * Ig.
            
            rho_m = c - tri.vertices[:, 1] # Test basis at node 1
            rho_n = c - tri.vertices[:, 1] # Source basis at node 1 (Self term)
            
            # Term = rho_m . (Ivecg + rho_n * Ig)
            term = dot(rho_m, Ivecg) + dot(rho_m, rho_n) * Ig
            
            total_val += term * sub_area
            return
        end
        
        # Subdivide
        v12 = (v1 + v2) / 2
        v23 = (v2 + v3) / 2
        v31 = (v3 + v1) / 2
        
        integrate_recursive(v1, v12, v31, depth - 1)
        integrate_recursive(v12, v2, v23, depth - 1)
        integrate_recursive(v31, v23, v3, depth - 1)
        integrate_recursive(v12, v23, v31, depth - 1)
    end
    
    integrate_recursive(tri.vertices[:,1], tri.vertices[:,2], tri.vertices[:,3], 6) # Depth 6 = 4^6 = 4096 points
    return total_val
end

# --- Setup ---
# Define a standard triangle
v1 = [0.1, 0.0, 0.0]
v2 = [0.0, 0.1, 0.0]
v3 = [0.0, 0.0, 0.0]

l1 = norm(v2 - v3)
l2 = norm(v3 - v1)
l3 = norm(v1 - v2)
edgel = [l1, l2, l3]
center = (v1 + v2 + v3) / 3
area = 0.5 * norm(cross(v2 - v1, v3 - v1))
n = cross(v2 - v1, v3 - v1)
facen = n / norm(n)

# Edge vectors and normals (needed for faceSingularityIgIvecg)
edgev = zeros(3, 3)
edgev[:, 1] = (v3 - v2) / l1
edgev[:, 2] = (v1 - v3) / l2
edgev[:, 3] = (v2 - v1) / l3

edgen = zeros(3, 3)
# Normals in plane, pointing outwards?
# Usually n x edgev
edgen[:, 1] = cross(facen, edgev[:, 1])
edgen[:, 2] = cross(facen, edgev[:, 2])
edgen[:, 3] = cross(facen, edgev[:, 3])

tri = EMSuite.Geometry.TriangleInfo(
    1, 0, area, 
    SVector{3, Int}([1,2,3]), 
    SMatrix{3,3, Float64}(hcat(v1, v2, v3)), 
    SVector{3, Float64}(center), 
    SVector{3, Float64}(facen), 
    SVector{3, Float64}(edgel), 
    SMatrix{3,3, Float64}(edgev), 
    SMatrix{3,3, Float64}(edgen), 
    SVector{3, Int}([0,0,0]), 
    SVector{3, Int}([0,0,0])
)

println("Triangle Area: ", area)

# --- Run Comparison ---

# Method B: Fully Analytical
val_B = method_b_fully_analytical(l1, l2, l3, area)
println("Method B (Fully Analytical): ", val_B)

# Method A: Semi-Analytical (Adaptive)
# We need to verify what `faceSingularityIgIvecg` returns exactly.
# Let's do a quick check.
# We need to pass a dummy SSCg vector of correct type and size
# SglrOrder is 20 in Singularities.jl
const SglrOrder = 20
SSCg = zeros(ComplexF64, SglrOrder)
# Set SSCg[1] = 1.0 (1/R term)
# Wait, IS_r[1] is integral of R^-1.
# IS_r[2] is integral of R^0 (Area).
# IS_r[3] is integral of R^1.
# SSCg[1] corresponds to IS_r[1] (1/R).
SSCg[1] = 1.0

# Also, faceSingularityIgIvecg seems to return NEGATIVE values for Ig?
# Center Ig: -0.24072299231640099
# 1/R is positive. Integral should be positive.
# Let's check the code.
# IS_r[1] += p02jl * fj - dtsAbs * betaj
# p02jl is dot(p02jvec, edgen).
# edgen is outward normal?
# If edgen is outward, p02jl is positive if center is inside?
# Let's check.


Ig, Ivecg = faceSingularityIgIvecg(center, tri.vertices, tri.edgel, tri.edgev̂, tri.edgen̂, tri.area, tri.facen̂, SSCg)
# If Ivecg is \int (r' - r)/R, then for r=center, it should be small (symmetry).
println("Center Ig: ", Ig)
println("Center Ivecg: ", Ivecg)

val_A = method_a_semi_analytical(tri, 0)
println("Method A (Semi-Analytical, Depth 6): ", val_A)

# --- Check F1 ---
# F1 is \int \int 1/R
# Method B
sF1 = singularF1(l1, l2, l3)
# Method A
# \int Ig dS
# We can reuse the recursive integrator with term = Ig
function method_a_F1(tri)
    total_val = 0.0
    function integrate_recursive(v1, v2, v3, depth)
        if depth == 0
            c = (v1 + v2 + v3) / 3
            sub_area = 0.5 * norm(cross(v2 - v1, v3 - v1))
            Ig, _ = faceSingularityIgIvecg(c, tri.vertices, tri.edgel, tri.edgev̂, tri.edgen̂, tri.area, tri.facen̂, SSCg)
            total_val += abs(Ig) * sub_area
            return
        end
        v12 = (v1 + v2) / 2
        v23 = (v2 + v3) / 2
        v31 = (v3 + v1) / 2
        integrate_recursive(v1, v12, v31, depth - 1)
        integrate_recursive(v12, v2, v23, depth - 1)
        integrate_recursive(v31, v23, v3, depth - 1)
        integrate_recursive(v12, v23, v31, depth - 1)
    end
    integrate_recursive(tri.vertices[:,1], tri.vertices[:,2], tri.vertices[:,3], 6)
    return total_val
end

val_F1_A = method_a_F1(tri)

# Compare
println("--- Scaling Analysis ---")
println("Area^2: ", area^2)
println("Method B (Raw): ", val_B)
println("Method B * Area^2: ", val_B * area^2)
println("Method A (Raw): ", val_A)
println("Ratio (B*Area^2 / A): ", (val_B * area^2) / val_A)

println("F1 Method B (Raw): ", sF1)
println("F1 Method B * Area^2: ", sF1 * area^2)
println("F1 Method A (Raw): ", val_F1_A)
println("Ratio (F1 B*Area^2 / A): ", (sF1 * area^2) / val_F1_A)

function method_a_F1(tri)
    total_val = 0.0
    function integrate_recursive(v1, v2, v3, depth)
        if depth == 0
            c = (v1 + v2 + v3) / 3
            sub_area = 0.5 * norm(cross(v2 - v1, v3 - v1))
            Ig, _ = faceSingularityIgIvecg(c, tri.vertices, tri.edgel, tri.edgev̂, tri.edgen̂, tri.area, tri.facen̂, SSCg)
            # Ig seems to be negative in implementation?
            # Let's take absolute value for now or check sign.
            # 1/R integral must be positive.
            # Also, F1 Method B is ~40. F1 Method A is ~0.001.
            # This is a HUGE difference.
            # F1 Method B is singularF1(l1, l2, l3).
            # singularF1 returns -4 * (...) / 3.
            # Let's check singularF1 formula again.
            # It returns integral of 1/R.
            # For a triangle of size 0.1, 1/R ~ 10. Area ~ 0.005. Integral ~ 0.05.
            # Why is singularF1 returning 40?
            # singularF1(l1, l2, l3)
            # l1 ~ 0.14, l2 ~ 0.1, l3 ~ 0.1.
            # s = (0.34)/2 = 0.17.
            # 1/a * log(1 - a/s) = 1/0.14 * log(1 - 0.14/0.17) = 7 * log(0.03/0.17) = 7 * -1.7 = -12.
            # Sum ~ -30.
            # -4 * (-30) / 3 = 40.
            # So singularF1 returns ~40.
            # But physically, \int 1/R dS should be small.
            # Wait, singularF1 formula might be normalized by Area?
            # "Integral of 1/R over the triangle."
            # If it returns 40, and Area is 0.005.
            # Maybe it returns Integral / Area^2 ?
            # 0.05 / 0.005^2 = 0.05 / 2.5e-5 = 2000.
            # Maybe Integral / Area ?
            # 0.05 / 0.005 = 10.
            
            # Let's check `singularF1` documentation or derivation.
            # It seems `singularF1` returns something very large.
            
            # Method A calculates \int Ig dS.
            # Ig is \int 1/R dS'.
            # So Method A calculates \int \int 1/R dS' dS.
            # Ig at center is ~0.24.
            # Area is 0.005.
            # \int Ig dS ~ 0.24 * 0.005 = 0.0012.
            # This matches Method A result (0.001).
            
            # So Method A result (~0.001) is physically reasonable for \int \int 1/R.
            # Method B result (~40) is NOT reasonable for the raw integral.
            # It must be scaled.
            
            # If singularF1 returns 40, and we expect 0.001.
            # Ratio is 40000.
            # Area is 0.005. Area^2 is 2.5e-5.
            # 1/Area^2 = 40000.
            # So singularF1 returns Integral / Area^2 ?
            # Or Integral * (something).
            
            # Let's check `singularF21` scaling too.
            # Method B: 0.24.
            # Method A: 7e-6.
            # Ratio: 0.24 / 7e-6 ~ 34000.
            # Again, close to 40000 (1/Area^2).
            
            # HYPOTHESIS: singularF1 and singularF21 return values normalized by Area^2 (or similar).
            # i.e. They return I / Area^2.
            # So to get the integral, we must multiply by Area^2.
            
            # Let's test this hypothesis.
            # Multiply Method B results by Area^2.
            
            total_val += abs(Ig) * sub_area
            return
        end
        v12 = (v1 + v2) / 2
        v23 = (v2 + v3) / 2
        v31 = (v3 + v1) / 2
        integrate_recursive(v1, v12, v31, depth - 1)
        integrate_recursive(v12, v2, v23, depth - 1)
        integrate_recursive(v31, v23, v3, depth - 1)
        integrate_recursive(v12, v23, v31, depth - 1)
    end
    integrate_recursive(tri.vertices[:,1], tri.vertices[:,2], tri.vertices[:,3], 6)
    return total_val
end

val_F1_A = method_a_F1(tri)
println("F1 Method B: ", sF1)
println("F1 Method A: ", val_F1_A)
println("F1 Diff: ", abs(sF1 - val_F1_A))

