using EMMoMSuite
using EMMoMSuite.IntegralEquations.EFIEModule.Singularities
using EMMoMSuite.Geometry
using LinearAlgebra
using StaticArrays

function verify_singular_near()
    println("=== Verifying faceSingularityIgIvecg (Near Field) ===")
    
    # Define Source Triangle (T_src)
    # Right isosceles in xy-plane
    v1 = [0.0, 0.0, 0.0]
    v2 = [1.0, 0.0, 0.0]
    v3 = [0.0, 1.0, 0.0]
    
    vertices = hcat(v1, v2, v3)
    
    # Edge lengths
    l1 = norm(v2 - v3)
    l2 = norm(v3 - v1)
    l3 = norm(v1 - v2)
    edgel = [l1, l2, l3]
    
    # Edge vectors (v_next - v_prev) ?
    # Need to match definition in Singularities.jl
    # Usually: e1 = v3 - v2 ? Or v2 - v3?
    # Let's check Singularities.jl or assume standard cyclic.
    # In EFIE.jl:
    # edgev[:, 1] = v3 - v2
    # edgev[:, 2] = v1 - v3
    # edgev[:, 3] = v2 - v1
    edgev = zeros(3, 3)
    edgev[:, 1] = v3 - v2
    edgev[:, 2] = v1 - v3
    edgev[:, 3] = v2 - v1
    
    # Normalize edgev
    edgev[:, 1] /= norm(edgev[:, 1])
    edgev[:, 2] /= norm(edgev[:, 2])
    edgev[:, 3] /= norm(edgev[:, 3])
    
    # Edge normals (in plane, outward?)
    # edgen = cross(edgev, normal) / length?
    # Let's compute normal first.
    normal = normalize(cross(v2 - v1, v3 - v1)) # +z
    
    edgen = zeros(3, 3)
    edgen[:, 1] = cross(edgev[:, 1], normal)
    edgen[:, 2] = cross(edgev[:, 2], normal)
    edgen[:, 3] = cross(edgev[:, 3], normal)
    # Normalize?
    edgen[:, 1] /= norm(edgen[:, 1])
    edgen[:, 2] /= norm(edgen[:, 2])
    edgen[:, 3] /= norm(edgen[:, 3])
    
    area = 0.5
    
    # Observation point (r_obs)
    # 1. Far enough (should match standard quadrature)
    r_obs_far = [0.3, 0.3, 1.0]
    
    # 2. Close (near field)
    r_obs_near = [0.3, 0.3, 0.1]
    
    # 3. Very close (near field)
    r_obs_very_near = [0.3, 0.3, 0.01]
    
    k = 2π / 1.0 # lambda = 1.0
    SSCg = compute_SSCg(k)
    
    # Test Points
    test_points = [r_obs_far, r_obs_near, r_obs_very_near]
    labels = ["Far (z=1.0)", "Near (z=0.1)", "Very Near (z=0.01)"]
    
    for (i, r_obs) in enumerate(test_points)
        println("\n--- Test Case: $(labels[i]) ---")
        println("Observation Point: $r_obs")
        
        # 1. Analytical / Semi-Analytical (faceSingularityIgIvecg)
        # We need to shift vertices relative to r_obs
        verts_shifted = vertices .- r_obs
        
        # Facen?
        # facen is vector from r_obs to plane? Or normal?
        # In Legacy: facen = dot(r_obs - v1, normal) * normal ?
        # No, facen is likely the normal vector scaled by distance?
        # Let's check how it's called in EFIERWGTri.jl
        # faceSingularityIgIvecg(rgt[:, gi], tris, abs(tris.area), tris.facen̂)
        # Wait, `tris.facen̂` is just the normal vector of the triangle.
        # But `faceSingularityIgIvecg` takes `facen`.
        # Let's check `Singularities.jl` signature.
        # faceSingularityIgIvecg(rgt, vertices, edgel, edgev, edgen, area, facen, SSCg)
        # `rgt` is the observation point? No, `rgt` is passed as first arg.
        # In Legacy: `faceSingularityIgIvecg(rgt[:, gi], ...)`
        # So first arg is observation point.
        # `vertices` should be absolute vertices?
        # Let's check `Singularities.jl` implementation.
        
        # It calculates `R` from `vertices`.
        # So `vertices` should be relative to `rgt`?
        # "rgt: Observation point"
        # "vertices: Triangle vertices"
        # Inside: `v1 = vertices[:, 1] - rgt` ?
        # Let's check the code.
        
        Ig, Ivecg = faceSingularityIgIvecg(r_obs, vertices, edgel, edgev, edgen, area, normal, SSCg)
        
        println("Analytical Ig: $Ig")
        println("Analytical Ivecg (Raw): $Ivecg")
        
        # Hypothesis: Ivecg is int (r' - r_obs) G dS'
        # So int r' G dS' = Ivecg + r_obs * Ig
        Ivecg_corrected = Ivecg + r_obs * Ig
        println("Analytical Ivecg (Corrected): $Ivecg_corrected")
        
        # 2. Numerical Reference (Adaptive Quadrature)
        # We'll use a high-order Gaussian quadrature (7 points) and subdivide if necessary?
        # For simplicity, let's use 7-point rule on the whole triangle.
        # For z=0.01, this might be inaccurate, but for z=1.0 it should match.
        
        Ig_num, Ivecg_num = numerical_integration(r_obs, vertices, k)
        
        println("Numerical Ig:  $Ig_num")
        println("Numerical Ivecg: $Ivecg_num")
        
        diff_Ig = abs(Ig - Ig_num)
        diff_Ivecg = norm(Ivecg - Ivecg_num)
        
        println("Diff Ig: $diff_Ig")
        println("Diff Ivecg: $diff_Ivecg")
        
        if i == 1 && diff_Ig > 1e-4
            println("WARNING: Far field mismatch!")
        end
    end
end

function numerical_integration(r_obs, vertices, k)
    # 7-point rule
    # Coordinates (u, v, w)
    # Weights
    # From GaussQuadrature4Geos.jl (copied values)
    a1 = 0.059715871790; b1 = 0.470142064105
    a2 = 0.797426985353; b2 = 0.101286507323
    
    coords = [
        (1/3, 1/3, 1/3),
        (a1, b1, b1), (b1, a1, b1), (b1, b1, a1),
        (a2, b2, b2), (b2, a2, b2), (b2, b2, a2)
    ]
    
    weights = [
        9/40,
        0.132394152788506, 0.132394152788506, 0.132394152788506,
        0.125939180544827, 0.125939180544827, 0.125939180544827
    ]
    
    Ig = 0.0im
    Ivecg = zeros(ComplexF64, 3)
    
    area = 0.5 # Hardcoded for this test
    
    for i in 1:7
        u, v, w = coords[i]
        weight = weights[i]
        
        r_src = u * vertices[:, 1] + v * vertices[:, 2] + w * vertices[:, 3]
        
        R_vec = r_obs - r_src
        R = norm(R_vec)
        
        # Green's function: exp(-jkR) / R
        # Note: 4pi is usually outside.
        # The singular integrals usually compute int exp(-jkR)/R dS.
        
        val = exp(-im * k * R) / R
        
        Ig += val * weight * area
        Ivecg += r_src * val * weight * area
    end
    
    return Ig, Ivecg
end

function verify_edge_integrals()
    println("\n=== Verifying Edge Integrals ===")
    # Edge from (-0.5, 0, 0) to (0.5, 0, 0)
    # Obs at (0, 0.5, 0) -> d = 0.5
    # Obs at (0, 0, 1) -> d = 1.0
    
    v1 = [-0.5, 0.0, 0.0]
    v2 = [0.5, 0.0, 0.0]
    r_obs = [0.0, 0.5, 0.0]
    
    # Parameters for formula
    # Projection of r_obs onto line containing edge
    # Line direction: (1, 0, 0)
    # r_obs projection: (0, 0, 0)
    # r0gi = (0, 0, 0)
    # dts (distance to plane) = 0 (in plane)
    # p02jvec = v1 - r0gi = (-0.5, 0, 0)
    # p02jl (dist to edge line) = 0.5?
    # Wait, p02jl is dot(p02jvec, edgen).
    # edgen is normal to edge in plane. (0, 1, 0).
    # p02jvec = (-0.5, -0.5, 0).
    # p02jl = -0.5.
    
    # Let's just use the formula directly
    # R^n integral along line segment
    
    # Numerical
    Il_r1_num = 0.0
    Il_r3_num = 0.0
    
    steps = 100
    dl = norm(v2 - v1) / steps
    for i in 1:steps
        t = (i - 0.5) / steps
        r = v1 + t * (v2 - v1)
        R = norm(r - r_obs)
        Il_r1_num += (1/R) * dl
        Il_r3_num += (R) * dl
    end
    
    println("Numerical Il_r[1] (1/R): $Il_r1_num")
    println("Numerical Il_r[3] (R):   $Il_r3_num")
    
    # Analytical
    # Using logic from Singularities.jl
    # R0^2 = 0.5^2 = 0.25
    # l_minus = -0.5
    # l_plus = 0.5
    # R_minus = sqrt(0.25 + 0.25) = sqrt(0.5)
    # R_plus = sqrt(0.25 + 0.25) = sqrt(0.5)
    
    R02 = 0.25
    lj_plus = 0.5
    lj_minus = -0.5
    R_plus = sqrt(0.5)
    R_minus = sqrt(0.5)
    
    fj = log((lj_plus + R_plus) / (lj_minus + R_minus))
    Il_r1_ana = fj
    
    Il_r3_ana = (lj_plus * R_plus - lj_minus * R_minus + 1 * R02 * Il_r1_ana) / 2
    
    println("Analytical Il_r[1]: $Il_r1_ana")
    println("Analytical Il_r[3]: $Il_r3_ana")
    
end

function verify_IS_r1()
    println("\n=== Verifying IS_r[1] (1/R Integral) ===")
    v1 = [0.0, 0.0, 0.0]
    v2 = [1.0, 0.0, 0.0]
    v3 = [0.0, 1.0, 0.0]
    r_obs = [0.3, 0.3, 1.0]
    
    # Numerical
    IS_r1_num = 0.0
    
    # 7-point rule
    a1 = 0.059715871790; b1 = 0.470142064105
    a2 = 0.797426985353; b2 = 0.101286507323
    coords = [(1/3, 1/3, 1/3), (a1, b1, b1), (b1, a1, b1), (b1, b1, a1), (a2, b2, b2), (b2, a2, b2), (b2, b2, a2)]
    weights = [9/40, 0.132394152788506, 0.132394152788506, 0.132394152788506, 0.125939180544827, 0.125939180544827, 0.125939180544827]
    area = 0.5
    
    for i in 1:7
        u, v, w = coords[i]
        r_src = u * v1 + v * v2 + w * v3
        R = norm(r_src - r_obs)
        IS_r1_num += (1/R) * weights[i] * area
    end
    println("Numerical IS_r[1]: $IS_r1_num")
    
    # Analytical (Manual implementation of logic)
    vertices = hcat(v1, v2, v3)
    edgev = zeros(3, 3)
    edgev[:, 1] = v3 - v2
    edgev[:, 2] = v1 - v3
    edgev[:, 3] = v2 - v1
    
    normal = [0.0, 0.0, 1.0]
    edgen = zeros(3, 3)
    edgen[:, 1] = cross(edgev[:, 1], normal); edgen[:, 1] /= norm(edgen[:, 1])
    edgen[:, 2] = cross(edgev[:, 2], normal); edgen[:, 2] /= norm(edgen[:, 2])
    edgen[:, 3] = cross(edgev[:, 3], normal); edgen[:, 3] /= norm(edgen[:, 3])
    
    edgel = [norm(v3-v2), norm(v1-v3), norm(v2-v1)]
    
    dts = dot(normal, r_obs - v1) # 1.0
    r0gi = r_obs - dts * normal # (0.3, 0.3, 0)
    dtsAbs = abs(dts)
    
    IS_r1_ana = 0.0
    
    edge_start_indices = [2, 3, 1]
    for edgej in 1:3
        edgeNodei_minus = vertices[:, edge_start_indices[edgej]]
        lj = edgel[edgej]
        
        lj_minus = dot(edgeNodei_minus - r0gi, edgev[:, edgej]) / norm(edgev[:, edgej]) # Normalize edgev for dot product?
        # In Singularities.jl, edgev is NOT normalized?
        # "edgev: Edge vectors"
        # Usually edge vectors are v_next - v_prev. Length is lj.
        # Code: lj_minus = dot(..., edgev)
        # If edgev has length lj, then dot product scales by lj.
        # But lj_minus should be length.
        # Let's check Singularities.jl again.
        
        # "lj_minus = dot(edgeNodei_minus - r0gi, edgev[:, edgej])"
        # If edgev is not normalized, this is Length * Length * cos(theta).
        # This is WRONG if edgev is not normalized.
        
        # Let's check if edgev is normalized in Singularities.jl or caller.
        # In verify_singular_near.jl:
        # edgev[:, 1] = v3 - v2. (Length sqrt(2)).
        # So edgev is NOT normalized.
        
        # If Singularities.jl expects normalized edgev, then my input is wrong.
        # If Singularities.jl expects unnormalized edgev, then the dot product is wrong.
        
        # Let's check `faceSingularityIgIvecg` implementation.
        # `lj_minus = dot(edgeNodei_minus - r0gi, edgev[:, edgej])`
        # This looks suspicious.
        
        # Let's check `p02jvec`.
        # `p02jvec = edgeNodei_minus - lj_minus * edgev[:, edgej] - r0gi`
        # If `lj_minus` is length along edge, then `lj_minus * edgev` (where edgev has length L) would be L^2.
        # This confirms `edgev` MUST be normalized for this logic to work.
        
        # BUT `edgev` passed to `faceSingularityIgIvecg` usually comes from `TriangleInfo`.
        # In `TriangleInfo`, `edgev` are usually the edge vectors (v_next - v_prev).
        
        # Let's check `MoM_Kernels` or `EMMoMSuite` `TriangleInfo` definition.
    end
end

verify_IS_r1()
verify_edge_integrals()
verify_singular_near()
