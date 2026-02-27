
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "LegacyBenchmark", "packages"))

using EMSuite
using EMSuite.IntegralEquations.EFIEModule.Singularities: singularF1, singularF21, singularF22, faceSingularityIgIvecg, compute_SSCg
using MoM_Kernels
using MoM_Basics
using StaticArrays

function debug_singularities()
    println("Debugging Singularities...")
    
    # Define a triangle (Right isosceles)
    # (0,0,0), (1,0,0), (0,1,0)
    vertices = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0] # 3x3 Matrix
    # Edges: 1, sqrt(2), 1
    a = 1.0
    b = sqrt(2.0)
    c = 1.0
    area = 0.5
    area2 = area^2
    
    # Edge lengths vector
    edgel = [a, b, c]
    
    # Edge vectors (normalized?)
    # MoM_Basics uses specific definitions.
    # Let's just use dummy values for now, or try to construct a proper triangle info.
    # But faceSingularityIgIvecg takes raw arrays.
    
    # We need edgev (vectors along edges) and edgen (normals to edges in plane)
    # This is complicated to setup manually.
    
    # Let's just check if the function exists and runs.
    # And compare output for dummy inputs.
    
    r = [0.5, 0.5, 0.1]
    edgev = zeros(3, 3)
    edgen = zeros(3, 3)
    facen = [0.0, 0.0, 1.0]
    
    # SSCg
    k = 6.28
    SSCg_em = compute_SSCg(k)
    # Legacy SSCg?
    # MoM_Kernels.Params.SSCg?
    # We can compute it using coeffgreen.
    
    println("Triangle: a=$a, b=$b, c=$c, area=$area, area2=$area2")
    
    # EMSuite
    sF1_em = singularF1(a, b, c)
    sF21_em = singularF21(a, b, c, area2)
    sF22_em = singularF22(a, b, c, area2)
    
    println("\n--- EMSuite ---")
    println("sF1: $sF1_em")
    println("sF21: $sF21_em")
    println("sF22: $sF22_em")
    
    # Legacy
    sF1_leg = MoM_Kernels.singularF1(a, b, c)
    sF21_leg = MoM_Kernels.singularF21(a, b, c, area2)
    sF22_leg = MoM_Kernels.singularF22(a, b, c, area2)
    
    println("\n--- Legacy ---")
    println("sF1: $sF1_leg")
    println("sF21: $sF21_leg")
    println("sF22: $sF22_leg")
    
    println("\n--- Comparison ---")
    println("Diff sF1: $(sF1_em - sF1_leg)")
    println("Diff sF21: $(sF21_em - sF21_leg)")
    println("Diff sF22: $(sF22_em - sF22_leg)")
    
    # Check faceSingularityIgIvecg
    println("\n--- Checking faceSingularityIgIvecg ---")
    # We need proper inputs.
    # r is point.
    # vertices is 3x3 matrix.
    # edgel is vector.
    # edgev is 3x3 matrix (vectors along edges).
    # edgen is 3x3 matrix (normals to edges).
    # area is scalar.
    # facen is vector (normal to face).
    # SSCg is vector.
    
    # Construct dummy inputs
    edgev = [1.0 0.0 0.0; -0.707 0.707 0.0; 0.0 -1.0 0.0]' # Transpose to match column vectors?
    # edgev columns are vectors along edges.
    # edge 1: (1,0,0) - (0,1,0) = (1,-1,0). Length sqrt(2).
    # Wait, vertices are columns?
    # vertices = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
    # v1=(0,0,0), v2=(1,0,0), v3=(0,1,0).
    # e1 = v2-v1 = (1,0,0).
    # e2 = v3-v2 = (-1,1,0).
    # e3 = v1-v3 = (0,-1,0).
    
    vertices = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
    edgev = [1.0 -0.707 0.0; 0.0 0.707 -1.0; 0.0 0.0 0.0] # Dummy
    edgen = [0.0 0.707 1.0; 1.0 0.707 0.0; 0.0 0.0 0.0] # Dummy
    
    # Just call it with zeros to see if it crashes or returns same
    edgev = zeros(3, 3)
    edgen = zeros(3, 3)
    
    Ig_em, Ivecg_em = faceSingularityIgIvecg(r, vertices, edgel, edgev, edgen, area, facen, SSCg_em)
    
    # Legacy
    # We need Legacy SSCg.
    # MoM_Kernels.Params.SSCg is not available directly?
    # We can use compute_SSCg from EMSuite (it matches Legacy formula).
    SSCg_leg = SSCg_em
    
    # MoM_Kernels.faceSingularityIgIvecg is not exported.
    # It is in ZmatAndVvec/Singularity/FaceSingularity.jl
    # We can access it via MoM_Kernels.faceSingularityIgIvecg if included.
    # It is included in ZmatVvec.jl.
    
    Ig_leg, Ivecg_leg = MoM_Kernels.faceSingularityIgIvecg(r, vertices, edgel, edgev, edgen, area, facen, SSCg_leg)
    
    println("Ig EM: $Ig_em")
    println("Ig Leg: $Ig_leg")
    println("Diff Ig: $(Ig_em - Ig_leg)")
    
    println("Ivecg EM: $Ivecg_em")
    println("Ivecg Leg: $Ivecg_leg")
    println("Diff Ivecg: $(Ivecg_em - Ivecg_leg)")
    
end

debug_singularities()
