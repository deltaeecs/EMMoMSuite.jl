
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
    
end

debug_singularities()
