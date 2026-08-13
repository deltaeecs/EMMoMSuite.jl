# test_cov_fa.jl — FastAlgorithms 结构层覆盖（Octree/Level/Lebedev LVI 路径）
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.FastAlgorithms
using LinearAlgebra

@testset "Octree/Level 结构路径" begin
    mesh = generate_sphere_mesh(0.3, 4, 8)
    basis = RWGBasis(mesh)
    centers = reduce(hcat, [bf.center for bf in basis.functions])
    for meth in (Val(:Lagrange2Step), Val(:FFTSpectral), Val(:LbTrained1Step))
        oct, sorted_ids = build_octree(centers, 0.15; λ = 1.0, interp_method = meth)
        @test oct isa OctreeInfo
        @test length(sorted_ids) == num_basis(basis)
        leaf = oct.levels[oct.nLevels]
        @test leaf isa LevelInfo
        @test leaf.cubeEdgel > 0
        @test length(leaf.cubes) > 0
        @test leaf.nCubes == length(leaf.cubes)
        @test length(leaf.poles.r̂sθsϕs) > 0
    end
end
