# test_cov_post2.jl — PostProcessing 电流/体积电流覆盖
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.PostProcessing
using LinearAlgebra, Random

@testset "表面/体积电流" begin
    mesh = generate_sphere_mesh(0.3, 4, 8)
    basis = RWGBasis(mesh)
    Random.seed!(6)
    I = randn(ComplexF64, num_basis(basis))
    J = geoElectricJCal(I, basis)
    @test size(J) == (3, num_elements(mesh))
    @test all(isfinite, J)

    vol = generate_box_tet_mesh(0.4, 0.4, 0.4, 1, 1, 1)
    perms = fill(ComplexF64(2.0), num_elements(vol))
    pwc = PWCBasis(vol)
    Ip = randn(ComplexF64, num_basis(pwc))
    Jv = geoVolumeCurrentCal(Ip, pwc, perms)
    @test length(Jv) == num_elements(vol)
    @test all(v -> length(v) == 3, Jv)
    @test all(v -> all(isfinite, v), Jv)
end
