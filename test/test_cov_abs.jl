# test_cov_abs.jl — Absorption（SAR/吸收功率）轻量覆盖
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.PostProcessing
using LinearAlgebra, Random

@testset "吸收功率与 SAR" begin
    vol = generate_box_tet_mesh(0.4, 0.4, 0.4, 1, 1, 1)
    basis = PWCBasis(vol)
    perms = fill(ComplexF64(2.0 - 0.05im), num_elements(vol))
    Random.seed!(4)
    I = randn(ComplexF64, num_basis(basis))
    r = absorbed_power(basis, I, perms)
    @test r.P_total >= 0
    @test length(r.P_density) == num_elements(vol)
    @test all(isfinite, r.P_density)
    rho = fill(1000.0, num_elements(vol))
    s = sar(basis, I, perms, rho)
    @test s.SAR_total >= 0
    @test length(s.SAR_per_element) == num_elements(vol)
    @test all(isfinite, s.SAR_per_element)
end
