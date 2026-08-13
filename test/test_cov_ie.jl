# test_cov_ie.jl — IntegralEquations 边界覆盖（CFIE/FastExp/Singularities）
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using LinearAlgebra

const _SING = EMMoMSuite.IntegralEquations.EFIEModule.Singularities
const _FEXP = EMMoMSuite.IntegralEquations.VEFIEModule.FastExpModule

@testset "CFIE 装配与结构" begin
    mesh = generate_sphere_mesh(0.3, 4, 8)
    basis = RWGBasis(mesh)
    for alpha in (0.5, 0.2)
        cfie = CFIE(300e6, alpha)
        Z = assemble_impedance_matrix(cfie, basis)
        @test size(Z) == (num_basis(basis), num_basis(basis))
        @test norm(Z) > 0
    end
end

@testset "FastExp 查表" begin
    t = _FEXP.FastExpTable(2π * 300e6 / 299792458.0, 1.0; n_entries = 1000)
    R = 0.3
    v = _FEXP.fast_exp_ikr(t, R)
    @test isfinite(v)
    @test abs(abs(v) - 1.0) < 0.1  # 查表近似（模应 ≈ 1）
    @test abs(v - exp(-im * 2π * 300e6 / 299792458.0 * R)) < 1e-2
    @test isfinite(_FEXP.fast_exp_ikr(t, Float32(0.2)))
end

@testset "奇异积分函数" begin
    # 单位直角三角形的边长 (1, 1, √2), 面积 0.5 → area2 = 1.0
    a, b, c = 1.0, 1.0, sqrt(2.0)
    @test isfinite(_SING.singularF1(a, b, c))
    @test isfinite(_SING.singularF21(a, b, c, 1.0))
    @test isfinite(_SING.singularF22(a, b, c, 1.0))
    # 退化：等边三角形
    e = 1.0
    @test isfinite(_SING.singularF1(e, e, e))
    @test isfinite(_SING.singularF21(e, e, e, sqrt(3.0) / 2))
end
