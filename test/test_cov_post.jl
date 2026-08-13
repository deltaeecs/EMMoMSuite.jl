# test_cov_post.jl — PostProcessing 模块覆盖率补测
using Test
using EMMoMSuite
using EMMoMSuite.PostProcessing
using LinearAlgebra, Random
using StaticArrays

@testset "FarFieldPattern" begin
    freqs = [1e9]
    theta = collect(0.0:0.1:pi)
    phi = collect(0.0:0.1:2pi)
    Nθ = length(theta); Nφ = length(phi)
    Eθ = zeros(ComplexF64, 1, Nθ, Nφ)
    Eφ = zeros(ComplexF64, 1, Nθ, Nφ)
    for j in 1:Nθ, k in 1:Nφ
        Eθ[1, j, k] = sin(theta[j])
        Eφ[1, j, k] = cos(theta[j])
    end
    ff = FarFieldPattern(freqs, theta, phi, Eθ, Eφ)
    @test size(ff.E_theta) == (1, Nθ, Nφ)
    g = gain(ff)
    @test all(isfinite, g)
    @test gain_db(ff, pi / 2, 0.0) isa Real
    @test hpbw(ff) isa Real
    @test side_lobe_level(ff) isa Real
    @test all(isfinite, axial_ratio(ff))
    cc = co_cross_decompose(ff)
    @test cc isa NamedTuple
    @test all(isfinite, xpd(ff))
    # 尺寸校验
    @test_throws ArgumentError FarFieldPattern([1e9], theta, phi, zeros(ComplexF64, 1, 2, 2), Eφ)
end

@testset "MLFMA 缓存与多右端求解" begin
    c = MLFMACache(; freq = 300e6, is_valid = true, metadata = Dict("n_basis" => 10))
    @test validate_cache(c, 300e6)
    @test !validate_cache(c, 500e6)
    invalidate_cache!(c)
    @test !c.is_valid
    Random.seed!(3)
    A = randn(ComplexF64, 5, 5) + 5I
    B = randn(ComplexF64, 5, 2)
    X = solve_multi_rhs(A, B)
    @test size(X) == (5, 2)
    @test norm(A * X - B) / norm(B) < 1e-10
end

@testset "SolverResult 与近场容器" begin
    sr = SolverResult(nothing, [1.0 + 0im, 2.0 + 0im], 1e9; metadata = Dict("eq" => "EFIE"))
    @test sr.freq == 1e9
    @test length(sr.I_coeffs) == 2
    @test_throws ArgumentError SolverResult(nothing, [1.0 + 0im], -1.0)
    # NearFieldGrid / NearFieldLine
    cu = collect(0.0:0.1:0.4)
    cv = collect(0.0:0.1:0.4)
    E = fill(SVector(0.0 + 0im, 0.0 + 0im, 0.0 + 0im), length(cu), length(cv))
    H = fill(SVector(0.0 + 0im, 0.0 + 0im, 0.0 + 0im), length(cu), length(cv))
    g = NearFieldGrid(cu, cv, E, H, 1e9, Dict("kind" => "cut"))
    @test g isa NearFieldGrid
    arc = collect(0.0:0.1:pi)
    pts = [SVector(0.0, 0.0, Float64(i)) for i in eachindex(arc)]
    Ez = [SVector(0.0 + 0im, 0.0 + 0im, 0.0 + 0im) for _ in eachindex(arc)]
    l = NearFieldLine(arc, pts, Ez, copy(Ez), 1e9)
    @test l isa NearFieldLine
end
