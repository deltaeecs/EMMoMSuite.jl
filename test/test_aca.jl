# test/test_aca.jl — M1 ACA 核心模块测试（TDD）
#
# 转置约定：Z ≈ U * transpose(V)（无共轭），已数值验证收敛。
# 参考：Gibson, The Method of Moments in Electromagnetics, 3rd, Ch9 Algorithm 6。

using Test
using EMMoMSuite.FastAlgorithms.ACA
using LinearAlgebra
using Random

@testset "ACA core" begin
    Random.seed!(42)

    @testset "exact low-rank reconstruction" begin
        m, n, k0 = 60, 50, 4
        U0 = randn(ComplexF64, m, k0) .+ im .* randn(ComplexF64, m, k0) ./ 10
        V0 = randn(ComplexF64, n, k0) .+ im .* randn(ComplexF64, n, k0) ./ 10
        Z = U0 * transpose(V0)
        B = aca(Z; tol = 1e-12)
        err = norm(Z - B.U * transpose(B.V)) / norm(Z)
        @test size(B.U, 2) <= 2 * k0
        @test err < 1e-10
    end

    @testset "Green's function well-separated block" begin
        m = n = 100
        rng = Random.MersenneTwister(7)
        src = 3 .* rand(rng, 3, n)
        dst = 3 .* rand(rng, 3, m) .+ 8.0
        k = 1.0
        Z = zeros(ComplexF64, m, n)
        for i in 1:m, j in 1:n
            R = norm(dst[:, i] - src[:, j])
            Z[i, j] = exp(-im * k * R) / (4π * R)
        end
        B = aca(Z; tol = 1e-4)
        err = norm(Z - B.U * transpose(B.V)) / norm(Z)
        @test err < 5e-4
        @test size(B.U, 2) < 40
    end

    @testset "symmetric block transpose application" begin
        m = n = 40
        U0 = randn(ComplexF64, m, 3)
        V0 = randn(ComplexF64, n, 3)
        Z = U0 * transpose(V0)
        B = aca(Z; tol = 1e-12)
        x = randn(ComplexF64, m)
        y1 = B.U * (transpose(B.V) * x)
        # 转置块应用：transpose(Z) * x == V * (transpose(U) * x)
        y2 = B.V * (transpose(B.U) * x)
        @test norm(y1 - Z * x) / norm(Z * x) < 1e-10
        @test norm(y2 - transpose(Z) * x) / norm(Z * x) < 1e-10
    end

    @testset "zero block early termination" begin
        Z = zeros(ComplexF64, 30, 30)
        B = aca(Z; tol = 1e-4)
        @test size(B.U, 2) == 0
    end

    @testset "QR/SVD recompression reduces rank" begin
        m, n = 80, 70
        U0 = randn(ComplexF64, m, 5)
        V0 = randn(ComplexF64, n, 5)
        Z = U0 * transpose(V0) .+ 1e-8 .* randn(ComplexF64, m, n)
        B0 = aca(Z; tol = 1e-6, recompress = false)
        B1 = recompress!(B0; tol = 1e-6)
        @test size(B1.U, 2) <= size(B0.U, 2)
        @test norm(Z - B1.U * transpose(B1.V)) / norm(Z) < 1e-6
        rank, ratio = compression_stats(B1, m, n)
        @test rank == size(B1.U, 2)
        @test 0 <= ratio <= 1
    end
end
