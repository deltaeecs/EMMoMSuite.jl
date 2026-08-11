# test/test_aca_operator.jl — M2 ACAOperator 集成测试（TDD）
#
# ACAOperator 复用 MLFMA 八叉树聚类与近场稀疏装配；叶层非邻盒子对经 ACA
# 压缩为低秩块。验收：MatVec 与稠密参照误差 < 1e-2，GMRES 解与稠密一致。

using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.ACA
using EMMoMSuite.FastAlgorithms.ACAOperatorModule: ACAOperator
using LinearAlgebra
using IterativeSolvers

@testset "ACAOperator" begin
    mesh = generate_sphere_mesh(0.5, 6, 10)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    λ = 299792458.0 / 300e6
    leaf = 0.25 * λ

    @testset "EFIE symmetric MatVec matches dense" begin
        efie = EFIE(300e6)
        Z = assemble_impedance_matrix(efie, basis)
        op = ACAOperator(efie, basis, leaf; tol = 1e-4, near_range = 1)

        @test size(op) == (N, N)
        @test !isempty(op.blocks)
        @test length(get_leaf_intervals(op)) >= 1

        x = randn(ComplexF64, N)
        y = op * x
        err = norm(y - Z * x) / norm(Z * x)
        @test err < 1e-2
    end

    @testset "EFIE GMRES solve matches direct" begin
        efie = EFIE(300e6)
        Z = assemble_impedance_matrix(efie, basis)
        op = ACAOperator(efie, basis, leaf; tol = 1e-4, near_range = 1)

        src = PlaneWave(300e6, 0.0, 0.0, [1.0, 0.0, 0.0])
        V = excitation_vector(efie, src, basis)
        I_aca = gmres(op, V; abstol = 1e-6, reltol = 1e-8, maxiter = 300, restart = 50)
        I_dir = Z \ V
        @test norm(I_aca - I_dir) / norm(I_dir) < 1e-2
    end

    @testset "CFIE non-symmetric MatVec matches dense" begin
        cfie = CFIE(300e6)
        Z = assemble_impedance_matrix(cfie, basis)
        op = ACAOperator(cfie, basis, leaf; tol = 1e-4, near_range = 1, symmetric = false)
        @test !isempty(op.blocks)
        x = randn(ComplexF64, N)
        y = op * x
        err = norm(y - Z * x) / norm(Z * x)
        @test err < 1e-2
    end

    @testset "preconditioner wiring" begin
        efie = EFIE(300e6)
        op = ACAOperator(efie, basis, leaf; tol = 1e-4, near_range = 1)
        x = randn(ComplexF64, N)

        P1 = ILUPreconditioner(op)
        P2 = ILUPreconditioner(op.Z_near; τ = 0.01)
        @test norm(P1 \ x - P2 \ x) / norm(P2 \ x) < 1e-8

        P3 = SPAIPreconditioner(op)
        @test size(P3.M) == (N, N)

        P4 = BlockJacobiPreconditioner(op)
        @test !isempty(P4.blocks)
        @test norm(P4 \ x) > 0
    end
end
