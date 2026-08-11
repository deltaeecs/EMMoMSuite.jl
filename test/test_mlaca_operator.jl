# test/test_mlaca_operator.jl — M3 MLACAOperator 测试（TDD）
#
# MLACA 在八叉树多层结构上做 H-矩阵风格递归块压缩：可容许远块 ACA 压缩、
# 近块下钻、叶层邻对由 Z_near 覆盖。验收：存在跨多个叶层盒子的多层块，
# MatVec 与 GMRES 解对照稠密误差 < 1e-2。

using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.ACA
using EMMoMSuite.FastAlgorithms.MLACAOperatorModule: MLACAOperator
using LinearAlgebra
using IterativeSolvers

@testset "MLACAOperator" begin
    mesh = generate_sphere_mesh(0.5, 8, 12)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    λ = 299792458.0 / 300e6

    efie = EFIE(300e6)
    Z = assemble_impedance_matrix(efie, basis)

    # 叶层 0.125λ → nLevels=4，level3 出现非邻对 → 真正的多层压缩
    op = MLACAOperator(efie, basis, 0.125 * λ; tol = 1e-4, near_range = 1)

    @test size(op) == (N, N)
    @test !isempty(op.blocks)
    @test length(get_leaf_intervals(op)) >= 1

    # 多层证据：至少一个低秩块的行覆盖超过单个叶层盒子的基函数数
    leaf_level = op.octree.levels[op.octree.nLevels]
    max_leaf_bs = maximum(length, [c.bfInterval for c in leaf_level.cubes])
    @test any(length(blk.rows) > max_leaf_bs for blk in op.blocks)

    x = randn(ComplexF64, N)
    y = op * x
    err = norm(y - Z * x) / norm(Z * x)
    @test err < 1e-2

    src = PlaneWave(300e6, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(efie, src, basis)
    I_mlaca = gmres(op, V; abstol = 1e-6, reltol = 1e-8, maxiter = 300, restart = 50)
    I_dir = Z \ V
    @test norm(I_mlaca - I_dir) / norm(I_dir) < 1e-2
end

@testset "MLACAOperator non-symmetric (CFIE)" begin
    mesh = generate_sphere_mesh(0.5, 8, 12)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    λ = 299792458.0 / 300e6

    cfie = CFIE(300e6)
    Z = assemble_impedance_matrix(cfie, basis)
    op = MLACAOperator(cfie, basis, 0.125 * λ; tol = 1e-4, near_range = 1, symmetric = false)

    @test size(op) == (N, N)
    @test !isempty(op.blocks)

    x = randn(ComplexF64, N)
    y = op * x
    err = norm(y - Z * x) / norm(Z * x)
    @test err < 1e-2

    src = PlaneWave(300e6, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(cfie, src, basis)
    I_mlaca = gmres(op, V; abstol = 1e-6, reltol = 1e-8, maxiter = 300, restart = 50)
    I_dir = Z \ V
    @test norm(I_mlaca - I_dir) / norm(I_dir) < 1e-2
end
