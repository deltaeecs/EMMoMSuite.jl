# test/test_block_evaluator.jl — M2 块求值器测试（TDD）
#
# BlockEvaluator 预计算三角形信息与基函数映射，eval_block 按全局基函数索引
# 求值 Z[rows, cols]，与稠密装配子矩阵逐元素对照。

using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.ACA
using LinearAlgebra

@testset "BlockEvaluator" begin
    mesh = generate_sphere_mesh(0.5, 8, 12)
    basis = RWGBasis(mesh)
    N = num_basis(basis)

    rows = [1, 3, 7, 15, N]
    cols = [2, 4, 8, N - 1]
    @test all(x -> x <= N, rows) && all(x -> x <= N, cols)

    @testset "EFIE block matches dense submatrix" begin
        efie = EFIE(300e6)
        Z = assemble_impedance_matrix(efie, basis)
        ev = BlockEvaluator(efie, basis)
        B = eval_block(ev, rows, cols)
        @test size(B) == (length(rows), length(cols))
        for (ii, i) in enumerate(rows), (jj, j) in enumerate(cols)
            @test isapprox(B[ii, jj], Z[i, j]; atol = 1e-10)
        end
    end

    @testset "CFIE block matches dense submatrix" begin
        cfie = CFIE(300e6)
        Z = assemble_impedance_matrix(cfie, basis)
        ev = BlockEvaluator(cfie, basis)
        B = eval_block(ev, rows, cols)
        @test size(B) == (length(rows), length(cols))
        for (ii, i) in enumerate(rows), (jj, j) in enumerate(cols)
            @test isapprox(B[ii, jj], Z[i, j]; atol = 1e-10)
        end
    end
end
