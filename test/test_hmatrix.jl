# test/test_hmatrix.jl - H1a HMatrix tree + builder from MLACA (TDD)
#
# hmatrix_from_mlaca 将 MLACA 算子的分层低秩结构重建为显式 H 矩阵树
# （:dense / :lowrank / :split），materialize 稠密化后应与算子的稠密参照一致。

using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.ACA
using EMMoMSuite.FastAlgorithms.MLACAOperatorModule: MLACAOperator
using EMMoMSuite.FastAlgorithms.ACAOperatorModule: ACAOperator
using EMMoMSuite.FastAlgorithms.ACAOperatorModule: ACAOperator as ACAOp
using EMMoMSuite.FastAlgorithms.BlockLUModule: extract_block
using EMMoMSuite.FastAlgorithms.HMatrixModule: hmatrix_from_mlaca, materialize
using EMMoMSuite.FastAlgorithms.HMatrixModule: h_lu!, h_lu_solve
using LinearAlgebra

function operator_dense(op, N)
    # 由近场 + 低秩块重建算子的稠密参照（叶层块对）
    leaf = op.octree.levels[op.octree.nLevels]
    idxs = [collect(op.sorted_ids[c.bfInterval]) for c in leaf.cubes]
    filter!(!isempty, idxs)
    Z = zeros(ComplexF64, N, N)
    for a in eachindex(idxs), b in eachindex(idxs)
        Z[idxs[a], idxs[b]] = extract_block(op, idxs[a], idxs[b])
    end
    return Z
end

function compare_hmatrix(H, Zop)
    # Zop 为全局索引稠密参照；H 的行/列为 DFS 顺序，需按全局索引对齐
    Zm = materialize(H)
    ref = Zop[H.rows, H.cols]
    return norm(Zm - ref) / norm(ref)
end

@testset "HMatrix" begin
    mesh = generate_sphere_mesh(0.5, 6, 10)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    efie = EFIE(300e6)
    lambda = 299792458.0 / 300e6

    @testset "build from MLACA multi-level reconstructs operator" begin
        op = MLACAOperator(efie, basis, 0.125 * lambda; tol = 1e-4, near_range = 1)
        H = hmatrix_from_mlaca(op)
        @test H.kind == :split
        Zop = operator_dense(op, N)
        @test size(materialize(H)) == (N, N)
        @test compare_hmatrix(H, Zop) < 1e-2
    end

    @testset "build from ACA single-level reconstructs operator" begin
        op = ACAOperator(efie, basis, 0.25 * lambda; tol = 1e-4, near_range = 1)
        H = hmatrix_from_mlaca(op)
        Zop = operator_dense(op, N)
        @test compare_hmatrix(H, Zop) < 1e-2
    end

    @testset "H-LU EFIE exact factorization + multi-RHS solve" begin
        op = MLACAOperator(efie, basis, 0.125 * lambda; tol = 1e-4, near_range = 1)
        H = hmatrix_from_mlaca(op)
        h_lu!(H)
        X_true = randn(ComplexF64, N, 2)
        B = reduce(hcat, [op * X_true[:, j] for j in 1:2])
        X = h_lu_solve(H, B)
        @test norm(X - X_true) / norm(X_true) < 1e-6

        b = randn(ComplexF64, N)
        x = H \ b
        @test norm(op * x - b) / norm(b) < 1e-6
    end

    @testset "H-LU ACA single-level solve" begin
        op = ACAOperator(efie, basis, 0.25 * lambda; tol = 1e-4, near_range = 1)
        H = hmatrix_from_mlaca(op)
        h_lu!(H)
        X_true = randn(ComplexF64, N, 2)
        B = reduce(hcat, [op * X_true[:, j] for j in 1:2])
        X = h_lu_solve(H, B)
        @test norm(X - X_true) / norm(X_true) < 1e-6
    end

    @testset "H-LU CFIE (non-symmetric) solve" begin
        cfie = CFIE(300e6)
        op = ACAOperator(cfie, basis, 0.25 * lambda; tol = 1e-4, near_range = 1, symmetric = false)
        H = hmatrix_from_mlaca(op)
        h_lu!(H)
        X_true = randn(ComplexF64, N, 2)
        B = reduce(hcat, [op * X_true[:, j] for j in 1:2])
        X = h_lu_solve(H, B)
        @test norm(X - X_true) / norm(X_true) < 1e-6
    end

    @testset "H-LU PMCHW (2N system) solve" begin
        pmchw = PMCHW(300e6, 4.0)
        op = ACAOperator(pmchw, basis, 0.25 * lambda; tol = 1e-4, near_range = 1)
        H = hmatrix_from_mlaca(op)
        @test length(H.rows) == 2N
        h_lu!(H)
        X_true = randn(ComplexF64, 2N, 2)
        B = reduce(hcat, [op * X_true[:, j] for j in 1:2])
        X = h_lu_solve(H, B)
        @test norm(X - X_true) / norm(X_true) < 1e-6
    end

    @testset "H-LU low-frequency EFIE solve" begin
        lf = EFIE(30e6)
        lambda_lf = 299792458.0 / 30e6
        op = ACAOperator(lf, basis, 0.03 * lambda_lf; tol = 1e-4, near_range = 1)
        H = hmatrix_from_mlaca(op)
        h_lu!(H)
        X_true = randn(ComplexF64, N, 2)
        B = reduce(hcat, [op * X_true[:, j] for j in 1:2])
        X = h_lu_solve(H, B)
        @test norm(X - X_true) / norm(X_true) < 1e-6
    end
end
