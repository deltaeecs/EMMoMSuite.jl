# test/test_block_lu.jl - F3 ACA/MLACA direct block LU (TDD)
#
# Block LU (Gibson Ch9 Algorithm 7): leaf-cube block partition, factor once,
# then multi-RHS direct solves.
# Gates: reconstructed L*U residual < 1e-2; multi-RHS and single-RHS solves
# match the dense LU within 1e-2.

using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.ACA
using EMMoMSuite.FastAlgorithms.ACAOperatorModule: ACAOperator
using EMMoMSuite.FastAlgorithms.MLACAOperatorModule: MLACAOperator
using EMMoMSuite.FastAlgorithms.BlockLUModule: block_lu, block_lu_solve
using LinearAlgebra

function densify_block(X)
    X isa LowRankBlock && return X.U * transpose(X.V)
    return Matrix(X)
end

function reconstruct_lu(F::BlockLUFactorization, N::Int)
    M = length(F.blocks)
    Zr = zeros(ComplexF64, N, N)
    for b in 1:M, s in 1:M
        Ib = F.blocks[b]
        Is = F.blocks[s]
        blk = zeros(ComplexF64, length(Ib), length(Is))
        for p in 1:min(b, s)
            Lbp = if b == p
                Matrix(F.diag[b].L)
            else
                densify_block(F.L[b][p])
            end
            Ups = if p == s
                Matrix(F.diag[s].U)
            else
                densify_block(F.U[p][s])
            end
            blk .+= Lbp * Ups
        end
        Zr[Ib, Is] = blk
    end
    return Zr
end

@testset "Block LU" begin
    mesh = generate_sphere_mesh(0.5, 6, 10)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    lambda = 299792458.0 / 300e6
    efie = EFIE(300e6)
    Z = assemble_impedance_matrix(efie, basis)

    @testset "ACA EFIE block LU factorization + multi-RHS solve" begin
        op = ACAOperator(efie, basis, 0.25 * lambda; tol = 1e-4, near_range = 1)
        F = block_lu(op; tol = 1e-4)
        Zr = reconstruct_lu(F, N)
        @test norm(Zr - Z) / norm(Z) < 1e-2

        B = randn(ComplexF64, N, 3)
        X = block_lu_solve(F, B)
        Xd = Z \ B
        @test norm(X - Xd) / norm(Xd) < 1e-2

        b = randn(ComplexF64, N)
        x = F \ b
        @test norm(x - Z \ b) / norm(Z \ b) < 1e-2
    end

    @testset "MLACA EFIE block LU multi-RHS solve" begin
        op = MLACAOperator(efie, basis, 0.125 * lambda; tol = 1e-4, near_range = 1)
        F = block_lu(op; tol = 1e-4)
        B = randn(ComplexF64, N, 2)
        X = block_lu_solve(F, B)
        Xd = Z \ B
        @test norm(X - Xd) / norm(Xd) < 1e-2
    end

    @testset "PMCHW block LU (2N system) multi-RHS solve" begin
        pmchw = PMCHW(300e6, 4.0)
        op = ACAOperator(pmchw, basis, 0.25 * lambda; tol = 1e-4, near_range = 1)
        S = 2N
        @test size(op, 1) == S

        # 分解+求解应精确反演压缩算子：X_true → B = op*X_true → block_lu_solve ≈ X_true
        X_true = randn(ComplexF64, S, 2)
        B = reduce(hcat, [op * X_true[:, j] for j in 1:2])
        # PMCHW 病态（cond≈4e6）：默认精确分解（recompress=false）须精确反演；
        # 再压缩在病态系统下不稳定（因子误差被块更新复合放大），不作为门控。
        F = block_lu(op; tol = 1e-4, recompress = false)
        @test sum(length, F.blocks) == S
        X = block_lu_solve(F, B)
        @test norm(X - X_true) / norm(X_true) < 1e-8
    end

    @testset "EFIE recompressed block LU (well-conditioned opt-in)" begin
        op = ACAOperator(efie, basis, 0.25 * lambda; tol = 1e-4, near_range = 1)
        F = block_lu(op; tol = 1e-4, recompress = true)
        B = randn(ComplexF64, N, 2)
        X = block_lu_solve(F, B)
        Xd = Z \ B
        @test norm(X - Xd) / norm(Xd) < 1e-2
    end

    @testset "CFIE block LU (non-symmetric) multi-RHS solve" begin
        cfie = CFIE(300e6)
        op = ACAOperator(cfie, basis, 0.25 * lambda; tol = 1e-4, near_range = 1, symmetric = false)
        Zc = assemble_impedance_matrix(cfie, basis)
        F = block_lu(op; tol = 1e-4)
        B = randn(ComplexF64, N, 2)
        X = block_lu_solve(F, B)
        Xd = Zc \ B
        @test norm(X - Xd) / norm(Xd) < 1e-2
    end
end
