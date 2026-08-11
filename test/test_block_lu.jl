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
end
