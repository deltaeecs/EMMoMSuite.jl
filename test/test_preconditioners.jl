# test_preconditioners.jl
# 覆盖率目标: Preconditioners.jl — 所有预条件器类型的构造与应用

using Test
using EMMoMSuite
using LinearAlgebra
using SparseArrays

# Preconditioners 通过 EMMoMSuite.Solvers 导出
using EMMoMSuite.Solvers:
    IdentityPreconditioner,
    DiagonalPreconditioner,
    ILUPreconditioner,
    SPAIPreconditioner,
    BlockJacobiPreconditioner
using EMMoMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator

@testset "Preconditioners" begin
    # 构造简单测试矩阵：5×5 对角占优复数稀疏矩阵
    n = 8
    A_dense = (10.0 + 0im) * I + 0.5 * (rand(ComplexF64, n, n) .- (0.5 + 0.5im))
    # 确保对角占优
    for i = 1:n
        A_dense[i, i] = 20.0 + 1im
    end
    A_sparse = sparse(A_dense)

    b = rand(ComplexF64, n)
    y = similar(b)

    # ── Identity ──
    @testset "IdentityPreconditioner" begin
        P = IdentityPreconditioner()
        @test (P \ b) ≈ b
        ldiv!(y, P, b)
        @test y ≈ b
        ldiv!(P, y)
        @test y ≈ b
    end

    # ── Diagonal ──
    @testset "DiagonalPreconditioner" begin
        P = DiagonalPreconditioner(A_dense)
        x = P \ b
        @test length(x) == n
        @test all(isfinite, x)
        ldiv!(y, P, b)
        @test y ≈ x
        ldiv!(P, y)     # in-place version
        @test all(isfinite, y)
    end

    # ── ILU ──
    @testset "ILUPreconditioner" begin
        # ILU on real sparse matrix (IncompleteLU expects real)
        A_real = sparse(20.0 * I + rand(n, n))
        P = ILUPreconditioner(A_real; τ = 0.01)
        b_real = rand(n)
        y_real = similar(b_real)
        ldiv!(y_real, P, b_real)
        @test all(isfinite, y_real)
        ldiv!(P, y_real)
        @test all(isfinite, y_real)
        x_bs = P \ b_real
        @test all(isfinite, x_bs)
    end

    # ── SPAI ──
    @testset "SPAIPreconditioner" begin
        A_real_sp = sparse(20.0 * I + 0.5 * rand(n, n))
        P = SPAIPreconditioner(A_real_sp)
        b_real = rand(n)
        y_real = similar(b_real)
        ldiv!(y_real, P, b_real)
        @test all(isfinite, y_real)
        ldiv!(P, copy(b_real))
        x_bs = P \ b_real
        @test all(isfinite, x_bs)
    end

    # ── BlockJacobi ──
    @testset "BlockJacobiPreconditioner" begin
        # BlockJacobi requires SparseMatrixCSC{Real} + block_intervals
        A_real_sp = sparse(Diagonal(fill(20.0, n)) + 0.3 * sprand(n, n, 0.3))
        # 4 blocks of size 2 each
        block_intervals = [UnitRange(1,2), UnitRange(3,4), UnitRange(5,6), UnitRange(7,8)]
        P = BlockJacobiPreconditioner(A_real_sp, block_intervals)
        b_real = rand(n)
        y_real = similar(b_real)
        ldiv!(y_real, P, b_real)
        @test all(isfinite, y_real)
        ldiv!(P, copy(b_real))
        x_bs = P \ b_real
        @test all(isfinite, x_bs)
    end

    @testset "BlockJacobiPreconditioner arbitrary index blocks" begin
        A_real_sp = sparse(Diagonal(fill(20.0, n)) + 0.3 * sprand(n, n, 0.3))
        block_indices = [[1, 3], [2, 4], [5, 7], [6, 8]]
        P = BlockJacobiPreconditioner(A_real_sp, block_indices)
        b_real = rand(n)
        y_real = similar(b_real)
        ldiv!(y_real, P, b_real)
        @test all(isfinite, y_real)
        x_bs = P \ b_real
        @test all(isfinite, x_bs)
    end

    @testset "BlockJacobiPreconditioner PMCHW MLFMA constructor" begin
        mesh = generate_sphere_mesh(0.5, 3, 6)
        basis = RWGBasis(mesh)
        pmchw = PMCHW(300e6, 4.0)
        op = PMCHWMLFMAOperator(pmchw, basis, 0.10)

        P = BlockJacobiPreconditioner(op)
        P_compat = BlockJacobiPreconditioner(op, basis)
        x = randn(ComplexF64, 2 * num_basis(basis))
        y = similar(x)
        ldiv!(y, P, x)
        y_compat = P_compat \ x

        @test length(y) == length(x)
        @test y_compat ≈ (P \ x)
        @test all(isfinite, real.(y))
        @test all(isfinite, imag.(y))
    end
end
