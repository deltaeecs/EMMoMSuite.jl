# test_cov_solvers.jl — Solvers 模块边界覆盖（初值/预条件/verbose/SPAI/BlockJacobi）
using Test
using EMMoMSuite
using EMMoMSuite.Solvers
using LinearAlgebra, SparseArrays, Random

@testset "求解器参数边界" begin
    Random.seed!(21)
    A = randn(ComplexF64, 8, 8) + 8I
    b = randn(ComplexF64, 8)
    # GMRES 带初值/预条件/verbose
    x0 = ones(ComplexF64, 8)
    x = solve!(
        GMRESSolver(restart = 4, maxiter = 50, tol = 1e-8, verbose = true),
        A, b, x0; Pl = Identity(), Pr = Identity(),
    )
    @test norm(A * x - b) / norm(b) < 1e-6
    # BiCGSTAB 带初值
    xb = solve!(BiCGSTABSolver(maxiter = 50, tol = 1e-8), A, b, x0)
    @test norm(A * xb - b) / norm(b) < 1e-6
    # LU 带初值（忽略）
    xl = solve!(LUSolver(), A, b, x0)
    @test norm(A * xl - b) / norm(b) < 1e-10
    # 实矩阵
    Ar = randn(8, 8) + 8I
    br = randn(8)
    xr = solve!(GMRESSolver(tol = 1e-10), Ar, br)
    @test norm(Ar * xr - br) / norm(br) < 1e-8
end

@testset "预条件构造与应用" begin
    n = 12
    Random.seed!(22)
    A_full = (20.0 + 0im) * I + 0.5 * randn(ComplexF64, n, n)
    A_sp = sparse(A_full)
    # ILU
    pilu = ILUPreconditioner(A_sp; τ = 0.01)
    x = randn(ComplexF64, n)
    y = pilu \ x
    @test all(isfinite, y)
    # SPAI（实矩阵）
    Ar_sp = sparse(20.0 * I + 0.5 * rand(n, n))
    psp = SPAIPreconditioner(Ar_sp)
    xr = rand(n)
    yr = psp \ xr
    @test all(isfinite, yr)
    # BlockJacobi：非连续块索引
    block_indices = [collect(1:4), collect(6:9), collect(10:12)]
    pbj = BlockJacobiPreconditioner(A_sp, block_indices)
    yb = pbj \ x
    @test all(isfinite, yb)
    # 空块过滤
    pbj2 = BlockJacobiPreconditioner(A_sp, [Int[], collect(1:3)])
    @test length(pbj2.blocks) == 1
end
