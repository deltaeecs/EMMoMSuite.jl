# test_distributed_gmres.jl
# Phase 9 测试覆盖率补全：验证 DistributedGMRES.jl 实现
# 覆盖目标:
#   - row_partition：边界条件（N 不能整除 nproc）
#   - distributed_gmres!：log=true 历史记录，verbose=true 路径
#   - mpi_gmres!：收敛历史 API
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using LinearAlgebra
using MPI

# ── row_partition 单元测试（内部函数，通过 Parallel 模块访问）─────────────────
@testset "row_partition" begin
    # 内部函数：通过 EMMoMSuite.Parallel 访问
    rp = EMMoMSuite.Parallel.row_partition

    # N=10, nproc=3：分区应为 [1:4], [5:7], [8:10]
    r0 = rp(10, 0, 3);  @test r0 == 1:4
    r1 = rp(10, 1, 3);  @test r1 == 5:7
    r2 = rp(10, 2, 3);  @test r2 == 8:10
    # 总行数必须等于 N
    @test length(r0) + length(r1) + length(r2) == 10

    # N=6, nproc=3：均匀分配
    s0 = rp(6, 0, 3);  @test length(s0) == 2
    s1 = rp(6, 1, 3);  @test length(s1) == 2
    s2 = rp(6, 2, 3);  @test length(s2) == 2
    @test length(s0) + length(s1) + length(s2) == 6

    # N=1, nproc=1
    t0 = rp(1, 0, 1);  @test t0 == 1:1

    # nproc > N：有些进程分配 0 行
    u0 = rp(3, 0, 5);  @test length(u0) == 1
    u4 = rp(3, 4, 5);  @test length(u4) == 0  # 最后两个进程无行

    # 验证所有分区不重叠且连续覆盖 1:N
    N_test = 17; P = 4
    ranges = [rp(N_test, r, P) for r in 0:P-1]
    all_idcs = vcat([collect(r) for r in ranges]...)
    sort!(all_idcs)
    @test all_idcs == collect(1:N_test)
end

# ── mpi_gmres! 收敛历史 (log=true) ──────────────────────────────────────────
@testset "mpi_gmres! with log history" begin
    if !MPI.Initialized()
        MPI.Init()
    end

    radius = 0.3; freq = 300e6
    mesh   = generate_sphere_mesh(radius, 5, 10)
    basis  = RWGBasis(mesh)
    N      = num_basis(basis)

    efie = EFIE(freq)
    Z    = assemble_impedance_matrix_parallel(efie, basis)

    import Random: seed!
    seed!(7)
    x_ref = randn(ComplexF64, N)
    b     = zeros(ComplexF64, N);  mul!(b, Z, x_ref)

    x0 = zeros(ComplexF64, N)
    x_sol, history = mpi_gmres!(x0, Z, b;
        restart  = min(30, N),
        maxiter  = 3 * N,
        reltol   = 1e-8,
        log      = true,
    )

    # history 不为 nothing
    @test history !== nothing
    # 最终残差满足收敛容差
    r = zeros(ComplexF64, N);  mul!(r, Z, x_sol);  r .-= b
    @test norm(r) / norm(b) < 1e-6
end

# ── mpi_gmres! verbose 路径 ──────────────────────────────────────────────────
@testset "mpi_gmres! verbose path" begin
    if !MPI.Initialized()
        MPI.Init()
    end

    radius = 0.25; freq = 300e6
    mesh   = generate_sphere_mesh(radius, 4, 8)
    basis  = RWGBasis(mesh)
    N      = num_basis(basis)

    efie = EFIE(freq)
    Z    = assemble_impedance_matrix_parallel(efie, basis)

    import Random: seed!
    seed!(99)
    x_ref = randn(ComplexF64, N)
    b     = zeros(ComplexF64, N);  mul!(b, Z, x_ref)

    # verbose=true 应打印迭代信息（不抛异常）
    x_sol = mpi_gmres(Z, b;
        restart  = min(20, N),
        maxiter  = 2 * N,
        reltol   = 1e-6,
        verbose  = true,
    )

    r = zeros(ComplexF64, N);  mul!(r, Z, x_sol);  r .-= b
    @test norm(r) / norm(b) < 1e-4
end
