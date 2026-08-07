# test_parallel_mfie_cfie.jl
# Phase 9 测试覆盖率补全：验证 Assembly.jl 中 MFIE 和 CFIE 并行装配路径
# 覆盖目标: _fill_local!(MFIE)、_fill_local!(CFIE) 两条 dispatch 分支
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.Parallel
using LinearAlgebra
using MPI

@testset "Parallel MFIE Assembly (single-process)" begin
    if !MPI.Initialized()
        MPI.Init()
    end

    radius = 0.3; freq = 300e6
    mesh  = generate_sphere_mesh(radius, 5, 10)
    basis = RWGBasis(mesh)
    N     = num_basis(basis)
    @test N > 0

    mfie = MFIE(freq)

    # 并行 MFIE 装配 (singe-process, P=1)
    Z_mpi = assemble_impedance_matrix_parallel(mfie, basis)
    Z_gathered = gather(Z_mpi)

    if Z_gathered !== nothing
        # 串行 MFIE 装配
        Z_serial = assemble_impedance_matrix(mfie, basis)

        # 矩阵对角元应与串行结果匹配（放宽至 1e-6 相对误差）
        @test size(Z_gathered) == (N, N)
        rel_err = norm(Z_gathered .- Z_serial) / norm(Z_serial)
        @test rel_err < 1e-6
    else
        # 非 rank-0 进程：仅验证本地列分区正确
        @test size(Z_mpi) == (N, N)
    end
end

@testset "Parallel CFIE Assembly (single-process)" begin
    if !MPI.Initialized()
        MPI.Init()
    end

    radius = 0.3; freq = 300e6
    mesh  = generate_sphere_mesh(radius, 5, 10)
    basis = RWGBasis(mesh)
    N     = num_basis(basis)
    @test N > 0

    cfie = CFIE(freq, 0.5)

    # 并行 CFIE 装配
    Z_mpi = assemble_impedance_matrix_parallel(cfie, basis)
    Z_gathered = gather(Z_mpi)

    if Z_gathered !== nothing
        Z_serial = assemble_impedance_matrix(cfie, basis)

        @test size(Z_gathered) == (N, N)
        rel_err = norm(Z_gathered .- Z_serial) / norm(Z_serial)
        @test rel_err < 1e-6
    else
        @test size(Z_mpi) == (N, N)
    end
end
