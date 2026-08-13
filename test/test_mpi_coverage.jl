# test_mpi_coverage.jl — 单进程（P=1）MPI 覆盖率补测
#
# 覆盖 src/Parallel/* 中未进入 light_cov 的模块：
#   - MPIUtils / Threading / Parallel 顶层（init_parallel!、@root 等）
#   - MPIArray 数据通路（构造/索引/fill!/getdata/sum/local_col_index）
#   - Assembly.jl（RWG MPI 列分区装配，P=1 与串行对照）
#   - DistributedGMRES.jl + Solver.jl（分布式 GMRES 小系统）
#   - Preconditioners.jl（分布式块 Jacobi / 对角，构造与施加）
#   - ScaLAPACKLU.jl 纯逻辑（numroc / distribute / gather_solution，不依赖外部库）
#
# 运行方式：单进程 julia（MPI.Init() 单进程即可），无需 mpiexec。
using Test
using EMMoMSuite
using EMMoMSuite.Parallel
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using LinearAlgebra
using SparseArrays
using Random
using MPI

if !MPI.Initialized()
    MPI.Init()
end
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

@testset "Parallel 顶层与工具" begin
    init_parallel!()
    @test mpi_rank() == rank
    @test nprocs() == MPI.Comm_size(comm)
    @test num_threads() >= 1
    @test thread_id() >= 1
    @root begin
        @test mpi_rank() == 0
    end
    barrier()
    @test true
end

@testset "MPIArray 数据通路" begin
    A = mpiarray(ComplexF64, (6, 6); comm = comm)
    @test size(A) == (6, 6)
    @test length(A) == 36
    fill!(A, 2.0 + 0im)
    @test all(EMMoMSuite.Parallel.getdata(A) .== 2.0)
    # 索引读写
    A[1, 1] = 3.0 + 0im
    @test A[1, 1] == 3.0 + 0im
    # 本地列索引
    local_cols = A.indices[2]
    @test local_col_index(A, first(local_cols)) == 1
    @test local_col_index(A, last(local_cols)) == length(local_cols)
    # 向量构造与 gather
    x = mpiarray(ComplexF64, 8; comm = comm)
    fill!(x, 0.5 + 0im)
    xg = gather(x)
    if rank == 0
        @test all(xg .≈ 0.5)
    else
        @test isnothing(xg)
    end
    # sum（分布式归约）
    @test real(sum(A)) == 2.0 * 35 + 3.0  # 35 个 2.0 + 1 个 3.0
    # similar / copyto!
    B = similar(A)
    copyto!(B, A)
    @test B[1, 1] == 3.0 + 0im
    sync!(A)
    @test true
end

@testset "MPI RWG 装配（P=1，与串行对照）" begin
    freq = 300e6
    mesh = generate_sphere_mesh(0.3, 4, 8)
    basis = RWGBasis(mesh)
    efie = EFIE(freq)
    Z = assemble_impedance_matrix_parallel(efie, basis)
    N = num_basis(basis)
    @test size(Z) == (N, N)
    Zs = assemble_impedance_matrix(efie, basis)
    Zd = EMMoMSuite.Parallel.getdata(Z)
    # P=1 时本地数据即完整列分区（全部列）
    @test size(Zd, 1) == N
    # 抽样对比对角
    for i in 1:min(4, N)
        @test abs(Zd[i, i] - Zs[i, i]) < 1e-10 * max(abs(Zs[i, i]), 1e-30)
    end
end

@testset "分布式 GMRES（P=1 小系统）" begin
    rp = EMMoMSuite.Parallel.row_partition
    @test rp(10, 0, 1) == 1:10
    @test rp(0, 0, 1) == 1:0
    N = 16
    A = mpiarray(ComplexF64, N, N; comm = comm)
    Random.seed!(3)
    Afull = randn(ComplexF64, N, N) + 4I
    for j in A.indices[2], i in 1:N
        A[i, j] = Afull[i, j]
    end
    sync!(A)
    b = randn(ComplexF64, N)
    x = zeros(ComplexF64, N)
    x, hist = mpi_gmres!(x, A, b; reltol = 1e-8, maxiter = 200, restart = 16, log = true)
    @test hist.isconverged
    @test norm(Afull \ b - x) / norm(b) < 1e-6
    # Solver.jl 包装 API
    x2 = mpi_gmres(A, b; reltol = 1e-8, maxiter = 200, restart = 16)
    @test norm(x2 - x) / norm(x) < 1e-6
end

@testset "分布式预条件（P=1）" begin
    n = 8
    Random.seed!(5)
    A_full = (10.0 + 0im) * I + 0.5 * randn(ComplexF64, n, n)
    Asp = sparse(A_full)
    block_rows = [collect(i:i) for i in 1:n]
    P = DistributedBlockJacobiPreconditioner(Asp, block_rows, comm)
    x = randn(ComplexF64, n)
    y = zeros(ComplexF64, n)
    apply_mpi_preconditioner!(y, P, x)
    # 1×1 块：M⁻¹ = diag(A)⁻¹
    @test norm(y - x ./ diag(A_full)) / norm(x) < 1e-10
    # 对角预条件
    Pd = DistributedDiagonalPreconditioner([inv(A_full[i, i]) for i in 1:n])
    y2 = zeros(ComplexF64, n)
    apply_mpi_preconditioner!(y2, Pd, x)
    @test norm(y2 - x ./ diag(A_full)) / norm(x) < 1e-10
    # 无预条件（恒等）
    y3 = zeros(ComplexF64, n)
    apply_mpi_preconditioner!(y3, nothing, x)
    @test y3 == x
end

@testset "ScaLAPACK 纯逻辑（numroc/distribute/gather，P=1）" begin
    @test EMMoMSuite.Parallel.numroc(200, 64, 0, 0, 1) == 200
    @test EMMoMSuite.Parallel.numroc(1, 64, 0, 0, 1) == 1
    @test EMMoMSuite.Parallel.numroc(0, 64, 0, 0, 1) == 0

    lib = EMMoMSuite.Parallel.SCALAPACK_LIB
    if lib === nothing
        @test_skip "ScaLAPACK 动态库未安装（无 SCALAPACK_LIB_PATH），跳过 BLACS 网格测试"
    else
        N = 8
        A = randn(ComplexF64, N, N) + N * I
        grid = EMMoMSuite.Parallel.init_grid(comm; MB = 4, NB = 4)
        grid = EMMoMSuite.Parallel.ScaLAPACKGrid(
            grid.ictxt, grid.nprow, grid.npcol, grid.myrow, grid.mycol, 4, 4, N,
        )
        Aloc = EMMoMSuite.Parallel.distribute(A, grid; MB = 4, NB = 4)
        @test size(Aloc, 1) == N
        @test size(Aloc, 2) == N
        # 分发往返：本地块重新映射回全局应等于 A
        Atest = zeros(ComplexF64, N, N)
        myrow = Int(grid.myrow); mycol = Int(grid.mycol)
        lr = 1
        for bi in 1:cld(N, 4)
            (bi - 1) % grid.nprow == myrow || continue
            i0 = (bi - 1) * 4 + 1
            nrows = min(4, N - i0 + 1)
            lc = 1
            for bj in 1:cld(N, 4)
                (bj - 1) % grid.npcol == mycol || continue
                j0 = (bj - 1) * 4 + 1
                ncols = min(4, N - j0 + 1)
                Atest[i0:i0+nrows-1, j0:j0+ncols-1] .= Aloc[lr:lr+nrows-1, lc:lc+ncols-1]
                lc += ncols
            end
            lr += nrows
        end
        @test norm(Atest - A) / norm(A) < 1e-12
        EMMoMSuite.Parallel._blacs_gridexit!(grid.ictxt)
    end
end
