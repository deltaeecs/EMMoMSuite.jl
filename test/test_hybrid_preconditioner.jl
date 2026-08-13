# test_hybrid_preconditioner.jl — MPI 分布式预条件精度/分布门
# 运行: mpiexec -n <P> julia -t <T> --project=. test/test_hybrid_preconditioner.jl
#
# 门：
#   1. 分布式 BlockJacobi 分块分布：各秩块数 ≈ 总块数/P（内存不再每秩全量复制）
#   2. 施加一致性：apply_mpi_preconditioner!(y, P_mpi, x) == 串行 BlockJacobi \ x
#   3. PMCHW 2N×2N 块（J+M 行）同样分布且施加一致
#   4. distributed_gmres! 接入 Pl（对角预条件）收敛到稠密 LU 参照
using Test
using MPI
using EMMoMSuite
using LinearAlgebra, Random, SparseArrays
using EMMoMSuite.Geometry, EMMoMSuite.BasisFunctions, EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.MLFMA: MLFMAOperatorMPI, MLFMAOperator, PMCHWMLFMAOperatorMPI, PMCHWMLFMAOperator
using EMMoMSuite.Solvers: BlockJacobiPreconditioner
using EMMoMSuite.Parallel:
    DistributedBlockJacobiPreconditioner,
    DistributedDiagonalPreconditioner,
    apply_mpi_preconditioner!,
    mpi_gmres!,
    mpiarray,
    sync!

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)

    @testset "MPI 分布式预条件 (P=$P)" begin
        freq = 300e6
        λ = 299792458.0 / freq

        # ── 1. MLFMA 算子钩子：分块分布 + 施加一致性 ─────────────────────
        mesh = generate_sphere_mesh(0.6, 6, 10)
        basis = RWGBasis(mesh)
        efie = EFIE(freq)
        N = num_basis(basis)
        leaf = 0.25 * λ

        op_mpi = MLFMAOperatorMPI(efie, basis, leaf; comm = comm)
        op_ser = rank == 0 ? MLFMAOperator(efie, basis, leaf) : nothing

        P_mpi = DistributedBlockJacobiPreconditioner(op_mpi)
        n_my = length(P_mpi.block_rows)
        n_total = MPI.Allreduce(n_my, +, comm)

        leaf_level = op_mpi.octree.levels[op_mpi.octree.nLevels]
        n_nonempty = count(c -> !isempty(c.bfInterval), leaf_level.cubes)
        @test n_total == n_nonempty
        @test n_my <= ceil(Int, n_nonempty / P)
        @info "BlockJacobi blocks: total=$n_total my=$n_my (P=$P)"

        Random.seed!(7)
        x = randn(ComplexF64, N)
        y_mpi = zeros(ComplexF64, N)
        apply_mpi_preconditioner!(y_mpi, P_mpi, x)
        MPI.Barrier(comm)
        if rank == 0
            P_ser = BlockJacobiPreconditioner(op_ser)
            y_ser = P_ser \ x
            rel = norm(y_mpi - y_ser) / norm(y_ser)
            @info "BlockJacobi apply rel diff (MPI vs serial)" rel
            @test rel < 1e-10
        end
        MPI.Barrier(comm)

        # ── 2. 对角预条件：施加一致性 ────────────────────────────────────
        P_diag = DistributedDiagonalPreconditioner(op_mpi)
        Random.seed!(8)
        x2 = randn(ComplexF64, N)
        y2 = zeros(ComplexF64, N)
        apply_mpi_preconditioner!(y2, P_diag, x2)
        MPI.Barrier(comm)
        if rank == 0
            d = [op_ser.Z_near[i, i] for i in 1:N]
            d = [abs(d[i]) < 1e-15 ? 1.0 : inv(d[i]) for i in 1:N]
            rel2 = norm(y2 - d .* x2) / norm(d .* x2)
            @info "Diagonal apply rel diff" rel2
            @test rel2 < 1e-10
        end
        MPI.Barrier(comm)

        # ── 3. PMCHW 2N×2N 块（J+M 行）分布 + 施加一致性 ─────────────────
        mesh2 = generate_sphere_mesh(0.5, 4, 6)
        basis2 = RWGBasis(mesh2)
        pmchw = PMCHW(freq, 4.0, 1.0)
        op_pmchw_mpi = PMCHWMLFMAOperatorMPI(pmchw, basis2, 0.10; comm = comm)
        op_pmchw_ser = rank == 0 ? PMCHWMLFMAOperator(pmchw, basis2, 0.10) : nothing

        Pp_mpi = DistributedBlockJacobiPreconditioner(op_pmchw_mpi)
        np_my = length(Pp_mpi.block_rows)
        np_total = MPI.Allreduce(np_my, +, comm)
        pleaf = op_pmchw_mpi.octree0.levels[op_pmchw_mpi.octree0.nLevels]
        np_nonempty = count(c -> !isempty(c.bfInterval), pleaf.cubes)
        @test np_total == np_nonempty

        N2 = size(op_pmchw_mpi, 1)
        Random.seed!(9)
        xp = randn(ComplexF64, N2)
        yp = zeros(ComplexF64, N2)
        apply_mpi_preconditioner!(yp, Pp_mpi, xp)
        MPI.Barrier(comm)
        if rank == 0
            Pp_ser = BlockJacobiPreconditioner(op_pmchw_ser)
            yp_ser = Pp_ser \ xp
            relp = norm(yp - yp_ser) / norm(yp_ser)
            @info "PMCHW BlockJacobi apply rel diff" relp
            @test relp < 1e-10
        end
        MPI.Barrier(comm)

        # ── 4. distributed_gmres! 接入 Pl（对角预条件，稠密参照） ────────
        N3 = 150
        A3 = mpiarray(ComplexF64, N3, N3; comm = comm)
        A_full = Matrix{ComplexF64}(undef, N3, N3)
        Random.seed!(11)
        A_full .= randn(ComplexF64, N3, N3)
        for i in 1:N3
            A_full[i, i] += (50.0 + 0im)
        end
        for j in A3.indices[2], i in 1:N3
            A3[i, j] = A_full[i, j]
        end
        sync!(A3)
        b3 = randn(ComplexF64, N3)
        P3 = DistributedDiagonalPreconditioner([inv(A_full[i, i]) for i in 1:N3])

        x3 = zeros(ComplexF64, N3)
        x3, hist = mpi_gmres!(
            x3, A3, b3; Pl = P3, reltol = 1e-10, maxiter = 400, restart = 60, log = true,
        )
        MPI.Barrier(comm)
        if rank == 0
            x_ref = A_full \ b3
            rel3 = norm(x3 - x_ref) / norm(x_ref)
            @info "distributed_gmres + Diagonal" rel3 iters = hist.mvps converged = hist.isconverged
            @test rel3 < 1e-8
        end
        MPI.Barrier(comm)
    end
    MPI.Finalize()
end

main()
