# test_mpi_sai.jl — MPI 分布式逐八叉树块 SAI 预条件
# 运行: mpiexec -n <P> julia --project=. test/test_mpi_sai.jl
#       （本地单进程也可直接 julia --project=. test/test_mpi_sai.jl，覆盖 P=1 路径）
#
# 门：
#   1. 块分布：各秩拥有 cube 数之和 == 非空叶 cube 总数（与 BlockJacobi 同一归属约定）
#   2. 施加一致性：apply_mpi_preconditioner!(y, P_mpi, x) == 串行 SPAIPreconditioner \ x
#   3. GMRES 迭代次数：块 SAI 显著优于无预条件（闭环 t3 的迭代数改善结论）
using Test
using MPI
using EMMoMSuite
using LinearAlgebra, Random
using IterativeSolvers
using EMMoMSuite.Geometry, EMMoMSuite.BasisFunctions, EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.MLFMA: MLFMAOperatorMPI, MLFMAOperator
using EMMoMSuite.Solvers: SPAIPreconditioner
using EMMoMSuite.Parallel:
    DistributedSPAIPreconditioner, apply_mpi_preconditioner!, mpi_gmres!

# 预条件包装：把 apply_mpi_preconditioner! 适配成 IterativeSolvers 的 ldiv! 接口
struct FuncPrecond{F}
    f::F
end
LinearAlgebra.ldiv!(y, P::FuncPrecond, x) = P.f(y, x)
function LinearAlgebra.ldiv!(P::FuncPrecond, x)
    y = similar(x)
    P.f(y, x)
    copyto!(x, y)
    return x
end

# 算子包装：mul! + size，满足 IterativeSolvers 接口
struct OpWrap{F}
    f::F
    n::Int
end
Base.size(A::OpWrap) = (A.n, A.n)
Base.size(A::OpWrap, d::Int) = d == 1 || d == 2 ? A.n : 1
LinearAlgebra.mul!(y, A::OpWrap, v) = y .= A.f(v)

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)

    @testset "MPI 分布式块 SAI (P=$P)" begin
        freq = 300e6
        λ = 299792458.0 / freq

        mesh = generate_sphere_mesh(0.5, 4, 6)
        basis = RWGBasis(mesh)
        efie = EFIE(freq)
        N = num_basis(basis)
        leaf = 0.25 * λ

        op_mpi = MLFMAOperatorMPI(efie, basis, leaf; comm = comm)
        op_ser = rank == 0 ? MLFMAOperator(efie, basis, leaf) : nothing

        P_mpi = DistributedSPAIPreconditioner(op_mpi)

        # ── 门 1：块归属分布 ─────────────────────────────────────────────
        leaf_level = op_mpi.octree.levels[op_mpi.octree.nLevels]
        n_nonempty = count(c -> !isempty(c.bfInterval), leaf_level.cubes)
        n_my_owned = count(((i, c),) -> !isempty(c.bfInterval) &&
                                                (i - 1) % P == rank,
                           enumerate(leaf_level.cubes))
        n_total_owned = MPI.Allreduce(n_my_owned, +, comm)
        @test n_total_owned == n_nonempty

        # ── 门 2：施加一致（与串行块 SAI 相同 M）─────────────────────────
        Random.seed!(7)
        x = randn(ComplexF64, N)
        y_mpi = zeros(ComplexF64, N)
        apply_mpi_preconditioner!(y_mpi, P_mpi, x)
        MPI.Barrier(comm)
        if rank == 0
            P_ser = SPAIPreconditioner(op_ser)
            y_ser = P_ser \ x
            rel = norm(y_mpi - y_ser) / norm(y_ser)
            @info "SAI apply rel diff (MPI vs serial)" rel
            @test rel < 1e-10
        end
        MPI.Barrier(comm)

        # ── 门 3：GMRES 迭代次数改善（P=1 时分布式施加与串行一致，做闭环验证） ──
        if P == 1 && rank == 0
            Afun = v -> x - op_ser * v   # I - Z（EFIE 无单位项写法）
            b = randn(ComplexF64, N)
            x0 = zeros(ComplexF64, N)
            Aw = OpWrap(Afun, N)
            xg0, hist0 = gmres!(copy(x0), Aw, b; restart = 40, maxiter = 300,
                                log = true, abstol = 1e-10, reltol = 1e-10,
                                initially_zero = true)
            pwrap = FuncPrecond((y, t) -> apply_mpi_preconditioner!(y, P_mpi, t))
            xg1, hist1 = gmres!(copy(x0), Aw, b; Pl = pwrap, restart = 40,
                                maxiter = 300, log = true, abstol = 1e-10,
                                reltol = 1e-10, initially_zero = true)
            @info "GMRES iters" none = hist0.iters sai = hist1.iters
            @test hist1.iters < hist0.iters
        end
    end
    MPI.Finalize()
end

main()
