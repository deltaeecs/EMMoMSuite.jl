# 混合并行（MPI × Julia 线程）PMCHW MLFMA 精度门
# 运行: mpiexec -n <P> julia -t <T> --project=. test/test_hybrid_pmchw.jl
# 门：混合 matvec 与串行 PMCHWMLFMAOperator 参照相对差 < 1e-10。
using Test
using MPI
using EMMoMSuite
using LinearAlgebra, Random
using EMMoMSuite.Geometry, EMMoMSuite.BasisFunctions, EMMoMSuite.IntegralEquations

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    T = Threads.nthreads()

    @testset "混合并行 PMCHW 精度门 (P=$P, T=$T)" begin
        freq = 300e6
        # 与 test_pmchw_mlfma_operator.jl 相同夹具（r=0.5, 4×6 → N=54）
        mesh = generate_sphere_mesh(0.5, 4, 6)
        basis = RWGBasis(mesh)
        pmchw = PMCHW(freq, 4.0, 1.0)
        leaf = 0.1

        op_mpi = PMCHWMLFMAOperatorMPI(pmchw, basis, leaf; comm = comm)
        op_ser = rank == 0 ? PMCHWMLFMAOperator(pmchw, basis, leaf) : nothing

        Random.seed!(789)
        x = randn(ComplexF64, 2 * num_basis(basis))
        y_mpi = op_mpi * x
        MPI.Barrier(comm)
        if rank == 0
            y_ser = op_ser * x
            rel = norm(y_mpi - y_ser) / norm(y_ser)
            @info "PMCHW matvec rel diff (MPI×threads vs serial)" rel P T
            @test rel < 1e-10
        end
        MPI.Barrier(comm)
    end
    MPI.Finalize()
end

main()
