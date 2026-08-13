# ScaLAPACK（MinGW/MSMPI）分布式稠密 LU 精度门
# 运行: mpiexec -n <P> julia -t <T> --project=. test/test_scalapack_lu.jl
# 门：ScaLAPACK 解 == 串行 LU 解（相对差 < 1e-8）。
using Test
using MPI
using EMMoMSuite
using LinearAlgebra, Random

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)

    @testset "ScaLAPACK 分布式 LU 精度门 (P=$P)" begin
        N = 200
        Random.seed!(17)
        A = randn(ComplexF64, N, N) + N * I
        b = randn(ComplexF64, N)
        x_ref = A \ b
        x = EMMoMSuite.Parallel.scalapack_lu_solve(A, b, comm)
        if rank == 0
            rel = norm(x - x_ref) / norm(x_ref)
            @info "ScaLAPACK LU solve rel err vs serial" rel
            @test rel < 1e-8
        end
        MPI.Barrier(comm)
    end
    MPI.Finalize()
end

main()
