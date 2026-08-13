# ScaLAPACK（MinGW/MSMPI）分布式稠密 LU 效率基准
# 用法: mpiexec -n <P> julia -t <T> --project=. benchmark/benchmark_scalapack_lu.jl [N [MB [NB]]]
using MPI, LinearAlgebra, Random
using EMMoMSuite.Parallel: scalapack_lu_solve

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 600
    MB = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 64
    NB = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : MB
    Random.seed!(5)
    A = randn(ComplexF64, N, N) + N * I
    b = randn(ComplexF64, N)

    MPI.Barrier(comm)
    t0 = MPI.Wtime()
    x = scalapack_lu_solve(A, b, comm; MB = MB, NB = NB)
    MPI.Barrier(comm)
    t_total = MPI.Wtime() - t0

    if rank == 0
        xr = A \ b
        println("SLU_BENCH P=", P, " N=", N, " MB=", MB, " NB=", NB, " total_s=", round(t_total, digits=3),
                " relerr=", round(norm(x - xr) / norm(xr), digits=2))
    end
    MPI.Finalize()
end

main()
