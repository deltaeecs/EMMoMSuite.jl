# 混合并行（MPI × Julia 线程）AIM/IE-FFT 精度门
# 运行: mpiexec -n <P> julia -t <T> --project=. test/test_hybrid_aim.jl
# 门：混合 matvec 与串行 AIM 参照相对差 < 1e-10（SPMD + 行分区应接近逐位一致）。
using Test
using MPI
using EMMoMSuite
using LinearAlgebra, Random
using FFTW
using EMMoMSuite.Geometry, EMMoMSuite.BasisFunctions, EMMoMSuite.IntegralEquations

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    T = Threads.nthreads()
    FFTW.set_num_threads(T)

    @testset "混合并行 AIM 精度门 (P=$P, T=$T)" begin
        freq = 300e6
        λ = 299792458.0 / freq
        mesh = generate_sphere_mesh(0.8, 14, 28)
        basis = RWGBasis(mesh)
        efie = EFIE(freq)
        op_mpi = AIMOperatorMPI(efie, basis; h_ratio = 0.1, near_radius = 0.35, comm = comm)
        op_ser = rank == 0 ? AIMOperator(efie, basis; h_ratio = 0.1, near_radius = 0.35) : nothing

        Random.seed!(321)
        x = randn(ComplexF64, num_basis(basis))
        y_mpi = op_mpi * x

        if rank == 0
            y_ser = op_ser * x
            rel = norm(y_mpi - y_ser) / norm(y_ser)
            @info "AIM matvec rel diff (MPI×threads vs serial)" rel P T
            @test rel < 1e-10
        end
        MPI.Barrier(comm)
    end
    MPI.Finalize()
end

main()
