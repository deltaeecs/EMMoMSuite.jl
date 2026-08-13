# 混合并行（MPI × Julia 线程）AIM/IE-FFT 效率基准
# 用法: mpiexec -n <P> julia -t <T> --project=. benchmark/benchmark_hybrid_aim.jl [radius] [nθ] [nφ] [maxiter]
# 输出一行 AIM_HYBRID_BENCH 记录：P、T、N、setup/matvec/solve 墙钟时间。
using MPI
using EMMoMSuite
using LinearAlgebra, Random, Printf
using FFTW
using EMMoMSuite.Geometry, EMMoMSuite.BasisFunctions, EMMoMSuite.IntegralEquations

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    T = Threads.nthreads()
    FFTW.set_num_threads(T)   # 秩内 FFTW 线程（与 Julia @threads 并行路径互补）
    radius = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1.5
    nθ = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20
    nφ = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 40
    maxiter = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 40

    freq = 300e6
    λ = 299792458.0 / freq
    mesh = generate_sphere_mesh(radius, nθ, nφ)
    basis = RWGBasis(mesh)
    efie = EFIE(freq)
    N = num_basis(basis)

    t0 = MPI.Wtime()
    op = AIMOperatorMPI(efie, basis; h_ratio = 0.1, near_radius = 0.35, comm = comm)
    MPI.Barrier(comm)
    t_setup = MPI.Wtime() - t0

    Random.seed!(7)
    x = randn(ComplexF64, N)
    op * x
    MPI.Barrier(comm)
    t0 = MPI.Wtime()
    nrep = 10
    for _ in 1:nrep
        op * x
    end
    MPI.Barrier(comm)
    t_matvec = (MPI.Wtime() - t0) / nrep

    src = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(efie, src, basis)
    sol = GMRESSolver(restart = 30, maxiter = maxiter, tol = 1e-5, verbose = false)
    t0 = MPI.Wtime()
    I = solve!(sol, op, V)
    MPI.Barrier(comm)
    t_solve = MPI.Wtime() - t0

    if rank == 0
        println(
            "AIM_HYBRID_BENCH P=", P, " T=", T, " N=", N,
            " setup_s=", round(t_setup, digits = 2),
            " matvec_ms=", round(t_matvec * 1000, digits = 2),
            " solve_s=", round(t_solve, digits = 2),
            " |I|=", round(norm(I), digits = 4),
        )
    end
    MPI.Finalize()
end

main()
