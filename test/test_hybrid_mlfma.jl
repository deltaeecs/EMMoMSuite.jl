# 混合并行（MPI × Julia 线程）MLFMA 精度门
# 运行: mpiexec -n <P> julia -t <T> --project=. test/test_hybrid_mlfma.jl
# 门：混合 matvec/求解与串行参照相对差在容差内（SPMD 设计应接近逐位一致）。
using Test
using MPI
using EMMoMSuite
using LinearAlgebra, Random
using EMMoMSuite.Geometry, EMMoMSuite.BasisFunctions, EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.MLFMA: MLFMAOperatorMPI

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    T = Threads.nthreads()

    @testset "混合并行 MLFMA 精度门 (P=$P, T=$T)" begin
        @test T >= 1
        freq = 300e6
        λ = 299792458.0 / freq
        mesh = generate_sphere_mesh(1.0, 14, 28)
        basis = RWGBasis(mesh)
        efie = EFIE(freq)
        N = num_basis(basis)
        leaf = 0.25 * λ

        op_mpi = MLFMAOperatorMPI(efie, basis, leaf; comm = comm)
        op_ser = rank == 0 ? MLFMAOperator(efie, basis, leaf) : nothing

        Random.seed!(123)
        x = randn(ComplexF64, N)
        y_mpi = op_mpi * x

        if rank == 0
            y_ser = op_ser * x
            rel = norm(y_mpi - y_ser) / norm(y_ser)
            @info "matvec rel diff (MPI×threads vs serial)" rel P T
            @test rel < 1e-10
        end
        MPI.Barrier(comm)

        # 求解（SPMD：所有秩以相同初值/容差迭代）
        src = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
        V = excitation_vector(efie, src, basis)
        sol = GMRESSolver(restart = 30, maxiter = 120, tol = 1e-6, verbose = false)
        I_mpi = solve!(sol, op_mpi, V)

        if rank == 0
            sol_ser = GMRESSolver(restart = 30, maxiter = 120, tol = 1e-6, verbose = false)
            I_ser = solve!(sol_ser, op_ser, V)
            rel = norm(I_mpi - I_ser) / norm(I_ser)
            @info "solve rel diff (MPI×threads vs serial)" rel P T
            @test rel < 1e-4
        end
        MPI.Barrier(comm)
    end

    @testset "混合并行 MLFMA 精度门 (FFTSpectral, P=$P, T=$T)" begin
        freq = 300e6
        λ = 299792458.0 / freq
        mesh = generate_sphere_mesh(0.8, 14, 28)
        basis = RWGBasis(mesh)
        efie = EFIE(freq)
        leaf = 0.25 * λ
        op_mpi = MLFMAOperatorMPI(
            efie, basis, leaf;
            comm = comm, interp_method = Val(:FFTSpectral),
        )
        op_ser = rank == 0 ?
            MLFMAOperator(efie, basis, leaf, Val(:FFTSpectral)) : nothing

        Random.seed!(456)
        x = randn(ComplexF64, num_basis(basis))
        y_mpi = op_mpi * x
        MPI.Barrier(comm)
        if rank == 0
            y_ser = op_ser * x
            rel = norm(y_mpi - y_ser) / norm(y_ser)
            @info "FFTSpectral matvec rel diff (MPI×threads vs serial)" rel P T
            @test rel < 1e-10
        end
        MPI.Barrier(comm)
    end

    @testset "混合并行 MLFMA 精度门 (Lebedev 1-step, P=$P, T=$T)" begin
        freq = 300e6
        λ = 299792458.0 / freq
        mesh = generate_sphere_mesh(0.8, 14, 28)
        basis = RWGBasis(mesh)
        efie = EFIE(freq)
        leaf = 0.25 * λ
        op_mpi = MLFMAOperatorMPI(
            efie, basis, leaf;
            comm = comm, interp_method = Val(:LbTrained1Step),
        )
        op_ser = rank == 0 ?
            MLFMAOperator(efie, basis, leaf, Val(:LbTrained1Step)) : nothing

        Random.seed!(987)
        x = randn(ComplexF64, num_basis(basis))
        y_mpi = op_mpi * x
        MPI.Barrier(comm)
        if rank == 0
            y_ser = op_ser * x
            rel = norm(y_mpi - y_ser) / norm(y_ser)
            @info "Lebedev matvec rel diff (MPI×threads vs serial)" rel P T
            @test rel < 1e-10
        end
        MPI.Barrier(comm)
    end
    MPI.Finalize()
end

main()
