using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using LinearAlgebra
using MPI

@testset "Parallel Computing" begin
    # Initialize MPI
    init_parallel!()
    
    r = mpi_rank()
    n = nprocs()
    
    @info "Running on rank $r of $n processes"
    
    @test r >= 0
    @test n >= 1
    
    # Test @root macro
    @root begin
        @info "This should only appear on rank 0"
        @test mpi_rank() == 0
    end
    
    # Test Threading
    nt = num_threads()
    tid = thread_id()
    
    @info "Running on thread $tid of $nt threads"
    
    @test nt >= 1
    @test tid >= 1
    @test tid <= nt
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 14.3: mpi_gmres! — MPIMatrix + 分布式 GMRES (P=1 单进程验证)
# ─────────────────────────────────────────────────────────────────────────────
@testset "Phase 14.3: mpi_gmres! (single-process)" begin
    if !MPI.Initialized()
        MPI.Init()
    end

    # Small EFIE problem: sphere at 300 MHz, N ≈ 100-200 unknowns
    radius = 0.3
    freq   = 300e6
    mesh   = generate_sphere_mesh(radius, 5, 10)
    basis  = RWGBasis(mesh)
    N      = num_basis(basis)
    @test N > 0

    efie = EFIE(freq)

    # Assemble distributed impedance matrix (MPIMatrix, column-partitioned)
    Z = assemble_impedance_matrix_parallel(efie, basis)

    # Build a known RHS via b = Z_full * x_ref using mul!
    import Random: seed!
    seed!(42)
    x_ref = randn(ComplexF64, N)
    b     = zeros(ComplexF64, N)
    mul!(b, Z, x_ref)

    # Solve A*x = b with mpi_gmres
    x_sol = mpi_gmres(Z, b; restart = min(50, N), maxiter = 3*N, reltol = 1e-8)

    # Verify residual via mul!
    r = zeros(ComplexF64, N)
    mul!(r, Z, x_sol)
    r .-= b
    rel_res = norm(r) / norm(b)
    @test rel_res < 1e-6

    # Solution should be close to x_ref
    @test norm(x_sol .- x_ref) / norm(x_ref) < 1e-4
end
