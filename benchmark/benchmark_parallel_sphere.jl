using MPI
using EMSuite
using LinearAlgebra
using Printf
using CSV
using DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
const REPORT_DIR = joinpath(ROOT, "test_results", "reports")
const PARALLEL_SAMPLE_CSV = joinpath(REPORT_DIR, "PARALLEL_MPI_SAMPLE.csv")

# Benchmark: Parallel PEC Sphere Scattering
# -----------------------------------------
# This script tests the MPI parallel execution capability.
# Run with: mpiexecjl -n <N> julia benchmark/benchmark_parallel_sphere.jl

function run_parallel_benchmark()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    comm_size = MPI.Comm_size(comm)
    
    if rank == 0
        println("==================================================")
        println("   Benchmark: Parallel PEC Sphere (MPI n=$comm_size)   ")
        println("==================================================")
    end
    
    # Initialize Parallel module in EMSuite
    init_parallel!()
    
    # 1. Parameters
    freq = 300e6
    radius = 1.0
    
    # 2. Mesh (Only on root or all? Usually all generate or root broadcasts)
    # For simplicity, all generate.
    mesh = generate_sphere_mesh(radius, 12, 24) # Coarser for quick test
    
    if rank == 0
        println("Mesh: $(num_vertices(mesh)) vertices, $(num_elements(mesh)) elements")
    end
    
    # 3. Basis
    basis = RWGBasis(mesh)
    if rank == 0
        println("Unknowns: $(num_basis(basis))")
    end
    
    # 4. Assembly (Parallel)
    MPI.Barrier(comm)
    t_start = MPI.Wtime()
    
    efie = EFIE(freq)
    
    if nprocs() > 1
        if rank == 0 println("Assembling in Parallel...") end
        Z = assemble_impedance_matrix_parallel(efie, basis)
    else
        if rank == 0 println("Assembling in Serial...") end
        Z = assemble_impedance_matrix(efie, basis)
    end
    
    MPI.Barrier(comm)
    t_end = MPI.Wtime()
    assembly_time = t_end - t_start
    
    if rank == 0
        println(@sprintf("Assembly Time: %.4f s", assembly_time))
        println("Matrix size: $(size(Z))")
    end
    
    # 5. Solve
    # Parallel solve?
    # If Z is a distributed matrix (e.g. from PartitionedArrays or similar), solve! needs to handle it.
    # If Z is local (replicated), then it's not true parallel solve.
    # Let's check if Z is distributed.
    
    # For this benchmark, we just run the solve and see if it works.
    source = PlaneWave(freq, pi, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(source, basis)
    
    solver = GMRESSolver(tol=1e-4, maxiter=500)
    I = solve!(solver, Z, V)
    
    if rank == 0
        current_norm = norm(I)
        println("Solved. Current magnitude: $(current_norm)")
        mkpath(REPORT_DIR)
        CSV.write(PARALLEL_SAMPLE_CSV, DataFrame(
            case_name = ["Parallel PEC Sphere"],
            mpi_procs = [comm_size],
            threads_per_rank = [Threads.nthreads()],
            unknowns = [num_basis(basis)],
            assembly_time_s = [assembly_time],
            current_norm = [current_norm],
        ))
        println("Parallel sample CSV written: $(PARALLEL_SAMPLE_CSV)")
    end
    
    MPI.Finalize()
end

run_parallel_benchmark()
