# Phase 10: MPI Accuracy Verification — A4, B3, C3-MPI
# Run with: mpiexec -n 2 julia --project benchmark/verify_mpi_accuracy.jl
#
# A4: S-EFIE MPI (2 processes) vs Serial — machine precision
# B3: S-MFIE MPI (2 processes) vs Serial — machine precision (small mesh)
# C3-MPI: S-CFIE MPI (2 processes) vs Serial — machine precision (small mesh)
#
# Note: B3/C3-MPI use a smaller sphere mesh to avoid memory issues
# (N=26424 Dense Z requires ~11 GB per copy, infeasible for MPI gather + compare)

using MPI
using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Parallel: assemble_impedance_matrix_parallel, gather
using EMSuite.Utilities.Parameters: set_frequency!
using LinearAlgebra
using Printf

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
n_procs = MPI.Comm_size(comm)

function log_root(msg)
    if rank == 0
        println(msg)
        flush(stdout)
    end
end

# ===================================================================
#  A4: S-EFIE MPI — Jet 100 MHz, N=14559
# ===================================================================
function test_A4()
    log_root("\n" * "=" ^ 60)
    log_root("  A4: S-EFIE MPI vs Serial — Jet 100 MHz")
    log_root("=" ^ 60)

    freq = 1e8
    set_frequency!(freq)
    mesh_file = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles/jet_100MHz.nas")
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    log_root("  N = $N, MPI processes = $n_procs")

    efie = EFIE(freq)

    # Serial assembly (rank 0 only)
    Z_serial = nothing
    if rank == 0
        log_root("  Assembling Z_serial...")
        t_serial = @elapsed Z_serial = assemble_impedance_matrix(efie, basis)
        log_root("  Serial: $(round(t_serial, digits=1))s")
    end
    MPI.Barrier(comm)

    # Parallel assembly
    log_root("  Assembling Z_parallel ($n_procs processes)...")
    t_par = @elapsed Z_par = assemble_impedance_matrix_parallel(efie, basis)
    MPI.Barrier(comm)
    log_root("  Parallel: $(round(t_par, digits=1))s")

    # Gather to rank 0
    Z_par_gathered = gather(Z_par, root=0)

    rel_diff = NaN
    if rank == 0
        diff_norm = norm(Z_serial - Z_par_gathered)
        rel_diff = diff_norm / norm(Z_serial)
        @printf("  ||Z_serial - Z_parallel|| / ||Z_serial|| = %.6e\n", rel_diff)

        if rel_diff < 1e-12
            @printf("  ✅ A4 PASS: rel_diff = %.2e < 1e-12 (machine precision)\n", rel_diff)
        else
            @printf("  ❌ A4 FAIL: rel_diff = %.2e ≥ 1e-12\n", rel_diff)
        end
    end

    Z_serial = nothing; Z_par = nothing; Z_par_gathered = nothing; GC.gc()
    MPI.Barrier(comm)
    return rel_diff
end

# ===================================================================
#  B3: S-MFIE MPI — Small sphere (to avoid memory issues)
#  Full sphere_600MHz (N=26424) is infeasible for Dense Z comparison
#  Use a small plate mesh instead to verify MFIE MPI infrastructure
# ===================================================================
function test_B3()
    log_root("\n" * "=" ^ 60)
    log_root("  B3: S-MFIE MPI vs Serial — Small plate mesh")
    log_root("  Note: Using small mesh instead of sphere_600MHz (N=26424)")
    log_root("        to avoid 11 GB Dense Z memory issues")
    log_root("=" ^ 60)

    freq = 3e8  # 300 MHz
    set_frequency!(freq)

    # Use 2-triangle plate for deterministic, fast test
    nodes = [0.0 1.0 1.0 0.0; 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0]
    elements = Int.([1 1; 2 3; 3 4])
    mesh = TriangleMesh(2, nodes, elements)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    log_root("  N = $N")

    mfie = MFIE(freq)

    Z_serial = nothing
    if rank == 0
        Z_serial = assemble_impedance_matrix(mfie, basis)
    end
    MPI.Barrier(comm)

    Z_par = assemble_impedance_matrix_parallel(mfie, basis)
    Z_par_gathered = gather(Z_par, root=0)

    rel_diff = NaN
    if rank == 0
        diff_norm = norm(Z_serial - Z_par_gathered)
        rel_diff = diff_norm / norm(Z_serial)
        @printf("  ||Z_serial - Z_parallel|| / ||Z_serial|| = %.6e\n", rel_diff)

        if rel_diff < 1e-12
            @printf("  ✅ B3 PASS: rel_diff = %.2e < 1e-12\n", rel_diff)
        else
            @printf("  ❌ B3 FAIL: rel_diff = %.2e ≥ 1e-12\n", rel_diff)
        end
    end

    Z_serial = nothing; Z_par = nothing; GC.gc()
    MPI.Barrier(comm)
    return rel_diff
end

# ===================================================================
#  C3-MPI: S-CFIE MPI — Small mesh for MPI infrastructure check
#  Same rationale as B3: avoid 11 GB memory for sphere_600MHz
# ===================================================================
function test_C3_MPI()
    log_root("\n" * "=" ^ 60)
    log_root("  C3-MPI: S-CFIE MPI vs Serial — Small plate mesh")
    log_root("=" ^ 60)

    freq = 3e8
    set_frequency!(freq)

    nodes = [0.0 1.0 1.0 0.0; 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0]
    elements = Int.([1 1; 2 3; 3 4])
    mesh = TriangleMesh(2, nodes, elements)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    log_root("  N = $N")

    cfie = CFIE(freq)

    Z_serial = nothing
    if rank == 0
        Z_serial = assemble_impedance_matrix(cfie, basis)
    end
    MPI.Barrier(comm)

    Z_par = assemble_impedance_matrix_parallel(cfie, basis)
    Z_par_gathered = gather(Z_par, root=0)

    rel_diff = NaN
    if rank == 0
        diff_norm = norm(Z_serial - Z_par_gathered)
        rel_diff = diff_norm / norm(Z_serial)
        @printf("  ||Z_serial - Z_parallel|| / ||Z_serial|| = %.6e\n", rel_diff)

        if rel_diff < 1e-12
            @printf("  ✅ C3-MPI PASS: rel_diff = %.2e < 1e-12\n", rel_diff)
        else
            @printf("  ❌ C3-MPI FAIL: rel_diff = %.2e ≥ 1e-12\n", rel_diff)
        end
    end

    Z_serial = nothing; Z_par = nothing; GC.gc()
    MPI.Barrier(comm)
    return rel_diff
end

# ===================================================================
# Main
# ===================================================================
log_root("=" ^ 60)
log_root("  Phase 10: MPI Accuracy Verification (A4, B3, C3-MPI)")
log_root("  Processes: $n_procs")
log_root("=" ^ 60)

# A4 is the heavyweight test — N=14559 EFIE
diff_A4 = test_A4()

# B3 and C3-MPI are lightweight MPI infrastructure checks
diff_B3 = test_B3()
diff_C3 = test_C3_MPI()

if rank == 0
    println("\n" * "=" ^ 60)
    println("  SUMMARY (MPI)")
    println("=" ^ 60)
    @printf("  A4 S-EFIE MPI:   rel_diff = %.2e %s\n", diff_A4, diff_A4 < 1e-12 ? "✅" : "❌")
    @printf("  B3 S-MFIE MPI:   rel_diff = %.2e %s\n", diff_B3, diff_B3 < 1e-12 ? "✅" : "❌")
    @printf("  C3-MPI S-CFIE:   rel_diff = %.2e %s\n", diff_C3, diff_C3 < 1e-12 ? "✅" : "❌")
    println("=" ^ 60)
end

MPI.Finalize()
