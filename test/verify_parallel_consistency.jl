using MPI
using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Parallel: assemble_impedance_matrix_parallel, gather
using Test
using LinearAlgebra

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)

if rank == 0
    println("Running parallel consistency check on $comm_size processes")
end

# 1. Setup Geometry (Simple Plate)
# 2 triangles forming a square
# Nodes: (0,0,0), (1,0,0), (1,1,0), (0,1,0)
nodes = [0.0 1.0 1.0 0.0; 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0] # 3x4
elements = [1 1; 2 3; 3 4] # 3x2
# Note: TriangleMesh constructor expects Int for elements
elements_int = Int.(elements)
mesh = TriangleMesh(2, nodes, elements_int)

# 2. Basis
basis = RWGBasis(mesh)
if rank == 0
    println("Number of basis functions: $(num_basis(basis))")
end

# 3. Operator
freq = 300e6
k = 2 * pi * freq / 299792458.0
efie = EFIE(k)

# 4. Serial Assembly (Reference)
Z_serial = nothing
if rank == 0
    println("Assembling Serial Matrix...")
    Z_serial = assemble_impedance_matrix(efie, basis)
end

# 5. Parallel Assembly
if rank == 0
    println("Assembling Parallel Matrix...")
end
# assemble_impedance_matrix_parallel returns an MPIArray
Z_parallel_mpi = assemble_impedance_matrix_parallel(efie, basis)

# 6. Gather
Z_parallel_gathered = gather(Z_parallel_mpi, root=0)

# 7. Compare
if rank == 0
    println("Comparing results...")
    # Z_parallel_gathered might be nothing on other ranks, but we are in rank 0 block
    
    # Check dimensions
    @test size(Z_serial) == size(Z_parallel_gathered)
    
    diff = norm(Z_serial - Z_parallel_gathered)
    rel_diff = diff / norm(Z_serial)
    println("Difference norm: $diff")
    println("Relative difference: $rel_diff")
    
    @test rel_diff < 1e-12
    
    if rel_diff < 1e-12
        println("SUCCESS: Parallel assembly matches Serial assembly.")
    else
        println("FAILURE: Parallel assembly mismatch.")
        exit(1)
    end
end

# MPI.Finalize() is usually handled by atexit in Julia MPI, but explicit is fine if not using interactive
