using EMSuite
using EMSuite.Parallel
using MPI
using Test
using LinearAlgebra

# Initialize MPI if not already initialized
if !MPI.Initialized()
    init_parallel!()
end

comm = MPI.COMM_WORLD
rank = mpi_rank()
np = nprocs()

@testset "MPIArray basics" begin
    # Basic construction
    A = mpiarray(ComplexF64, (10, 10); partition = (1, np), buffersize = 3)
    @test size(A) == (10, 10)
    @test length(A) == 100

    # 3D array
    if np >= 2
        A = mpiarray(ComplexF64, 10, 10, 10; partition = (2, 1, np÷2), buffersize = 3)
        
        fill!(A, rank)
        sync!(A)
        
        # Check ghost data
        # Accessing internal fields for testing might require importing them or using property access if allowed
        # Assuming we can access fields:
        for (gr, gidcs) in A.grank2ghostindices
            @test all(A.ghostdata[gidcs...] .== gr)
        end
    end

    # Vector
    x = mpiarray(ComplexF64, np + 1; partition = (np, ), buffersize = 3)
    fill!(x, 0.1)
    xlc = gather(x)
    if rank == 0
        @test all(xlc .≈ 0.1)
    else
        @test isnothing(xlc)
    end
end

@testset "MPIMatrix local column indexing" begin
    A = mpiarray(Float64, (8, 8); partition = (1, np), buffersize = 0)
    local_cols = A.indices[2]

    @test local_col_index(A, first(local_cols)) == 1
    @test local_col_index(A, last(local_cols)) == length(local_cols)

    if first(local_cols) > 1
        @test_throws BoundsError local_col_index(A, first(local_cols) - 1)
    end
    if last(local_cols) < size(A, 2)
        @test_throws BoundsError local_col_index(A, last(local_cols) + 1)
    end
end
