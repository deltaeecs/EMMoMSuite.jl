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

@testset "MPIArray LinearAlgebra - rmul!/axpy!/norm/dot" begin
    n = 4
    # rmul!: scale MPIVector in place
    v1 = mpiarray(Float64, n; partition=(1,), buffersize=0)
    fill!(v1, 3.0)
    LinearAlgebra.rmul!(v1, 2.0)
    @test all(v1.data .≈ 6.0)

    # axpy!: y = a*x + y
    v1 = mpiarray(Float64, n; partition=(1,), buffersize=0)
    v2 = mpiarray(Float64, n; partition=(1,), buffersize=0)
    fill!(v1, 1.0)
    fill!(v2, 0.5)
    LinearAlgebra.axpy!(2.0, v1, v2)
    @test all(v2.data .≈ 2.5)  # 2*1.0 + 0.5 = 2.5

    # norm (Allreduce, root=-1 default branch)
    fill!(v1, 1.0)
    n_val = LinearAlgebra.norm(v1)
    @test isfinite(n_val) && n_val > 0

    # dot (Allreduce, root=-1 default branch)
    fill!(v1, 1.0); fill!(v2, 2.0)
    d = LinearAlgebra.dot(v1, v2)
    @test isfinite(real(d)) && real(d) > 0

    # mul! for SubOrMPIArray (view of MPIMatrix — not MPIMatrix path)
    mat = mpiarray(Float64, (4, 4); partition=(1, 1), buffersize=0)
    for i = 1:4; mat[i, i] = Float64(i); end   # diagonal
    sub_mat = view(mat, 1:2, 1:2)               # SubMPIMatrix
    y_vec = zeros(Float64, 2)
    x_vec = [1.0, 0.0]
    LinearAlgebra.mul!(y_vec, sub_mat, x_vec)
    @test isapprox(y_vec[1], 1.0; atol=1e-10)  # [1,0] · [1,0] via 1x1 block
    @test isapprox(y_vec[2], 0.0; atol=1e-10)

    # Base.:* error throw for MPIMatrix × AbstractVector (throws a String, not Exception)
    threw = false
    try
        mat * [1.0, 2.0, 3.0, 4.0]
    catch
        threw = true
    end
    @test threw
end
