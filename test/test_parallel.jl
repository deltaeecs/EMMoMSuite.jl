using Test
using EMSuite
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
