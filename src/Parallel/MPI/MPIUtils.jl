using MPI

"""
    init_parallel!()

Initialize the parallel environment (MPI).
"""
function init_parallel!()
    if !MPI.Initialized()
        MPI.Init()
    end
end

"""
    mpi_rank()

Return the rank of the current process.
"""
mpi_rank() = MPI.Comm_rank(MPI.COMM_WORLD)

"""
    nprocs()

Return the number of processes.
"""
nprocs() = MPI.Comm_size(MPI.COMM_WORLD)

"""
    barrier()

Block until all processes have reached this routine.
"""
barrier() = MPI.Barrier(MPI.COMM_WORLD)

"""
    @root expr

Execute the expression only on the root process (rank 0).
"""
macro root(expr)
    quote
        if mpi_rank() == 0
            $(esc(expr))
        end
    end
end
