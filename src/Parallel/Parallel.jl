module Parallel

using ..CoreModule
using MPI

include("Threading.jl")
include("MPI/MPIUtils.jl")
include("MPI/MPIArray/MPIArrays.jl")
using .MPIArrays
include("MPI/Assembly.jl")
using .Assembly

export init_parallel!, mpi_rank, nprocs, barrier, @root
export MPIArray, mpiarray, sync!, gather
export assemble_impedance_matrix_parallel

end
