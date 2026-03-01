module Parallel

using ..CoreModule
using MPI

include("Threading.jl")
include("MPI/MPIUtils.jl")
include("MPI/MPIArray/MPIArrays.jl")
using .MPIArrays
include("MPI/Assembly.jl")
using .Assembly
# VolumeAssembly extends Assembly.assemble_impedance_matrix_parallel in-place
include("MPI/VolumeAssembly.jl")
# Phase 14.3: 分布式 GMRES (SPMD)
include("MPI/Solver.jl")

export init_parallel!, mpi_rank, nprocs, barrier, @root
export MPIArray, mpiarray, sync!, gather
export assemble_impedance_matrix_parallel
export mpi_gmres!, mpi_gmres

end
