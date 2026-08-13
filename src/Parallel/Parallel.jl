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
# Phase 15.4 Method C: 分布式 Krylov GMRES（消除 I-3）
include("MPI/DistributedGMRES.jl")
# ScaLAPACK（本机 MinGW，MSMPI 版）分布式稠密 LU：pzgesv + BLACS
include("MPI/ScaLAPACKLU.jl")
# MPI 分布式预条件（块按 cube 归属分秩；接入 distributed_gmres! 的 Pl）
include("MPI/Preconditioners.jl")
# Public solver API (mpi_gmres! / mpi_gmres) — 委托给 distributed_gmres!
include("MPI/Solver.jl")

export init_parallel!, mpi_rank, nprocs, barrier, @root
export MPIArray, mpiarray, sync!, gather, local_col_index
export assemble_impedance_matrix_parallel
export mpi_gmres!, mpi_gmres
export scalapack_lu_solve
export DistributedBlockJacobiPreconditioner, DistributedDiagonalPreconditioner
export apply_mpi_preconditioner!

end
