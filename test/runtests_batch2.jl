# 分批覆盖率测试脚本 - 第2批: MLFMA + 并行（中速）
using Test
using EMMoMSuite

include("test_mlfma.jl")
include("test_parallel.jl")
include("test_mpiarray.jl")
include("test_parallel_mfie_cfie.jl")
include("test_distributed_gmres.jl")
include("test_mpi_array_utils.jl")
