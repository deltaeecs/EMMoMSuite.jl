#!/usr/bin/env julia
"""
benchmark_mpi_performance.jl

MPI 并行性能基准 — V-EFIE 强扩展性测试.

用法（需从 MoM 根目录运行）:
  1 进程: julia --project=EMSuite -t 4 EMSuite/benchmark/benchmark_mpi_performance.jl
  N 进程: <mpiexec> -n <N> julia --project=EMSuite -t 4 EMSuite/benchmark/benchmark_mpi_performance.jl
"""

using MPI
MPI.Init()

using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.IntegralEquations.VEFIEModule: VEFIE
using EMSuite.Parallel
using LinearAlgebra

comm    = MPI.COMM_WORLD
rank    = MPI.Comm_rank(comm)
n_procs = MPI.Comm_size(comm)

init_parallel!()

# ══════════════════════════════════════════════════════════════════════════════
#  测试配置
# ══════════════════════════════════════════════════════════════════════════════
MOM_DIR    = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles")
mesh_file  = joinpath(MOM_DIR, "plate_and_metal_1dot2GHz.nas")
freq       = 1.2e9
eps_r_diel = complex(2.56, -0.02)   # 电介质相对介电常数
eps_r_bg   = complex(1.0,  0.0)     # 背景（金属内部不参与, 自动跳过）

rank == 0 && println("="^70)
rank == 0 && println("  MPI 并行性能基准 — V-EFIE 强扩展性")
rank == 0 && println("="^70)
rank == 0 && println("  MPI processes : $n_procs")
rank == 0 && println("  Threads/proc  : $(Threads.nthreads())")
rank == 0 && println("="^70)

if !isfile(mesh_file)
    rank == 0 && @warn "Mesh not found: $mesh_file\nSkipping V-EFIE benchmark."
    MPI.Finalize()
    exit(1)
end

# ── 加载网格 ─────────────────────────────────────────────────────────────────
mesh = read_nas_mesh(mesh_file; scale=1.0)
basis = SWGBasis(mesh)
N = num_basis(basis)

rank == 0 && println("  Mesh  : $(basename(mesh_file))")
rank == 0 && println("  N_SWG : $N")
rank == 0 && println()

# 材料属性 (dielectric slab; background = vacuum)
n_tets = length(mesh.tetras)
permittivities = fill(eps_r_bg, n_tets)
# 对所有四面体赋予介质参数 (简化: 整体介质体)
fill!(permittivities, eps_r_diel)

# ── 并行组装 ─────────────────────────────────────────────────────────────────
vefie = VEFIE(freq, permittivities)

MPI.Barrier(comm)
t_start = MPI.Wtime()

Z = assemble_impedance_matrix_parallel(vefie, basis, permittivities)

MPI.Barrier(comm)
t_assembly = MPI.Wtime() - t_start

# ── 输出结果 ─────────────────────────────────────────────────────────────────
if rank == 0
    println("  Assembly time : $(round(t_assembly; digits=3)) s")
    println("  Matrix size   : $(size(Z))")
    println()

    # 保存结果供自动化脚本读取
    results_file = joinpath(@__DIR__, "mpi_result_n$(n_procs)_t$(Threads.nthreads()).txt")
    open(results_file, "w") do io
        println(io, "Test: VEFIE_$(basename(mesh_file))")
        println(io, "N: $N")
        println(io, "MPI_processes: $n_procs")
        println(io, "Threads: $(Threads.nthreads())")
        println(io, "Assembly_time: $t_assembly")
    end
    println("  Result saved → $(basename(results_file))")
    println("="^70)
end

MPI.Finalize()

