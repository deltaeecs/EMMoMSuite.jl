#!/usr/bin/env julia
"""
benchmark_mpi_performance.jl

系统测试 MPI 并行性能 — 强扩展性（固定问题规模，增加进程数）。

用法:
  单次运行: mpiexecjl -n <N> julia --project=EMSuite EMSuite/benchmark/benchmark_mpi_performance.jl
  自动化扫描: julia EMSuite/benchmark/run_mpi_sweep.jl
"""

using MPI
using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Parallel
using LinearAlgebra
using Printf

function benchmark_mpi_efies_efie()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    n_procs = MPI.Comm_size(comm)
    
    # 初始化并行模块
    init_parallel!()
    
    # ══════════════════════════════════════════════════════════════════════
    #  测试配置
    # ══════════════════════════════════════════════════════════════════════
    const MOM_DIR = joinpath(@__DIR__, "../../MoM_AllinOne/meshfiles")
    
    # 测试用例1: Plate 300MHz (N=2640) — 适合快速测试
    test_cases = [
        (name="Plate_300MHz", file="plate_standard_0p3GHz.nas", freq=300e6, scale=1.0),
    ]
    
    rank == 0 && println("="^70)
    rank == 0 && println("  MPI 并行性能基准 — S-EFIE 强扩展性")
    rank == 0 && println("="^70)
    rank == 0 && @printf("  MPI processes: %d\n", n_procs)
    rank == 0 && @printf("  Threads per process: %d\n", Threads.nthreads())
    rank == 0 && println("="^70)
    
    for (idx, tc) in enumerate(test_cases)
        mesh_file = joinpath(MOM_DIR, tc.file)
        
        if !isfile(mesh_file)
            rank == 0 && @warn "Mesh not found: $(tc.file), skipping..."
            continue
        end
        
        rank == 0 && println("\n[$idx/$(length(test_cases))] $(tc.name)")
        
        # ── 加载网格 ────────────────────────────────────────────────────
        mesh = read_nas_mesh(mesh_file; scale=tc.scale)
        basis = RWGBasis(mesh)
        N = num_basis(basis)
        
        rank == 0 && @printf("  N = %d basis functions\n", N)
        
        # ── 装配 Z 矩阵（并行）──────────────────────────────────────────
        efie = EFIE(tc.freq)
        
        MPI.Barrier(comm)
        t_assembly_start = MPI.Wtime()
        
        Z = assemble_impedance_matrix_parallel(efie, basis)
        
        MPI.Barrier(comm)
        t_assembly_end = MPI.Wtime()
        t_assembly = t_assembly_end - t_assembly_start
        
        rank == 0 && @printf("  Assembly time: %.4f s\n", t_assembly)
        rank == 0 && @printf("  Matrix size: %s\n", size(Z))
        
        # ══════════════════════════════════════════════════════════════════
        #  输出结果（仅 rank 0）
        # ══════════════════════════════════════════════════════════════════
        if rank == 0
            # 保存到文件供外部脚本分析
            results_file = "benchmark_mpi_n$(n_procs)_$(tc.name).txt"
            open(results_file, "w") do io
                println(io, "Test: $(tc.name)")
                println(io, "N: $N")
                println(io, "MPI_processes: $n_procs")
                println(io, "Threads: $(Threads.nthreads())")
                println(io, "Assembly_time: $t_assembly")
            end
            println("  Result saved to: $results_file")
        end
    end
    
    rank == 0 && println("\n" * "="^70)
    rank == 0 && println("  Benchmark complete.")
    rank == 0 && println("="^70)
    
    MPI.Finalize()
end

benchmark_mpi_efies_efie()
