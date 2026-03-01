#!/usr/bin/env julia
"""
run_mpi_sweep.jl

自动化运行 MPI 性能扫描：测试 1, 2, 4, 8 进程的加速效果。

用法:
  julia --project=EMSuite EMSuite/benchmark/run_mpi_sweep.jl
"""

using Printf

println("="^70)
println("  MPI 并行性能扫描 — 自动化测试")
println("="^70)

# 测试进程数列表
proc_counts = [1, 2, 4]

# 检查 mpiexecjl 是否可用
mpiexec_cmd = "mpiexecjl"
try
    run(`$mpiexec_cmd --version`)
    println("✓ mpiexecjl found")
catch
    @warn "mpiexecjl not found. Trying system mpiexec..."
    mpiexec_cmd = "mpiexec"
    try
        run(`$mpiexec_cmd --version`)
        println("✓ mpiexec found")
    catch
        error("Neither mpiexecjl nor mpiexec found. Please install MPI.jl and configure.")
    end
end

println("\nRunning benchmarks for process counts: $proc_counts")
println("="^70)

for n in proc_counts
    println("\n[$(findfirst(==(n), proc_counts))/$(length(proc_counts))] Testing with $n processes...")
    
    cmd = `$mpiexec_cmd -n $n julia --project=EMSuite EMSuite/benchmark/benchmark_mpi_performance.jl`
    
    try
        run(cmd)
        println("✓ Completed n=$n")
    catch e
        @warn "Failed for n=$n: $e"
    end
    
    sleep(2)  # 短暂延迟避免资源竞争
end

println("\n" * "="^70)
println("  扫描完成。分析结果...")
println("="^70)

# ══════════════════════════════════════════════════════════════════════════
#  读取并分析结果
# ══════════════════════════════════════════════════════════════════════════
result_files = filter(f -> startswith(f, "benchmark_mpi_n"), readdir("."))

if isempty(result_files)
    @warn "No result files found."
    exit(1)
end

# 解析结果
results = Dict{Int, Dict{String, Float64}}()

for file in result_files
    lines = readlines(file)
    n_procs = 0
    t_assembly = 0.0
    test_name = ""
    
    for line in lines
        if startswith(line, "MPI_processes:")
            n_procs = parse(Int, split(line, ":")[2])
        elseif startswith(line, "Assembly_time:")
            t_assembly = parse(Float64, split(line, ":")[2])
        elseif startswith(line, "Test:")
            test_name = strip(split(line, ":")[2])
        end
    end
    
    if n_procs > 0 && t_assembly > 0.0
        if !haskey(results, n_procs)
            results[n_procs] = Dict{String, Float64}()
        end
        results[n_procs]["Assembly"] = t_assembly
    end
end

# 排序并打印
sorted_procs = sort(collect(keys(results)))

println("\n" * "="^70)
println("  MPI 并行加速结果")
println("="^70)
println()
@printf("%-12s | %-15s | %-10s | %-12s\n", "Processes", "Assembly (s)", "Speedup", "Efficiency")
println("-"^70)

t_serial = get(results[sorted_procs[1]], "Assembly", 0.0)

for n in sorted_procs
    t = get(results[n], "Assembly", 0.0)
    speedup = t_serial / t
    efficiency = speedup / n * 100
    
    @printf("%-12d | %-15.4f | %-10.2f× | %-11.1f%%\n", n, t, speedup, efficiency)
end

println("="^70)
println("\n关键指标:")
if length(sorted_procs) >= 2
    n_max = maximum(sorted_procs)
    t_max = results[n_max]["Assembly"]
    speedup_max = t_serial / t_max
    
    if speedup_max >= 0.8 * n_max
        status = "✅ GOOD"
    elseif speedup_max >= 0.5 * n_max
        status = "⚠️ FAIR"
    else
        status = "❌ POOR"
    end
    
    @printf("  最大进程数 (%d): 加速比 %.2f× / 理想 %d× = %.1f%% 效率 %s\n",
            n_max, speedup_max, n_max, speedup_max/n_max*100, status)
end

println("\n完成！结果文件: $result_files")
