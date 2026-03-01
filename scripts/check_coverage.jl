#!/usr/bin/env julia
# scripts/check_coverage.jl
# 运行测试覆盖率分析
#
# 用法:
#   julia --project=. scripts/check_coverage.jl
#
# 前置条件: Coverage.jl 需已全局安装
#   julia -e 'using Pkg; Pkg.add("Coverage")'
#
# 流程:
#   1. 使用 --code-coverage=user 重新运行测试
#   2. 用 Coverage.jl 处理 .cov 文件
#   3. 打印覆盖率摘要

project_dir = joinpath(@__DIR__, "..")

println("=" ^ 60)
println("  EMSuite.jl — 测试覆盖率检查")
println("=" ^ 60)

# Step 1: 运行测试（带 coverage 标志）
println("\n[1/3] 运行测试 (--code-coverage=user) ...")
test_cmd = `julia --project=$project_dir --threads=auto --code-coverage=user
             --startup-file=no $(joinpath(project_dir, "test", "runtests.jl"))`
run(test_cmd)

# Step 2: 处理覆盖率数据
println("\n[2/3] 处理 .cov 文件 ...")
try
    using Coverage
catch e
    println("Coverage.jl 未安装，正在安装 ...")
    import Pkg
    Pkg.add("Coverage")
    using Coverage
end

coverage = process_folder(joinpath(project_dir, "src"))

# Step 3: 打印摘要
covered_lines, total_lines = get_summary(coverage)
pct = 100.0 * covered_lines / total_lines

println("\n[3/3] 覆盖率摘要")
println("-" ^ 40)
@printf("  总行数:   %6d\n", total_lines)
@printf("  已覆盖:   %6d\n", covered_lines)
@printf("  未覆盖:   %6d\n", total_lines - covered_lines)
@printf("  覆盖率:   %6.2f%%\n", pct)
println("-" ^ 40)

if pct >= 80.0
    println("  ✅ 覆盖率目标达成 (≥ 80%)")
else
    println("  ⚠️  覆盖率不足 (< 80%)，需要添加更多测试")
end

# 打印覆盖率最低的 10 个文件
println("\n📋 覆盖率最低的 10 个文件:")
file_coverages = [(f.filename,
    let (c, t) = get_summary([f]); t > 0 ? 100.0 * c / t : 100.0 end)
                  for f in coverage]
sort!(file_coverages, by = x -> x[2])
for (fname, fc) in file_coverages[1:min(10, end)]
    rel = relpath(fname, project_dir)
    @printf("  %5.1f%%  %s\n", fc, rel)
end

# 清理 .cov 文件（可选）
# clean_folder(joinpath(project_dir, "src"))
println("\n完成。.cov 文件保留在 src/ 中 (可手动删除)。")
