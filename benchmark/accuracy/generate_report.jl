"""
generate_report.jl — Phase 14 全量精度报告汇总

此脚本汇总所有 CSV 结果文件（由 run_F*.jl 和 run_P*.jl 生成），
合并生成 ACCURACY_REPORT.md。

运行方式:
  julia --project=. benchmark/accuracy/generate_report.jl
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Printf, Dates

const ROOT_DIR   = joinpath(@__DIR__, "..", "..")
const RESULT_DIR = joinpath(ROOT_DIR, "test_results", "accuracy")

function read_result_csv(csv_path)
    isfile(csv_path) || return nothing
    lines = readlines(csv_path)
    length(lines) < 2 && return nothing
    header = split(lines[1], ",")
    θ = Float64[]; ems = Float64[]; ref = Float64[]; diff = Float64[]
    for l in lines[2:end]
        isempty(strip(l)) && continue
        parts = split(l, ",")
        length(parts) >= 4 || continue
        push!(θ,   parse(Float64, parts[1]))
        push!(ems, parse(Float64, parts[2]))
        push!(ref, parse(Float64, parts[3]))
        push!(diff, parse(Float64, parts[4]))
    end
    isempty(θ) && return nothing
    return (theta=θ, rcs_ems=ems, rcs_ref=ref, diff=diff)
end

function compute_stats(data)
    isnothing(data) && return nothing
    d = data.diff
    rmse = sqrt(sum(d .^ 2) / length(d))
    maxe = maximum(abs.(d))
    bias = sum(d) / length(d)
    # 后向散射（θ≈180°）
    _, ib = findmin(abs.(data.theta .- 180.0))
    bs = abs(d[ib])
    return (rmse=rmse, max_err=maxe, bias=bias, bs_err=bs, n=length(d))
end

# ─── 用例注册表 ─────────────────────────────────────────────────────────────
# (label_prefix, description, threshold)
const EXPECTED_CASES = [
    # F1
    ("F1_SEFIE_Jet_Direct_phi0_vs_Feko",   "F1 S-EFIE Direct φ=0°",   2.0),
    ("F1_SEFIE_Jet_Direct_phi90_vs_Feko",  "F1 S-EFIE Direct φ=90°",  2.0),
    # F2
    ("F2_CFIE_Jet_Direct_phi0_vs_Feko",    "F2 S-CFIE Direct φ=0°",   2.0),
    ("F2_CFIE_Jet_Direct_phi90_vs_Feko",   "F2 S-CFIE Direct φ=90°",  2.0),
    # F3
    ("F3_SEFIE_Jet_MLFMA_phi0_vs_Feko",   "F3 S-EFIE MLFMA φ=0°",    3.0),
    ("F3_SEFIE_Jet_MLFMA_phi90_vs_Feko",  "F3 S-EFIE MLFMA φ=90°",   3.0),
    # F4
    ("F4_CFIE_Jet_MLFMA_phi0_vs_Feko",    "F4 S-CFIE MLFMA φ=0°",    3.0),
    ("F4_CFIE_Jet_MLFMA_phi90_vs_Feko",   "F4 S-CFIE MLFMA φ=90°",   3.0),
    # F5
    ("F5_CFIE_Sphere_Direct_phi0_vs_Feko", "F5 S-CFIE Sphere Direct φ=0°", 2.0),
    ("F5_CFIE_Sphere_Direct_phi90_vs_Feko","F5 S-CFIE Sphere Direct φ=90°",2.0),
    ("F5_CFIE_Sphere_Direct_phi0_vs_Mie",  "F5 S-CFIE Sphere Direct vs Mie", 2.0),
    # F6
    ("F6_CFIE_Sphere_MLFMA_phi0_vs_Feko",  "F6 S-CFIE Sphere MLFMA φ=0°",  3.0),
    ("F6_CFIE_Sphere_MLFMA_phi90_vs_Feko", "F6 S-CFIE Sphere MLFMA φ=90°", 3.0),
    ("F6_CFIE_Sphere_MLFMA_phi0_vs_Mie",   "F6 S-CFIE Sphere MLFMA vs Mie",3.0),
    # F7
    ("F7_VEFIE_Plate_Direct_phi0_vs_Feko",  "F7 V-EFIE Plate φ=0°",  2.0),
    ("F7_VEFIE_Plate_Direct_phi90_vs_Feko", "F7 V-EFIE Plate φ=90°", 2.0),
    # F9
    ("F9_SCFIE_PlateMetal_Direct_phi0_vs_Feko",  "F9 SCFIE Plate+Metal φ=0°",  2.0),
    ("F9_SCFIE_PlateMetal_Direct_phi90_vs_Feko", "F9 SCFIE Plate+Metal φ=90°", 2.0),
    # P1
    ("P1_PMCHW_Sphere_Direct_phi0_vs_Mie",  "P1 PMCHW εᵣ=4 φ=0°",  2.0),
    ("P1_PMCHW_Sphere_Direct_phi90_vs_Mie", "P1 PMCHW εᵣ=4 φ=90°", 2.0),
    # P3
    ("P3_PMCHW_LossySphere_Direct_phi0_vs_Mie",  "P3 PMCHW (有损) φ=0°",  2.0),
    ("P3_PMCHW_LossySphere_Direct_phi90_vs_Mie", "P3 PMCHW (有损) φ=90°", 2.0),
    # X1
    ("X1_SEFIE_Sphere_Direct_phi0_vs_Mie",  "X1 S-EFIE Sphere vs Mie φ=0°",  2.0),
    ("X1_SEFIE_Sphere_Direct_phi90_vs_Mie", "X1 S-EFIE Sphere vs Mie φ=90°", 2.0),
]

mkpath(RESULT_DIR)
report_path = joinpath(RESULT_DIR, "ACCURACY_REPORT.md")

open(report_path, "w") do io
    ts = Dates.format(now(), "yyyy-mm-dd HH:MM")
    println(io, "# EMSuite Phase 14 全量精度对比报告")
    println(io, "")
    println(io, "- **生成时间**: $ts")
    println(io, "- **基准**: Feko 商业软件 / Mie 解析级数")
    println(io, "- **门限**: Direct ≤ 2 dB RMSE；MLFMA ≤ 3 dB RMSE")
    println(io, "")
    println(io, "---")
    println(io, "")

    # 生成汇总表
    println(io, "## 精度汇总表\n")
    println(io, "| 用例 | 对比基准 | N点 | RMSE(dB) | MaxErr(dB) | Bias(dB) | BS_Err(dB) | 门限 | 结论 |")
    println(io, "|------|---------|-----|---------|----------|--------|----------|------|------|")

    n_pass = 0; n_total = 0; n_missing = 0

    for (fname, desc, threshold) in EXPECTED_CASES
        csv_path = joinpath(RESULT_DIR, fname * ".csv")
        data = read_result_csv(csv_path)
        stats = compute_stats(data)
        n_total += 1

        if isnothing(stats)
            n_missing += 1
            println(io, "| $desc | — | — | — | — | — | — | $(threshold) dB | ⚠ 未运行 |")
            continue
        end

        pass = stats.rmse < threshold
        pass && (n_pass += 1)
        status = pass ? "✓ PASS" : "✗ FAIL"
        @printf(io, "| %s | Feko/Mie | %d | %.3f | %.3f | %+.3f | %.3f | %.1f dB | %s |\n",
            desc, stats.n, stats.rmse, stats.max_err, stats.bias, stats.bs_err,
            threshold, status)
    end

    println(io, "")
    println(io, "**总通过率**: $n_pass / $(n_total - n_missing) ($(n_missing) 项未运行)")
    println(io, "")
    println(io, "---")
    println(io, "")
    println(io, "## 备注")
    println(io, "")
    println(io, "- **F8**: 纯介质板 SCFIE，网格不含 CTRIA3，暂跳过")
    println(io, "- **P2**: PMCHW MLFMA，`PMCHWMLFMAOperator` 待实现")
    println(io, "- **A1-A4**: 偶极子天线端口，运行 `run_A1_A4_antenna.jl`")
end

println("报告已保存: $report_path")

# 控制台摘要
let results_found = 0, results_pass = 0
    for (fname, desc, threshold) in EXPECTED_CASES
        csv_path = joinpath(RESULT_DIR, fname * ".csv")
        data  = read_result_csv(csv_path)
        stats = compute_stats(data)
        isnothing(stats) && continue
        results_found += 1
        stats.rmse < threshold && (results_pass += 1)
        status = stats.rmse < threshold ? "✓" : "✗"
        @printf("  %s %-46s RMSE=%.3f dB\n", status, desc, stats.rmse)
    end
    println("\n已有结果: $results_found / $(length(EXPECTED_CASES))  通过: $results_pass")
end
