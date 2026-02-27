"""
    precision_comparison.jl

Phase 10 精度对齐: 全角度 RCS 定量对比分析。
读取现有 Legacy 基线和 EMSuite 验证结果，计算全面统计指标。

用法: julia --project=. scripts/precision_comparison.jl
"""

using DelimitedFiles
using Printf
using Statistics
using LinearAlgebra

# ============================================================
#  工具函数
# ============================================================

"""读取 Legacy 基线 CSV (theta_rad, theta_deg, phi0_dB, phi90_dB)"""
function read_legacy_csv(filepath::String)
    lines = readlines(filepath)
    header = split(lines[1], ',')
    n = length(lines) - 1
    θ_rad = Vector{Float64}(undef, n)
    θ_deg = Vector{Float64}(undef, n)
    rcs_phi0 = Vector{Float64}(undef, n)
    rcs_phi90 = Vector{Float64}(undef, n)
    for i in 1:n
        fields = split(lines[i+1], ',')
        θ_rad[i] = parse(Float64, fields[1])
        θ_deg[i] = parse(Float64, fields[2])
        rcs_phi0[i] = parse(Float64, fields[3])
        rcs_phi90[i] = parse(Float64, fields[4])
    end
    return θ_deg, rcs_phi0, rcs_phi90
end

"""读取 EMSuite 验证 CSV (theta_rad, phi0_dB, phi90_dB)"""
function read_emsuite_csv(filepath::String)
    lines = readlines(filepath)
    header = split(lines[1], ',')
    n = length(lines) - 1
    
    # 自动检测列: 有些文件第一列是 Theta_Deg, 有些是 Theta_Rad
    has_deg = any(occursin("Deg", h) for h in header)
    
    θ_deg = Vector{Float64}(undef, n)
    rcs_phi0 = Vector{Float64}(undef, n)
    rcs_phi90 = Vector{Float64}(undef, n)
    for i in 1:n
        fields = split(lines[i+1], ',')
        θ_val = parse(Float64, fields[1])
        θ_deg[i] = has_deg ? θ_val : rad2deg(θ_val)
        rcs_phi0[i] = parse(Float64, fields[2])
        rcs_phi90[i] = parse(Float64, fields[3])
    end
    return θ_deg, rcs_phi0, rcs_phi90
end

# ============================================================
#  统计指标
# ============================================================

"""计算全面统计指标 (dB 域比较)"""
function compute_metrics(legacy_dB::Vector{Float64}, emsuite_dB::Vector{Float64}, label::String)
    diff = emsuite_dB .- legacy_dB
    abs_diff = abs.(diff)
    
    metrics = Dict{String, Float64}()
    metrics["mean_diff"]   = mean(diff)           # 平均偏差 (系统性偏移)
    metrics["std_diff"]    = std(diff)             # 偏差标准差
    metrics["rmse"]        = sqrt(mean(diff.^2))   # 均方根误差
    metrics["max_abs"]     = maximum(abs_diff)     # 最大绝对偏差
    metrics["median_abs"]  = median(abs_diff)      # 中位数绝对偏差
    metrics["pct95"]       = quantile(abs_diff, 0.95) # 95% 分位数
    metrics["pct99"]       = quantile(abs_diff, 0.99) # 99% 分位数
    
    # 合格率 (误差 < 阈值的比例)
    metrics["pass_0.5dB"] = count(d -> d < 0.5, abs_diff) / length(abs_diff) * 100
    metrics["pass_1.0dB"] = count(d -> d < 1.0, abs_diff) / length(abs_diff) * 100
    metrics["pass_2.0dB"] = count(d -> d < 2.0, abs_diff) / length(abs_diff) * 100
    
    return metrics
end

"""打印统计报告"""
function print_report(metrics::Dict, label::String)
    println("  ┌─ $label")
    @printf("  │  平均偏差 (dB):    %+.4f  (系统偏移)\n", metrics["mean_diff"])
    @printf("  │  标准差 (dB):      %.4f\n", metrics["std_diff"])
    @printf("  │  RMSE (dB):        %.4f\n", metrics["rmse"])
    @printf("  │  最大绝对差 (dB):  %.4f\n", metrics["max_abs"])
    @printf("  │  中位数差 (dB):    %.4f\n", metrics["median_abs"])
    @printf("  │  P95 (dB):         %.4f\n", metrics["pct95"])
    @printf("  │  P99 (dB):         %.4f\n", metrics["pct99"])
    @printf("  │  < 0.5 dB 比例:    %.1f%%\n", metrics["pass_0.5dB"])
    @printf("  │  < 1.0 dB 比例:    %.1f%%\n", metrics["pass_1.0dB"])
    @printf("  │  < 2.0 dB 比例:    %.1f%%\n", metrics["pass_2.0dB"])
    println("  └─")
end

"""分段统计 (前向/侧向/后向)"""
function segmented_analysis(θ_deg, diff, label)
    # 前向: θ ∈ [-30°, 30°] 
    # 侧向: θ ∈ [-90°,-30°] ∪ [30°,90°]
    # 后向: θ ∈ [-180°,-90°) ∪ (90°,180°]
    
    forward_mask = abs.(θ_deg) .<= 30.0
    side_mask = (abs.(θ_deg) .> 30.0) .& (abs.(θ_deg) .<= 90.0)
    back_mask = abs.(θ_deg) .> 90.0
    
    function seg_stats(mask, name)
        d = diff[mask]
        if isempty(d)
            println("      $name: 无数据")
            return
        end
        @printf("      %-8s  N=%3d  RMSE=%.3f  Max=%.3f  Mean=%+.3f dB\n",
                name, count(mask), sqrt(mean(d.^2)), maximum(abs.(d)), mean(d))
    end
    
    println("    分段分析 ($label):")
    seg_stats(forward_mask, "前向")
    seg_stats(side_mask,    "侧向")
    seg_stats(back_mask,    "后向")
end

# ============================================================
#  判定标准
# ============================================================

"""
    pass_criteria(metrics) -> (passed::Bool, reason::String)

Phase 10 通过标准:
  - RMSE < 1.0 dB
  - 最大绝对差 < 3.0 dB  
  - >90% 采样点误差 < 1.0 dB
  - 系统偏移 |mean_diff| < 0.5 dB
"""
function pass_criteria(metrics::Dict)
    reasons = String[]
    passed = true
    
    if metrics["rmse"] >= 1.0
        push!(reasons, @sprintf("RMSE %.3f ≥ 1.0 dB", metrics["rmse"]))
        passed = false
    end
    if metrics["max_abs"] >= 3.0
        push!(reasons, @sprintf("Max %.3f ≥ 3.0 dB", metrics["max_abs"]))
        passed = false
    end
    if metrics["pass_1.0dB"] < 90.0
        push!(reasons, @sprintf("<1dB合格率 %.1f%% < 90%%", metrics["pass_1.0dB"]))
        passed = false
    end
    if abs(metrics["mean_diff"]) >= 0.5
        push!(reasons, @sprintf("|MeanDiff| %.3f ≥ 0.5 dB", abs(metrics["mean_diff"])))
        passed = false
    end
    
    return passed, isempty(reasons) ? "ALL PASS" : join(reasons, "; ")
end

# ============================================================
#  A1: S-EFIE Direct Jet
# ============================================================

function compare_A1_sefie_direct()
    println("\n" * "="^72)
    println("  A1: S-EFIE Direct — Jet 100MHz PEC (开体)")
    println("="^72)
    
    legacy_file  = joinpath(@__DIR__, "..", "test_results", "legacy_baseline", "SEFIE_Direct_Jet.csv")
    emsuite_file = joinpath(@__DIR__, "..", "test_results", "emsuite_verification", "SEFIE_Direct_Jet_PostFix.csv")
    
    if !isfile(legacy_file)
        println("  ❌ Legacy baseline not found: $legacy_file")
        return nothing
    end
    if !isfile(emsuite_file)
        println("  ❌ EMSuite results not found: $emsuite_file")
        return nothing
    end
    
    θ_legacy, legacy_phi0, legacy_phi90 = read_legacy_csv(legacy_file)
    θ_emsuite, emsuite_phi0, emsuite_phi90 = read_emsuite_csv(emsuite_file)
    
    @assert length(θ_legacy) == length(θ_emsuite) "采样点数不一致: Legacy=$(length(θ_legacy)), EMSuite=$(length(θ_emsuite))"
    
    println("\n  采样点: $(length(θ_legacy)), θ ∈ [$(minimum(θ_legacy))°, $(maximum(θ_legacy))°]")
    
    println("\n  --- ϕ = 0° 切面 ---")
    m0 = compute_metrics(legacy_phi0, emsuite_phi0, "Phi=0°")
    print_report(m0, "Phi=0°")
    segmented_analysis(θ_legacy, emsuite_phi0 .- legacy_phi0, "Phi=0°")
    
    println("\n  --- ϕ = 90° 切面 ---")
    m90 = compute_metrics(legacy_phi90, emsuite_phi90, "Phi=90°")
    print_report(m90, "Phi=90°")
    segmented_analysis(θ_legacy, emsuite_phi90 .- legacy_phi90, "Phi=90°")
    
    # 合并两个切面
    all_legacy = vcat(legacy_phi0, legacy_phi90)
    all_emsuite = vcat(emsuite_phi0, emsuite_phi90)
    println("\n  --- 综合 (两个切面合并) ---")
    m_all = compute_metrics(all_legacy, all_emsuite, "综合")
    print_report(m_all, "综合")
    
    passed, reason = pass_criteria(m_all)
    println("  ★ 判定: $(passed ? "✅ PASS" : "❌ FAIL") — $reason")
    
    # 找到最大误差点
    diff_all = abs.(all_emsuite .- all_legacy)
    idx_max = argmax(diff_all)
    if idx_max <= length(θ_legacy)
        @printf("  ★ 最大误差点: θ=%.1f° (ϕ=0°), Δ=%.3f dB\n", θ_legacy[idx_max], diff_all[idx_max])
    else
        j = idx_max - length(θ_legacy)
        @printf("  ★ 最大误差点: θ=%.1f° (ϕ=90°), Δ=%.3f dB\n", θ_legacy[j], diff_all[idx_max])
    end
    
    return Dict("phi0" => m0, "phi90" => m90, "combined" => m_all, "passed" => passed)
end

# ============================================================
#  A3: S-EFIE MLFMA Jet (对比 Direct baseline)
# ============================================================

function compare_A3_sefie_mlfma()
    println("\n" * "="^72)
    println("  A3: S-EFIE MLFMA — Jet 100MHz (对比 Direct Legacy 基线)")
    println("="^72)
    
    legacy_file  = joinpath(@__DIR__, "..", "test_results", "legacy_baseline", "SEFIE_Direct_Jet.csv")
    emsuite_file = joinpath(@__DIR__, "..", "test_results", "emsuite_verification", "SEFIE_MLFMA_Jet_PostFix.csv")
    
    if !isfile(legacy_file)
        println("  ❌ Legacy baseline not found")
        return nothing
    end
    if !isfile(emsuite_file)
        println("  ❌ EMSuite MLFMA results not found")
        return nothing
    end
    
    θ_legacy, legacy_phi0, legacy_phi90 = read_legacy_csv(legacy_file)
    θ_emsuite, emsuite_phi0, emsuite_phi90 = read_emsuite_csv(emsuite_file)
    
    println("\n  采样点: $(length(θ_legacy))")
    
    println("\n  --- ϕ = 0° 切面 ---")
    m0 = compute_metrics(legacy_phi0, emsuite_phi0, "Phi=0°")
    print_report(m0, "Phi=0°")
    
    println("\n  --- ϕ = 90° 切面 ---")
    m90 = compute_metrics(legacy_phi90, emsuite_phi90, "Phi=90°")
    print_report(m90, "Phi=90°")
    
    all_legacy = vcat(legacy_phi0, legacy_phi90)
    all_emsuite = vcat(emsuite_phi0, emsuite_phi90)
    println("\n  --- 综合 ---")
    m_all = compute_metrics(all_legacy, all_emsuite, "综合")
    print_report(m_all, "综合")
    
    passed, reason = pass_criteria(m_all)
    println("  ★ 判定: $(passed ? "✅ PASS" : "❌ FAIL") — $reason")
    
    return Dict("phi0" => m0, "phi90" => m90, "combined" => m_all, "passed" => passed)
end

# ============================================================
#  C1: S-CFIE Direct (对比 SCFIE_Direct_Sphere Legacy 和 EMSuite SCFIE_Direct_Jet)
# ============================================================

function compare_C_scfie()
    println("\n" * "="^72)
    println("  C: S-CFIE — 可用数据对比")
    println("="^72)
    
    # 检查 EMSuite SCFIE Direct Jet 数据
    emsuite_jet_file = joinpath(@__DIR__, "..", "test_results", "emsuite_verification", "SCFIE_Direct_Jet.csv")
    emsuite_mlfma_file = joinpath(@__DIR__, "..", "test_results", "emsuite_verification", "SCFIE_MLFMA_Jet.csv")
    
    # 对比 SCFIE Direct vs MLFMA (EMSuite 内部一致性)
    if isfile(emsuite_jet_file) && isfile(emsuite_mlfma_file)
        println("\n  --- SCFIE Direct vs MLFMA (EMSuite 内部一致性) ---")
        θ_d, d0, d90 = read_emsuite_csv(emsuite_jet_file)
        θ_m, m0, m90 = read_emsuite_csv(emsuite_mlfma_file)
        
        if length(θ_d) == length(θ_m)
            println("  采样点: $(length(θ_d))")
            all_d = vcat(d0, d90)
            all_m = vcat(m0, m90)
            m_dm = compute_metrics(all_d, all_m, "Direct vs MLFMA")
            print_report(m_dm, "SCFIE Direct vs MLFMA")
        end
    end
    
    println("\n  注: SCFIE 完整精度验证需在 Sphere 网格上运行 EMSuite，生成与 Legacy SCFIE_Direct_Sphere 可比数据。")
end

# ============================================================
#  主入口
# ============================================================

function main()
    println("╔════════════════════════════════════════════════════════════════════════╗")
    println("║     Phase 10: 全角度 RCS 精度对齐报告                                 ║")
    println("║     EMSuite vs Legacy (MoM_AllinOne)                                  ║")
    println("╚════════════════════════════════════════════════════════════════════════╝")
    
    results = Dict{String, Any}()
    
    # A1: S-EFIE Direct
    results["A1"] = compare_A1_sefie_direct()
    
    # A3: S-EFIE MLFMA
    results["A3"] = compare_A3_sefie_mlfma()
    
    # C: S-CFIE
    compare_C_scfie()
    
    # 汇总
    println("\n" * "="^72)
    println("  汇总")
    println("="^72)
    for key in sort(collect(keys(results)))
        val = results[key]
        if val !== nothing && haskey(val, "passed")
            status = val["passed"] ? "✅ PASS" : "❌ FAIL"
            rmse = val["combined"]["rmse"]
            maxd = val["combined"]["max_abs"]
            @printf("    %s  %s  RMSE=%.3f dB  MaxDiff=%.3f dB\n", key, status, rmse, maxd)
        else
            println("    $key  ⚠️  数据不足")
        end
    end
end

main()
