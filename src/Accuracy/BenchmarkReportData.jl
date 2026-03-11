using CSV
using DataFrames
using Statistics

"""
    AccuracyCurve

统一表示一条远场对比曲线。
"""
struct AccuracyCurve
    label::String
    theta_deg::Vector{Float64}
    rcs_model_dBsm::Vector{Float64}
    rcs_reference_dBsm::Vector{Float64}
    diff_dB::Vector{Float64}
    model_column::String
    reference_column::String
    source_path::String
end

"""
    AccuracyCurveSummary

远场曲线的摘要统计，用于总报告汇总表。
"""
struct AccuracyCurveSummary
    label::String
    rmse_dB::Float64
    max_err_dB::Float64
    mean_bias_dB::Float64
    n_points::Int
end

"""
    PerformanceBenchmarkResult

统一表示一条性能基线结果记录。
"""
struct PerformanceBenchmarkResult
    case_name::String
    equation::String
    solver::String
    N::Int
    t_mesh::Float64
    t_assembly::Float64
    t_precond::Float64
    t_solve::Float64
    t_rcs::Float64
    t_total::Float64
    notes::String
end

function _require_column(df, candidates::Vector{String}, kind::String)
    for candidate in candidates
        if candidate in names(df)
            return candidate
        end
    end
    throw(ArgumentError("缺少$(kind)列，候选列: $(join(candidates, ", "))"))
end

"""
    load_accuracy_curve(csv_path; label=nothing)

读取精度对比 CSV，兼容当前报告脚本中的列命名差异。
"""
function load_accuracy_curve(csv_path::String; label::Union{Nothing, String} = nothing)
    isfile(csv_path) || throw(ArgumentError("精度 CSV 不存在: $csv_path"))

    df = DataFrames.DataFrame(CSV.File(csv_path))
    theta_col = _require_column(df, ["theta_deg"], "theta")
    model_col = _require_column(df, ["rcs_emsuite_dBsm", "rcs_ems_dBsm"], "EMSuite 曲线")
    ref_col = _require_column(df, ["rcs_feko_dBsm", "rcs_mie_dBsm", "rcs_ref_dBsm"], "参考曲线")
    diff_col = _require_column(df, ["diff_dB"], "误差")

    curve_label = isnothing(label) ? splitext(basename(csv_path))[1] : label
    return AccuracyCurve(
        curve_label,
        Float64.(df[!, theta_col]),
        Float64.(df[!, model_col]),
        Float64.(df[!, ref_col]),
        Float64.(df[!, diff_col]),
        model_col,
        ref_col,
        csv_path,
    )
end

"""
    summarize_accuracy_curve(curve)

对单条远场曲线计算 RMSE、最大误差和均值偏差。
"""
function summarize_accuracy_curve(curve::AccuracyCurve)
    rmse = sqrt(sum(abs2, curve.diff_dB) / length(curve.diff_dB))
    return AccuracyCurveSummary(
        curve.label,
        rmse,
        maximum(abs.(curve.diff_dB)),
        mean(curve.diff_dB),
        length(curve.diff_dB),
    )
end

"""
    accuracy_curve_group(label)

返回图表分组名。`*_phi0_*` 和 `*_phi90_*` 会归入同一组。
"""
function accuracy_curve_group(label::String)
    return replace(label, r"_phi(?:0|90)_vs_.*$" => "")
end

"""
    accuracy_curve_cut(label)

从曲线标签提取切面名，当前支持 `phi0` / `phi90`。
"""
function accuracy_curve_cut(label::String)
    match_obj = match(r"_(phi(?:0|90))_vs_", label)
    return isnothing(match_obj) ? "full" : match_obj.captures[1]
end

"""
    load_performance_results(csv_path)

读取性能基线 CSV。
"""
function load_performance_results(csv_path::String)
    isfile(csv_path) || throw(ArgumentError("性能 CSV 不存在: $csv_path"))

    df = DataFrames.DataFrame(CSV.File(csv_path))
    required = [
        "case_name", "equation", "solver", "N",
        "t_mesh", "t_assembly", "t_precond", "t_solve", "t_rcs", "t_total", "notes",
    ]
    for name in required
        name in names(df) || throw(ArgumentError("性能 CSV 缺少列: $name"))
    end

    return [
        PerformanceBenchmarkResult(
            String(row.case_name),
            String(row.equation),
            String(row.solver),
            Int(row.N),
            Float64(row.t_mesh),
            Float64(row.t_assembly),
            Float64(row.t_precond),
            Float64(row.t_solve),
            Float64(row.t_rcs),
            Float64(row.t_total),
            ismissing(row.notes) ? "" : String(row.notes),
        ) for row in eachrow(df)
    ]
end