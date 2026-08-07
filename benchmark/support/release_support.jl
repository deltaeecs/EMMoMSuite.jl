module ReleaseSupport

using Dates
using Statistics
using TOML

export AccuracyCurve,
    AccuracyCurveSummary,
    PerformanceBenchmarkResult,
    KnownException,
    RunArtifactPaths,
    ReleaseCaseStatus,
    parse_csv_rows,
    csv_row_has,
    csv_row_string,
    csv_row_float,
    csv_row_int,
    csv_row_bool,
    load_accuracy_curve,
    summarize_accuracy_curve,
    accuracy_curve_group,
    accuracy_curve_cut,
    load_performance_results,
    generate_run_id,
    create_run_artifact_dirs,
    load_release_profile,
    load_known_exceptions,
    match_known_exception,
    default_accuracy_threshold,
    accuracy_case_status,
    performance_case_status,
    build_case_statuses,
    write_run_manifest,
    write_case_status_csv,
    write_artifact_index_csv,
    mirror_tree

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

struct AccuracyCurveSummary
    label::String
    rmse_dB::Float64
    max_err_dB::Float64
    mean_bias_dB::Float64
    n_points::Int
end

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

struct KnownException
    id::String
    label_pattern::String
    scope::String
    rationale::String
    disposition::String
    tracking::String
end

struct RunArtifactPaths
    run_id::String
    root_dir::String
    manifest_dir::String
    logs_dir::String
    metrics_dir::String
    metrics_accuracy_dir::String
    metrics_performance_dir::String
    metrics_summary_dir::String
    plots_dir::String
    plots_accuracy_dir::String
    plots_accuracy_polar_dir::String
    plots_performance_dir::String
    report_dir::String
end

struct ReleaseCaseStatus
    case_name::String
    category::String
    status::String
    metric_name::String
    metric_value::Float64
    threshold::Float64
    note::String
end

_split_csv_line(line::String) = split(chomp(line), ','; keepempty = true)

function parse_csv_rows(csv_path::String)
    isfile(csv_path) || throw(ArgumentError("CSV 不存在: $csv_path"))
    lines = readlines(csv_path)
    isempty(lines) && return Dict{String, String}[]
    header = replace.(String.(_split_csv_line(first(lines))), '\ufeff' => "")
    rows = Dict{String, String}[]
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        values = _split_csv_line(line)
        if length(values) < length(header)
            append!(values, fill("", length(header) - length(values)))
        end
        push!(rows, Dict(header[i] => values[i] for i in eachindex(header)))
    end
    return rows
end

csv_row_has(row::AbstractDict, key::AbstractString) = haskey(row, String(key))
csv_row_string(row::AbstractDict, key::AbstractString, default::String = "") = get(row, String(key), default)

function csv_row_float(row::AbstractDict, key::AbstractString, default::Float64 = NaN)
    value = strip(csv_row_string(row, key, ""))
    return isempty(value) ? default : parse(Float64, value)
end

function csv_row_int(row::AbstractDict, key::AbstractString, default::Int = 0)
    value = strip(csv_row_string(row, key, ""))
    return isempty(value) ? default : parse(Int, value)
end

function csv_row_bool(row::AbstractDict, key::AbstractString, default::Bool = false)
    value = lowercase(strip(csv_row_string(row, key, "")))
    return isempty(value) ? default : value in ("true", "1", "yes")
end

function _require_column(row::AbstractDict, candidates::Vector{String}, kind::String)
    for candidate in candidates
        csv_row_has(row, candidate) && return candidate
    end
    throw(ArgumentError("缺少$(kind)列，候选列: $(join(candidates, ", "))"))
end

function load_accuracy_curve(csv_path::String; label::Union{Nothing, String} = nothing)
    rows = parse_csv_rows(csv_path)
    isempty(rows) && throw(ArgumentError("精度 CSV 为空: $csv_path"))

    sample = first(rows)
    theta_col = _require_column(sample, ["theta_deg"], "theta")
    model_col = _require_column(sample, ["rcs_emsuite_dBsm", "rcs_ems_dBsm"], "EMMoMSuite 曲线")
    ref_col = _require_column(sample, ["rcs_feko_dBsm", "rcs_mie_dBsm", "rcs_ref_dBsm"], "参考曲线")
    diff_col = _require_column(sample, ["diff_dB"], "误差")

    curve_label = isnothing(label) ? splitext(basename(csv_path))[1] : label
    return AccuracyCurve(
        curve_label,
        [csv_row_float(row, theta_col) for row in rows],
        [csv_row_float(row, model_col) for row in rows],
        [csv_row_float(row, ref_col) for row in rows],
        [csv_row_float(row, diff_col) for row in rows],
        model_col,
        ref_col,
        csv_path,
    )
end

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

accuracy_curve_group(label::String) = replace(label, r"_phi(?:0|90)_vs_.*$" => "")

function accuracy_curve_cut(label::String)
    match_obj = match(r"_(phi(?:0|90))_vs_", label)
    return isnothing(match_obj) ? "full" : match_obj.captures[1]
end

function load_performance_results(csv_path::String)
    rows = parse_csv_rows(csv_path)
    isempty(rows) && return PerformanceBenchmarkResult[]

    required = [
        "case_name", "equation", "solver", "N",
        "t_mesh", "t_assembly", "t_precond", "t_solve", "t_rcs", "t_total", "notes",
    ]
    for name in required
        csv_row_has(first(rows), name) || throw(ArgumentError("性能 CSV 缺少列: $name"))
    end

    return [
        PerformanceBenchmarkResult(
            csv_row_string(row, "case_name"),
            csv_row_string(row, "equation"),
            csv_row_string(row, "solver"),
            csv_row_int(row, "N"),
            csv_row_float(row, "t_mesh"),
            csv_row_float(row, "t_assembly"),
            csv_row_float(row, "t_precond"),
            csv_row_float(row, "t_solve"),
            csv_row_float(row, "t_rcs"),
            csv_row_float(row, "t_total"),
            csv_row_string(row, "notes"),
        ) for row in rows
    ]
end

generate_run_id(now_utc::DateTime = now(UTC)) = Dates.format(now_utc, dateformat"yyyymmdd_HHMMSS")

function create_run_artifact_dirs(runs_root::String; run_id::Union{Nothing, String} = nothing)
    actual_run_id = isnothing(run_id) ? generate_run_id() : run_id
    root_dir = joinpath(runs_root, actual_run_id)
    paths = RunArtifactPaths(
        actual_run_id,
        root_dir,
        joinpath(root_dir, "manifest"),
        joinpath(root_dir, "logs"),
        joinpath(root_dir, "metrics"),
        joinpath(root_dir, "metrics", "accuracy"),
        joinpath(root_dir, "metrics", "performance"),
        joinpath(root_dir, "metrics", "summary"),
        joinpath(root_dir, "plots"),
        joinpath(root_dir, "plots", "accuracy"),
        joinpath(root_dir, "plots", "accuracy_polar"),
        joinpath(root_dir, "plots", "performance"),
        joinpath(root_dir, "report"),
    )

    for dir in (
        paths.root_dir,
        paths.manifest_dir,
        paths.logs_dir,
        paths.metrics_dir,
        paths.metrics_accuracy_dir,
        paths.metrics_performance_dir,
        paths.metrics_summary_dir,
        paths.plots_dir,
        paths.plots_accuracy_dir,
        paths.plots_accuracy_polar_dir,
        paths.plots_performance_dir,
        paths.report_dir,
    )
        mkpath(dir)
    end

    return paths
end

function load_release_profile(config_path::String)
    isfile(config_path) || throw(ArgumentError("workflow 配置不存在: $config_path"))
    return TOML.parsefile(config_path)
end

function load_known_exceptions(config_path::String)
    isfile(config_path) || return KnownException[]

    parsed = TOML.parsefile(config_path)
    raw_entries = get(parsed, "exceptions", Any[])
    return [
        KnownException(
            String(get(entry, "id", "")),
            String(get(entry, "label_pattern", "")),
            String(get(entry, "scope", "")),
            String(get(entry, "rationale", "")),
            String(get(entry, "disposition", "")),
            String(get(entry, "tracking", "")),
        ) for entry in raw_entries
    ]
end

function match_known_exception(label::String, exceptions::Vector{KnownException})
    for exception in exceptions
        isempty(exception.label_pattern) && continue
        occursin(exception.label_pattern, label) && return exception
    end
    return nothing
end

default_accuracy_threshold(label::String) = occursin("MLFMA", label) ? 3.0 : 2.5

function accuracy_case_status(curve::AccuracyCurve, exceptions::Vector{KnownException})
    summary = summarize_accuracy_curve(curve)
    threshold = default_accuracy_threshold(curve.label)
    matched = match_known_exception(curve.label, exceptions)
    status = summary.rmse_dB <= threshold ? "PASS" : isnothing(matched) ? "FAIL" : "KNOWN_EXCEPTION"
    note = isnothing(matched) ? "" : matched.rationale
    return ReleaseCaseStatus(curve.label, "accuracy", status, "rmse_dB", summary.rmse_dB, threshold, note)
end

performance_case_status(row::PerformanceBenchmarkResult) = ReleaseCaseStatus(row.case_name, "performance", "RECORDED", "t_total_s", row.t_total, 0.0, row.notes)

function build_case_statuses(curves::Vector{AccuracyCurve}, rows::Vector{PerformanceBenchmarkResult}, exceptions::Vector{KnownException})
    statuses = ReleaseCaseStatus[]
    append!(statuses, accuracy_case_status.(curves, Ref(exceptions)))
    append!(statuses, performance_case_status.(rows))
    return statuses
end

function write_run_manifest(
    manifest_path::String;
    run_id::String,
    workspace_root::String,
    profile_path::String,
    plot_style_path::String,
    thresholds_path::String,
    known_exceptions_path::String,
    profile::Dict,
)
    manifest = Dict(
        "run" => Dict(
            "id" => run_id,
            "generated_at_utc" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
            "workspace_root" => workspace_root,
        ),
        "config" => Dict(
            "profile_path" => profile_path,
            "plot_style_path" => plot_style_path,
            "thresholds_path" => thresholds_path,
            "known_exceptions_path" => known_exceptions_path,
        ),
        "profile" => profile,
        "environment" => Dict(
            "julia_version" => string(VERSION),
            "threads" => Threads.nthreads(),
            "os" => Sys.iswindows() ? "Windows" : Sys.islinux() ? "Linux" : Sys.isapple() ? "macOS" : "unknown",
        ),
    )

    open(manifest_path, "w") do io
        TOML.print(io, manifest)
    end
    return manifest_path
end

function write_case_status_csv(csv_path::String, statuses::Vector{ReleaseCaseStatus})
    open(csv_path, "w") do io
        println(io, "case_name,category,status,metric_name,metric_value,threshold,note")
        for status in statuses
            note = replace(status.note, ',' => ';')
            println(io,
                string(
                    status.case_name, ",",
                    status.category, ",",
                    status.status, ",",
                    status.metric_name, ",",
                    status.metric_value, ",",
                    status.threshold, ",",
                    note,
                ),
            )
        end
    end
    return csv_path
end

function write_artifact_index_csv(csv_path::String, root_dir::String)
    files = String[]
    for (dir, _, names) in walkdir(root_dir)
        for name in names
            path = joinpath(dir, name)
            path == csv_path && continue
            push!(files, path)
        end
    end
    sort!(files)

    open(csv_path, "w") do io
        println(io, "relative_path,bytes")
        for path in files
            rel = replace(relpath(path, root_dir), '\\' => '/')
            println(io, string(rel, ",", filesize(path)))
        end
    end
    return csv_path
end

function mirror_tree(source_dir::String, destination_dir::String)
    isdir(source_dir) || return destination_dir
    mkpath(destination_dir)
    for entry in readdir(source_dir; join = true)
        target = joinpath(destination_dir, basename(entry))
        if isdir(entry)
            mirror_tree(entry, target)
        else
            cp(entry, target; force = true)
        end
    end
    return destination_dir
end

end