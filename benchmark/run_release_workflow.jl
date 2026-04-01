"""
run_release_workflow.jl

Phase 17 统一发布编排入口：
  - 读取 TOML profile
  - 建立标准 run artifact 目录
  - 调度 accuracy / performance / report 三条链
  - 生成 manifest / run_status / artifact_index
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using Printf

include(joinpath(@__DIR__, "support", "release_support.jl"))
using .ReleaseSupport

const ROOT = normpath(joinpath(@__DIR__, ".."))
const BENCHMARK_ENV_ROOT = normpath(@__DIR__)
const DEFAULT_PROFILE = joinpath(ROOT, "benchmark", "configs", "release_default.toml")
const RUNS_ROOT = joinpath(ROOT, "test_results", "runs")
const STABLE_ACCURACY_DIR = joinpath(ROOT, "test_results", "accuracy")
const STABLE_REPORT_DIR = joinpath(ROOT, "test_results", "reports")

function quote_cmd(cmd::Cmd)
    return join(cmd.exec, ' ')
end

function run_logged_step(name::String, cmd::Cmd, log_path::String)
    open(log_path, "w") do io
        println(io, "[$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))] STEP=$name")
        println(io, "COMMAND=$(quote_cmd(cmd))")
        println(io)
        run(pipeline(cmd, stdout = io, stderr = io))
    end
end

function copy_file_if_exists(source_path::String, destination_path::String)
    isfile(source_path) || return false
    mkpath(dirname(destination_path))
    cp(source_path, destination_path; force = true)
    return true
end

function sync_run_outputs(paths::RunArtifactPaths)
    mirror_tree(STABLE_ACCURACY_DIR, paths.metrics_accuracy_dir)

    copy_file_if_exists(joinpath(STABLE_REPORT_DIR, "PERFORMANCE_BASELINE.csv"), joinpath(paths.metrics_performance_dir, "PERFORMANCE_BASELINE.csv"))
    copy_file_if_exists(joinpath(STABLE_REPORT_DIR, "PERFORMANCE_BASELINE.md"), joinpath(paths.report_dir, "PERFORMANCE_BASELINE.md"))
    copy_file_if_exists(joinpath(STABLE_REPORT_DIR, "PARALLEL_MPI_SAMPLE.csv"), joinpath(paths.metrics_performance_dir, "PARALLEL_MPI_SAMPLE.csv"))
    copy_file_if_exists(joinpath(STABLE_REPORT_DIR, "RELEASE_VALIDATION_REPORT.md"), joinpath(paths.report_dir, "RELEASE_VALIDATION_REPORT.md"))

    mirror_tree(joinpath(STABLE_REPORT_DIR, "assets", "accuracy"), paths.plots_accuracy_dir)
    mirror_tree(joinpath(STABLE_REPORT_DIR, "assets", "accuracy_polar"), paths.plots_accuracy_polar_dir)
    mirror_tree(joinpath(STABLE_REPORT_DIR, "assets", "performance"), paths.plots_performance_dir)
end

function promote_run_outputs(paths::RunArtifactPaths)
    mirror_tree(paths.metrics_accuracy_dir, STABLE_ACCURACY_DIR)

    copy_file_if_exists(joinpath(paths.metrics_performance_dir, "PERFORMANCE_BASELINE.csv"), joinpath(STABLE_REPORT_DIR, "PERFORMANCE_BASELINE.csv"))
    copy_file_if_exists(joinpath(paths.metrics_performance_dir, "PARALLEL_MPI_SAMPLE.csv"), joinpath(STABLE_REPORT_DIR, "PARALLEL_MPI_SAMPLE.csv"))
    copy_file_if_exists(joinpath(paths.report_dir, "PERFORMANCE_BASELINE.md"), joinpath(STABLE_REPORT_DIR, "PERFORMANCE_BASELINE.md"))
    copy_file_if_exists(joinpath(paths.report_dir, "RELEASE_VALIDATION_REPORT.md"), joinpath(STABLE_REPORT_DIR, "RELEASE_VALIDATION_REPORT.md"))

    mirror_tree(paths.plots_accuracy_dir, joinpath(STABLE_REPORT_DIR, "assets", "accuracy"))
    mirror_tree(paths.plots_accuracy_polar_dir, joinpath(STABLE_REPORT_DIR, "assets", "accuracy_polar"))
    mirror_tree(paths.plots_performance_dir, joinpath(STABLE_REPORT_DIR, "assets", "performance"))
end

function collect_curves_from_dir(accuracy_dir::String)
    isdir(accuracy_dir) || return AccuracyCurve[]
    curve_paths = sort(filter(path -> occursin(r"_phi(?:0|90)_vs_.*\.csv$", basename(path)), readdir(accuracy_dir; join = true)))
    return [load_accuracy_curve(path) for path in curve_paths]
end

function main()
    profile_path = isempty(ARGS) ? DEFAULT_PROFILE : normpath(joinpath(ROOT, ARGS[1]))
    profile = load_release_profile(profile_path)

    relpath_or(path) = haskey(profile, "paths") ? normpath(joinpath(ROOT, profile["paths"][path])) : ""
    plot_style_path = relpath_or("plot_style")
    thresholds_path = relpath_or("thresholds")
    known_exceptions_path = relpath_or("known_exceptions")
    known_exceptions = load_known_exceptions(known_exceptions_path)

    paths = create_run_artifact_dirs(RUNS_ROOT)
    manifest_path = joinpath(paths.manifest_dir, "run_manifest.toml")
    write_run_manifest(
        manifest_path;
        run_id = paths.run_id,
        workspace_root = ROOT,
        profile_path = profile_path,
        plot_style_path = plot_style_path,
        thresholds_path = thresholds_path,
        known_exceptions_path = known_exceptions_path,
        profile = profile,
    )

    copy_file_if_exists(profile_path, joinpath(paths.manifest_dir, basename(profile_path)))
    copy_file_if_exists(plot_style_path, joinpath(paths.manifest_dir, basename(plot_style_path)))
    copy_file_if_exists(thresholds_path, joinpath(paths.manifest_dir, basename(thresholds_path)))
    copy_file_if_exists(known_exceptions_path, joinpath(paths.manifest_dir, basename(known_exceptions_path)))

    steps = get(profile, "steps", Dict{String, Any}())
    workflow = get(profile, "workflow", Dict{String, Any}())
    scripts = get(profile, "scripts", Dict{String, Any}())
    accuracy = get(profile, "accuracy", Dict{String, Any}())

    run_accuracy = get(steps, "run_accuracy", false)
    run_performance = get(steps, "run_performance", true)
    run_report = get(steps, "run_report", true)
    reuse_existing_accuracy = get(workflow, "reuse_existing_accuracy", true)
    reuse_existing_performance = get(workflow, "reuse_existing_performance", true)

    if run_accuracy && !reuse_existing_accuracy
        for (index, script_rel) in enumerate(get(accuracy, "scripts", String[]))
            script_abs = normpath(joinpath(ROOT, script_rel))
            log_path = joinpath(paths.logs_dir, @sprintf("accuracy_%02d.log", index))
            cmd = `$(Base.julia_cmd()) --project=$(ROOT) $script_abs`
            run_logged_step("accuracy:$script_rel", cmd, log_path)
        end
    end

    if run_performance && !reuse_existing_performance
        perf_script = normpath(joinpath(ROOT, scripts["performance"]))
        perf_log = joinpath(paths.logs_dir, "performance.log")
        perf_cmd = `$(Base.julia_cmd()) -t auto --project=$(ROOT) $perf_script`
        run_logged_step("performance", perf_cmd, perf_log)
    end

    if run_report
        report_script = normpath(joinpath(ROOT, scripts["report"]))
        report_log = joinpath(paths.logs_dir, "report.log")
        report_cmd = `$(Base.julia_cmd()) --project=$(BENCHMARK_ENV_ROOT) $report_script`
        run_logged_step("report", report_cmd, report_log)
    end

    sync_run_outputs(paths)
    promote_run_outputs(paths)

    curves = collect_curves_from_dir(paths.metrics_accuracy_dir)
    perf_csv = joinpath(paths.metrics_performance_dir, "PERFORMANCE_BASELINE.csv")
    performance_rows = isfile(perf_csv) ? load_performance_results(perf_csv) : PerformanceBenchmarkResult[]
    statuses = build_case_statuses(curves, performance_rows, known_exceptions)

    status_csv = joinpath(paths.metrics_summary_dir, "run_status.csv")
    artifact_index_csv = joinpath(paths.metrics_summary_dir, "artifact_index.csv")
    write_case_status_csv(status_csv, statuses)
    write_artifact_index_csv(artifact_index_csv, paths.root_dir)

    open(joinpath(RUNS_ROOT, "LATEST_RUN.txt"), "w") do io
        println(io, paths.run_id)
    end

    println("Release workflow completed.")
    println("  Run ID: $(paths.run_id)")
    println("  Run root: $(relpath(paths.root_dir, ROOT))")
    println("  Manifest: $(relpath(manifest_path, ROOT))")
    println("  Status CSV: $(relpath(status_csv, ROOT))")
    println("  Artifact index: $(relpath(artifact_index_csv, ROOT))")
end

main()