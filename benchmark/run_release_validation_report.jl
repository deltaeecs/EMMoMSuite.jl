"""
run_release_validation_report.jl

Purpose:
    - Generate a unified release validation report from existing accuracy and performance artifacts.
    - Produce far-field comparison plots and performance comparison plots alongside markdown summary.
    - Auto-activate the benchmark environment when launched from the root project for a smoother CLI entry.

Usage:
    julia --project=benchmark benchmark/run_release_validation_report.jl
    julia --project=. benchmark/run_release_validation_report.jl
"""

using Pkg

const BENCHMARK_ENV_ROOT = normpath(@__DIR__)
const ACTIVE_PROJECT = Base.active_project()

if isnothing(ACTIVE_PROJECT) || normpath(dirname(ACTIVE_PROJECT)) != BENCHMARK_ENV_ROOT
    @info "Activating benchmark environment for release report" benchmark_env = BENCHMARK_ENV_ROOT active_project = ACTIVE_PROJECT
    Pkg.activate(BENCHMARK_ENV_ROOT; io = devnull)
end

using Dates
using Printf
using TOML
ENV["GKSwstype"] = "100"
using Plots

gr()

const ROOT = normpath(joinpath(@__DIR__, ".."))
const ACCURACY_DIR = joinpath(ROOT, "test_results", "accuracy")
const REPORT_DIR = joinpath(ROOT, "test_results", "reports")
const ASSET_DIR = joinpath(REPORT_DIR, "assets")
const ACCURACY_ASSET_DIR = joinpath(ASSET_DIR, "accuracy")
const ACCURACY_POLAR_ASSET_DIR = joinpath(ASSET_DIR, "accuracy_polar")
const PERFORMANCE_ASSET_DIR = joinpath(ASSET_DIR, "performance")
const REPORT_PATH = joinpath(REPORT_DIR, "RELEASE_VALIDATION_REPORT.md")
const PERFORMANCE_CSV_PATH = joinpath(REPORT_DIR, "PERFORMANCE_BASELINE.csv")
const PARALLEL_SAMPLE_CSV_PATH = joinpath(REPORT_DIR, "PARALLEL_MPI_SAMPLE.csv")
const ACCURACY_REPORT_PATH = joinpath(ACCURACY_DIR, "ACCURACY_REPORT.md")
const CONFIG_DIR = joinpath(ROOT, "benchmark", "configs")
const PLOT_STYLE_PATH = joinpath(CONFIG_DIR, "plot_style.toml")
const KNOWN_EXCEPTIONS_PATH = joinpath(CONFIG_DIR, "known_exceptions.toml")

include(joinpath(@__DIR__, "support", "release_support.jl"))
using .ReleaseSupport
include(joinpath(@__DIR__, "reporting", "collector.jl"))
include(joinpath(@__DIR__, "reporting", "plotting.jl"))
include(joinpath(@__DIR__, "reporting", "writer.jl"))

build_report()
