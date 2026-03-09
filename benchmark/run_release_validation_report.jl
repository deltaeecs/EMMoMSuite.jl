"""
run_release_validation_report.jl

Purpose:
- Build a release validation report skeleton from existing benchmark outputs.
- Keep a single markdown entrypoint for pre-release quality review.

Usage:
  julia --project=. benchmark/run_release_validation_report.jl
"""

using Dates
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const REPORT_DIR = joinpath(ROOT, "test_results", "reports")
const REPORT_PATH = joinpath(REPORT_DIR, "RELEASE_VALIDATION_REPORT.md")
const TEMPLATE_PATH = joinpath(REPORT_DIR, "RELEASE_VALIDATION_REPORT_TEMPLATE.md")

function read_text_or_empty(path::String)
    return isfile(path) ? read(path, String) : ""
end

function find_existing_outputs()
    return Dict(
        "coverage" => isfile(joinpath(ROOT, "test_results", "COVERAGE_REPORT.md")),
        "accuracy_report" => isfile(joinpath(ROOT, "test_results", "accuracy", "ACCURACY_REPORT.md")),
        "budget_medium" => isfile(joinpath(ROOT, "test_results", "accuracy", "PMCHW_MLFMA_budget_sweep_medium.csv")),
        "budget_krylov_medium" => isfile(joinpath(ROOT, "test_results", "accuracy", "PMCHW_MLFMA_budget_krylov_medium.csv")),
        "block_fidelity_medium" => isfile(joinpath(ROOT, "test_results", "accuracy", "PMCHW_block_fidelity_medium.csv")),
        "nmuller_medium" => isfile(joinpath(ROOT, "test_results", "accuracy", "PMCHW_NMuller_gmres_trajectory_medium.csv")),
    )
end

function build_report()
    mkpath(REPORT_DIR)
    template = read_text_or_empty(TEMPLATE_PATH)
    outputs = find_existing_outputs()

    open(REPORT_PATH, "w") do io
        println(io, "# RELEASE_VALIDATION_REPORT")
        println(io)
        println(io, "- Generated at: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "- Generator: `benchmark/run_release_validation_report.jl`")
        println(io)

        println(io, "## Collected Artifacts")
        println(io)
        for (k, v) in outputs
            status = v ? "FOUND" : "MISSING"
            println(io, "- `$(k)`: $(status)")
        end

        println(io)
        println(io, "## Template")
        println(io)
        if isempty(template)
            println(io, "Template file not found: `test_results/reports/RELEASE_VALIDATION_REPORT_TEMPLATE.md`")
        else
            println(io, template)
        end
    end

    println("Release validation report generated: ", relpath(REPORT_PATH, ROOT))
end

build_report()
