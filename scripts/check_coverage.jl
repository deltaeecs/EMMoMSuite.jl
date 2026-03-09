using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC_ROOT = joinpath(ROOT, "src")
const REPORT_DIR = joinpath(ROOT, "test_results")
const REPORT_PATH = joinpath(REPORT_DIR, "COVERAGE_REPORT.md")

struct FileCoverage
    path::String
    covered::Int
    executable::Int
end

"""
    coverage_token(line) -> Union{Nothing, Int}

Parse first token in a Julia `.cov` line.
- `nothing`: non-executable line (`-` token)
- `Int`: execution count for executable line
"""
function coverage_token(line::AbstractString)
    s = lstrip(line)
    isempty(s) && return nothing
    first = split(s, ';'; limit = 2)[1]
    tok = split(first)[1]
    tok == "-" && return nothing
    try
        return parse(Int, tok)
    catch
        return nothing
    end
end

"""
    source_path_from_cov(cov_path) -> String

Convert `foo.jl.<pid>.cov` to `foo.jl`.
"""
function source_path_from_cov(cov_path::String)
    return replace(cov_path, r"\.\d+\.cov$" => "")
end

function collect_cov_files(src_root::String)
    cov_files = String[]
    for (dir, _, files) in walkdir(src_root)
        for file in files
            endswith(file, ".cov") || continue
            push!(cov_files, joinpath(dir, file))
        end
    end
    return cov_files
end

function merge_line_counts(cov_paths::Vector{String})
    max_lines = 0
    for path in cov_paths
        n = countlines(path)
        max_lines = max(max_lines, n)
    end

    # `nothing` => never observed executable on this line across shards.
    # Int => total execution count summed across shards.
    merged = Vector{Union{Nothing, Int}}(undef, max_lines)
    fill!(merged, nothing)

    for path in cov_paths
        line_no = 0
        for line in eachline(path)
            line_no += 1
            tok = coverage_token(line)
            tok === nothing && continue
            if merged[line_no] === nothing
                merged[line_no] = tok
            else
                merged[line_no] = merged[line_no] + tok
            end
        end
    end

    return merged
end

function evaluate_file_coverage(src_path::String, cov_paths::Vector{String})
    merged = merge_line_counts(cov_paths)
    executable = 0
    covered = 0

    for tok in merged
        tok === nothing && continue
        executable += 1
        tok > 0 && (covered += 1)
    end

    rel = relpath(src_path, ROOT)
    return FileCoverage(rel, covered, executable)
end

function summarize_coverage(file_cov::Vector{FileCoverage})
    total_cov = sum(fc.covered for fc in file_cov)
    total_exec = sum(fc.executable for fc in file_cov)
    pct = total_exec == 0 ? 0.0 : 100 * total_cov / total_exec
    return total_cov, total_exec, pct
end

function write_report(file_cov::Vector{FileCoverage}, total_cov::Int, total_exec::Int, pct::Float64)
    mkpath(REPORT_DIR)
    open(REPORT_PATH, "w") do io
        println(io, "# Coverage Report")
        println(io)
        println(io, "- Scope: `src/**/*.jl.*.cov`")
        println(io, "- Covered executable lines: **$(total_cov)**")
        println(io, "- Total executable lines: **$(total_exec)**")
        @printf(io, "- Coverage: **%.2f%%**\n", pct)
        println(io)
        println(io, "## Lowest Coverage Files (Top 20)")
        println(io)
        println(io, "| File | Covered | Executable | Coverage |")
        println(io, "|------|---------|------------|----------|")

        ranked = sort(file_cov; by = fc -> fc.executable == 0 ? 1.0 : fc.covered / fc.executable)
        for fc in first(ranked, min(20, length(ranked)))
            fpct = fc.executable == 0 ? 0.0 : 100 * fc.covered / fc.executable
            @printf(io, "| `%s` | %d | %d | %.2f%% |\n", fc.path, fc.covered, fc.executable, fpct)
        end
    end
end

function main()
    threshold = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 80.0

    cov_files = collect_cov_files(SRC_ROOT)
    isempty(cov_files) && error("No .cov files found under src/. Run tests with --code-coverage=user first.")

    grouped = Dict{String, Vector{String}}()
    for cov in cov_files
        src = source_path_from_cov(cov)
        if !haskey(grouped, src)
            grouped[src] = String[]
        end
        push!(grouped[src], cov)
    end

    file_cov = FileCoverage[]
    for (src, paths) in grouped
        push!(file_cov, evaluate_file_coverage(src, paths))
    end

    total_cov, total_exec, pct = summarize_coverage(file_cov)
    write_report(file_cov, total_cov, total_exec, pct)

    @printf("Coverage summary: %d/%d executable lines => %.2f%%\n", total_cov, total_exec, pct)
    println("Report written to: ", relpath(REPORT_PATH, ROOT))

    if pct + 1e-9 < threshold
        @printf("Coverage below threshold %.2f%%\n", threshold)
        exit(1)
    end
end

main()
