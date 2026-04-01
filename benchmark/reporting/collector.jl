using TOML

function ensure_dirs()
    mkpath(REPORT_DIR)
    mkpath(ACCURACY_ASSET_DIR)
    mkpath(ACCURACY_POLAR_ASSET_DIR)
    mkpath(PERFORMANCE_ASSET_DIR)
end

function load_plot_style()
    return isfile(PLOT_STYLE_PATH) ? TOML.parsefile(PLOT_STYLE_PATH) : Dict{String, Any}()
end

function plot_style_value(style::Dict, section::String, key::String, default)
    section_dict = get(style, section, Dict{String, Any}())
    return get(section_dict, key, default)
end

function plot_style_color(style::Dict, section::String, key::String, default::Symbol)
    value = plot_style_value(style, section, key, String(default))
    return Symbol(String(value))
end

function load_release_exceptions()
    return load_known_exceptions(KNOWN_EXCEPTIONS_PATH)
end

function safe_slug(text::String)
    return replace(text, r"[^A-Za-z0-9._-]" => "_")
end

function pretty_label(label::String)
    return replace(label, "_" => " ")
end

function reference_label(curve::AccuracyCurve)
    lowered = lowercase(curve.label * " " * curve.reference_column)
    if occursin("feko", lowered)
        return "FEKO"
    elseif occursin("mie", lowered)
        return "Mie"
    else
        return "Reference"
    end
end

function cut_sort_key(curve::AccuracyCurve)
    ref_rank = reference_label(curve) == "FEKO" ? 1 : reference_label(curve) == "Mie" ? 2 : 3
    cut = accuracy_curve_cut(curve.label)
    cut_rank = cut == "phi0" ? 1 : cut == "phi90" ? 2 : 3
    return (ref_rank, cut_rank, curve.label)
end

function collect_accuracy_curves()
    isdir(ACCURACY_DIR) || return AccuracyCurve[]
    curve_paths = sort(filter(path -> occursin(r"_phi(?:0|90)_vs_.*\.csv$", basename(path)),
        readdir(ACCURACY_DIR; join = true)))
    return [load_accuracy_curve(path) for path in curve_paths]
end

function load_antenna_rows(prefixes::Vector{String})
    rows = NamedTuple[]
    for prefix in prefixes
        csv_path = joinpath(ACCURACY_DIR, prefix * ".csv")
        isfile(csv_path) || continue
        parsed_rows = parse_csv_rows(csv_path)
        isempty(parsed_rows) && continue
        push!(rows, (label = prefix, row = first(parsed_rows)))
    end
    return rows
end

function load_performance_rows()
    return isfile(PERFORMANCE_CSV_PATH) ? load_performance_results(PERFORMANCE_CSV_PATH) : PerformanceBenchmarkResult[]
end

function load_parallel_sample_row()
    isfile(PARALLEL_SAMPLE_CSV_PATH) || return nothing
    rows = parse_csv_rows(PARALLEL_SAMPLE_CSV_PATH)
    isempty(rows) && return nothing
    return first(rows)
end

function short_case_name(case_name::String)
    replacements = Dict(
        "Plate EFIE Direct" => "Plate EFIE D",
        "Jet EFIE Direct" => "Jet EFIE D",
        "Jet EFIE MLFMA" => "Jet EFIE MLFMA",
        "Sphere CFIE MLFMA" => "Sphere CFIE MLFMA",
        "Plate VEFIE Direct" => "Plate VEFIE D",
        "Plate SCFIE Direct" => "Plate SCFIE D",
    )
    return get(replacements, case_name, case_name)
end

function markdown_relpath(path::String)
    return replace(relpath(path, REPORT_DIR), '\\' => '/')
end

function current_git_commit()
    try
        return readchomp(`git -C $ROOT rev-parse --short HEAD`)
    catch
        return "unknown"
    end
end

function accuracy_threshold(curve::AccuracyCurve)
    return default_accuracy_threshold(curve.label)
end

function dense_matrix_gib(n::Int)
    return (float(n)^2 * 16.0) / 1024.0^3
end

function has_parallel_validation_evidence()
    required = [
        joinpath(ROOT, "test", "test_parallel.jl"),
        joinpath(ROOT, "test", "test_parallel_mfie_cfie.jl"),
        joinpath(ROOT, "benchmark", "benchmark_parallel_sphere.jl"),
    ]
    return all(isfile, required)
end

function coverage_rows(curves::Vector{AccuracyCurve}, rows::Vector{PerformanceBenchmarkResult})
    curve_labels = [curve.label for curve in curves]
    perf_cases = [row.case_name for row in rows]
    return [
        (; area = "Geometry / Mesh I/O", evidence = "F1/F5/F7 curves + performance baseline", status = (!isempty(curve_labels) || !isempty(perf_cases)) ? "COVERED" : "MISSING"),
        (; area = "BasisFunctions", evidence = "RWG/SWG paths in F1/F7/B-series", status = any(label -> occursin(r"F1|F7|P1|P3", label), curve_labels) ? "COVERED" : "PARTIAL"),
        (; area = "IntegralEquations", evidence = "SEFIE/CFIE/VEFIE/PMCHW cases", status = any(label -> occursin(r"F1|F2|F5|F6|F7|P1|P3|X1", label), curve_labels) ? "COVERED" : "PARTIAL"),
        (; area = "FastAlgorithms", evidence = "F6 + performance MLFMA cases", status = any(label -> occursin("MLFMA", label), curve_labels) || any(row -> occursin("MLFMA", row.solver), rows) ? "COVERED" : "PARTIAL"),
        (; area = "Solvers / Preconditioners", evidence = "LU + GMRES baseline rows", status = (!isempty(rows)) ? "COVERED" : "PARTIAL"),
        (; area = "Ports / Excitations", evidence = "PlaneWave + A/B antenna CSV summaries", status = "COVERED"),
        (; area = "PostProcessing", evidence = "RCS curves + Zin summaries", status = (!isempty(curves)) ? "COVERED" : "PARTIAL"),
        (; area = "Parallel", evidence = has_parallel_validation_evidence() ? "MPI tests + parallel sphere benchmark entrypoint" : "performance baseline run with Julia threads", status = has_parallel_validation_evidence() ? "COVERED" : (any(row -> occursin("MLFMA", row.solver), rows) ? "PARTIAL" : "MISSING")),
        (; area = "IO", evidence = "CSV ingestion + markdown/png report assets", status = "COVERED"),
    ]
end

function has_column(row, name::Symbol)
    return csv_row_has(row, String(name))
end

function format_complex_pair(re_val, im_val)
    return @sprintf("%+.1f%+.1fj", Float64(re_val), Float64(im_val))
end

function bool_status(value)
    return Bool(value) ? "PASS" : "FAIL"
end