function write_execution_summary(io, curves::Vector{AccuracyCurve}, rows::Vector{PerformanceBenchmarkResult}, exceptions::Vector{KnownException})
    statuses = accuracy_case_status.(curves, Ref(exceptions))
    passed_curves = count(status -> status.status == "PASS", statuses)
    known_exception_curves = count(status -> status.status == "KNOWN_EXCEPTION", statuses)
    failed_curves = count(status -> status.status == "FAIL", statuses)
    worst_curve = isempty(curves) ? nothing : begin
        paired = [(curve, summarize_accuracy_curve(curve)) for curve in curves]
        sort!(paired; by = item -> item[2].rmse_dB, rev = true)
        first(paired)
    end
    slowest_case = isempty(rows) ? nothing : argmax(getfield.(rows, :t_total))
    println(io, "## Executive Summary")
    println(io)
    println(io, "- Accuracy curves within threshold: $(passed_curves) / $(length(curves))")
    println(io, "- Accuracy curves accepted as known exceptions: $(known_exception_curves)")
    println(io, "- Accuracy curves above threshold: $(failed_curves)")
    if !isnothing(worst_curve)
        curve, summary = worst_curve
        @printf(io, "- Worst far-field gap: %s (RMSE %.3f dB, threshold %.1f dB)\n", curve.label, summary.rmse_dB, accuracy_threshold(curve))
    end
    if !isnothing(slowest_case)
        row = rows[slowest_case]
        @printf(io, "- Slowest performance baseline case: %s (total %.2f s, assembly %.2f s, solve %.2f s)\n", row.case_name, row.t_total, row.t_assembly, row.t_solve)
    end
    if any(status -> status.status == "KNOWN_EXCEPTION", statuses)
        println(io, "- Known exception policy active: waived release exceptions are tracked separately and excluded from hard blocker count.")
    end
    println(io)
end

function write_coverage_section(io, curves::Vector{AccuracyCurve}, rows::Vector{PerformanceBenchmarkResult})
    println(io, "## Coverage Matrix")
    println(io)
    println(io, "| Module | Evidence | Status |")
    println(io, "|--------|----------|--------|")
    for item in coverage_rows(curves, rows)
        println(io, "| $(item.area) | $(item.evidence) | $(item.status) |")
    end
    println(io)
end

function write_efficiency_diagnostics(io, rows::Vector{PerformanceBenchmarkResult}, parallel_sample)
    println(io, "## Efficiency and Memory Diagnostics")
    println(io)
    if isempty(rows)
        println(io, "- Performance CSV missing, diagnostics deferred.")
        println(io)
        return
    end
    direct_rows = filter(row -> row.solver == "LU", rows)
    if !isempty(direct_rows)
        densest = direct_rows[argmax(getfield.(direct_rows, :N))]
        @printf(io, "- Largest dense direct case in the current baseline: %s, N=%d, single dense complex matrix footprint is about %.2f GiB.\n", densest.case_name, densest.N, dense_matrix_gib(densest.N))
    end
    slowest_assembly = rows[argmax(getfield.(rows, :t_assembly))]
    slowest_total = rows[argmax(getfield.(rows, :t_total))]
    @printf(io, "- Dominant setup hotspot: %s assembly/setup %.2f s.\n", slowest_assembly.case_name, slowest_assembly.t_assembly)
    @printf(io, "- Dominant end-to-end runtime hotspot: %s total %.2f s.\n", slowest_total.case_name, slowest_total.t_total)
    println(io, "- Threaded Phase 16 baseline was executed with `julia -t auto --project=.`, so the plotted MLFMA timings reflect multi-thread runtime rather than single-thread debug timing.")
    if !isnothing(parallel_sample)
        @printf(io, "- Fresh MPI sample: %s used %d ranks × %d threads, N=%d, assembly %.4f s.\n", csv_row_string(parallel_sample, "case_name"), csv_row_int(parallel_sample, "mpi_procs"), csv_row_int(parallel_sample, "threads_per_rank"), csv_row_int(parallel_sample, "unknowns"), csv_row_float(parallel_sample, "assembly_time_s"))
    end
    println(io)
end

function write_parallel_section(io, parallel_sample)
    println(io)
    println(io, "## Parallel Sample")
    println(io)
    if isnothing(parallel_sample)
        println(io, "- No fresh multi-process runtime sample recorded yet.")
        println(io)
        return
    end
    println(io, "| Case | MPI Ranks | Threads/Rank | N | Assembly (s) | Current Norm |")
    println(io, "|------|-----------|--------------|---|--------------|--------------|")
    @printf(io, "| %s | %d | %d | %d | %.4f | %.6f |\n", csv_row_string(parallel_sample, "case_name"), csv_row_int(parallel_sample, "mpi_procs"), csv_row_int(parallel_sample, "threads_per_rank"), csv_row_int(parallel_sample, "unknowns"), csv_row_float(parallel_sample, "assembly_time_s"), csv_row_float(parallel_sample, "current_norm"))
    println(io)
end

function write_release_risks(io, curves::Vector{AccuracyCurve}, parallel_sample, exceptions::Vector{KnownException})
    println(io, "## Risks and Release Notes")
    println(io)
    flagged = [(curve, summarize_accuracy_curve(curve)) for curve in curves if summarize_accuracy_curve(curve).rmse_dB > accuracy_threshold(curve)]
    hard_blockers = filter(item -> isnothing(match_known_exception(item[1].label, exceptions)), flagged)
    waived = filter(item -> !isnothing(match_known_exception(item[1].label, exceptions)), flagged)
    if isempty(flagged)
        println(io, "- No accuracy curve exceeds the current release threshold.")
    else
        for (curve, summary) in sort(hard_blockers; by = item -> item[2].rmse_dB, rev = true)
            @printf(io, "- Accuracy exception: %s RMSE %.3f dB exceeds threshold %.1f dB.\n", curve.label, summary.rmse_dB, accuracy_threshold(curve))
        end
        for (curve, summary) in sort(waived; by = item -> item[2].rmse_dB, rev = true)
            @printf(io, "- Waived known exception: %s RMSE %.3f dB exceeds threshold %.1f dB but is covered by registry.\n", curve.label, summary.rmse_dB, accuracy_threshold(curve))
        end
    end
    if isnothing(parallel_sample)
        println(io, "- Parallel coverage is backed by dedicated MPI tests and a benchmark entrypoint, but the release report still does not record a fresh multi-process runtime sample.")
    end
    println(io)
end

function write_release_recommendation(io, curves::Vector{AccuracyCurve}, exceptions::Vector{KnownException})
    println(io, "## Release Recommendation")
    println(io)
    flagged = [(curve, summarize_accuracy_curve(curve)) for curve in curves if summarize_accuracy_curve(curve).rmse_dB > accuracy_threshold(curve)]
    hard_blockers = filter(item -> isnothing(match_known_exception(item[1].label, exceptions)), flagged)
    waived = filter(item -> !isnothing(match_known_exception(item[1].label, exceptions)), flagged)
    if isempty(flagged)
        println(io, "- Final decision: GO")
        println(io, "- Rationale: no accuracy or engineering blockers remain in the current release validation report.")
    elseif isempty(hard_blockers) && !isempty(waived)
        println(io, "- Final decision: GO WITH KNOWN LEGACY EXCEPTION")
        println(io, "- Exception set kept open: all curves above threshold are covered by the known-exception registry and therefore tracked as accepted release waivers.")
        println(io, "- Engineering blockers: closed.")
    else
        println(io, "- Final decision: NO-GO")
        println(io, "- Rationale: non-legacy blockers still remain in the current validation report.")
    end
    println(io)
end

function write_accuracy_section(io, curves::Vector{AccuracyCurve}, groups, plot_paths, polar_plot_paths)
    println(io, "## Accuracy Summary")
    println(io)
    if isempty(curves)
        println(io, "No far-field accuracy CSV artifacts found under test_results/accuracy.")
        println(io)
        return
    end
    println(io, "| Curve | Reference | Points | RMSE (dB) | MaxErr (dB) | Bias (dB) |")
    println(io, "|------|-----------|--------|-----------|-------------|-----------|")
    for curve in sort(curves; by = c -> c.label)
        summary = summarize_accuracy_curve(curve)
        @printf(io, "| %s | %s | %d | %.3f | %.3f | %+.3f |\n", curve.label, reference_label(curve), summary.n_points, summary.rmse_dB, summary.max_err_dB, summary.mean_bias_dB)
    end
    println(io)
    println(io, "## Far-Field Curve Plots")
    println(io)
    for group_name in sort(collect(keys(groups)))
        println(io, "### $(pretty_label(group_name))")
        println(io)
        println(io, "Cartesian view:")
        println(io)
        println(io, "![$(group_name)]($(markdown_relpath(plot_paths[group_name])))")
        println(io)
        println(io, "Polar view:")
        println(io)
        println(io, "![$(group_name)_polar]($(markdown_relpath(polar_plot_paths[group_name])))")
        println(io)
    end
end

function write_case_status_section(io, curves::Vector{AccuracyCurve}, rows::Vector{PerformanceBenchmarkResult}, exceptions::Vector{KnownException})
    println(io)
    println(io, "## Case Status Matrix")
    println(io)
    statuses = build_case_statuses(curves, rows, exceptions)
    println(io, "| Case | Category | Status | Metric | Value | Threshold | Note |")
    println(io, "|------|----------|--------|--------|-------|-----------|------|")
    for status in statuses
        @printf(io, "| %s | %s | %s | %s | %.3f | %.3f | %s |\n", status.case_name, status.category, status.status, status.metric_name, status.metric_value, status.threshold, isempty(status.note) ? "-" : replace(status.note, "|" => "/"))
    end
    println(io)
end

function write_known_exceptions_section(io, exceptions::Vector{KnownException}, curves::Vector{AccuracyCurve})
    println(io)
    println(io, "## Known Exceptions")
    println(io)
    if isempty(exceptions)
        println(io, "- No known exceptions registry configured.")
        println(io)
        return
    end
    println(io, "| ID | Pattern | Disposition | Tracking | Matched Curves | Rationale |")
    println(io, "|----|---------|-------------|----------|----------------|-----------|")
    for exception in exceptions
        matched = join([curve.label for curve in curves if !isnothing(match_known_exception(curve.label, [exception]))], "; ")
        println(io, "| $(exception.id) | $(exception.label_pattern) | $(exception.disposition) | $(exception.tracking) | $(isempty(matched) ? "-" : matched) | $(replace(exception.rationale, "|" => "/")) |")
    end
    println(io)
end

function write_antenna_section(io)
    a_rows = load_antenna_rows(["A1_halfwave_direct", "A2_halfwave_mlfma", "A3_resonant_direct", "A4_50ohm_s11"])
    b_rows = load_antenna_rows(["B1_PMCHW_sphere_eps4_Zin", "B1_PMCHW_eps1_vs_2xEFIE_Zin", "B2_PMCHW_MLFMA_Zin", "B3_VEFIE_TriTetra_direct_Zin", "B4_VCFIE_TriTetra_direct_Zin", "B5_VCFIE_TriTetra_MLFMA_Zin"])
    println(io)
    println(io, "## Antenna and Port Metrics")
    println(io)
    if !isempty(a_rows)
        println(io, "### A-Series")
        println(io)
        println(io, "| Case | Z_in | Z_ref | Status |")
        println(io, "|------|------|-------|--------|")
        for item in a_rows
            row = item.row
            zin = format_complex_pair(csv_row_float(row, "Z_in_re"), csv_row_float(row, "Z_in_im"))
            zref = format_complex_pair(csv_row_float(row, "Z_ref_re"), csv_row_float(row, "Z_ref_im"))
            status = bool_status(csv_row_bool(row, "passed"))
            println(io, "| $(item.label) | $(zin) | $(zref) | $(status) |")
        end
        println(io)
    end
    if !isempty(b_rows)
        println(io, "### B-Series")
        println(io)
        println(io, "| Case | Z_in | Z_ref | Status |")
        println(io, "|------|------|-------|--------|")
        for item in b_rows
            row = item.row
            zin = format_complex_pair(csv_row_float(row, "Zin_re"), csv_row_float(row, "Zin_im"))
            zref = has_column(row, :Zref_re) ? format_complex_pair(csv_row_float(row, "Zref_re"), csv_row_float(row, "Zref_im")) : "n/a"
            status = bool_status(csv_row_bool(row, "passed"))
            println(io, "| $(item.label) | $(zin) | $(zref) | $(status) |")
        end
        println(io)
    end
end

function write_performance_section(io, rows::Vector{PerformanceBenchmarkResult})
    println(io)
    println(io, "## Performance Summary")
    println(io)
    if isempty(rows)
        println(io, "性能基线 CSV 缺失：请先运行 `julia -t auto --project=. benchmark/performance_baseline.jl`。")
        println(io)
        return
    end
    println(io, "| Case | Equation | Solver | N | Assembly/Setup (s) | Solve (s) | Total (s) |")
    println(io, "|------|----------|--------|---|--------------------|-----------|-----------|")
    for row in rows
        @printf(io, "| %s | %s | %s | %d | %.2f | %.2f | %.2f |\n", row.case_name, row.equation, row.solver, row.N, row.t_assembly, row.t_solve, row.t_total)
    end
    total_plot = generate_total_runtime_plot(rows)
    breakdown_plot = generate_breakdown_plot(rows)
    println(io)
    println(io, "## Performance Plots")
    println(io)
    println(io, "### Total Runtime")
    println(io)
    println(io, "![performance_total_runtime]($(markdown_relpath(total_plot)))")
    println(io)
    println(io, "### Breakdown")
    println(io)
    println(io, "![performance_breakdown]($(markdown_relpath(breakdown_plot)))")
    println(io)
end

function build_report()
    ensure_dirs()
    style = load_plot_style()
    exceptions = load_release_exceptions()
    curves = collect_accuracy_curves()
    groups, plot_paths, polar_plot_paths = generate_accuracy_plots(curves, style)
    performance_rows = load_performance_rows()
    parallel_sample = load_parallel_sample_row()
    open(REPORT_PATH, "w") do io
        println(io, "# RELEASE_VALIDATION_REPORT")
        println(io)
        println(io, "- Generated at: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "- Generator: `benchmark/run_release_validation_report.jl`")
        println(io, "- Julia: $(VERSION)")
        println(io, "- Threads: $(Threads.nthreads())")
        println(io, "- Git commit: $(current_git_commit())")
        println(io)
        println(io, "## Artifacts")
        println(io)
        println(io, "- Accuracy report: $(isfile(ACCURACY_REPORT_PATH) ? "FOUND" : "MISSING")")
        println(io, "- Accuracy curves: $(length(curves))")
        println(io, "- Performance CSV: $(isfile(PERFORMANCE_CSV_PATH) ? "FOUND" : "MISSING")")
        println(io, "- Parallel sample CSV: $(isfile(PARALLEL_SAMPLE_CSV_PATH) ? "FOUND" : "MISSING")")
        println(io)
        write_execution_summary(io, curves, performance_rows, exceptions)
        write_accuracy_section(io, curves, groups, plot_paths, polar_plot_paths)
        write_antenna_section(io)
        write_performance_section(io, performance_rows)
        write_parallel_section(io, parallel_sample)
        write_case_status_section(io, curves, performance_rows, exceptions)
        write_known_exceptions_section(io, exceptions, curves)
        write_coverage_section(io, curves, performance_rows)
        write_efficiency_diagnostics(io, performance_rows, parallel_sample)
        write_release_risks(io, curves, parallel_sample, exceptions)
        write_release_recommendation(io, curves, exceptions)
    end
    println("Release validation report generated: ", relpath(REPORT_PATH, ROOT))
end