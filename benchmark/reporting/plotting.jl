function generate_accuracy_group_plot(group_name::String, curves::Vector{AccuracyCurve}, style::Dict)
    ordered = sort(curves; by = cut_sort_key)
    nrows = length(ordered)
    width = Int(plot_style_value(style, "cartesian", "width", 1200))
    height_per_panel = Int(plot_style_value(style, "cartesian", "height_per_panel", 320))
    dpi = Int(plot_style_value(style, "cartesian", "dpi", 150))
    lw_main = Float64(plot_style_value(style, "cartesian", "line_width_main", 2.2))
    lw_ref = Float64(plot_style_value(style, "cartesian", "line_width_ref", 1.8))
    main_color = plot_style_color(style, "cartesian", "main_color", :steelblue)
    ref_color = plot_style_color(style, "cartesian", "reference_color", :darkorange)
    plt = plot(layout = (nrows, 1), size = (width, height_per_panel * nrows), dpi = dpi)

    for (idx, curve) in enumerate(ordered)
        summary = summarize_accuracy_curve(curve)
        cut_name = uppercase(accuracy_curve_cut(curve.label))
        ref_name = reference_label(curve)
        title = "$(pretty_label(group_name)) | $(cut_name) | $(ref_name) | RMSE=$(round(summary.rmse_dB, digits = 3)) dB"
        y_all = vcat(curve.rcs_model_dBsm, curve.rcs_reference_dBsm)
        x_annot = maximum(curve.theta_deg) - 0.02 * (maximum(curve.theta_deg) - minimum(curve.theta_deg))
        y_annot = minimum(y_all) + 0.08 * (maximum(y_all) - minimum(y_all) + eps())
        plot!(plt[idx], curve.theta_deg, curve.rcs_model_dBsm; label = "EMSuite", lw = lw_main, color = main_color)
        plot!(plt[idx], curve.theta_deg, curve.rcs_reference_dBsm; label = ref_name, lw = lw_ref, ls = :dash, color = ref_color)
        xlabel!(plt[idx], "theta (deg)")
        ylabel!(plt[idx], "RCS (dBsm)")
        title!(plt[idx], title)
        annotate!(plt[idx], x_annot, y_annot, text("Max=$(round(summary.max_err_dB, digits = 3)) dB  Bias=$(round(summary.mean_bias_dB, digits = 3)) dB", 9, :right))
    end

    outpath = joinpath(ACCURACY_ASSET_DIR, safe_slug(group_name) * ".png")
    savefig(plt, outpath)
    return outpath
end

function generate_accuracy_group_polar_plot(group_name::String, curves::Vector{AccuracyCurve}, style::Dict)
    ordered = sort(curves; by = cut_sort_key)
    nrows = length(ordered)
    width = Int(plot_style_value(style, "polar", "width", 1200))
    height_per_panel = Int(plot_style_value(style, "polar", "height_per_panel", 320))
    dpi = Int(plot_style_value(style, "polar", "dpi", 150))
    lw_main = Float64(plot_style_value(style, "polar", "line_width_main", 2.0))
    lw_ref = Float64(plot_style_value(style, "polar", "line_width_ref", 1.6))
    padding = Float64(plot_style_value(style, "polar", "radial_padding", 1.0))
    main_color = plot_style_color(style, "polar", "main_color", :steelblue)
    ref_color = plot_style_color(style, "polar", "reference_color", :darkorange)
    plt = plot(layout = (nrows, 1), size = (width, height_per_panel * nrows), dpi = dpi)

    for (idx, curve) in enumerate(ordered)
        summary = summarize_accuracy_curve(curve)
        cut_name = uppercase(accuracy_curve_cut(curve.label))
        ref_name = reference_label(curve)
        angles = deg2rad.(mod.(curve.theta_deg .+ 360.0, 360.0))
        order = sortperm(angles)
        y_all = vcat(curve.rcs_model_dBsm, curve.rcs_reference_dBsm)
        radial_offset = abs(minimum(y_all)) + padding
        plot!(plt[idx], angles[order], curve.rcs_model_dBsm[order] .+ radial_offset; proj = :polar, label = "EMSuite", lw = lw_main, color = main_color)
        plot!(plt[idx], angles[order], curve.rcs_reference_dBsm[order] .+ radial_offset; proj = :polar, label = ref_name, lw = lw_ref, ls = :dash, color = ref_color)
        title!(plt[idx], "$(pretty_label(group_name)) | $(cut_name) | $(ref_name) | RMSE=$(round(summary.rmse_dB, digits = 3)) dB | offset=$(round(radial_offset, digits = 2))")
    end

    outpath = joinpath(ACCURACY_POLAR_ASSET_DIR, safe_slug(group_name) * "_polar.png")
    savefig(plt, outpath)
    return outpath
end

function generate_accuracy_plots(curves::Vector{AccuracyCurve}, style::Dict)
    groups = Dict{String, Vector{AccuracyCurve}}()
    for curve in curves
        push!(get!(groups, accuracy_curve_group(curve.label), AccuracyCurve[]), curve)
    end
    plot_paths = Dict{String, String}()
    polar_plot_paths = Dict{String, String}()
    for group_name in sort(collect(keys(groups)))
        plot_paths[group_name] = generate_accuracy_group_plot(group_name, groups[group_name], style)
        polar_plot_paths[group_name] = generate_accuracy_group_polar_plot(group_name, groups[group_name], style)
    end
    return groups, plot_paths, polar_plot_paths
end

function generate_total_runtime_plot(rows::Vector{PerformanceBenchmarkResult})
    labels = short_case_name.(getfield.(rows, :case_name))
    totals = getfield.(rows, :t_total)
    x = collect(eachindex(labels))
    plt = bar(x, totals; legend = false, xlabel = "case", ylabel = "time (s)", title = "Performance Total Runtime", color = :teal, size = (1200, 500), xticks = (x, labels), xrotation = 20, dpi = 150)
    outpath = joinpath(PERFORMANCE_ASSET_DIR, "performance_total_runtime.png")
    savefig(plt, outpath)
    return outpath
end

function generate_breakdown_plot(rows::Vector{PerformanceBenchmarkResult})
    labels = short_case_name.(getfield.(rows, :case_name))
    x = collect(eachindex(labels))
    mesh = getfield.(rows, :t_mesh)
    assembly = getfield.(rows, :t_assembly)
    precond = getfield.(rows, :t_precond)
    solve = getfield.(rows, :t_solve)
    rcs = getfield.(rows, :t_rcs)
    values = hcat(mesh, assembly, precond, solve, rcs)
    plt = bar(x, values; label = ["mesh+basis" "assembly/setup" "preconditioner" "solve" "postprocess"], xlabel = "case", ylabel = "time (s)", title = "Performance Breakdown by Stage", size = (1280, 560), xticks = (x, labels), xrotation = 20, dpi = 150)
    outpath = joinpath(PERFORMANCE_ASSET_DIR, "performance_breakdown.png")
    savefig(plt, outpath)
    return outpath
end