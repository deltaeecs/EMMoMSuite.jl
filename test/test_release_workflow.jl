using Test

include(joinpath(@__DIR__, "..", "benchmark", "support", "release_support.jl"))
using .ReleaseSupport

@testset "Phase 17 ReleaseWorkflow" begin
    mktempdir() do tmpdir
        profile_path = joinpath(tmpdir, "release_default.toml")
        open(profile_path, "w") do io
            println(io, "[workflow]")
            println(io, "name = \"release-default\"")
            println(io, "[steps]")
            println(io, "run_accuracy = true")
            println(io, "reuse_accuracy = false")
        end

        known_ex_path = joinpath(tmpdir, "known_exceptions.toml")
        open(known_ex_path, "w") do io
            println(io, "[[exceptions]]")
            println(io, "id = \"legacy-f2\"")
            println(io, "label_pattern = \"F2_CFIE_Jet_Direct\"")
            println(io, "scope = \"release\"")
            println(io, "rationale = \"legacy inherited\"")
            println(io, "disposition = \"waived\"")
            println(io, "tracking = \"phase16\"")
        end

        profile = load_release_profile(profile_path)
        exceptions = load_known_exceptions(known_ex_path)
        paths = create_run_artifact_dirs(joinpath(tmpdir, "runs"); run_id = "20260311_120000")

        @test profile["workflow"]["name"] == "release-default"
        @test length(exceptions) == 1
        @test exceptions[1].id == "legacy-f2"
        @test paths.run_id == "20260311_120000"
        @test isdir(paths.report_dir)
        @test isdir(paths.plots_accuracy_polar_dir)

        manifest_path = joinpath(paths.manifest_dir, "run_manifest.toml")
        write_run_manifest(
            manifest_path;
            run_id = paths.run_id,
            workspace_root = tmpdir,
            profile_path = profile_path,
            plot_style_path = joinpath(tmpdir, "plot_style.toml"),
            thresholds_path = joinpath(tmpdir, "thresholds.toml"),
            known_exceptions_path = known_ex_path,
            profile = profile,
        )
        @test isfile(manifest_path)

        curve = AccuracyCurve(
            "F2_CFIE_Jet_Direct_phi0_vs_Feko",
            [-180.0, 0.0, 180.0],
            [1.0, 2.0, 3.0],
            [3.0, 4.0, 5.0],
            [-2.0, -2.0, -2.0],
            "rcs_emsuite_dBsm",
            "rcs_feko_dBsm",
            "curve.csv",
        )
        status = accuracy_case_status(curve, exceptions)
        @test status.status == "PASS"

        failing_curve = AccuracyCurve(
            "F2_CFIE_Jet_Direct_phi90_vs_Feko",
            [-180.0, 0.0, 180.0],
            [1.0, 2.0, 3.0],
            [6.0, 7.0, 8.0],
            [-5.0, -5.0, -5.0],
            "rcs_emsuite_dBsm",
            "rcs_feko_dBsm",
            "curve_fail.csv",
        )
        known_exception_status = accuracy_case_status(failing_curve, exceptions)
        @test known_exception_status.status == "KNOWN_EXCEPTION"

        perf_row = PerformanceBenchmarkResult("Jet EFIE Direct", "EFIE", "LU", 14559, 1.0, 2.0, 0.0, 3.0, 0.5, 6.5, "")
        statuses = build_case_statuses([curve, failing_curve], [perf_row], exceptions)
        @test length(statuses) == 3

        status_csv = joinpath(paths.metrics_summary_dir, "run_status.csv")
        artifact_csv = joinpath(paths.metrics_summary_dir, "artifact_index.csv")
        write_case_status_csv(status_csv, statuses)
        @test isfile(status_csv)

        sample_plot = joinpath(paths.plots_accuracy_dir, "plot.png")
        write(sample_plot, "png")
        write_artifact_index_csv(artifact_csv, paths.root_dir)
        @test isfile(artifact_csv)

        latest_src = joinpath(tmpdir, "source_tree")
        mkpath(latest_src)
        write(joinpath(latest_src, "example.txt"), "hello")
        latest_dst = joinpath(tmpdir, "dest_tree")
        mirror_tree(latest_src, latest_dst)
        @test read(joinpath(latest_dst, "example.txt"), String) == "hello"
    end
end