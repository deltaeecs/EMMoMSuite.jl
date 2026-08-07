using Test

if !isdefined(Main, :ReleaseSupport)
    include(joinpath(@__DIR__, "..", "benchmark", "support", "release_support.jl"))
end
using .ReleaseSupport

@testset "Phase 16 BenchmarkReportData" begin
    mktempdir() do tmpdir
        accuracy_csv = joinpath(tmpdir, "accuracy.csv")
        open(accuracy_csv, "w") do io
            println(io, "theta_deg,rcs_emsuite_dBsm,rcs_feko_dBsm,diff_dB")
            println(io, "-180.0,1.0,1.5,-0.5")
            println(io, "0.0,2.0,1.0,1.0")
            println(io, "180.0,3.0,2.0,1.0")
        end

        curve = load_accuracy_curve(accuracy_csv; label = "F1_SEFIE_Jet_Direct_phi0_vs_Feko")
        summary = summarize_accuracy_curve(curve)

        @test curve.label == "F1_SEFIE_Jet_Direct_phi0_vs_Feko"
        @test curve.theta_deg == [-180.0, 0.0, 180.0]
        @test curve.rcs_model_dBsm == [1.0, 2.0, 3.0]
        @test curve.rcs_reference_dBsm == [1.5, 1.0, 2.0]
        @test curve.diff_dB == [-0.5, 1.0, 1.0]
        @test summary.n_points == 3
        @test summary.rmse_dB ≈ sqrt((0.25 + 1.0 + 1.0) / 3) atol = 1e-12
        @test accuracy_curve_group(curve.label) == "F1_SEFIE_Jet_Direct"
        @test accuracy_curve_cut(curve.label) == "phi0"

        alt_accuracy_csv = joinpath(tmpdir, "accuracy_alt.csv")
        open(alt_accuracy_csv, "w") do io
            println(io, "theta_deg,rcs_ems_dBsm,rcs_ref_dBsm,diff_dB")
            println(io, "0.0,4.0,3.0,1.0")
        end
        alt_curve = load_accuracy_curve(alt_accuracy_csv)
        @test alt_curve.rcs_model_dBsm == [4.0]
        @test alt_curve.rcs_reference_dBsm == [3.0]
        @test accuracy_curve_group(alt_curve.label) == "accuracy_alt"
        @test accuracy_curve_cut(alt_curve.label) == "full"

        perf_csv = joinpath(tmpdir, "performance.csv")
        open(perf_csv, "w") do io
            println(io, "case_name,equation,solver,N,t_mesh,t_assembly,t_precond,t_solve,t_rcs,t_total,notes")
            println(io, "Jet EFIE Direct,EFIE,LU,14559,1.0,2.0,0.0,3.0,0.5,6.5,")
            println(io, "Sphere CFIE MLFMA,CFIE,MLFMA+GMRES,26424,2.0,8.0,0.8,5.0,0.7,16.5,threaded")
        end

        perf = load_performance_results(perf_csv)
        @test length(perf) == 2
        @test perf[1].case_name == "Jet EFIE Direct"
        @test perf[1].t_total ≈ 6.5 atol = 1e-12
        @test perf[1].notes == ""
        @test perf[2].solver == "MLFMA+GMRES"
        @test perf[2].N == 26424
    end
end
