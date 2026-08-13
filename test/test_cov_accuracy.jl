# test_cov_accuracy.jl — Accuracy 模块轻量覆盖
#
# 覆盖 src/Accuracy/：
#   FekoReader.jl（read_feko_rcs / split_phi_cuts）
#   ReferenceData.jl（Mie PEC/介质球 RCS、偶极子解析、球半径提取）
#   AccuracyMetrics.jl（RCS/天线精度指标、文本报告）
using Test
using EMMoMSuite
using EMMoMSuite.Accuracy
using EMMoMSuite.Accuracy.FekoReader: read_feko_rcs, split_phi_cuts
using EMMoMSuite.Accuracy.ReferenceData:
    mie_pec_rcs_dBsm, mie_pec_bistatic_rcs_dBsm,
    mie_dielectric_rcs_dBsm, mie_dielectric_bistatic_rcs_dBsm,
    dipole_halfwave_Zin_analytic, dipole_resonant_Zin_analytic,
    dipole_halfwave_farfield_analytic, extract_sphere_radius
using EMMoMSuite.Accuracy.AccuracyMetrics:
    compute_rcs_accuracy, compute_antenna_accuracy, print_accuracy_report
using LinearAlgebra

@testset "Accuracy: FekoReader 解析" begin
    fpath = tempname() * ".csv"
    open(fpath, "w") do io
        println(io, "THETA PHI FREQ RCS 0.0 0.0 RCS_LIN")
        for th in 0.0:10.0:180.0
            println(io, "$th 0.0 3.0e8 0.0 0.0 0.0 1.0e-3")
        end
    end
    try
        theta, phi, rcs_sqm, rcs_dBsm = read_feko_rcs(fpath)
        @test length(theta) == 19
        @test all(theta .≈ collect(0.0:10.0:180.0))
        @test all(phi .== 0.0)
        @test all(rcs_sqm .≈ 1.0e-3)
        @test all(rcs_dBsm .≈ -30.0)
        cuts = split_phi_cuts(theta, phi, rcs_dBsm)
        @test length(cuts) == 1
        key = first(keys(cuts))
        @test length(cuts[key].theta) == 19
        @test issorted(cuts[key].theta)
    finally
        rm(fpath; force = true)
    end
    @test_throws ArgumentError read_feko_rcs(tempname() * "_missing.csv")
end

@testset "Accuracy: Mie 参考与偶极子解析" begin
    r = 0.1
    f = 300e6
    th = collect(range(0.0, π; length = 37))

    rcs_pec = mie_pec_rcs_dBsm(r, f, th)
    @test length(rcs_pec) == 37
    @test all(isfinite, rcs_pec)

    tho = collect(range(0.0, π; length = 19))
    rcs_bi = mie_pec_bistatic_rcs_dBsm(r, f, tho, 0.0, π, 0.0, [1.0, 0.0, 0.0])
    @test length(rcs_bi) == 19
    @test all(isfinite, rcs_bi)

    rcs_diel = mie_dielectric_rcs_dBsm(r, f, 2.0, 1.0, th)
    @test length(rcs_diel) == 37
    @test all(isfinite, rcs_diel)

    rcs_diel_bi = mie_dielectric_bistatic_rcs_dBsm(r, f, 2.0, 1.0, tho, 0.0, π, 0.0, [1.0, 0.0, 0.0])
    @test length(rcs_diel_bi) == 19
    @test all(isfinite, rcs_diel_bi)

    @test dipole_halfwave_Zin_analytic() == ComplexF64(73.1, 42.5)
    @test abs(imag(dipole_resonant_Zin_analytic())) < 1e-12
    fpat = dipole_halfwave_farfield_analytic(th)
    @test length(fpat) == 37
    @test maximum(fpat) ≈ 1.0
end

@testset "Accuracy: 球半径提取" begin
    mesh = generate_sphere_mesh(0.1, 4, 8)
    nas = tempname() * ".nas"
    write_nas_mesh(nas, mesh)
    try
        @test extract_sphere_radius(nas) ≈ 0.1 atol = 1e-3
    finally
        rm(nas; force = true)
    end
end

@testset "Accuracy: 精度指标与报告" begin
    theta_deg = collect(0.0:5.0:180.0)
    ref = -30.0 .+ 0.01 .* sin.(deg2rad.(theta_deg))

    res = compute_rcs_accuracy(ref, ref, theta_deg, "identical")
    @test res.pass
    @test res.rmse_dB < 1e-12

    res2 = compute_rcs_accuracy(ref .+ 1.0, ref, theta_deg, "bias1"; threshold = 0.5)
    @test !res2.pass
    @test res2.mean_bias_dB ≈ 1.0 atol = 1e-12

    th_rad = collect(range(0.0, π; length = 37))
    ares = compute_antenna_accuracy(
        73.1 + 42.5im, 73.1 + 42.5im, 1.64,
        fill(1.0, 37), fill(1.0, 37), th_rad, "perfect",
    )
    @test ares.pass
    @test ares.Zin_rel_err < 1e-12

    io = IOBuffer()
    print_accuracy_report(io, [res, ares]; title = "coverage-test")
    text = String(take!(io))
    @test occursin("coverage-test", text)
    @test occursin("identical", text)
end
