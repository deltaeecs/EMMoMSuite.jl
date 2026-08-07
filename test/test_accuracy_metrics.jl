"""
test_accuracy_metrics.jl — Phase 14 TDD: AccuracyMetrics 单元测试

测试覆盖：
  14.1a  compute_rcs_accuracy — 完全相同输入 → rmse=0, pass=true
  14.1b  compute_rcs_accuracy — 常数偏移 → bias/rmse 精确
  14.1c  compute_rcs_accuracy — 门限通过/不通过
  14.1d  compute_rcs_accuracy — 后向散射误差定位
  14.1e  compute_rcs_accuracy — 长度不匹配抛出 ArgumentError
  14.1f  compute_antenna_accuracy — 理想输入全部通过
  14.1g  compute_antenna_accuracy — 阻抗偏差大时不通过
  14.1h  print_accuracy_report — 混合列表不报错，含 PASS/FAIL 字样
"""

using Test
using EMMoMSuite
using EMMoMSuite.Accuracy: AccuracyResult, AntennaAccuracyResult,
                        compute_rcs_accuracy, compute_antenna_accuracy,
                        print_accuracy_report, extract_sphere_radius as extract_sphere_radius_public,
                        mie_dielectric_bistatic_rcs_dBsm
using EMMoMSuite: calculate_mie_rcs_dielectric_sphere
using EMMoMSuite.Accuracy.ReferenceData: extract_sphere_radius,
                                      mie_pec_bistatic_rcs_dBsm,
                                      mie_pec_rcs_dBsm

@testset "14.1 AccuracyMetrics" begin

    # 使用确定性测试数据（无需仿真）
    theta = collect(-180.0:0.5:180.0)   # 721 pts
    ref   = fill(-10.0, length(theta))  # 均匀 RCS 参考场景

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1a  完全相同输入 → 误差为 0
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1a 零误差（相同输入）" begin
        r = compute_rcs_accuracy(ref, ref, theta, "test_zero")
        @test r.rmse_dB     ≈ 0.0  atol = 1e-12
        @test r.max_err_dB  ≈ 0.0  atol = 1e-12
        @test r.mean_bias_dB ≈ 0.0 atol = 1e-12
        @test r.pass == true
        @test r.n_points == length(theta)
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1b  常数偏移 +1 dB → rmse = 1.0, bias = +1.0
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1b 常数偏移 +1 dB" begin
        ems = ref .+ 1.0
        r = compute_rcs_accuracy(ems, ref, theta, "test_offset1")
        @test r.rmse_dB      ≈ 1.0  atol = 1e-10
        @test r.max_err_dB   ≈ 1.0  atol = 1e-10
        @test r.mean_bias_dB ≈ 1.0  atol = 1e-10
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1c  门限通过/不通过
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1c 门限判断" begin
        # +1 dB 偏移，门限 2.0 → 通过
        r_pass = compute_rcs_accuracy(ref .+ 1.0, ref, theta, "pass"; threshold = 2.0)
        @test r_pass.pass == true
        @test r_pass.threshold_dB ≈ 2.0

        # +3 dB 偏移，门限 2.0 → 不通过
        r_fail = compute_rcs_accuracy(ref .+ 3.0, ref, theta, "fail"; threshold = 2.0)
        @test r_fail.pass == false
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1d  后向散射误差：在 θ = 180° 处设置大误差
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1d 后向散射误差定位" begin
        ems = copy(ref)
        # θ=180° 位于 theta 末尾（idx = 721）
        ems[end] = ref[end] + 5.0
        r = compute_rcs_accuracy(ems, ref, theta, "bs_test"; backscatter_theta = 180.0)
        @test r.backscatter_err_dB ≈ 5.0 atol = 1e-10
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1e  长度不匹配 → ArgumentError
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1e 长度不匹配异常" begin
        @test_throws ArgumentError compute_rcs_accuracy(
            rand(100), rand(99), rand(100), "bad"
        )
        @test_throws ArgumentError compute_rcs_accuracy(
            rand(99), rand(99), rand(100), "bad2"
        )
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1f  compute_antenna_accuracy — 理想输入
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1f 天线精度 — 理想输入通过" begin
        # 半波偶极子解析值
        Zin_analytic  = 73.1 + 42.5im
        Zin_emsuite   = 73.1 + 42.5im   # 完全对齐
        D_max         = 1.64             # 半波偶极子

        # 解析方向图 f(θ) = cos(π/2·cosθ)/sinθ
        theta_r = range(0.01, π - 0.01, length = 181)
        pattern = [cos(π/2 * cos(t)) / sin(t) for t in theta_r]

        r = compute_antenna_accuracy(
            Zin_emsuite, Zin_analytic,
            D_max, pattern, pattern,  # EMMoMSuite = Analytic
            collect(theta_r), "ideal_dipole";
            Zin_threshold    = 0.05,
            D_max_threshold  = 1.0,
            pattern_threshold = 1.0,
        )
        @test r.Zin_rel_err    ≈ 0.0  atol = 1e-10
        @test r.D_max_err_dBi  ≈ 0.0  atol = 1e-6
        @test r.pattern_rmse_dB ≈ 0.0 atol = 1e-10
        @test r.pass == true
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1g  compute_antenna_accuracy — 阻抗偏差大时不通过
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1g 天线精度 — 阻抗偏差大不通过" begin
        Zin_analytic = 73.0 + 0im
        Zin_bad      = 110.0 + 30im   # 相对误差 >> 5%

        theta_r = range(0.01, π - 0.01, length = 91)
        pattern = ones(length(theta_r))

        r = compute_antenna_accuracy(
            Zin_bad, Zin_analytic,
            1.64, pattern, pattern,
            collect(theta_r), "bad_impedance";
            Zin_threshold = 0.05,
        )
        @test r.Zin_rel_err > 0.05
        @test r.pass == false
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1h  print_accuracy_report — 包含 ✓ 和 ✗ 不报错
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1h print_accuracy_report 不崩溃" begin
        rcs1 = compute_rcs_accuracy(ref .+ 1.0, ref, theta, "rcs_pass"; threshold = 2.0)
        rcs2 = compute_rcs_accuracy(ref .+ 5.0, ref, theta, "rcs_fail"; threshold = 2.0)

        Zin_an = 73.0 + 0im
        theta_r = range(0.01, π - 0.01, length = 91)
        pattern = ones(length(theta_r))
        ant1 = compute_antenna_accuracy(
            73.0 + 0im, Zin_an, 1.64, pattern, pattern,
            collect(theta_r), "ant_pass")

        buf = IOBuffer()
        print_accuracy_report(buf, [rcs1, rcs2, ant1])
        report_str = String(take!(buf))

        @test occursin("✓", report_str)
        @test occursin("✗", report_str)
        @test occursin("rcs_pass", report_str)
        @test occursin("rcs_fail", report_str)
        @test occursin("ant_pass", report_str)
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1i  extract_sphere_radius — Nastran 球面网格半径提取
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1i 球半径提取" begin
        mesh_file = joinpath(@__DIR__, "..", "deps", "fixtures", "AllinOne", "meshfiles", "sphere_600MHz.nas")
        if isfile(mesh_file)
            r = extract_sphere_radius(mesh_file)
            @test 0.8 < r < 1.2
            @test extract_sphere_radius_public(mesh_file) == r
        else
            @test_skip "外部基线网格不存在，跳过半径提取测试"
        end
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1j  bistatic Mie — 规范 +z/x 情形应退化为原始散射角定义
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1j 双站 Mie 规范情形退化" begin
        radius = 0.3
        freq = 6.0e8
        theta_obs = collect(range(0.0, π, length = 37))

        ref_phi0 = mie_pec_rcs_dBsm(radius, freq, theta_obs)
        bistatic_phi0 = mie_pec_bistatic_rcs_dBsm(
            radius, freq, theta_obs, 0.0, 0.0, 0.0, [1.0, 0.0, 0.0])

        @test bistatic_phi0 ≈ ref_phi0 atol = 1e-10 rtol = 1e-10
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1k  bistatic Mie — 非 +z 入射时，全局 φ 切面参考不能强行视为相同
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1k 双站 Mie 非规范入射区分全局切面" begin
        radius = 0.3
        freq = 6.0e8
        theta_obs = collect(range(-π, π, length = 181))

        mie_phi0 = mie_pec_bistatic_rcs_dBsm(
            radius, freq, theta_obs, 0.0, π / 2, π, [0.0, 0.0, 1.0])
        mie_phi90 = mie_pec_bistatic_rcs_dBsm(
            radius, freq, theta_obs, π / 2, π / 2, π, [0.0, 0.0, 1.0])

        @test maximum(abs.(mie_phi0 .- mie_phi90)) > 1e-3
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1l  dielectric bistatic Mie — 规范 +z/x 情形退化到原始 E/H 切面
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1l 介质球双站 Mie 规范情形退化" begin
        radius = 0.3
        freq = 6.0e8
        eps_r = 4.0
        mu_r = 1.0
        theta_obs = collect(range(0.0, π, length = 37))

        rcs_e, rcs_h, _ = calculate_mie_rcs_dielectric_sphere(radius, freq, theta_obs, eps_r, mu_r)
        ref_phi0 = 10.0 .* log10.(max.(rcs_e, 1e-100))
        ref_phi90 = 10.0 .* log10.(max.(rcs_h, 1e-100))

        bistatic_phi0 = mie_dielectric_bistatic_rcs_dBsm(
            radius, freq, eps_r, mu_r, theta_obs, 0.0, 0.0, 0.0, [1.0, 0.0, 0.0])
        bistatic_phi90 = mie_dielectric_bistatic_rcs_dBsm(
            radius, freq, eps_r, mu_r, theta_obs, π / 2, 0.0, 0.0, [1.0, 0.0, 0.0])

        @test bistatic_phi0 ≈ ref_phi0 atol = 1e-10 rtol = 1e-10
        @test bistatic_phi90 ≈ ref_phi90 atol = 1e-10 rtol = 1e-10
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.1m  lossy dielectric Mie — EMMoMSuite 的 e^{-jωt} 负虚部损耗约定
    #         应映射到 B&H 原始实现所需的共轭材料参数
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.1m 介质球 Mie 有损符号约定" begin
        radius = 0.15
        freq = 3.0e8
        eps_r = 2.2 - 0.1im
        mu_r = 1.0 + 0im
        theta_obs = collect(range(0.0, π, length = 37))

        rcs_e_public, rcs_h_public, rcs_u_public = calculate_mie_rcs_dielectric_sphere(
            radius, freq, theta_obs, eps_r, mu_r)
        rcs_e_bh, rcs_h_bh, rcs_u_bh = EMMoMSuite.Utilities.MieSeries._calculate_mie_rcs_dielectric_sphere_bh(
            radius, freq, theta_obs, conj(eps_r), conj(mu_r))

        @test rcs_e_public ≈ rcs_e_bh atol = 1e-12 rtol = 1e-12
        @test rcs_h_public ≈ rcs_h_bh atol = 1e-12 rtol = 1e-12
        @test rcs_u_public ≈ rcs_u_bh atol = 1e-12 rtol = 1e-12
    end

end  # @testset "14.1 AccuracyMetrics"
