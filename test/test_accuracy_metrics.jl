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
using EMSuite.Accuracy: AccuracyResult, AntennaAccuracyResult,
                        compute_rcs_accuracy, compute_antenna_accuracy,
                        print_accuracy_report

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
            D_max, pattern, pattern,  # EMSuite = Analytic
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

end  # @testset "14.1 AccuracyMetrics"
