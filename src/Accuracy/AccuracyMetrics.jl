"""
    AccuracyMetrics — Phase 14 精度指标计算

定义 `AccuracyResult` 和 `AntennaAccuracyResult` 结构体，
提供 `compute_rcs_accuracy`、`compute_antenna_accuracy` 计算函数，
以及 `print_accuracy_report` 文本报告函数。
"""
module AccuracyMetrics

using LinearAlgebra
using Printf
using Statistics: mean

export AccuracyResult, AntennaAccuracyResult,
       compute_rcs_accuracy, compute_antenna_accuracy,
       print_accuracy_report

# ─────────────────────────────────────────────────────────────────────────────
# RCS 精度结构体
# ─────────────────────────────────────────────────────────────────────────────

"""
    AccuracyResult

RCS 精度对比结果（单个 φ 切面或两切面组合）。

# 字段
- `label`:               测试标识（如 "F1_SEFIE_Jet_Direct_phi0"）
- `n_points`:            对比点数
- `rmse_dB`:             均方根误差（dB）= √(mean((EMMoMSuite - Ref)²))
- `max_err_dB`:          最大绝对误差（dB）= max(|EMMoMSuite - Ref|)
- `mean_bias_dB`:        均值偏差（dB）= mean(EMMoMSuite - Ref)
- `backscatter_err_dB`:  后向散射处（θ=180°）的单点误差（dB）
- `pass`:                是否通过精度门限
- `threshold_dB`:        精度门限（dB）
"""
struct AccuracyResult
    label             :: String
    n_points          :: Int
    rmse_dB           :: Float64
    max_err_dB        :: Float64
    mean_bias_dB      :: Float64
    backscatter_err_dB :: Float64
    pass              :: Bool
    threshold_dB      :: Float64
end

# ─────────────────────────────────────────────────────────────────────────────
# 天线精度结构体
# ─────────────────────────────────────────────────────────────────────────────

"""
    AntennaAccuracyResult

天线辐射精度对比结果。

# 字段
- `label`:           测试标识（如 "A1_Dipole_LumpedPort_Direct"）
- `Zin_rel_err`:     输入阻抗相对误差 = |ΔZ_in| / |Z_analytic|（无量纲）
- `D_max_err_dBi`:   最大方向性误差（dBi）= |D_max_EMSuite - 1.64| in (10log10)
- `pattern_rmse_dB`: 方向图 RMSE（dB，E 面归一化对数坐标）
- `S11_err_dB`:      S11 误差（dB）
- `Im_Zin`:          Im(Z_in) 绝对值（Ω），谐振测试专用
- `pass`:            是否通过全部精度门限
- `detail`:          详细描述字典（可选）
"""
struct AntennaAccuracyResult
    label          :: String
    Zin_rel_err    :: Float64
    D_max_err_dBi  :: Float64
    pattern_rmse_dB :: Float64
    S11_err_dB     :: Float64
    Im_Zin         :: Float64
    pass           :: Bool
    detail         :: Dict{String,Any}
end

AntennaAccuracyResult(label, Zin_r, D_err, pat_rmse, s11_e, im_z, pass) =
    AntennaAccuracyResult(label, Zin_r, D_err, pat_rmse, s11_e, im_z, pass, Dict{String,Any}())

# ─────────────────────────────────────────────────────────────────────────────
# RCS 精度计算
# ─────────────────────────────────────────────────────────────────────────────

"""
    compute_rcs_accuracy(
        rcs_emsuite_dBsm, rcs_ref_dBsm, theta_deg, label;
        threshold = 2.0,
        backscatter_theta = 180.0,
    ) -> AccuracyResult

计算 EMMoMSuite 仿真 RCS 与参考 RCS 之间的精度指标。

# 参数
- `rcs_emsuite_dBsm`: EMMoMSuite 计算结果（dBsm）
- `rcs_ref_dBsm`:     参考值（Feko 或 Mie，dBsm）
- `theta_deg`:        对应的 θ 角度值（度）
- `label`:            测试标识字符串
- `threshold`:        RMSE 通过门限（dB），默认 2.0 dB
- `backscatter_theta`:後向散射角度（度），默认 180.0°

# 返回
`AccuracyResult`

# 注意
两个向量必须长度相同且角度对齐。
"""
function compute_rcs_accuracy(
    rcs_emsuite_dBsm::AbstractVector{<:Real},
    rcs_ref_dBsm::AbstractVector{<:Real},
    theta_deg::AbstractVector{<:Real},
    label::AbstractString;
    threshold::Float64 = 2.0,
    backscatter_theta::Float64 = 180.0,
)
    length(rcs_emsuite_dBsm) == length(rcs_ref_dBsm) ||
        throw(ArgumentError("EMMoMSuite 与参考 RCS 长度不同"))
    length(rcs_emsuite_dBsm) == length(theta_deg) ||
        throw(ArgumentError("RCS 与 theta_deg 长度不同"))

    diff = rcs_emsuite_dBsm .- rcs_ref_dBsm

    rmse      = sqrt(mean(diff .^ 2))
    max_err   = maximum(abs.(diff))
    bias      = mean(diff)

    # 后向散射误差：找最近的点
    _, idx_bs = findmin(abs.(theta_deg .- backscatter_theta))
    bs_err = abs(diff[idx_bs])

    pass = rmse < threshold

    return AccuracyResult(
        string(label), length(diff),
        rmse, max_err, bias, bs_err, pass, threshold,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# 天线精度计算
# ─────────────────────────────────────────────────────────────────────────────

"""
    compute_antenna_accuracy(
        Zin_emsuite, Zin_analytic,
        D_max_emsuite,
        pattern_emsuite, pattern_analytic, theta_rad_vec,
        label;
        S11_emsuite = nothing,
        S11_analytic = nothing,
        Zin_threshold = 0.05,
        D_max_threshold = 1.0,
        pattern_threshold = 1.0,
    ) -> AntennaAccuracyResult

计算天线仿真精度，综合输入阻抗、方向性和方向图三项指标。
"""
function compute_antenna_accuracy(
    Zin_emsuite::Number,
    Zin_analytic::Number,
    D_max_emsuite::Real,
    pattern_emsuite::AbstractVector{<:Real},
    pattern_analytic::AbstractVector{<:Real},
    theta_rad_vec::AbstractVector{<:Real},
    label::AbstractString;
    S11_emsuite::Union{Nothing,Number} = nothing,
    S11_analytic::Union{Nothing,Number} = nothing,
    Zin_threshold::Float64 = 0.05,
    D_max_threshold::Float64 = 1.0,    # dBi
    pattern_threshold::Float64 = 1.0,  # dB RMSE
)
    # 输入阻抗相对误差
    Zin_rel_err = abs(Zin_emsuite - Zin_analytic) / max(abs(Zin_analytic), 1e-10)

    # 方向性误差（dBi）：理论最大方向性 1.64（半波偶极子）
    D_max_analytic_dBi = 10.0 * log10(1.64)
    D_max_emsuite_dBi  = 10.0 * log10(max(D_max_emsuite, 1e-10))
    D_max_err_dBi      = abs(D_max_emsuite_dBi - D_max_analytic_dBi)

    # 方向图 RMSE（对数坐标，dB，经归一化后）
    norm_e = maximum(abs.(pattern_emsuite))
    norm_a = maximum(abs.(pattern_analytic))
    p_e = norm_e > 0 ? 20.0 .* log10.(max.(abs.(pattern_emsuite) ./ norm_e, 1e-6)) : zeros(length(pattern_emsuite))
    p_a = norm_a > 0 ? 20.0 .* log10.(max.(abs.(pattern_analytic) ./ norm_a, 1e-6)) : zeros(length(pattern_analytic))

    # 只比较主波束区域（|p_a| > -30 dB，避免零点噪声）
    mask = p_a .> -30.0
    pattern_rmse_dB = if sum(mask) >= 2
        sqrt(mean((p_e[mask] .- p_a[mask]) .^ 2))
    else
        sqrt(mean((p_e .- p_a) .^ 2))
    end

    # S11（可选）
    S11_err_dB = if !isnothing(S11_emsuite) && !isnothing(S11_analytic)
        abs(20.0 * log10(max(abs(S11_emsuite), 1e-10)) -
            20.0 * log10(max(abs(S11_analytic), 1e-10)))
    else
        NaN
    end

    Im_Zin = abs(imag(Zin_emsuite))

    pass = (Zin_rel_err  < Zin_threshold) &&
           (D_max_err_dBi < D_max_threshold) &&
           (pattern_rmse_dB < pattern_threshold)

    detail = Dict{String,Any}(
        "Zin_emsuite"   => Zin_emsuite,
        "Zin_analytic"  => Zin_analytic,
        "D_max_emsuite" => D_max_emsuite,
        "S11_emsuite"   => S11_emsuite,
        "thresholds"    => (Zin_threshold, D_max_threshold, pattern_threshold),
    )

    return AntennaAccuracyResult(
        string(label), Zin_rel_err, D_max_err_dBi,
        pattern_rmse_dB, isnan(S11_err_dB) ? 0.0 : S11_err_dB,
        Im_Zin, pass, detail,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# 文本报告
# ─────────────────────────────────────────────────────────────────────────────

"""
    print_accuracy_report([io], results; title="精度对比报告")

打印或写入精度报告。接受 `AccuracyResult` 和 `AntennaAccuracyResult` 的混合向量。
"""
function print_accuracy_report(
    io::IO,
    results::AbstractVector;
    title::AbstractString = "Phase 14 精度对比报告",
)
    println(io, "\n" * "="^72)
    println(io, "  $title")
    println(io, "="^72)

    rcs_results    = filter(x -> x isa AccuracyResult,    results)
    ant_results    = filter(x -> x isa AntennaAccuracyResult, results)

    if !isempty(rcs_results)
        println(io, "\n── RCS 散射精度 ──────────────────────────────────────────────────────")
        @printf(io, "  %-38s %8s %8s %8s %6s  %s\n",
            "标签", "RMSE(dB)", "MaxErr", "Bias", "BS_Err", "PASS?")
        println(io, "  " * "-"^70)
        for r in rcs_results
            status = r.pass ? "✓" : "✗"
            @printf(io, "  %-38s %8.3f %8.3f %+8.3f %6.3f  %s\n",
                r.label, r.rmse_dB, r.max_err_dB, r.mean_bias_dB,
                r.backscatter_err_dB, status)
        end
        n_pass = count(r -> r.pass, rcs_results)
        println(io, "\n  RCS 通过: $n_pass / $(length(rcs_results))")
    end

    if !isempty(ant_results)
        println(io, "\n── 天线端口精度 ──────────────────────────────────────────────────────")
        @printf(io, "  %-32s %10s %10s %10s %8s  %s\n",
            "标签", "ΔZin(%)", "D_err(dBi)", "Pat_RMSE", "Im(Zin)", "PASS?")
        println(io, "  " * "-"^72)
        for r in ant_results
            status = r.pass ? "✓" : "✗"
            @printf(io, "  %-32s %10.2f %10.3f %10.3f %8.2f  %s\n",
                r.label, r.Zin_rel_err * 100, r.D_max_err_dBi,
                r.pattern_rmse_dB, r.Im_Zin, status)
        end
        n_pass = count(r -> r.pass, ant_results)
        println(io, "\n  天线通过: $n_pass / $(length(ant_results))")
    end

    n_total = length(results)
    n_pass  = count(results) do r
        r isa AccuracyResult ? r.pass : (r isa AntennaAccuracyResult ? r.pass : false)
    end
    println(io, "\n  总通过率: $n_pass / $n_total")
    println(io, "="^72 * "\n")
end

print_accuracy_report(results::AbstractVector; kwargs...) =
    print_accuracy_report(stdout, results; kwargs...)

end # module AccuracyMetrics
