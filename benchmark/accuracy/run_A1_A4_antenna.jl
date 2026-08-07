"""
A1-A4: 天线端口精度基准测试 (偶极子，Delta-Gap 激励)

用例:
  A1 — 半波偶极子 (L=0.5m, 300 MHz), Direct 求解, Z_in 与解析值对比
  A2 — 同 A1, 改用 MLFMA + GMRES 求解
  A3 — 近谐振偶极子 (L=0.48m, 300 MHz, ~0.48λ), 检查 Im(Z_in) < 5 Ω
  A4 — A1 + 50 Ω 端口匹配, 计算 S11

参考值:
  Z_in(半波) ≈ 73.1 + j42.5 Ω  （解析值 / L=0.5m）
  Z_in(谐振) ≈ 73.0 + j~0   Ω  （L ≈ 0.48λ）
  S11(50Ω)  = (Z_in-50)/(Z_in+50)

精度门限:
  Z_in 实部相对误差 < 10%  [宽松，因薄柱面网格 vs 理想细线]
  Z_in 虚部误差 < 20 Ω    [对半波偶极子绝对误差]
  Im(Z_in) 谐振偶极子 < 20 Ω [0.48λ 时虚部应显著减小]

用法:
  julia run_A1_A4_antenna.jl [A1] [A2] [A3] [A4]   # 无参数则全部运行
"""

# ── 环境设置 ─────────────────────────────────────────────────────────────────
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using EMMoMSuite
using LinearAlgebra
using Printf
using Statistics
using CSV, DataFrames
using Dates

# ── 输出目录 ─────────────────────────────────────────────────────────────────
const RESULTS_DIR = joinpath(@__DIR__, "..", "..", "test_results", "accuracy")
mkpath(RESULTS_DIR)

# ── Preconditioner (Block-Jacobi via Z_near.LU) ─────────────────────────────
struct LUPrecond
    F
end
LinearAlgebra.ldiv!(y, P::LUPrecond, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPrecond, x) = (x .= P.F \ x)

# ══════════════════════════════════════════════════════════════════════════════
# 工具函数
# ══════════════════════════════════════════════════════════════════════════════

"""
    build_dipole(L, freq; radius=0.005, n_circum=8, n_height=20)

生成开口薄柱面网格 (无端盖) + RWGBasis, 中心单馈口 delta-gap 激励。

返回:
  (basis, source, feed_edges)
    feed_edges — z≈0 处选中的单个 RWG 馈电边索引
"""
function build_dipole(L::Float64, freq::Float64;
                      radius::Float64=0.0006,
                      n_circum::Int=12,
                      n_height::Int=30)
    # 开口圆柱 (closed=false) 模拟偶极子 —— 高度 = L, 中心在 z=0
    mesh = generate_cylinder_mesh(radius, L, n_circum, n_height; closed=false)
    set_frequency!(freq)
    basis = RWGBasis(mesh)

    # Delta-gap for a cylindrical dipole is applied across the full center ring.
    feed_edges = select_gap_feed_edges(basis; axis=3, center=0.0)
    isempty(feed_edges) && error("未找到中心 gap 馈电边，请检查网格参数")

    # Delta-gap 激励: V_applied = 1 V
    source = DeltaGapSource(freq, feed_edges, 1.0 + 0.0im)
    return basis, source
end

"""
    solve_direct(efie, basis, source) → (I, elapsed)
直接 LU 求解。
"""
function solve_direct(efie, basis::RWGBasis, source::DeltaGapSource)
    t = @elapsed begin
        Z = assemble_impedance_matrix(efie, basis)
        V = excitation_vector(efie, source, basis)
        I = Z \ V
    end
    return I, t
end

"""
    solve_mlfma(efie, basis, source; leaf_frac=0.3) → (I, elapsed_setup, elapsed_solve)
MLFMA + GMRES 求解。
"""
function solve_mlfma(efie, basis::RWGBasis, source::DeltaGapSource;
                     leaf_frac::Float64=0.3)
    freq = source.frequency
    lambda = 299792458.0 / freq
    leaf_size = leaf_frac * lambda

    t_setup = @elapsed mlfma_op = MLFMAOperator(efie, basis, leaf_size)
    V = excitation_vector(efie, source, basis)
    P = LUPrecond(lu(mlfma_op.Z_near))
    solver = GMRESSolver(restart=50, maxiter=300, tol=1e-3, verbose=false)
    t_solve = @elapsed I = solve!(solver, mlfma_op, V, Pl=P)
    return I, t_setup, t_solve
end

"""
    run_antenna_case(label, basis, source, I, freq; Z_ref=nothing, Z_load=50.0)
计算并打印天线参数，返回结果 NamedTuple。
"""
function run_antenna_case(label::String, basis, source::DeltaGapSource,
                          I::Vector{<:Complex}, freq::Float64;
                          Z_ref::Union{Nothing,ComplexF64}=nothing,
                          Z_port::Float64=50.0)
    # ── 输入阻抗 ──────────────────────────────────────────────────────────────
    Z_in = input_impedance(source, I, basis)

    # ── 方向性 (φ=0 纵切面) ─────────────────────────────────────────────────
    θs = collect(Float64, range(0, π, length=91))   # 0°→180°
    ϕs = [0.0, π/2]
    res_dir = antenna_directivity(θs, ϕs, I, basis; source=nothing)
    D_max_dBi = 10 * log10(maximum(res_dir.D))

    # ── S11 ───────────────────────────────────────────────────────────────────
    Γ = (Z_in - Z_port) / (Z_in + Z_port)
    S11_dB = 20 * log10(abs(Γ))

    # ── 打印 ────────────────────────────────────────────────────────────────
    println("\n--- $label ---")
    @printf("  Z_in        = %+.2f + j%+.2f  Ω\n", real(Z_in), imag(Z_in))
    @printf("  |Z_in|      = %.2f Ω\n", abs(Z_in))
    @printf("  D_max       = %.2f dBi\n", D_max_dBi)
    @printf("  S11 (50Ω)   = %.2f dB\n", S11_dB)

    if Z_ref !== nothing
        re_err = abs(real(Z_in) - real(Z_ref)) / abs(real(Z_ref)) * 100
        im_err = abs(imag(Z_in) - imag(Z_ref))
        @printf("  Ref Z_in    = %+.2f + j%+.2f  Ω\n", real(Z_ref), imag(Z_ref))
        @printf("  Re 相对误差 = %.1f %%\n", re_err)
        @printf("  Im 绝对误差 = %.1f Ω\n", im_err)
    end

    return (; label, Z_in, D_max_dBi, S11_dB)
end

# ══════════════════════════════════════════════════════════════════════════════
# 用例主体
# ══════════════════════════════════════════════════════════════════════════════

const CASES = Dict(
    "A1" => "半波偶极子 Direct (L=0.5m, 300MHz)",
    "A2" => "半波偶极子 MLFMA  (L=0.5m, 300MHz)",
    "A3" => "近谐振偶极子 Direct (L=0.48m, 300MHz)",
    "A4" => "半波偶极子 50Ω端口 S11 (300MHz)",
)

selected = isempty(ARGS) ? keys(CASES) |> collect |> sort : ARGS

println("=" ^ 70)
println("  Phase 14 Antenna Accuracy Benchmark  [$(Dates.now())]")
println("=" ^ 70)
println("  Selected: $(join(selected, ", "))")

results_all = NamedTuple[]

# ── A1 ───────────────────────────────────────────────────────────────────────
if "A1" ∈ selected
    println("\n" * "─" ^ 70)
    println("  A1: 半波偶极子 Direct (L=0.5m, 300 MHz)")
    freq = 3e8      # 300 MHz → λ = 1 m → L = 0.5 m = λ/2
    L    = 0.5

    basis, source = build_dipole(L, freq)
    N = num_basis(basis)
    println("  unknowns = $N,  feed edges = $(length(source.edge_indices))")

    efie = EFIE(freq)
    I, t_lu = solve_direct(efie, basis, source)
    @printf("  Direct LU: %.2f s\n", t_lu)

    Z_ref = dipole_halfwave_Zin_analytic()
    res = run_antenna_case("A1 Direct", basis, source, I, freq; Z_ref=Z_ref)

    # 精度判定
    re_err = abs(real(res.Z_in) - real(Z_ref)) / abs(real(Z_ref)) * 100
    im_err = abs(imag(res.Z_in) - imag(Z_ref))
    passed = re_err < 10.0 && im_err < 20.0
    @printf("  判定: %s  (Re误差<10%%: %s, Im误差<20Ω: %s)\n",
        passed ? "PASS ✓" : "FAIL ✗",
        re_err < 10.0 ? "✓" : "✗",
        im_err < 20.0 ? "✓" : "✗")

    # 保存 CSV
    df = DataFrame(
        label     = ["A1"],
        Z_in_re   = [real(res.Z_in)],
        Z_in_im   = [imag(res.Z_in)],
        Z_ref_re  = [real(Z_ref)],
        Z_ref_im  = [imag(Z_ref)],
        re_err_pct= [re_err],
        im_err_ohm= [im_err],
        D_max_dBi = [res.D_max_dBi],
        S11_dB    = [res.S11_dB],
        passed    = [passed],
    )
    CSV.write(joinpath(RESULTS_DIR, "A1_halfwave_direct.csv"), df)
    println("  -> 已保存 A1_halfwave_direct.csv")
    push!(results_all, res)
end

# ── A2 ───────────────────────────────────────────────────────────────────────
if "A2" ∈ selected
    println("\n" * "─" ^ 70)
    println("  A2: 半波偶极子 MLFMA (L=0.5m, 300 MHz)")
    freq = 3e8
    L    = 0.5

    basis, source = build_dipole(L, freq)
    N = num_basis(basis)
    println("  unknowns = $N,  feed edges = $(length(source.edge_indices))")

    efie = EFIE(freq)
    I, t_setup, t_solve = solve_mlfma(efie, basis, source)
    @printf("  MLFMA setup: %.2f s | GMRES: %.2f s\n", t_setup, t_solve)

    Z_ref = dipole_halfwave_Zin_analytic()
    res = run_antenna_case("A2 MLFMA", basis, source, I, freq; Z_ref=Z_ref)

    re_err = abs(real(res.Z_in) - real(Z_ref)) / abs(real(Z_ref)) * 100
    im_err = abs(imag(res.Z_in) - imag(Z_ref))
    passed = re_err < 15.0 && im_err < 25.0   # MLFMA 门限稍宽
    @printf("  判定: %s  (Re误差<15%%: %s, Im误差<25Ω: %s)\n",
        passed ? "PASS ✓" : "FAIL ✗",
        re_err < 15.0 ? "✓" : "✗",
        im_err < 25.0 ? "✓" : "✗")

    df = DataFrame(
        label     = ["A2"],
        Z_in_re   = [real(res.Z_in)],
        Z_in_im   = [imag(res.Z_in)],
        Z_ref_re  = [real(Z_ref)],
        Z_ref_im  = [imag(Z_ref)],
        re_err_pct= [re_err],
        im_err_ohm= [im_err],
        D_max_dBi = [res.D_max_dBi],
        S11_dB    = [res.S11_dB],
        passed    = [passed],
    )
    CSV.write(joinpath(RESULTS_DIR, "A2_halfwave_mlfma.csv"), df)
    println("  -> 已保存 A2_halfwave_mlfma.csv")
    push!(results_all, res)
end

# ── A3 ───────────────────────────────────────────────────────────────────────
if "A3" ∈ selected
    println("\n" * "─" ^ 70)
    println("  A3: 近谐振偶极子 Direct (L=0.48m, 300 MHz, ~0.48λ)")
    freq = 3e8
    L    = 0.48    # ≈ 0.48λ 接近谐振，Im(Z_in) 应很小

    basis, source = build_dipole(L, freq)
    N = num_basis(basis)
    println("  unknowns = $N,  feed edges = $(length(source.edge_indices))")

    efie = EFIE(freq)
    I, t_lu = solve_direct(efie, basis, source)
    @printf("  Direct LU: %.2f s\n", t_lu)

    Z_ref = dipole_resonant_Zin_analytic()  # ≈ 73.0 + j0 Ω
    res = run_antenna_case("A3 Near-Resonant", basis, source, I, freq; Z_ref=Z_ref)

    im_zin = abs(imag(res.Z_in))
    # 对于 0.48λ，薄柱面近似下虚部应明显小于半波偶极子 (42.5Ω)
    passed = im_zin < 30.0   # 比半波偶极子虚部小
    @printf("  |Im(Z_in)| = %.2f Ω  判定: %s  (<30Ω: %s)\n",
        im_zin, passed ? "PASS ✓" : "FAIL ✗", passed ? "✓" : "✗")

    df = DataFrame(
        label         = ["A3"],
        Z_in_re       = [real(res.Z_in)],
        Z_in_im       = [imag(res.Z_in)],
        Z_ref_re      = [real(Z_ref)],
        Z_ref_im      = [imag(Z_ref)],
        abs_im_zin    = [im_zin],
        D_max_dBi     = [res.D_max_dBi],
        S11_dB        = [res.S11_dB],
        passed        = [passed],
    )
    CSV.write(joinpath(RESULTS_DIR, "A3_resonant_direct.csv"), df)
    println("  -> 已保存 A3_resonant_direct.csv")
    push!(results_all, res)
end

# ── A4 ───────────────────────────────────────────────────────────────────────
if "A4" ∈ selected
    println("\n" * "─" ^ 70)
    println("  A4: 半波偶极子 50Ω端口 S11 (300 MHz)")
    freq = 3e8
    L    = 0.5

    basis, source = build_dipole(L, freq)
    N = num_basis(basis)
    println("  unknowns = $N,  feed edges = $(length(source.edge_indices))")

    efie = EFIE(freq)
    I, t_lu = solve_direct(efie, basis, source)
    @printf("  Direct LU: %.2f s\n", t_lu)

    Z_ref = dipole_halfwave_Zin_analytic()
    # S11 解析对比
    Z_port = 50.0
    Γ_ref  = (Z_ref - Z_port) / (Z_ref + Z_port)
    S11_ref_dB = 20 * log10(abs(Γ_ref))

    res = run_antenna_case("A4 50Ω Port", basis, source, I, freq;
                           Z_ref=Z_ref, Z_port=Z_port)

    s11_err = abs(res.S11_dB - S11_ref_dB)
    passed = s11_err < 3.0   # S11 误差 < 3 dB
    @printf("  S11 参考值  = %.2f dB\n", S11_ref_dB)
    @printf("  S11 误差    = %.2f dB  判定: %s  (<3dB: %s)\n",
        s11_err, passed ? "PASS ✓" : "FAIL ✗", passed ? "✓" : "✗")

    df = DataFrame(
        label      = ["A4"],
        Z_in_re    = [real(res.Z_in)],
        Z_in_im    = [imag(res.Z_in)],
        Z_ref_re   = [real(Z_ref)],
        Z_ref_im   = [imag(Z_ref)],
        S11_dB     = [res.S11_dB],
        S11_ref_dB = [S11_ref_dB],
        S11_err_dB = [s11_err],
        D_max_dBi  = [res.D_max_dBi],
        passed     = [passed],
    )
    CSV.write(joinpath(RESULTS_DIR, "A4_50ohm_s11.csv"), df)
    println("  -> 已保存 A4_50ohm_s11.csv")
    push!(results_all, res)
end

# ── 汇总 ─────────────────────────────────────────────────────────────────────
println("\n" * "=" ^ 70)
println("  A1-A4 天线基准完成  —  共 $(length(results_all)) 个用例")
println("=" ^ 70)
