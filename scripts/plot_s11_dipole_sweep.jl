"""
plot_s11_dipole_sweep.jl — Dipole Antenna S11 Frequency Sweep

Simulation:
  - Dipole length L = 0.5 m, frequency sweep 200~500 MHz
  - Assemble EFIE + DeltaGapSource per frequency, solve, extract Z_in -> S11
  - Annotate with theoretical resonance f_res ~= c/(2L)

Outputs:
  - Plot 1: S11 magnitude (dB) vs frequency
  - Plot 2: Input impedance Re/Im vs frequency (theory annotation)
  - Plot 3: Smith chart (normalized to Z0=50 Ohm)

Usage:
  julia --project scripts/plot_s11_dipole_sweep.jl
"""

using EMSuite
using Plots
using LinearAlgebra, Printf, Statistics

gr()

# ─────────────────────────────────────────────────────────────────
# 参数
# ─────────────────────────────────────────────────────────────────
L_dip  = 0.5       # dipole length (m)
a_wire = 0.001     # thin wire radius (m)
Nz     = 60        # axial segments (increased from 20 for accuracy)
Nphi   = 8         # circumferential segments
Z0     = 50.0      # reference impedance (Ohm)
c0     = 3e8       # speed of light

f_res_theory = c0 / (2 * L_dip)   # 理论谐振频率
freqs_sweep  = LinRange(200e6, 500e6, 31)

println("=" ^ 60)
@printf "偶极子 S11 频率扫描  L=%.3f m  f_res(theory)=%.0f MHz\n" L_dip f_res_theory/1e6
println("=" ^ 60)

# ─────────────────────────────────────────────────────────────────
# 1. 建立固定网格（频率不影响网格拓扑）
# ─────────────────────────────────────────────────────────────────
println("\n[1] 生成偶极子网格 ...")
mesh  = generate_cylinder_mesh(a_wire, L_dip, Nphi, Nz)
basis = RWGBasis(mesh)
N     = num_basis(basis)
@printf "  RWG基函数数=%d\n" N

# 找中心馈电边
z_tol = L_dip * 0.02   # tighter tolerance: +/-1% length from center
feed_edges = [n for n in 1:N if abs(basis.functions[n].center[3]) < z_tol]
@printf "  馈电边: %d 条\n" length(feed_edges)
isempty(feed_edges) && error("未找到馈电边")

# ─────────────────────────────────────────────────────────────────
# 2. 频率扫描
# ─────────────────────────────────────────────────────────────────
println("\n[2] 频率扫描 ($(length(freqs_sweep)) 点) ...")
Z_in_vec = ComplexF64[]
S11_vec  = Float64[]

for (i, freq) in enumerate(freqs_sweep)
    set_frequency!(freq)
    efie = EFIE(freq)
    Z    = assemble_impedance_matrix(efie, basis)
    src  = DeltaGapSource(freq, feed_edges, 1.0 + 0.0im)
    V    = excitation_vector(src, basis)
    I    = Z \ V

    Zin = input_impedance(src, I, basis)
    s11 = (Zin - Z0) / (Zin + Z0)

    push!(Z_in_vec, Zin)
    push!(S11_vec,  20 * log10(abs(s11) + 1e-30))

    if i % 5 == 0 || i == 1
        @printf "  f=%.0f MHz  Z_in=%6.1f+j%6.1f Ohm  S11=%.1f dB\n" freq/1e6 real(Zin) imag(Zin) last(S11_vec)
    end
end

freqs_MHz  = freqs_sweep ./ 1e6
Rin_vec    = real.(Z_in_vec)
Xin_vec    = imag.(Z_in_vec)

# 找仿真谐振频率 (|S11| 最小)
idx_res = argmin(S11_vec)
@printf "\n  仿真谐振: f=%.1f MHz  Z_in=%.1f+j%.1f Ω  S11=%.1f dB\n" freqs_MHz[idx_res] Rin_vec[idx_res] Xin_vec[idx_res] S11_vec[idx_res]
@printf "  理论谐振: f=%.1f MHz (λ/2)\n" f_res_theory/1e6

# ─────────────────────────────────────────────────────────────────
# 图 1: S11 vs 频率
# ─────────────────────────────────────────────────────────────────
p1 = plot(freqs_MHz, S11_vec;
    label  = "S11 (Z0=$Z0 Ohm)",
    lc     = :steelblue, lw = 2,
    xlabel = "Frequency (MHz)",
    ylabel = "S11 (dB)",
    title  = "Dipole S11 vs Frequency  L=$(L_dip) m",
    legend = :topright,
    ylims  = (-40, 2),
)
vline!(p1, [f_res_theory/1e6];
    label = "Theory resonance $(round(Int, f_res_theory/1e6)) MHz",
    lc = :red, ls = :dash, lw = 1.5,
)
hline!(p1, [-10.0]; label = "-10 dB reference", lc = :gray, ls = :dot)
scatter!(p1, [freqs_MHz[idx_res]], [S11_vec[idx_res]];
    label = "Sim. min $(@sprintf(\"%.1f\",S11_vec[idx_res])) dB @ $(@sprintf(\"%.0f\",freqs_MHz[idx_res])) MHz",
    mc = :red, ms = 6, msw = 0
)

# ─────────────────────────────────────────────────────────────────
# 图 2: 输入阻抗实部/虚部 vs 频率
# ─────────────────────────────────────────────────────────────────
p2 = plot(freqs_MHz, Rin_vec;
    label  = "R_in",
    lc     = :steelblue, lw = 2,
    xlabel = "Frequency (MHz)",
    ylabel = "Impedance (Ohm)",
    title  = "Dipole Input Impedance  L=$(L_dip) m",
    legend = :topleft,
    ylims  = (-200, 300),
)
plot!(p2, freqs_MHz, Xin_vec;
    label = "X_in",
    lc    = :tomato, lw = 2,
)
hline!(p2, [73.0];
    label = "Theory R_in~73 Ohm (resonance)",
    lc = :steelblue, ls = :dash, lw = 1,
)
hline!(p2, [0.0]; lc = :gray, ls = :dot, label = "X_in=0")
vline!(p2, [f_res_theory/1e6]; lc = :red, ls = :dash, lw = 1.5, label = "")

# ─────────────────────────────────────────────────────────────────
# 图 3: Smith 圆图
# ─────────────────────────────────────────────────────────────────
# 归一化 Γ (反射系数) 坐标 → ΓRe, ΓIm
Γ_vec   = (Z_in_vec .- Z0) ./ (Z_in_vec .+ Z0)
Γre     = real.(Γ_vec)
Γim     = imag.(Γ_vec)

# 画 Smith 圆图背景圆环 (r=const, x=const)
theta_bg = range(0, 2π, 200)
p3 = plot(cos.(theta_bg), sin.(theta_bg);
    lc = :lightgray, lw = 1, label = "", aspect_ratio = 1,
    xlims = (-1.1, 1.1), ylims = (-1.1, 1.1),
    title = "Smith Chart (Z0=$Z0 Ohm)",
    xlabel = "ΓRe", ylabel = "ΓIm",
    legend = :topright,
)
hline!(p3, [0.0]; lc = :lightgray, lw = 0.5, label = "")
vline!(p3, [0.0]; lc = :lightgray, lw = 0.5, label = "")

# 绘制轨迹（按频率从低到高着色）
scatter!(p3, Γre, Γim;
    zcolor  = freqs_MHz,
    m       = :circle, ms = 5, msw = 0,
    label   = "Z_in trajectory",
    colorbar_title = "Frequency (MHz)",
)
# 标记谐振点
scatter!(p3, [Γre[idx_res]], [Γim[idx_res]];
    mc = :red, ms = 8, msw = 0, label = "Resonance",
)

combined = plot(p1, p2, p3, layout = (1, 3), size = (1500, 470), dpi = 110)
out = joinpath(dirname(@__DIR__), "docs", "images", "dipole_s11_sweep.png")
savefig(combined, out)
println("\nImage saved: $out")
