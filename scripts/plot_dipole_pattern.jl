"""
plot_dipole_pattern.jl — Half-wave Dipole Antenna Far-field Pattern

Comparison:
  - MoM (EFIE + RWG + DeltaGapSource) simulation
  - Analytical: D(theta) = [cos(pi/2 * cos(theta))/sin(theta)]^2 (normalized)

Also outputs:
  - Terminal: Input impedance Z_in, radiation efficiency eta_rad
  - Plot 1: E-plane normalized pattern (dB), polar coordinates
  - Plot 2: E-plane pattern + analytical comparison (line plot)

Usage:
  julia --project scripts/plot_dipole_pattern.jl
"""

using EMSuite
using Plots
using LinearAlgebra, Statistics, Printf

gr()

# ─────────────────────────────────────────────────────────────────
# Common Parameters: 300 MHz half-wave dipole
freq    = 300e6           # Hz
lambda  = 3e8 / freq      # 1 m
L_dip   = lambda / 2      # dipole length = lambda/2 = 0.5 m
a_wire  = 0.001           # thin wire radius 0.1 cm
Nz      = 60              # axial triangle segments (must be even for center feed)
Nphi    = 8               # circumferential segments (thin-wire approximation)
Z0      = 50.0            # reference impedance

set_frequency!(freq)

println("=" ^ 60)
@printf "半波偶极子天线  f=%.0f MHz  L=λ/2=%.3f m\n" freq/1e6 L_dip
println("=" ^ 60)

# ─────────────────────────────────────────────────────────────────
# 1. 生成圆柱形细线网格（模拟偶极子）
# ─────────────────────────────────────────────────────────────────
println("\n[1] 生成细线偶极子网格 ...")
mesh  = generate_cylinder_mesh(a_wire, L_dip, Nphi, Nz)
basis = RWGBasis(mesh)
N     = num_basis(basis)
@printf "  节点数=%d  三角形数=%d  RWG基函数数=%d\n" size(mesh.node, 2) mesh.trinum N

# ─────────────────────────────────────────────────────────────────
# 2. 查找中心 feed 边 (z ≈ 0)
#    利用 RWG.center z 坐标最小的若干条边作为馈电边
# ─────────────────────────────────────────────────────────────────
println("\n[2] 定位馈电边 ...")
feed_edges = Int[]
z_min_tol  = L_dip * 0.02   # within 1% of total length from center (tight feed)
for n in 1:N
    rwg = basis.functions[n]
    if abs(rwg.center[3]) < z_min_tol
        push!(feed_edges, n)
    end
end
@printf "  找到 %d 条馈电边\n" length(feed_edges)
isempty(feed_edges) && error("未找到馈电边，请调整 Nz/Nphi")

# ─────────────────────────────────────────────────────────────────
# 3. 组装 EFIE 矩阵 + DeltaGap 激励
# ─────────────────────────────────────────────────────────────────
println("\n[3] 组装阻抗矩阵 ...")
efie = EFIE(freq)
Z    = assemble_impedance_matrix(efie, basis)

src  = DeltaGapSource(freq, feed_edges, 1.0 + 0.0im)   # V = 1 V
V    = excitation_vector(src, basis)

println("\n[4] 直接求解 ...")
I = Z \ V

# Input impedance + S11
Z_in = input_impedance(src, I, basis)
S11  = (Z_in - Z0) / (Z_in + Z0)
S11_dB = 20 * log10(abs(S11))
I_feed = sum(I[n] * basis.functions[n].edge_length for n in feed_edges)
P_in   = 0.5 * real(conj(I_feed))
@printf "\n  Z_in = %.2f + j%.2f Ohm\n" real(Z_in) imag(Z_in)
@printf "  S11  = %.2f dB (Z0=%g Ohm)\n" S11_dB Z0
@printf "  (Theory half-wave dipole: Z_in approx 73 + j42.5 Ohm)\n"

# Far-field pattern + directivity
println("\n[5] Computing far-field pattern ...")
θs = collect(range(1e-3, π - 1e-3, 181))   # 避开奇点
ϕs = collect(range(0.0, 2π, 73))
result = antenna_directivity(θs, ϕs, I, basis; P_input = P_in)
D_dBi  = 10 .* log10.(result.D .+ 1e-30)

# E 平面 (phi=0): 取 phi 索引最近 0 的列
D_Eplane = D_dBi[:, 1]

@printf "\n  最大方向性 = %.2f dBi\n" maximum(D_dBi)
@printf "  辐射功率   = %.4e W\n" result.P_rad
@printf "  辐射效率   = %.2f %%\n" result.η_eff * 100

# ─────────────────────────────────────────────────────────────────
# 解析方向图: F(θ) = [cos(π/2·cosθ)/sinθ]²  (半波偶极子)
# ─────────────────────────────────────────────────────────────────
theta_a  = collect(range(1e-3, π - 1e-3, 360))
F_norm   = [abs(cos(π/2 * cos(t)) / sin(t))^2 for t in theta_a]
D_max_theory = 1.64     # 半波偶极子最大方向性 (理论)
D_theory = D_max_theory .* F_norm ./ maximum(F_norm)
D_theory_dBi = 10 .* log10.(D_theory .+ 1e-30)

theta_d = rad2deg.(θs)

# ─────────────────────────────────────────────────────────────────
# 图 1: 极坐标方向图
# ─────────────────────────────────────────────────────────────────
p1 = plot(θs, max.(D_Eplane, -30);
    proj        = :polar,
    label       = "MoM (EFIE)",
    lc          = :steelblue, lw = 2,
    title       = "Half-wave Dipole Pattern (E-plane, polar)\nf=300 MHz",
    legend      = :topright,
    ylims       = (-30, 4),
)
plot!(p1, theta_a, max.(D_theory_dBi, -30);
    label = "Analytical",
    lc = :black, lw = 1.5, ls = :dash,
)

# ─────────────────────────────────────────────────────────────────
# 图 2: 线形对比
# ─────────────────────────────────────────────────────────────────
p2 = plot(theta_d, D_Eplane;
    label  = "MoM (EFIE)",
    lc     = :steelblue, lw = 2,
    xlabel = "theta (deg)",
    ylabel = "Directivity (dBi)",
    title  = "Half-wave Dipole E-plane Pattern vs Analytical",
    legend = :topright,
    xticks = 0:30:180,
    ylims  = (-15, 4),
)
plot!(p2, rad2deg.(theta_a), D_theory_dBi;
    label = "Analytical [cos(pi/2*cos(theta))/sin(theta)]^2  (D_max=2.15 dBi)",
    lc = :black, lw = 1.5, ls = :dash,
)
annotate!(p2, 90, -12,
    text("Z_in=$(@sprintf(\"%.0f\",real(Z_in)))+j$(@sprintf(\"%.0f\",imag(Z_in))) Ohm  " *
         "S11=$(@sprintf(\"%.1f\",S11_dB)) dB", 9, :gray))

combined = plot(p1, p2, layout = (1, 2), size = (1100, 470), dpi = 120)
out = joinpath(dirname(@__DIR__), "docs", "images", "dipole_pattern.png")
savefig(combined, out)
println("\nImage saved: $out")
