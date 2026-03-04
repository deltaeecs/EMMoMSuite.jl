"""
plot_dipole_pattern.jl — 半波偶极子天线远场方向图

对比：
  - MoM (EFIE + RWG + DeltaGapSource) 仿真结果
  - 解析公式: D(θ) = [cos(π/2 · cosθ)/sinθ]²  (归一化)

还输出：
  - 终端打印: 输入阻抗 Z_in、辐射效率 η_rad
  - 图1: E 面归一化方向图 (dB), 极坐标
  - 图2: E 面方向图 + 解析比较 (线图)

用法:
  julia --project scripts/plot_dipole_pattern.jl
"""

using EMSuite
using Plots
using LinearAlgebra, Statistics, Printf

gr()

# ─────────────────────────────────────────────────────────────────
# 参数: 300 MHz 半波偶极子
# ─────────────────────────────────────────────────────────────────
freq    = 300e6           # Hz
lambda  = 3e8 / freq      # 1 m
L_dip   = lambda / 2      # 偶极子长度 = λ/2 = 0.5 m
a_wire  = 0.001           # 细线半径 0.1 cm
Nz      = 20              # 轴向三角形分段数 (需为偶数，中心处有 feed 边)
Nphi    = 6               # 周向分段数（细线近似）
Z0      = 50.0            # 参考阻抗

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
@printf "  节点数=%d  三角形数=%d  RWG基函数数=%d\n" mesh.nodenum mesh.trinum N

# ─────────────────────────────────────────────────────────────────
# 2. 查找中心 feed 边 (z ≈ 0)
#    利用 RWG.center z 坐标最小的若干条边作为馈电边
# ─────────────────────────────────────────────────────────────────
println("\n[2] 定位馈电边 ...")
feed_edges = Int[]
z_min_tol  = L_dip * 0.06   # 偶极子总长的 6% 范围内
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

# ─────────────────────────────────────────────────────────────────
# 4. 输入阻抗 + S11
# ─────────────────────────────────────────────────────────────────
Z_in = input_impedance(src, I, basis)
S11  = (Z_in - Z0) / (Z_in + Z0)
S11_dB = 20 * log10(abs(S11))
# P_in = 0.5 * Re(V · I_feed*)  (V=1, 馈电电流 I_feed = V/Z_in)
I_feed = sum(I[n] for n in feed_edges)
P_in   = 0.5 * real(conj(I_feed))   # = 0.5 Re(V · I_feed*), V=1
@printf "\n  Z_in = %.2f + j%.2f Ω\n" real(Z_in) imag(Z_in)
@printf "  S11  = %.2f dB (Z0=%g Ω)\n" S11_dB Z0
@printf "  （理论半波偶极子: Z_in ≈ 73 + j42.5 Ω）\n"

# ─────────────────────────────────────────────────────────────────
# 5. 远场方向图 + 方向性
# ─────────────────────────────────────────────────────────────────
println("\n[5] 计算远场方向图 ...")
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
    title       = "半波偶极子方向图 (极坐标E面)\nf=300 MHz",
    legend      = :topright,
    ylims       = (-30, 4),
)
plot!(p1, theta_a, max.(D_theory_dBi, -30);
    label = "解析公式",
    lc = :black, lw = 1.5, ls = :dash,
)

# ─────────────────────────────────────────────────────────────────
# 图 2: 线形对比
# ─────────────────────────────────────────────────────────────────
p2 = plot(theta_d, D_Eplane;
    label  = "MoM (EFIE)",
    lc     = :steelblue, lw = 2,
    xlabel = "θ (°)",
    ylabel = "方向性 (dBi)",
    title  = "半波偶极子 E 面方向图 vs 解析解",
    legend = :topright,
    xticks = 0:30:180,
    ylims  = (-15, 4),
)
plot!(p2, rad2deg.(theta_a), D_theory_dBi;
    label = "解析 [cos(π/2·cosθ)/sinθ]²  (D_max=2.15 dBi)",
    lc = :black, lw = 1.5, ls = :dash,
)
annotate!(p2, 90, -12,
    text("Z_in=$(@sprintf("%.0f",real(Z_in)))+j$(@sprintf("%.0f",imag(Z_in))) Ω  " *
         "S11=$(@sprintf("%.1f",S11_dB)) dB", 9, :gray))

combined = plot(p1, p2, layout = (1, 2), size = (1100, 470), dpi = 120)
out = joinpath(@__DIR__, "dipole_pattern.png")
savefig(combined, out)
println("\n图像已保存: $out")
