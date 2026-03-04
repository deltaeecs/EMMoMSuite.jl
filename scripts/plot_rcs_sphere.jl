"""
plot_rcs_sphere.jl — 球体双站 RCS 可视化

对比：
  - PEC 球 (EFIE + RWG) vs Mie 解析解
  - 均匀介质球 (PMCHW + RWG) vs Mie 解析解 (E 面 & H 面)

输出：
  - 图1: PEC 球 RCS (E 面, theta = 0→180°)
  - 图2: 介质球 PMCHW RCS (E 面 & H 面)

用法:
  julia --project scripts/plot_rcs_sphere.jl
"""

using EMSuite
using Plots
using LinearAlgebra, Statistics, Printf

gr()   # GR 后端，无需 GUI，可输出 PNG

# ─────────────────────────────────────────────────────────────────
# 公共参数
# ─────────────────────────────────────────────────────────────────
radius  = 0.5          # m  (电尺寸 ka ≈ π)
freq    = 300e6        # Hz
theta_v = collect(range(0.0, π, 181))   # 0°→180°
theta_d = rad2deg.(theta_v)

# ─────────────────────────────────────────────────────────────────
# 图 1: PEC 球 — EFIE vs Mie
# ─────────────────────────────────────────────────────────────────
println("=" ^ 60)
println("图 1: PEC 球 EFIE vs Mie")

mesh_pec   = generate_sphere_mesh(radius, 16, 32)
basis_pec  = RWGBasis(mesh_pec)
efie       = EFIE(freq)
Z_pec      = assemble_impedance_matrix(efie, basis_pec)
pw_pec     = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])   # +x 入射, z 极化
V_pec      = excitation_vector(efie, pw_pec, basis_pec)
I_pec      = Z_pec \ V_pec

set_frequency!(freq)
rcs_pec, _, rcs_pec_db = radarCrossSection(theta_v, [0.0], I_pec, basis_pec)
rcs_pec_Eplane = 10 .* log10.(rcs_pec[1, :, 1] .+ 1e-30)

mie_pec     = calculate_mie_rcs_pec_sphere(radius, freq, theta_v)
mie_pec_db  = 10 .* log10.(mie_pec .+ 1e-30)

err_rms = sqrt(mean((rcs_pec_Eplane .- mie_pec_db).^2))
@printf "  EFIE N=%d  RMS误差=%.2f dB\n" num_basis(basis_pec) err_rms

p1 = plot(theta_d, mie_pec_db,
    label = "Mie 解析解",
    lw = 2, lc = :black, ls = :dash,
    xlabel = "观测角 θ (°)",
    ylabel = "RCS (dBsm)",
    title  = "PEC 球双站 RCS (EFIE vs Mie)\nradius=0.5 m, f=300 MHz",
    legend = :bottomright,
    xticks = 0:30:180,
    ylims  = (-30, 5),
)
plot!(p1, theta_d, rcs_pec_Eplane,
    label = "EFIE (N=$(num_basis(basis_pec)))",
    lw = 1.5, lc = :steelblue,
)
annotate!(p1, 90, -27, text("RMS 误差 = $(@sprintf("%.2f", err_rms)) dB", 9, :gray))

# ─────────────────────────────────────────────────────────────────
# 图 2: 介质球 — PMCHW vs Mie (E 面 & H 面)
# ─────────────────────────────────────────────────────────────────
println("=" ^ 60)
println("图 2: 介质球 PMCHW vs Mie (E 面 & H 面)")

eps_r = 4.0;  mu_r = 1.0
mesh_die  = generate_sphere_mesh(radius, 16, 32)
basis_die = RWGBasis(mesh_die)
pmchw     = PMCHW(freq, eps_r, mu_r)
Z_die     = assemble_impedance_matrix(pmchw, basis_die)
pw_die    = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
V_die     = excitation_vector(pmchw, pw_die, basis_die)
I_die     = Z_die \ V_die

rcs_Ep, _, _ = radarCrossSection(theta_v, [0.0],  I_die, basis_die,
                                  pmchw.k0, pmchw.eta0)
rcs_Hp, _, _ = radarCrossSection(theta_v, [π/2], I_die, basis_die,
                                  pmchw.k0, pmchw.eta0)
mom_Eplane = 10 .* log10.(rcs_Ep[1, :, 1] .+ 1e-30)
mom_Hplane = 10 .* log10.(rcs_Hp[2, :, 1] .+ 1e-30)

mie_E, mie_H, _ = calculate_mie_rcs_dielectric_sphere(radius, freq, theta_v, eps_r, mu_r)
mie_Edb = 10 .* log10.(mie_E .+ 1e-30)
mie_Hdb = 10 .* log10.(mie_H .+ 1e-30)

err_E = sqrt(mean((mom_Eplane .- mie_Edb).^2))
err_H = sqrt(mean((mom_Hplane .- mie_Hdb).^2))
@printf "  PMCHW N=%d  E面RMS=%.2f dB  H面RMS=%.2f dB\n" num_basis(basis_die) err_E err_H

p2 = plot(theta_d, mie_Edb,
    label = "Mie E面",
    lw = 2, lc = :black, ls = :dash,
    xlabel = "观测角 θ (°)",
    ylabel = "RCS (dBsm)",
    title  = "均匀介质球双站 RCS (PMCHW vs Mie)\nradius=0.5 m, f=300 MHz, ε_r=4",
    legend = :bottomright,
    xticks = 0:30:180,
    ylims  = (-30, 5),
)
plot!(p2, theta_d, mie_Hdb,
    label = "Mie H面",
    lw = 2, lc = :gray, ls = :dash,
)
plot!(p2, theta_d, mom_Eplane,
    label = "PMCHW E面 (N=$(num_basis(basis_die)))",
    lw = 1.5, lc = :steelblue,
)
plot!(p2, theta_d, mom_Hplane,
    label = "PMCHW H面",
    lw = 1.5, lc = :tomato,
)
annotate!(p2, 90, -27, text("E面RMS=$(@sprintf("%.2f",err_E)) dB  H面RMS=$(@sprintf("%.2f",err_H)) dB", 9, :gray))

# ─────────────────────────────────────────────────────────────────
# 合并保存
# ─────────────────────────────────────────────────────────────────
combined = plot(p1, p2, layout = (1, 2), size = (1100, 450), dpi = 120)
out_path = joinpath(dirname(@__DIR__), "docs", "images", "rcs_sphere_comparison.png")
savefig(combined, out_path)
println("\n图像已保存: $out_path")
