"""
verify_pmchw_sphere.jl — Phase 22.5 验证脚本

用 PMCHW 公式求解均匀介质球 RCS，与 Mie 级数对比。

约定 (与 PlaneWave/RadiationIntegral 一致):
  入射波: PlaneWave(freq, 0, 0, [1,0,0]) → +z 方向传播，x 极化
  观测角: θ ∈ [0,π], φ=0 (E-plane), φ=π/2 (H-plane)
  Mie E-plane ↔ RCSθsϕs[1,:,1] (θ 分量，φ_obs=0)
  Mie H-plane ↔ RCSθsϕs[2,:,2] (ϕ 分量，φ_obs=π/2)
"""

using EMSuite
using LinearAlgebra, Statistics

# ─────────────────────────────────────────────
# 参数
# ─────────────────────────────────────────────
radius   = 0.5        # m
freq     = 300e6      # Hz (λ₀ = 1 m, ka ≈ π)
eps_r    = 4.0        # 相对介电常数 (无损)
mu_r     = 1.0
mesh_res = 8          # 网格精度 (球面细分数; 越大越精确)

θ_range  = range(0, π, 181)  # 0 → 180 度
φ_E      = [0.0]              # E-plane
φ_H      = [π/2]              # H-plane

println("="^60)
println("Phase 22.5 PMCHW + Mie 验证")
println("  球半径 = $radius m, freq = $(freq/1e6) MHz")
println("  eps_r = $eps_r, mu_r = $mu_r")
println("="^60)

# ─────────────────────────────────────────────
# 1. Mie 解析解
# ─────────────────────────────────────────────
println("\n[1] 计算 Mie 级数参考值 ...")
mie_E, mie_H, mie_U = calculate_mie_rcs_dielectric_sphere(
    radius, freq, collect(θ_range), eps_r, mu_r,
)
println("  Mie OK: max E-plane RCS = $(maximum(mie_E)*1e4:.2f) cm²")

# ─────────────────────────────────────────────
# 2. 生成网格
# ─────────────────────────────────────────────
println("\n[2] 建立球面 RWG 网格 (mesh_res=$mesh_res) ...")
mesh = sphere_mesh(radius, mesh_res)
basis = RWGBasis(mesh)
println("  RWG 基函数数: $(num_basis(basis))")

# ─────────────────────────────────────────────
# 3. 组装 PMCHW 矩阵并求解
# ─────────────────────────────────────────────
println("\n[3] 组装 PMCHW 矩阵 ...")
pmchw = PMCHWSetup(freq, eps_r, mu_r)
Z = assemble_pmchw(pmchw, basis)
println("  矩阵大小: $(size(Z))")

println("\n[4] 组装激励向量 ...")
pw = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
b  = excitation_vector(pw, pmchw, basis)

println("\n[5] 直接求解 ...")
I = Z \ b
println("  求解 OK, |I|_max = $(maximum(abs.(I)):.4e)")

# ─────────────────────────────────────────────
# 4. 计算双站 RCS
# ─────────────────────────────────────────────
println("\n[6] 计算双站 RCS ...")
θv = collect(θ_range)
RCS_E, _, _ = radarCrossSection(θv, φ_E, I, basis, pmchw.k0, pmchw.eta0)
RCS_H, _, _ = radarCrossSection(θv, φ_H, I, basis, pmchw.k0, pmchw.eta0)

# RCS_E[1,:,1] = θ 分量 at φ=0  →  Mie E-plane (S₂)
# RCS_H[2,:,1] = ϕ 分量 at φ=π/2 → Mie H-plane (S₁)
num_rcs_E = RCS_E[1, :, 1]
num_rcs_H = RCS_H[2, :, 1]

# ─────────────────────────────────────────────
# 5. 误差统计 (dB)
# ─────────────────────────────────────────────
ε  = 1e-20  # 避免 log(0)
dB_mie_E   = 10 .* log10.(mie_E .+ ε)
dB_num_E   = 10 .* log10.(num_rcs_E .+ ε)
dB_mie_H   = 10 .* log10.(mie_H .+ ε)
dB_num_H   = 10 .* log10.(num_rcs_H .+ ε)

diff_E = abs.(dB_num_E .- dB_mie_E)
diff_H = abs.(dB_num_H .- dB_mie_H)

println("\n[RESULT] dB 误差 (PMCHW vs Mie)")
println("  E-plane: mean = $(mean(diff_E):.2f) dB, max = $(maximum(diff_E):.2f) dB")
println("  H-plane: mean = $(mean(diff_H):.2f) dB, max = $(maximum(diff_H):.2f) dB")

tol_dB = 2.0
pass_E = maximum(diff_E) < tol_dB
pass_H = maximum(diff_H) < tol_dB
println("\n  E-plane pass (<$(tol_dB) dB): $pass_E")
println("  H-plane pass (<$(tol_dB) dB): $pass_H")
println(pass_E && pass_H ? "\n✓ 验证通过" : "\n✗ 验证未通过 — 检查网格精度或公式符号")

# ─────────────────────────────────────────────
# 6. 输出关键角度数值 (便于肉眼检查)
# ─────────────────────────────────────────────
println("\n[角度采样对比] E-plane RCS (m²)")
@printf "%8s  %12s  %12s\n" "θ(deg)" "Mie" "PMCHW"
for deg in [0, 30, 60, 90, 120, 150, 180]
    i    = round(Int, deg/180 * (length(θ_range)-1)) + 1
    @printf "%8.1f  %12.4e  %12.4e\n" deg mie_E[i] num_rcs_E[i]
end
