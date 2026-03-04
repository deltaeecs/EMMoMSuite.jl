"""
verify_pmchw_sphere.jl — Phase 22.5 验证脚本

用 PMCHW 公式求解均匀介质球 RCS，与 Mie 级数对比。

约定 (与 PlaneWave/RadiationIntegral 一致):
  入射波: PlaneWave(freq, 0.0, 0.0, [1,0,0]) → +z 方向传播，x 极化
  观测角: theta ∈ [0,π], phi=0 (E-plane), phi=π/2 (H-plane)
  Mie E-plane ↔ RCSθsϕs[1,:,1] (θ 分量，phi_obs=0)
  Mie H-plane ↔ RCSθsϕs[2,:,1] (ϕ 分量，phi_obs=π/2)
"""

using EMSuite
using LinearAlgebra, Statistics, Printf

# ─────────────────────────────────────────────
# 参数
# ─────────────────────────────────────────────
radius  = 0.5      # m
freq    = 300e6    # Hz  (lambda_0 = 1 m, ka ≈ π)
eps_r   = 4.0      # 相对介电常数 (无损)
mu_r    = 1.0
n_theta = 16       # 球面网格纬向数（nt=16 使 H-plane max < 2 dB）
n_phi   = 32       #             经向数

theta_range = range(0, pi, 181)  # 0 → 180 度
phi_E = [0.0]                    # E-plane
phi_H = [pi/2]                   # H-plane

println("="^60)
println("Phase 22.5 PMCHW + Mie 验证")
@printf "  球半径 = %.3f m, freq = %.1f MHz\n" radius freq/1e6
@printf "  eps_r = %.1f, mu_r = %.1f\n" eps_r mu_r
println("="^60)

# ─────────────────────────────────────────────
# 1. Mie 解析解
# ─────────────────────────────────────────────
println("\n[1] 计算 Mie 级数参考值 ...")
mie_E, mie_H, _ = calculate_mie_rcs_dielectric_sphere(
    radius, freq, collect(theta_range), eps_r, mu_r,
)
@printf "  Mie OK: max E-plane RCS = %.2f cm²\n" maximum(mie_E)*1e4

# ─────────────────────────────────────────────
# 2. 生成网格
# ─────────────────────────────────────────────
@printf "\n[2] 建立球面 RWG 网格 (n_theta=%d, n_phi=%d) ...\n" n_theta n_phi
mesh  = generate_sphere_mesh(radius, n_theta, n_phi)
basis = RWGBasis(mesh)
@printf "  RWG 基函数数: %d\n" num_basis(basis)

# ─────────────────────────────────────────────
# 3. 组装 PMCHW 矩阵
# ─────────────────────────────────────────────
println("\n[3] 组装 PMCHW 矩阵 ...")
pmchw = PMCHW(freq, eps_r, mu_r)
Z = assemble_impedance_matrix(pmchw, basis)
@printf "  矩阵大小: %d × %d\n" size(Z,1) size(Z,2)

# ─────────────────────────────────────────────
# 4. 组装激励向量并求解
# ─────────────────────────────────────────────
println("\n[4] 组装激励向量 ...")
pw = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
b  = excitation_vector(pmchw, pw, basis)

println("\n[5] 直接求解 ...")
I = Z \ b
@printf "  求解 OK, |I|_max = %.4e\n" maximum(abs.(I))

# ─────────────────────────────────────────────
# 5. 计算双站 RCS
# ─────────────────────────────────────────────
println("\n[6] 计算双站 RCS ...")
set_frequency!(freq)   # 确保 radiation_integral_rwg 中的 get_k0() 正确
theta_v = collect(theta_range)
RCSe, _, _ = radarCrossSection(theta_v, phi_E, I, basis, pmchw.k0, pmchw.eta0)
RCSh, _, _ = radarCrossSection(theta_v, phi_H, I, basis, pmchw.k0, pmchw.eta0)

# RCSe[1,:,1] = θ 分量 at phi=0   → Mie E-plane (S₂)
# RCSh[2,:,1] = ϕ 分量 at phi=π/2 → Mie H-plane (S₁)
num_rcs_E = RCSe[1, :, 1]
num_rcs_H = RCSh[2, :, 1]

# ─────────────────────────────────────────────
# 6. 误差统计 (dB)
# ─────────────────────────────────────────────
eps_guard = 1e-20
dB_mie_E = 10 .* log10.(mie_E .+ eps_guard)
dB_num_E = 10 .* log10.(num_rcs_E .+ eps_guard)
dB_mie_H = 10 .* log10.(mie_H .+ eps_guard)
dB_num_H = 10 .* log10.(num_rcs_H .+ eps_guard)

diff_E = abs.(dB_num_E .- dB_mie_E)
diff_H = abs.(dB_num_H .- dB_mie_H)

println("\n[RESULT] dB 误差 (PMCHW vs Mie)")
@printf "  E-plane: mean = %.2f dB,  max = %.2f dB\n" mean(diff_E) maximum(diff_E)
@printf "  H-plane: mean = %.2f dB,  max = %.2f dB\n" mean(diff_H) maximum(diff_H)

tol_dB = 2.0
pass_E = maximum(diff_E) < tol_dB
pass_H = maximum(diff_H) < tol_dB
@printf "\n  E-plane pass (< %.1f dB): %s\n" tol_dB string(pass_E)
@printf "  H-plane pass (< %.1f dB): %s\n" tol_dB string(pass_H)
println(pass_E && pass_H ? "\n✓ Phase 22.5 验证通过" : "\n✗ 验证未通过 — 检查网格精度或公式符号")

# ─────────────────────────────────────────────
# 7. 角度采样对比表
# ─────────────────────────────────────────────
println("\n[角度采样] E-plane RCS (m²):")
@printf "%8s  %12s  %12s  %8s\n" "θ(deg)" "Mie" "PMCHW" "|dB err|"
for deg in [0, 30, 60, 90, 120, 150, 180]
    i = round(Int, deg / 180 * (length(theta_range) - 1)) + 1
    @printf "%8.1f  %12.4e  %12.4e  %8.2f\n" deg mie_E[i] num_rcs_E[i] diff_E[i]
end
