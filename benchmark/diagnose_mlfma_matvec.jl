"""
    diagnose_mlfma_matvec.jl

诊断 MLFMA MatVec vs Direct MatVec 差异。
将 Direct Z*x 拆分为近场 + 远场，与 MLFMA 逐项对比。

用法: julia --project=. benchmark/diagnose_mlfma_matvec.jl
"""

using EMSuite
using LinearAlgebra
using Printf
using Statistics
using Random

println("=" ^ 72)
println("  MLFMA MatVec 诊断")
println("=" ^ 72)

mesh_file = joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles", "jet_100MHz.nas")
mesh = read_nas_mesh(mesh_file, scale=1.0)
freq = 1e8
set_frequency!(freq)

basis = RWGBasis(mesh)
N = num_basis(basis)
println("Unknowns: $N")

efie = EFIE(freq)

# ============================================================
#  1. Direct Z matrix (gold standard)
# ============================================================
println("\n--- 组装 Direct Z ---")
t1 = @elapsed Z = assemble_impedance_matrix(efie, basis)
@printf("  %d×%d, %.1f s\n", size(Z)..., t1)

# ============================================================
#  2. MLFMA Operator
# ============================================================
println("\n--- 构建 MLFMA ---")
λ = 299792458.0 / freq
basis_mlfma = RWGBasis(mesh)
t2 = @elapsed mlfma_op = MLFMAOperator(efie, basis_mlfma, 0.35λ)
@printf("  Setup: %.1f s\n", t2)

# ============================================================
#  3. MatVec 对比 (使用 x = 电流解作为测试向量)
# ============================================================
println("\n--- 求解 Direct ---")
source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
V = excitation_vector(efie, source, basis)
I_direct = Z \ V

# 检查 sorted_ids — MLFMA 重排了基函数顺序
sorted_ids = mlfma_op.sorted_ids
inv_sorted = zeros(Int, N)
for (new, old) in enumerate(sorted_ids)
    inv_sorted[old] = new
end

# 将 I_direct 转换到 MLFMA 排序
I_sorted = I_direct[sorted_ids]

println("\n--- MatVec 对比 (x = I_direct) ---")

# Direct: y_direct = Z * I_direct
y_direct = Z * I_direct

# MLFMA: 需要排列 input 为 MLFMA 顺序，output 再排回来
y_mlfma_sorted = zeros(ComplexF64, N)
mul!(y_mlfma_sorted, mlfma_op, I_sorted)

# 将 MLFMA 结果转换回原始顺序
y_mlfma = y_mlfma_sorted[inv_sorted]

diff_matvec = norm(y_direct - y_mlfma) / norm(y_direct)
@printf("  ||Z*x - MLFMA*x|| / ||Z*x|| = %.4f (%.2f%%)\n", diff_matvec, diff_matvec*100)

# ============================================================
#  4. 拆分贡献: 近场 vs 远场
# ============================================================
println("\n--- 近场/远场拆分 ---")

# MLFMA 近场: Z_near * x_sorted
y_near_sorted = mlfma_op.Z_near * I_sorted
y_near = y_near_sorted[inv_sorted]

# MLFMA 远场: y_mlfma - y_near (远场 = 总 - 近场)
y_far = y_mlfma - y_near

println("  ||y_near|| = $(norm(y_near))")
println("  ||y_far||  = $(norm(y_far))")
println("  ||y_direct|| = $(norm(y_direct))")
println("  ||y_mlfma||  = $(norm(y_mlfma))")

# 理论上 Z = Z_near + Z_far (远场隐含在 MLFMA 中)
# 所以 y_near 应该覆盖所有近邻交互
# y_far 应该覆盖所有远场交互

# ============================================================
#  5. 用 Direct Z 构建参考远场
# ============================================================
println("\n--- 构建 Direct 参考远场 ---")

# 从 Z_near 构建近场 mask
# Z_near 是 sorted 顺序的稀疏矩阵
# 需要把 Direct Z 也排成 sorted 顺序来比较

Z_sorted = Z[sorted_ids, sorted_ids]
Z_near_dense = Matrix(mlfma_op.Z_near)

# 近场 entries mask (Z_near != 0)
near_mask = abs.(Z_near_dense) .> 0

# Direct 中近场部分
Z_direct_near = Z_sorted .* near_mask

# Direct 中远场部分
Z_direct_far = Z_sorted .* .!near_mask

# 对比
y_direct_near = Z_direct_near * I_sorted
y_direct_far = Z_direct_far * I_sorted

println("  ||Z_direct_near * x|| = $(norm(y_direct_near))")
println("  ||Z_direct_far * x||  = $(norm(y_direct_far))")

# MLFMA 近场 vs Direct 近场
diff_near = norm(y_direct_near - y_near_sorted) / norm(y_direct_near)
@printf("  近场差异: ||Z_direct_near*x - Z_mlfma_near*x|| / ||Z_direct_near*x|| = %.6f (%.4f%%)\n",
        diff_near, diff_near*100)

# 近场矩阵 entry-by-entry 对比
Z_near_diff = Z_direct_near - Z_near_dense
near_entries = Z_near_dense[near_mask]
direct_near_entries = Z_sorted[near_mask]
@printf("  近场矩阵 Frobenius 比: ||Z_direct[near] - Z_mlfma_near|| / ||Z_direct[near]|| = %.6f\n",
        norm(Z_near_diff) / norm(Z_direct_near))

# 远场对比
diff_far = norm(y_direct_far - (y_mlfma_sorted - y_near_sorted)) / norm(y_direct_far)
@printf("  远场差异: ||Z_direct_far*x - MLFMA_far*x|| / ||Z_direct_far*x|| = %.6f (%.4f%%)\n",
        diff_far, diff_far*100)

# ============================================================
#  6. 元素级别的比率分析
# ============================================================
println("\n--- 远场比率分析 ---")
y_mlfma_far = y_mlfma_sorted - y_near_sorted
# 逐元素比率 (仅在 y_direct_far 不太小的位置)
threshold = 1e-10 * maximum(abs.(y_direct_far))
big_enough = abs.(y_direct_far) .> threshold
ratios = y_mlfma_far[big_enough] ./ y_direct_far[big_enough]
real_ratios = real.(ratios)
imag_ratios = imag.(ratios)
abs_ratios = abs.(ratios)

@printf("  元素比率 |MLFMA_far / Direct_far|:\n")
@printf("    Mean: %.4f\n", mean(abs_ratios))
@printf("    Std:  %.4f\n", std(abs_ratios))
@printf("    Min:  %.4f\n", minimum(abs_ratios))
@printf("    Max:  %.4f\n", maximum(abs_ratios))
@printf("    Mean Real: %.4f\n", mean(real_ratios))
@printf("    Mean Imag: %.4f\n", mean(imag_ratios))

# 检查是否有常数比率 (暗示缺少/多余的常数因子)
phase_ratios = angle.(ratios)
@printf("  相位偏差:\n")
@printf("    Mean phase: %.4f rad (%.1f deg)\n", mean(phase_ratios), rad2deg(mean(phase_ratios)))
@printf("    Std phase:  %.4f rad\n", std(phase_ratios))

# ============================================================
#  7. EFIE factor 信息
# ============================================================
println("\n--- EFIE factor ---")
println("  efie.factor = $(efie.factor)")
println("  |efie.factor| = $(abs(efie.factor))")
println("  phase(efie.factor) = $(angle(efie.factor)) rad = $(rad2deg(angle(efie.factor))) deg")
println("  jkη/(16π) = $(im * efie.k * efie.eta / (16π))")
println("  jkη/(4π)  = $(im * efie.k * efie.eta / (4π))")
println("  hasfield factor: $(hasfield(typeof(efie), :factor))")

# ============================================================
#  8. 对比 Z_near 对角线
# ============================================================
println("\n--- 对角线对比 ---")
diag_direct = diag(Z_sorted)
diag_near = diag(Z_near_dense)
diag_ratio = diag_near ./ diag_direct
@printf("  对角线比率: Mean=%.6f, Std=%.6f\n", mean(abs.(diag_ratio)), std(abs.(diag_ratio)))
@printf("  前5个: %s\n", join([@sprintf("%.4f", abs(r)) for r in diag_ratio[1:min(5,N)]], ", "))
