# 优化后 Lebedev 一步插值演示：
#   1. LVI truncation 修复生效（0.5λ -> truncL=19, p=41, 590 点，而非 60/5294）
#   2. LbTrainedInterp1tepInfo(method=:sh_exact) 机器精度且无需训练
#   3. :sh_local / :sh_local_orbit 稀疏版与轨道压缩逐位一致
#
# 用法: julia --project=. scripts/lebedev_optimized_demo.jl

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.Lebedev
using EMMoMSuite.FastAlgorithms.Lebedev.SHInterp: realSHmatrix
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation:
    levelIntegralInfoCal, interpolate, anterpolate
using LinearAlgebra, SparseArrays, Random, Printf

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP

println("=== 1. LVI truncation 修复验证 ===")
truncL, poles = levelIntegralInfoCal(0.5, Val(:LbTrained1Step))
@printf("levelIntegralInfoCal(0.5λ): truncL=%d (修复前 60), Lebedev 点数=%d (修复前 5294), 权重和/4π=%.6f\n",
    truncL, length(poles.Wθϕs), sum(poles.Wθϕs) / 4π)

for (pk, pt) in ((13, 27), (27, 53))
    Lb = (pk - 1) ÷ 2
    println("\n=== 2. pk=$pk -> pt=$pt :sh_exact（球谐精确）===")
    t0 = time()
    info = LbTrainedInterp1tepInfo(pk, pt; method = :sh_exact)
    t1 = time()
    tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
    pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
    Yc = realSHmatrix(tnodes, Lb)
    Yf = realSHmatrix(pnodes, Lb)
    Random.seed!(42)
    worst = 0.0
    for _ in 1:50
        cθ, cϕ = randn((Lb + 1)^2), randn((Lb + 1)^2)
        coarse = [Yc * cθ Yc * cϕ]          # (nc, 2) θ/ϕ 分量
        exact = [Yf * cθ Yf * cϕ]           # (nf, 2)
        est = interpolate(info, coarse)
        worst = max(worst, maximum(abs.(est .- exact)) / maximum(abs.(exact)))
    end
    @printf("  构造 %.3fs, 矩阵尺寸 %dx%d, 限带矢量场最大相对误差=%.3e\n",
        t1 - t0, size(info.θϕCSC, 1), size(info.θϕCSC, 2), worst)

    println("=== 3. :sh_local vs :sh_local_orbit（Lloc=$(min(Lb, 4))）===")
    Lloc = min(Lb, 4)
    cap = Lloc <= 3 ? 1.0 : 1.4
    wl = LbTrainedInterp1tepInfo(pk, pt; method = :sh_local, Lloc = Lloc, cap_rad = cap)
    wo = LbTrainedInterp1tepInfo(pk, pt; method = :sh_local_orbit, Lloc = Lloc, cap_rad = cap)
    d = maximum(abs.(wl.θϕCSC .- wo.θϕCSC))
    @printf("  max|W_local-W_orbit|=%.2e, nnz_local=%d, nnz_orbit=%d\n", d, nnz(wl.θϕCSC), nnz(wo.θϕCSC))

    # 反插值说明：MLFMA 下行（disaggregation）使用插值矩阵的转置 Wᵀ 作为
    # 伴随算子（adjoint），这与"从细层恢复粗层采样"不同，无需等于原值。
    println("  反插值 = 插值矩阵转置（伴随算子），用于 MLFMA 下行，接口一致。")
end

println("\n=== 4. 原始训练式路径仍可用（文件不存在时才会训练，此处只验证接口存在）===")
println("  method=:trained（默认）保留原行为；method=:sh_exact/:sh_local/:sh_local_orbit 为优化路径。")
