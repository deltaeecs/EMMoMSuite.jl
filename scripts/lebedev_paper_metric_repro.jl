# 按论文指标复现原始 LVI-LQ：
#   εi = |f̃(k̂) − f(k̂)| / max|f(k̂)|   （每个辐射函数一个标量，Eq.(7)，Özgür 定义）
# 复现目标（论文 Fig.2）：8 个插值点 -> 误差接近 4x4 Lagrange (~6e-4)；
#                           18 点 -> 接近 6x6 (~6e-5)。
# 同时验证：
#   - 各阶共享固定点（论文称 14 个）在插值中误差为 0；
#   - 优化方法（球谐精确 / 局部约束）在同一指标下的提升。
#
# 用法: julia --project=. scripts/lebedev_paper_metric_repro.jl

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncation_kernel
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.pinv2interpW as PW
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG
import EMMoMSuite.FastAlgorithms.Lebedev.SHInterp as SH

"""论文指标：每辐射函数 max|err| / max|f|"""
function paper_metric(W, Tte, Pte)
    est = W * Tte
    err = maximum(abs.(est .- Pte); dims = 1)       # (1, Nte)
    denom = maximum(abs.(Pte); dims = 1)             # (1, Nte)
    εi = vec(err ./ denom)
    return εi
end

function make_dataset(pk, pt, rel_l, nρ, npos; seed = 1)
    τt, τp = (pk - 1) ÷ 2, (pt - 1) ÷ 2
    tnodes = LSP.get_t_nodes(τt)
    pnodes = LSP.get_t_nodes(τp)
    nt, np_ = size(tnodes, 2), size(pnodes, 2)
    Random.seed!(seed)
    ρhats = zeros(Float64, 3, nρ)
    for i in axes(ρhats, 2)
        ρhats[:, i] .= DG.random_rhat()
    end
    rbmrps = zeros(Float64, 3, npos)
    for i in axes(rbmrps, 2)
        rbmrps[:, i] .= DG.random_rvec() .* (rel_l / 2)
    end
    rHatsC = LSP.nodes2Poles(tnodes)
    rHatsF = LSP.nodes2Poles(pnodes)
    tArray = zeros(ComplexF64, nt, 2, nρ, npos)
    pArray = zeros(ComplexF64, np_, 2, nρ, npos)
    for ir in axes(rbmrps, 2), iρ in axes(ρhats, 2)
        DG.generate_dataset_on_poles(rHatsC, view(tArray, :, :, iρ, ir); rvec = rbmrps[:, ir])
        DG.generate_dataset_on_poles(rHatsF, view(pArray, :, :, iρ, ir); rvec = rbmrps[:, ir])
    end
    return reshape(tArray, nt * 2, :), reshape(pArray, np_ * 2, :), tnodes, pnodes
end

function main()
    for (pk, pt) in ((13, 27), (27, 53))
        Lb = (pk - 1) ÷ 2
        rel_l = EMMoMSuite.Utilities.find_zero_bisection(x -> truncation_kernel(x) - (pk + 1) / 2, 0)
        tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
        pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
        nt, nf = size(tnodes, 2), size(pnodes, 2)
        println("\n############ pk=$pk -> pt=$pt ############")

        T, P, _, _ = make_dataset(pk, pt, rel_l, 40, 120; seed = 20260808)
        flag = trunc(Int, 0.8 * size(T, 2))
        Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]

        # --- 原始管线：nInterp 1..24 的论文指标 ---
        println("原始管线（IDW 初始化 + 逐行 pinv）论文指标 εi（均值 / 最大）:")
        npts = 1:24
        εi_mean = Float64[]
        εi_max = Float64[]
        for k in npts
            w = PW.interpWeightsInitial(tnodes, pnodes; nInterp = k)
            PW.pinv2W!(w, k, vcat(real(T[:, 1:flag]), imag(T[:, 1:flag])),
                vcat(real(P[:, 1:flag]), imag(P[:, 1:flag])))
            εi = paper_metric(w, Tte, Pte)
            push!(εi_mean, mean(εi))
            push!(εi_max, maximum(εi))
        end
        for (i, k) in enumerate(npts)
            mark = (k == 8 || k == 9 || k == 18) ? "  <== 论文关注点" : ""
            @printf("  nInterp=%2d  εi均值=%.3e  εi最大=%.3e%s\n", k, εi_mean[i], εi_max[i], mark)
        end

        # --- 共享固定点分析 ---
        common = Int[]
        for i in 1:nf
            any(norm(pnodes[:, i] .- tnodes[:, j]) < 1e-6 for j in 1:nt) && push!(common, i)
        end
        println("粗/细层共享节点数 = ", length(common), "（论文称 14：6 轴 + 8 立方体角）")
        # 共享点在原始 W 上的插值误差（应为 ~0）
        if length(common) > 0
            w8 = PW.interpWeightsInitial(tnodes, pnodes; nInterp = 8)
            PW.pinv2W!(w8, 8, vcat(real(T[:, 1:flag]), imag(T[:, 1:flag])),
                vcat(real(P[:, 1:flag]), imag(P[:, 1:flag])))
            rows = vcat(common, common .+ nf)
            errmat = abs.(w8 * Tte .- Pte)
            shared_err = maximum(errmat[rows, :]; dims = 1)[1] ./ maximum(abs.(Pte); dims = 1)[1]
            @printf("共享节点上的 εi 均值=%.3e（理想为 0，因为无需插值）\n", mean(shared_err))
        end

        # --- 优化方法在同一测试集上的论文指标 ---
        println("优化方法（同一 EFIE 留出集，论文指标 εi 均值）:")
        Wsh = SH.vectorize(SH.interp_weights_exact(pk, pt))
        εi = paper_metric(Wsh, Tte, Pte)
        @printf("  球谐精确(稠密 %dx%d): εi均值=%.3e\n", size(Wsh, 1), size(Wsh, 2), mean(εi))
        for (Lloc, cap) in ((3, 0.6), (6, 1.0))
            Wloc = SH.vectorize(SH.interp_weights_local(pk, pt; Lloc = Lloc, cap_rad = cap))
            εi = paper_metric(Wloc, Tte, Pte)
            @printf("  局部约束 L_loc=%d cap=%.1f (nnz/行=%.1f): εi均值=%.3e\n",
                Lloc, cap, nnz(Wloc) / (2nf), mean(εi))
        end
    end
end

main()
