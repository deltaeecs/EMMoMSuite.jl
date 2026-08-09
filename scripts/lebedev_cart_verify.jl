# 笛卡尔标量 SH 矢量插值（新算法）验证：
#   1. EFIE 矢量类 εi：笛卡尔(稠密/局部稀疏) vs VSH vs 训练式 k=8
#   2. 自旋/标量限带随机矢量场机器精度
#   3. 高阶 (65->131) 局部稀疏版构造与精度
#   4. 稀疏度-精度权衡（Lloc / cap 扫描）
#
# 用法: julia --project=. scripts/lebedev_cart_verify.jl

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncation_kernel
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG
import EMMoMSuite.FastAlgorithms.Lebedev.SHInterp as SH

paper_metric(W, Tte, Pte) =
    vec(maximum(abs.(W * Tte .- Pte); dims = 1) ./ maximum(abs.(Pte); dims = 1))

rel_l_for(pk) = begin
    f(x) = truncation_kernel(x) - (pk + 1) / 2
    lo, hi = 1e-4, 20.0
    for _ in 1:80
        mid = (lo + hi) / 2
        f(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    (lo + hi) / 2
end

"""修正几何数据集（源整体在盒内，与修复后生成器一致）"""
function make_corrected(pk, pt; nρ = 40, npos = 120, seed = 20260808)
    rel_l = rel_l_for(pk)
    arm = min(0.12, 0.5 * rel_l)
    rscale = max(rel_l / 2 - arm, 1e-4)
    tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
    pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    Random.seed!(seed)
    rbmrps = zeros(Float64, 3, npos)
    for i in axes(rbmrps, 2)
        rbmrps[:, i] .= DG.random_rvec() .* rscale
    end
    rC = LSP.nodes2Poles(tnodes)
    rF = LSP.nodes2Poles(pnodes)
    T = zeros(ComplexF64, nt, 2, nρ, npos)
    P = zeros(ComplexF64, nf, 2, nρ, npos)
    for ir in axes(rbmrps, 2), iρ in 1:nρ
        geom = DG.random_source_geometry(rbmrps[:, ir]; arm_max = arm, off_max = 0.125 * arm)
        DG.evaluate_poles!(rC, view(T, :, :, iρ, ir), geom)
        DG.evaluate_poles!(rF, view(P, :, :, iρ, ir), geom)
    end
    return reshape(T, nt * 2, :), reshape(P, nf * 2, :), tnodes, pnodes
end

"""训练式 k 点（修复后：共享几何 + hcat 复约束）"""
function trained_weights(pk, pt, k; nρ = 40, npos = 120, seed = 20260808)
    tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
    pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    T, P, _, _ = make_corrected(pk, pt; nρ = nρ, npos = npos, seed = seed)
    flag = trunc(Int, 0.8 * size(T, 2))
    Aall = hcat(real(T[:, 1:flag]), imag(T[:, 1:flag]))
    Ball = hcat(real(P[:, 1:flag]), imag(P[:, 1:flag]))
    idxs = [partialsortperm([norm(tnodes[:, j] .- pnodes[:, i]) for j in 1:nt], 1:k) for i in 1:nf]
    W = spzeros(nf * 2, nt * 2)
    for i in 1:nf
        S = idxs[i]
        cols = vcat(S, S .+ nt)
        A = Aall[cols, :]
        for comp in (0, 1)
            row = i + comp * nf
            W[row, cols] = reshape(Ball[row, :], 1, :) * pinv(A; rtol = 1e-12)
        end
    end
    return W
end

function main()
    for (pk, pt) in ((13, 27), (27, 53))
        Lb = (pk - 1) ÷ 2
        T, P, tnodes, pnodes = make_corrected(pk, pt)
        flag = trunc(Int, 0.8 * size(T, 2))
        Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]
        nt, nf = size(tnodes, 2), size(pnodes, 2)
        println("\n############ pk=$pk -> pt=$pt （Lb=$Lb）############")

        Wt = trained_weights(pk, pt, 8)
        @printf("训练式 k=8:        εi均值=%.3e\n", mean(paper_metric(Wt, Tte, Pte)))

        Wv = SH.interp_weights_vsh(pk, pt)
        @printf("VSH 精确(Lb+1):    εi均值=%.3e（极点 θ/ϕ 基奇异残留）\n", mean(paper_metric(Wv, Tte, Pte)))

        for deg in (Lb, Lb + 1)
            deg > isqrt(nt) - 1 && continue
            Wc = SH.interp_weights_cart(pk, pt; degree = deg)
            @printf("笛卡尔精确 degree=%d: εi均值=%.3e 最大=%.3e\n",
                deg, mean(paper_metric(Wc, Tte, Pte)), maximum(paper_metric(Wc, Tte, Pte)))
        end

        println("笛卡尔局部稀疏（Lloc, cap 扫描）:")
        for (Lloc, cap) in ((3, 0.6), (4, 0.8), (6, 1.0), (min(Lb, 8), 1.2), (Lb, 1.4))
            Lloc <= 0 && continue
            m = (Lloc + 1)^2
            m >= nt && continue
            Wl = SH.interp_weights_cart_local(pk, pt; Lloc = Lloc, cap_rad = cap)
            εi = paper_metric(Wl, Tte, Pte)
            @printf("  L_loc=%2d cap=%.1f: εi均值=%.3e 最大=%.3e nnz/行=%.1f\n",
                Lloc, cap, mean(εi), maximum(εi), nnz(Wl) / (2nf))
        end
    end

    println("\n############ 高阶 (65->131) 笛卡尔局部稀疏 ############")
    pk, pt = 65, 131
    Lb = (pk - 1) ÷ 2
    T, P, tnodes, pnodes = make_corrected(pk, pt; nρ = 20, npos = 60)
    flag = trunc(Int, 0.8 * size(T, 2))
    Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]
    for (Lloc, cap) in ((3, 0.35), (6, 0.5), (10, 0.7))
        t0 = time()
        Wl = SH.interp_weights_cart_local_orbit(pk, pt; Lloc = Lloc, cap_rad = cap)
        εi = paper_metric(Wl, Tte, Pte)
        @printf("  L_loc=%2d cap=%.2f(轨道压缩): εi均值=%.3e nnz/行=%.1f 构造 %.1fs\n",
            Lloc, cap, mean(εi), nnz(Wl) / (2 * 5810), time() - t0)
    end
end

main()
