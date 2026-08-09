# 权重矩阵修复验证：
#   1. 修复后的数据生成器（源几何按层盒子缩放）+ hcat 复约束逐行 LS
#      -> 论文 Fig.2 指标 εi 在 k=8 时应为 ~1e-4..1e-3（论文 ~6e-4）
#   2. 球谐精确 W（默认 degree=Lb+1）：限带函数机器精度、确定性、行和=1
#   3. 局部约束/轨道压缩 W：L_loc 精确性、行和=1、orbit==naive
#   4. 包内全链路 LbTrainedInterp1tepInfo(method=:sh_auto)
#
# 用法: julia --project=. scripts/lebedev_weights_fix_verify.jl

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.Lebedev
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation:
    truncation_kernel, interpolate, levelIntegralInfoCal
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.pinv2interpW as PW
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG
import EMMoMSuite.FastAlgorithms.Lebedev.SHInterp as SH

paper_metric(W, Tte, Pte) =
    vec(maximum(abs.(W * Tte .- Pte); dims = 1) ./ maximum(abs.(Pte); dims = 1))

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
        geom = DG.random_source_geometry(rbmrps[:, ir];
            arm_max = 0.5 * rel_l, off_max = 0.125 * rel_l)
        DG.evaluate_poles!(rHatsC, view(tArray, :, :, iρ, ir), geom)
        DG.evaluate_poles!(rHatsF, view(pArray, :, :, iρ, ir), geom)
    end
    return reshape(tArray, nt * 2, :), reshape(pArray, np_ * 2, :), tnodes, pnodes
end

"""修复后的逐行 LS：hcat(real, imag) 全复数约束（与论文 Algorithm 2 等价）"""
function fit_fixed(tnodes, pnodes, T, P, k; flag)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
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
    println("===== 0. 确认修复已生效（LVI truncation / 默认 :sh_auto）=====")
    truncL, poles = levelIntegralInfoCal(0.5, Val(:LbTrained1Step))
    @printf("levelIntegralInfoCal(0.5λ): truncL=%d, 点数=%d\n", truncL, length(poles.Wθϕs))

    for (pk, pt) in ((13, 27), (27, 53))
        Lb = (pk - 1) ÷ 2
        rel_l = EMMoMSuite.Utilities.find_zero_bisection(x -> truncation_kernel(x) - (pk + 1) / 2, 0)
        tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
        pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
        nt, nf = size(tnodes, 2), size(pnodes, 2)
        println("\n############ pk=$pk -> pt=$pt （Lb=$Lb, rel_l=$(round(rel_l, digits = 4))）############")

        # 1. 修复后管线：论文 Fig.2 复现
        T, P, _, _ = make_dataset(pk, pt, rel_l, 40, 120; seed = 20260808)
        flag = trunc(Int, 0.8 * size(T, 2))
        Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]
        println("1. 修复后管线（缩放生成器 + hcat 复约束），论文指标 εi:")
        for k in (4, 8, 9, 12, 18)
            W = fit_fixed(tnodes, pnodes, T, P, k; flag = flag)
            εi = paper_metric(W, Tte, Pte)
            @printf("   k=%2d: εi均值=%.3e εi最大=%.3e  nnz/行=%.1f\n", k, mean(εi), maximum(εi), nnz(W) / (2nf))
        end

        # 2. 球谐精确 W：确定性、限带函数机器精度、行和=1、EFIE 数据 εi
        println("2. 球谐精确 W:")
        t0 = time()
        W1 = SH.vectorize(SH.interp_weights_exact(pk, pt))
        W2 = SH.vectorize(SH.interp_weights_exact(pk, pt))
        @printf("   确定性(W1==W2)=%s, 构造 %.3fs\n", W1 == W2, time() - t0)
        @printf("   行和 max|W·1−1|=%.2e\n", maximum(abs.(W1 * ones(2nt) .- 1)))
        deg = Lb
        Yc = SH.realSHmatrix(tnodes, deg)
        Yf = SH.realSHmatrix(pnodes, deg)
        Random.seed!(7)
        worst = 0.0
        for _ in 1:100
            c = randn((deg + 1)^2)
            coarse = vcat(Yc * c, Yc * c)
            exact = vcat(Yf * c, Yf * c)
            worst = max(worst, maximum(abs.(W1 * coarse .- exact)) / maximum(abs.(exact)))
        end
        @printf("   限带矢量场(degree<=%d)最大相对误差=%.3e\n", deg, worst)
        εi = paper_metric(W1, Tte, Pte)
        @printf("   同一 EFIE 留出集 εi均值=%.3e（修复后数据可表示，应远小于坏权重时代的 O(1)）\n", mean(εi))

        # 3. 局部约束 + 轨道压缩
        Lloc = min(4, Lb)
        cap = 1.4
        Wloc = SH.vectorize(SH.interp_weights_local(pk, pt; Lloc = Lloc, cap_rad = cap))
        Worb = SH.vectorize(SH.interp_weights_local_orbit(pk, pt; Lloc = Lloc, cap_rad = cap))
        @printf("3. 局部约束(L_loc=%d): 行和=%.1e, max|W_local-W_orbit|=%.1e, nnz/行=%.1f\n",
            Lloc, maximum(abs.(Wloc * ones(2nt) .- 1)), maximum(abs.(Wloc .- Worb)), nnz(Wloc) / (2nf))

        # 4. 包内全链路 :sh_auto
        t0 = time()
        info = LbTrainedInterp1tepInfo(pk, pt; method = :sh_auto)
        @printf("4. LbTrainedInterp1tepInfo(:sh_auto): %.3fs, W 尺寸 %dx%d, nnz=%d\n",
            time() - t0, size(info.θϕCSC, 1), size(info.θϕCSC, 2), nnz(info.θϕCSC))
    end
end

main()
