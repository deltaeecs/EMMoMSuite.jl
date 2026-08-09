# 判别实验：为什么当前代码复现不出论文 Fig.2 的 ~6e-4？
# 假设1: vcat(real,imag) 布局导致只用了实部约束（论文伪代码是对复数数据直接 pinv）
# 假设2: θ/ϕ 分量非纯标量限带 -> 矢量耦合需要额外自由度
# 假设3: 球谐基在高阶有 bug
#
# 验证：
#   A. 实球谐 Gram 正交性（Lebedev 权重下）到 Lb=13/26
#   B. 修正版逐行 LS（hcat(real,imag)，全复数约束，与论文 Algorithm 2 等价）
#      在 nInterp=4/8/9/12/18 下的论文指标 εi
#   C. 全支撑矢量 LS（全局伪逆，含 θ↔ϕ 耦合）
#   D. 标量球谐 W 取 degree = Lb 与 Lb+1 的 εi（检验矢量投影的 L+1 分量）
#
# 用法: julia --project=. scripts/lebedev_paper_repro_v2.jl

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncation_kernel
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
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
        DG.generate_dataset_on_poles(rHatsC, view(tArray, :, :, iρ, ir); rvec = rbmrps[:, ir])
        DG.generate_dataset_on_poles(rHatsF, view(pArray, :, :, iρ, ir); rvec = rbmrps[:, ir])
    end
    return reshape(tArray, nt * 2, :), reshape(pArray, np_ * 2, :), tnodes, pnodes
end

"""修正版逐行 LS：支撑 = k 最近粗点（θ+ϕ 列），全复数约束（与论文 Algorithm 2 等价）"""
function fit_complex_local(tnodes, pnodes, T, P, k; flag)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    Aall = hcat(real(T[:, 1:flag]), imag(T[:, 1:flag]))    # (2nt, 2N)
    Ball = hcat(real(P[:, 1:flag]), imag(P[:, 1:flag]))    # (2nf, 2N)
    # 最近 k 个粗点（按弦距离）
    idxs = [partialsortperm([norm(tnodes[:, j] .- pnodes[:, i]) for j in 1:nt], 1:k) for i in 1:nf]
    W = spzeros(nf * 2, nt * 2)
    for i in 1:nf
        S = idxs[i]
        cols = vcat(S, S .+ nt)                              # θ+ϕ 列
        A = Aall[cols, :]
        for comp in (0, 1)
            row = i + comp * nf
            b = Ball[row, :]
            W[row, cols] = reshape(b, 1, :) * pinv(A; rtol = 1e-12)
        end
    end
    return W
end

"""全支撑矢量 LS：全局伪逆，含 θ↔ϕ 耦合，复数全约束"""
function fit_complex_global(T, P; flag)
    A = hcat(real(T[:, 1:flag]), imag(T[:, 1:flag]))       # (2nt, 2N)
    B = hcat(real(P[:, 1:flag]), imag(P[:, 1:flag]))       # (2nf, 2N)
    return Matrix(B * pinv(A; rtol = 1e-12))               # (2nf, 2nt) 实权重
end

function main()
    for (pk, pt) in ((13, 27), (27, 53))
        Lb = (pk - 1) ÷ 2
        rel_l = EMMoMSuite.Utilities.find_zero_bisection(x -> truncation_kernel(x) - (pk + 1) / 2, 0)
        tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
        pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
        nt, nf = size(tnodes, 2), size(pnodes, 2)
        println("\n############ pk=$pk -> pt=$pt （Lb=$Lb）############")

        # A. SH 基 Gram 检查
        for (p, Lg) in ((pk, Lb), (pt, min(Lb + 1, (pt - 1) ÷ 2)))
            nd, wts = LSP.getlbSortedData(p)
            Y = SH.realSHmatrix(nd, Lg)
            G = Y' * (wts .* Y)
            @printf("A. p=%d degree=%d: max|G−I|=%.2e\n", p, Lg, maximum(abs.(G - Matrix{Float64}(I, size(G, 1), size(G, 1)))))
        end

        T, P, _, _ = make_dataset(pk, pt, rel_l, 40, 120; seed = 20260808)
        flag = trunc(Int, 0.8 * size(T, 2))
        Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]

        # B. 修正版逐行 LS（复数全约束）
        println("B. 修正版逐行 LS（复数全约束，论文 Algorithm 2 语义）:")
        for k in (4, 8, 9, 12, 18)
            W = fit_complex_local(tnodes, pnodes, T, P, k; flag = flag)
            εi = paper_metric(W, Tte, Pte)
            @printf("   k=%2d: εi均值=%.3e εi最大=%.3e nnz/行=%.1f\n", k, mean(εi), maximum(εi), nnz(W) / (2nf))
        end

        # C. 全支撑矢量 LS
        Wg = fit_complex_global(T, P; flag = flag)
        εi = paper_metric(Wg, Tte, Pte)
        @printf("C. 全支撑矢量 LS: εi均值=%.3e εi最大=%.3e\n", mean(εi), maximum(εi))

        # D. 标量球谐 W degree = Lb / Lb+1
        for deg in (Lb, Lb + 1)
            Yc = SH.realSHmatrix(tnodes, deg)
            Yf = SH.realSHmatrix(pnodes, deg)
            Ws = Yf * pinv(Yc; rtol = 1e-13)
            Wv = [Ws zeros(nf, nt); zeros(nf, nt) Ws]
            εi = paper_metric(Wv, Tte, Pte)
            @printf("D. 标量球谐 degree=%d: εi均值=%.3e εi最大=%.3e\n", deg, mean(εi), maximum(εi))
        end
    end
end

main()
