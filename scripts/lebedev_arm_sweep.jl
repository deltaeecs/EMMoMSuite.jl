# 生成器臂长缩放扫描：找出 k=8 局部 LS（论文方法）能达到 ~1e-3 的源带宽范围，
# 用于确定 dataset_generator 修复后的默认几何缩放。
#
# 用法: julia --project=. scripts/lebedev_arm_sweep.jl

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncation_kernel
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG

paper_metric(W, Tte, Pte) =
    vec(maximum(abs.(W * Tte .- Pte); dims = 1) ./ maximum(abs.(Pte); dims = 1))

function make_dataset(pk, pt, rel_l, nρ, npos; seed = 1, arm_scale, off_scale = 0.25)
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
        DG.generate_dataset_on_poles(rHatsC, view(tArray, :, :, iρ, ir);
            rvec = rbmrps[:, ir], arm_max = arm_scale * rel_l, off_max = off_scale * arm_scale * rel_l)
        DG.generate_dataset_on_poles(rHatsF, view(pArray, :, :, iρ, ir);
            rvec = rbmrps[:, ir], arm_max = arm_scale * rel_l, off_max = off_scale * arm_scale * rel_l)
    end
    return reshape(tArray, nt * 2, :), reshape(pArray, np_ * 2, :), tnodes, pnodes
end

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
    for (pk, pt) in ((13, 27), (27, 53))
        Lb = (pk - 1) ÷ 2
        rel_l = EMMoMSuite.Utilities.find_zero_bisection(x -> truncation_kernel(x) - (pk + 1) / 2, 0)
        println("\n##### pk=$pk -> pt=$pt (Lb=$Lb, rel_l=$(round(rel_l, digits = 4))) #####")
        for asc in (0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
            T, P, tnodes, pnodes = make_dataset(pk, pt, rel_l, 40, 120; seed = 20260808, arm_scale = asc)
            flag = trunc(Int, 0.8 * size(T, 2))
            Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]
            W = fit_fixed(tnodes, pnodes, T, P, 8; flag = flag)
            εi = paper_metric(W, Tte, Pte)
            @printf("  arm=%.2f*rel_l (≈%.4fλ): k=8 εi均值=%.3e εi最大=%.3e\n",
                asc, asc * rel_l, mean(εi), maximum(εi))
        end
    end
end

main()
