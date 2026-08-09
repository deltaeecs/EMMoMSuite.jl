# 判别实验：合成 EFIE 数据的实际球谐带宽 vs 粗网格可表示带宽
#   A. 在细网格上用 LS 展开 F_θ/F_ϕ 到 degree<=20，看能量随 degree 的衰减
#   B. 按盒子尺寸缩放源几何（臂长 ~ 0.25*盒子边长）后，修正版局部 LS 的 εi
#
# 用法: julia --project=. scripts/lebedev_data_spectrum.jl

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncation_kernel
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG
import EMMoMSuite.FastAlgorithms.Lebedev.SHInterp as SH

paper_metric(W, Tte, Pte) =
    vec(maximum(abs.(W * Tte .- Pte); dims = 1) ./ maximum(abs.(Pte); dims = 1))

"""原版生成器：固定臂长 0.12λ（当前仓库代码）"""
function dataset_v1(pk, pt, rel_l, nρ, npos; seed = 1)
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

"""缩放版生成器：臂长/偏移按 rel_l 缩放（与盒子带宽一致）"""
function dataset_v2(pk, pt, rel_l, nρ, npos; seed = 1, arm_scale = 0.5)
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
    # 与 generate_dataset_on_poles 相同的结构，但几何按 rel_l 缩放
    JK_0 = im * 2π
    ws = [-4 / 5, 9 / 20, 9 / 20, 9 / 20, 9 / 20]
    for ir in axes(rbmrps, 2), iρ in axes(ρhats, 2)
        rvec = rbmrps[:, ir]
        rvecp = rvec
        rvecm = rvec .+ DG.random_rhat() .* (rand() * 0.12 * rel_l)
        offsets = [[DG.random_rhat()...] .* (rand() * 0.03 * rel_l) for _ in eachindex(ws)]
        offsets[1] .= 0
        r0p = rvecp .+ rvecp .- rvecm .+ DG.random_rhat() .* 0.005 * rel_l
        r0m = rvecm .+ rvecm .- rvecp .+ DG.random_rhat() .* 0.005 * rel_l
        for (iHats, arr) in ((rHatsC, view(tArray, :, :, iρ, ir)), (rHatsF, view(pArray, :, :, iρ, ir)))
            for iPole in eachindex(iHats)
                poler̂θϕ = iHats[iPole]
                for iw in eachindex(ws)
                    rp = rvecp .+ offsets[iw]
                    rm = rvecm .+ offsets[iw]
                    ρhatp_iw = rp .- r0p
                    ρhatm_iw = rm .- r0m
                    wpexptemp = ws[iw] * exp(JK_0 * dot(poler̂θϕ.r̂, rp))
                    wmexptemp = ws[iw] * exp(JK_0 * dot(poler̂θϕ.r̂, rm))
                    arr[iPole, 1] += dot(poler̂θϕ.θhat, ρhatp_iw) * wpexptemp
                    arr[iPole, 1] -= dot(poler̂θϕ.θhat, ρhatm_iw) * wmexptemp
                    arr[iPole, 2] += dot(poler̂θϕ.ϕhat, ρhatp_iw) * wpexptemp
                    arr[iPole, 2] -= dot(poler̂θϕ.ϕhat, ρhatm_iw) * wmexptemp
                end
            end
        end
    end
    return reshape(tArray, nt * 2, :), reshape(pArray, np_ * 2, :), tnodes, pnodes
end

"""细网格上 F_θ/F_ϕ 的 SH 频谱（单列辐射函数，累计到 99.9% 能量）"""
function sh_spectrum(pnodes, fine_data, nf, comp; jcol = 1, Lmax = 22)
    Y = SH.realSHmatrix(pnodes, Lmax)
    rows = comp == 1 ? (1:nf) : (nf+1:2nf)
    f = real(fine_data[rows, jcol])                        # 实部近似谱（带宽判断足够）
    c = pinv(Y; rtol = 1e-12) * f
    en = zeros(Lmax + 1)
    idx = 0
    for l in 0:Lmax
        for m in -l:l
            idx += 1
            en[l+1] += c[idx]^2
        end
    end
    tot = sum(en)
    cum = 0.0
    L99 = 0
    for l in 0:Lmax
        cum += en[l+1]
        if cum >= 0.999 * tot
            L99 = l
            break
        end
    end
    return en, tot, L99
end

function fit_complex_local(tnodes, pnodes, T, P, k; flag)
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
        tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
        pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
        nt, nf = size(tnodes, 2), size(pnodes, 2)
        println("\n############ pk=$pk -> pt=$pt （Lb=$Lb, rel_l=$(round(rel_l, digits=4))）############")

        T1, P1, _, _ = dataset_v1(pk, pt, rel_l, 40, 120; seed = 20260808)
        enθ, _, L99θ = sh_spectrum(pnodes, P1, nf, 1)
        enϕ, _, L99ϕ = sh_spectrum(pnodes, P1, nf, 2)
        @printf("A. 原版生成器(固定臂0.12λ): F_θ 99.9%%能量阶=%d, F_ϕ 99.9%%能量阶=%d（粗网格容量 degree<=%d）\n",
            L99θ, L99ϕ, Lb)
        @printf("   各阶能量占比: %s\n",
            join([@sprintf("l=%d:%.1e", l, enθ[l+1] / sum(enθ)) for l in 0:min(20, Lb + 3)], ", "))

        flag = trunc(Int, 0.8 * size(T1, 2))
        Tte1, Pte1 = T1[:, (flag+1):end], P1[:, (flag+1):end]
        W1 = fit_complex_local(tnodes, pnodes, T1, P1, 8; flag = flag)
        @printf("   原版数据 k=8 修正 LS: εi均值=%.3e\n", mean(paper_metric(W1, Tte1, Pte1)))

        T2, P2, _, _ = dataset_v2(pk, pt, rel_l, 40, 120; seed = 20260808)
        flag2 = trunc(Int, 0.8 * size(T2, 2))
        Tte2, Pte2 = T2[:, (flag2+1):end], P2[:, (flag2+1):end]
        enθ2, _, L99θ2 = sh_spectrum(pnodes, P2, nf, 1)
        enϕ2, _, L99ϕ2 = sh_spectrum(pnodes, P2, nf, 2)
        @printf("B. 缩放生成器(臂~rel_l): F_θ 99.9%%能量阶=%d, F_ϕ 99.9%%能量阶=%d\n", L99θ2, L99ϕ2)
        W2 = fit_complex_local(tnodes, pnodes, T2, P2, 8; flag = flag2)
        εi2 = paper_metric(W2, Tte2, Pte2)
        @printf("   缩放数据 k=8 修正 LS: εi均值=%.3e εi最大=%.3e\n", mean(εi2), maximum(εi2))
    end
end

main()
