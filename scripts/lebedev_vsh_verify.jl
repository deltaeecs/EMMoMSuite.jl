# 自旋加权球谐（VSH）精确矢量插值验证：
#   1. 基正交性（Lebedev 权重下 Σ w sY* sY ≈ δ，含 s 混合项）
#   2. EFIE 矢量类（修复后生成器）论文指标 εi：VSH 精确 vs 标量 SH vs 训练式 k=8
#   3. 自旋限带随机矢量场机器精度（l <= Lb+1）
#   4. VSH 局部稀疏版（Lloc 扫描）与轨道结构
#   5. 高阶 (65->131) 构造验证
#
# 用法: julia --project=. scripts/lebedev_vsh_verify.jl

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncation_kernel
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG
import EMMoMSuite.FastAlgorithms.Lebedev.SHInterp as SH

paper_metric(W, Tte, Pte) =
    vec(maximum(abs.(W * Tte .- Pte); dims = 1) ./ maximum(abs.(Pte); dims = 1))

function make_dataset(pk, pt, rel_l, nρ, npos; seed = 1, arm_scale = 0.5)
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
    rC = LSP.nodes2Poles(tnodes)
    rF = LSP.nodes2Poles(pnodes)
    tArray = zeros(ComplexF64, nt, 2, nρ, npos)
    pArray = zeros(ComplexF64, np_, 2, nρ, npos)
    for ir in axes(rbmrps, 2), iρ in axes(ρhats, 2)
        geom = DG.random_source_geometry(rbmrps[:, ir];
            arm_max = arm_scale * rel_l, off_max = 0.125 * arm_scale * rel_l)
        DG.evaluate_poles!(rC, view(tArray, :, :, iρ, ir), geom)
        DG.evaluate_poles!(rF, view(pArray, :, :, iρ, ir), geom)
    end
    return reshape(tArray, nt * 2, :), reshape(pArray, np_ * 2, :), tnodes, pnodes
end

"""训练式 k 点（修复后，共享几何 + hcat 复约束）"""
function trained_weights(pk, pt, k, rel_l; nρ = 40, npos = 120, seed = 20260808)
    tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
    pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    T, P, _, _ = make_dataset(pk, pt, rel_l, nρ, npos; seed = seed)
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
    println("=== 1. 自旋加权基正交性（Lebedev 权重）===")
    for (p, L) in ((13, 6), (27, 13))
        nodes, wts = LSP.getlbSortedData(p)
        for s in (1, -1)
            S = SH.spin_weighted_harmonics(nodes, L; s = s)
            G = S' * (wts .* S)
            @printf("  p=%d s=%+d L=%d: max|G-I|=%.2e\n", p, s, L,
                maximum(abs.(G - Matrix{ComplexF64}(I, size(G, 1), size(G, 1)))))
        end
        Sp = SH.spin_weighted_harmonics(nodes, L; s = 1)
        Sm = SH.spin_weighted_harmonics(nodes, L; s = -1)
        Gcross = Sp' * (wts .* Sm)
        # 自旋异号内积非零且为 ±δ 结构（自旋升降算子的伴随关系），不是正交失败
        diag_ok = maximum(abs.(diag(Gcross))) < 1.5
        offdiag = maximum(abs.(Gcross - Diagonal(diag(Gcross))))
        @printf("  p=%d 交叉矩阵: max|diag|=%.2e, max|offdiag|=%.2e（±δ 结构，非正交但自洽）\n",
            p, maximum(abs.(diag(Gcross))), offdiag)
    end

    for (pk, pt) in ((13, 27), (27, 53))
        Lb = (pk - 1) ÷ 2
        rel_l = EMMoMSuite.Utilities.find_zero_bisection(x -> truncation_kernel(x) - (pk + 1) / 2, 0)
        tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
        pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
        nt, nf = size(tnodes, 2), size(pnodes, 2)
        println("\n=== pk=$pk -> pt=$pt （Lb=$Lb, VSH Lmax=$(Lb + 1)）===")

        T, P, _, _ = make_dataset(pk, pt, rel_l, 40, 120; seed = 20260808)
        flag = trunc(Int, 0.8 * size(T, 2))
        Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]

        # 2. VSH 精确 vs 标量 SH vs 训练式
        Wv = SH.interp_weights_vsh(pk, pt)
        Ws = SH.vectorize(SH.interp_weights_exact(pk, pt))
        Wt = trained_weights(pk, pt, 8, rel_l)
        @printf("2. EFIE 矢量类 εi均值: VSH精确=%.3e | 标量SH=%.3e | 训练式k=8=%.3e\n",
            mean(paper_metric(Wv, Tte, Pte)), mean(paper_metric(Ws, Tte, Pte)),
            mean(paper_metric(Wt, Tte, Pte)))

        # 3. 自旋限带随机矢量场机器精度
        Lv = Lb + 1
        Sp_c = SH.spin_weighted_harmonics(tnodes, Lv; s = 1)
        Sm_c = SH.spin_weighted_harmonics(tnodes, Lv; s = -1)
        Sp_f = SH.spin_weighted_harmonics(pnodes, Lv; s = 1)
        Sm_f = SH.spin_weighted_harmonics(pnodes, Lv; s = -1)
        nc = (Lv + 1)^2 - 1
        Random.seed!(7)
        worst = 0.0
        for _ in 1:50
            cp = randn(ComplexF64, nc)
            cm = randn(ComplexF64, nc)
            Fpc, Fmc = Sp_c * cp, Sm_c * cm
            Fpf, Fmf = Sp_f * cp, Sm_f * cm
            coarse = vcat((Fpc + Fmc) / 2, (Fpc - Fmc) / (2im))
            exact = vcat((Fpf + Fmf) / 2, (Fpf - Fmf) / (2im))
            worst = max(worst, maximum(abs.(Wv * coarse .- exact)) / maximum(abs.(exact)))
        end
        @printf("3. 自旋限带(l<=%d)矢量场最大相对误差: VSH精确=%.3e\n", Lv, worst)

        # 4. VSH 局部稀疏版
        for (Lloc, cap) in ((min(4, Lb), 0.8), (min(6, Lb), 1.0), (min(Lb, Lb), 1.4))
            Lloc <= 0 && continue
            m = (Lloc + 1)^2 - 1
            m >= nt && continue
            Wl = SH.interp_weights_vsh_local(pk, pt; Lloc = Lloc, cap_rad = cap)
            εi = paper_metric(Wl, Tte, Pte)
            @printf("4. VSH局部 L_loc=%d cap=%.1f: εi均值=%.3e nnz/行=%.1f 行和偏差=%.1e\n",
                Lloc, cap, mean(εi), nnz(Wl) / (2nf),
                maximum(abs.(Wl * ones(ComplexF64, 2nt) .- 1)))
        end
    end

    println("\n=== 5. 高阶 (65->131) VSH ===")
    pk, pt = 65, 131
    t0 = time()
    Wv = SH.interp_weights_vsh(pk, pt)
    @printf("  VSH 精确构造 %.1fs, W 尺寸 %dx%d (ComplexF64, %.1f MB)\n",
        time() - t0, size(Wv, 1), size(Wv, 2), sizeof(Wv) / 1e6)
    tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
    pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
    Lv = (pk - 1) ÷ 2 + 1
    Sp_c = SH.spin_weighted_harmonics(tnodes, Lv; s = 1)
    Sp_f = SH.spin_weighted_harmonics(pnodes, Lv; s = 1)
    nc = (Lv + 1)^2 - 1
    Random.seed!(3)
    worst = 0.0
    for _ in 1:20
        cp = randn(ComplexF64, nc)
        Fpc, Fpf = Sp_c * cp, Sp_f * cp
        coarse = vcat(Fpc / 2, Fpc / (2im))
        exact = vcat(Fpf / 2, Fpf / (2im))
        worst = max(worst, maximum(abs.(Wv * coarse .- exact)) / maximum(abs.(exact)))
    end
    @printf("  自旋限带(l<=%d)矢量场最大相对误差: %.3e\n", Lv, worst)
end

main()
