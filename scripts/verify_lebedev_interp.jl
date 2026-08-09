# 验证 MoM_Lebedev(现 FastAlgorithms.Lebedev) 的关键疑点：
#   A. LVI.levelIntegralInfoCal 把 levelCubeEdgel*2π/λ (ka) 传给 truncation_kernel(期望 a/λ)
#   B. p=131 (5810 点) 是否因 off-by-one 成为死代码
#   C. 当前 IDW+逐行 pinv 训练插值权重 vs 球谐(real-SH)精确插值权重 的精度/稠密度对比
#   D. 训练数据只使用复数的实部（vcat(real,imag) 与 2np 行不匹配）
#   E. 训练权重不可复现（无 RNG 种子）
#
# 用法: julia --project=. scripts/verify_lebedev_interp.jl

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.Lebedev
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation:
    truncation_kernel, levelIntegralInfoCal
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.pinv2interpW as PW
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG

"""
实球谐基（包含 Condon-Shortley 相位），节点 (3,n) -> Y (n, (Lmax+1)^2)。
正交性由 Lebedev 权重近似验证：Σ_i w_i Y_k(r̂_i) Y_j(r̂_i) ≈ δ_kj。
"""
function realSHmatrix(nodes::AbstractMatrix, Lmax::Int)
    n = size(nodes, 2)
    Y = zeros(Float64, n, (Lmax + 1)^2)
    for i in 1:n
        x = Float64(nodes[3, i])
        s = sqrt(max(0.0, 1.0 - x * x))
        cphi = nodes[1, i] / max(s, 1e-300)
        sphi = nodes[2, i] / max(s, 1e-300)
        P = zeros(Lmax + 1, Lmax + 1)   # P[l+1, m+1] = P_l^m, m<=l
        P[1, 1] = 1.0
        for l in 1:Lmax
            P[l+1, l+1] = -(2l - 1) * s * P[l, l]                     # P_l^l
            for m in 0:(l-1)                                           # 一般递推
                pprev = l >= 2 ? P[l-1, m+1] : 0.0                     # P_{-1}^m := 0
                P[l+1, m+1] = ((2l - 1) * x * P[l, m+1] - (l - 1 + m) * pprev) / (l - m)
            end
        end
        col = 0
        for l in 0:Lmax, m in -l:l
            col += 1
            am = abs(m)
            nrm = sqrt((2l + 1) / (4π) * factorial(l - am) / factorial(l + am))
            val = nrm * P[l+1, am+1]
            if m == 0
                Y[i, col] = val
            elseif m > 0
                Y[i, col] = sqrt(2) * val * cos(am * atan(sphi, cphi))
            else
                Y[i, col] = sqrt(2) * val * sin(am * atan(sphi, cphi))
            end
        end
    end
    return Y
end

"""按现行 runpinvCal 的数据布局生成 (pk, pt) 训练数据集（可缩小规模）"""
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

"""现行算法：IDW 初始化 + 逐行 pinv（仅实部约束，忠实复现 pinv2W! 的调用方式）"""
function fit_shipped(tnodes, pnodes, T, P, nInterp; flag)
    xx2D = vcat(real(T[:, 1:flag]), imag(T[:, 1:flag]))
    yy2D = vcat(real(P[:, 1:flag]), imag(P[:, 1:flag]))
    w = PW.interpWeightsInitial(tnodes, pnodes; nInterp = nInterp)
    PW.pinv2W!(w, nInterp, xx2D, yy2D)
    return w
end

"""修正版：同一数据量，同时用实部+虚部做约束"""
function fit_both_parts(tnodes, pnodes, T, P, nInterp; flag)
    xxF = hcat(real(T[:, 1:flag]), imag(T[:, 1:flag]))
    yyF = hcat(real(P[:, 1:flag]), imag(P[:, 1:flag]))
    w = PW.interpWeightsInitial(tnodes, pnodes; nInterp = nInterp)
    for irow in axes(w, 1)
        wrow = w[irow, :]
        nnz(wrow) == 1 && continue
        nzind = wrow.nzind
        wi = view(yyF, irow:irow, :) * pinv(@view xxF[nzind, :])
        w[irow, nzind] .= reshape(wi, :)
    end
    return w
end

function main()
    println("=== A. truncation_kernel 参数不一致（LVI 传 ka，基实现传 a/λ）===")
    for a in (0.125, 0.25, 0.5)
        Lc = truncation_kernel(a)          # 正确: a 为以 λ 计的盒子边长
        Lbug = truncation_kernel(a * 2π)   # LVI.jl:39 实际传入 levelCubeEdgel*2π/λ
        @printf("a/λ=%.3f   L_correct(a/λ)=%8.3f   L_LVI(ka)=%8.3f   放大倍数=%.3f\n",
            a, Lc, Lbug, Lbug / Lc)
    end

    println("\n=== B. p=131 (5810 点) 死代码检查（2L+1 < max(p) 的 off-by-one）===")
    p2n = LSP.p2nDict
    pmax = maximum(keys(p2n))
    @printf("nodesSorted 最大阶 p=%d (n=%d)\n", pmax, p2n[pmax])
    for L in (64, 65, 66)
        cond = 2L + 1 < pmax
        @printf("  truncL=%3d: 2L+1=%3d, 条件 2L+1<%d = %s -> %s\n",
            L, 2L + 1, pmax, cond, cond ? "使用 Lebedev" : "触发 Lagrange 回退")
    end

    println("\n=== C. levelIntegralInfoCal(0.5, Val(:LbTrained1Step)) 实际选择 ===")
    truncL, poles = levelIntegralInfoCal(0.5, Val(:LbTrained1Step))
    @printf("返回 truncL=%d（正确公式应为 %d），Lebedev 点数=%d，权重和=%.6f (4π=%.6f)\n",
        truncL, floor(Int, truncation_kernel(0.5)), length(poles.Wθϕs), sum(poles.Wθϕs), 4π)

    println("\n=== D. pk=13 -> pt=27: 训练式权重 vs 球谐精确权重 ===")
    pk, pt = 13, 27
    nInterp = pk < 20 ? 9 : 8
    rel_l = EMMoMSuite.Utilities.find_zero_bisection(x -> truncation_kernel(x) - (pk + 1) / 2, 0)
    nρ, npos = 40, 120
    T, P, tnodes, pnodes = make_dataset(pk, pt, rel_l, nρ, npos; seed = 20260808)
    nt, np_ = size(tnodes, 2), size(pnodes, 2)
    @printf("粗层 p=%d (n=%d) -> 细层 p=%d (n=%d), 训练样本=%d\n", pk, nt, pt, np_, size(T, 2))
    N = size(T, 2)
    flag = trunc(Int, 0.8 * N)
    Tte = T[:, (flag+1):end]
    Pte = P[:, (flag+1):end]

    w0 = PW.interpWeightsInitial(tnodes, pnodes; nInterp = nInterp)     # 仅 IDW 初始化
    w = fit_shipped(tnodes, pnodes, T, P, nInterp; flag = flag)          # 现行（实部约束）
    wfix = fit_both_parts(tnodes, pnodes, T, P, nInterp; flag = flag)    # 修正（实+虚约束）

    @printf("IDW 初始化     EFIE 留出集平均相对误差: %.3e\n", mean(PW.acc(w0, Tte, Pte)))
    @printf("现行管线(实部) EFIE 留出集平均相对误差: %.3e\n", mean(PW.acc(w, Tte, Pte)))
    @printf("修正版(实+虚)  EFIE 留出集平均相对误差: %.3e\n", mean(PW.acc(wfix, Tte, Pte)))

    # 球谐精确插值：f 限带于 degree<=τk=6 (49 个系数)
    Lb = (pk - 1) ÷ 2
    Yc = realSHmatrix(tnodes, Lb)
    Yf = realSHmatrix(pnodes, Lb)
    _, wts = LSP.getlbSortedData(pk)
    Gram = Yc' * (wts .* Yc)
    @printf("SH 基正交性: 对角 max|G-I|=%.2e, 非对角 max|G|=%.2e\n",
        maximum(abs.(diag(Gram) .- 1)), maximum(abs.(Gram - Matrix{Float64}(I, size(Gram, 1), size(Gram, 1)))))
    Wsh = Yf * pinv(Yc; rtol = 1e-12)      # (np x nt) 稠密

    function maxrel(Wscalar, fc, ff)
        maximum(abs.(Wscalar * fc .- ff)) / maximum(abs.(ff))
    end

    Random.seed!(99)
    maxrel_shipped = maxrel_fixed = maxrel_sh = 0.0
    for trial in 1:200
        c = randn((Lb + 1)^2)
        fc, ff = Yc * c, Yf * c
        A = w[1:np_, 1:nt]      # 标量映射（ϕ 分量为 0 时交叉块不参与）
        maxrel_shipped = max(maxrel_shipped, maxrel(A, fc, ff))
        maxrel_fixed = max(maxrel_fixed, maxrel(wfix[1:np_, 1:nt], fc, ff))
        maxrel_sh = max(maxrel_sh, maxrel(Wsh, fc, ff))
    end
    @printf("限带函数(degree<=%d) 最大相对误差: 现行=%.3e, 修正版=%.3e, 球谐=%.3e\n",
        Lb, maxrel_shipped, maxrel_fixed, maxrel_sh)

    # 向量场（θ、ϕ 独立限带分量）: 现行 2x2 交叉耦合 vs 分块对角 SH
    Random.seed!(7)
    maxrel_v_shipped = maxrel_v_sh = 0.0
    Wshv = [Wsh zeros(np_, nt); zeros(np_, nt) Wsh]
    for trial in 1:100
        cθ, cϕ = randn((Lb + 1)^2), randn((Lb + 1)^2)
        fc = vcat(Yc * cθ, Yc * cϕ)
        ff = vcat(Yf * cθ, Yf * cϕ)
        maxrel_v_shipped = max(maxrel_v_shipped, maximum(abs.(w * fc .- ff)) / maximum(abs.(ff)))
        maxrel_v_sh = max(maxrel_v_sh, maximum(abs.(Wshv * fc .- ff)) / maximum(abs.(ff)))
    end
    @printf("向量场最大相对误差: 现行=%.3e, 分块对角球谐=%.3e\n", maxrel_v_shipped, maxrel_v_sh)

    @printf("现行 W: nnz=%d (每行 %.1f); 行和 max|A*1-1|=%.3e\n",
        nnz(w), nnz(w) / size(w, 1), maximum(abs.(w[1:np_, 1:nt] * ones(nt) .- 1)))
    @printf("球谐 W: 非零元=%d (稠密 %dx%d); 行和 max|Wsh*1-1|=%.3e\n",
        count(!iszero, Wsh), np_, nt, maximum(abs.(Wsh * ones(nt) .- 1)))

    println("\n=== E. 大阶数训练内存量级（generate_dataset_on_pkpt 预分配）===")
    for (pkE, ptE) in ((41, 83), (65, 131))
        nk = LSP.p2nDict[pkE]
        nf = LSP.p2nDict[ptE]
        @printf("pk=%d (n=%d) -> pt=%d (n=%d): pArray=%.2f GB, tArray=%.2f GB (各 50ρ x 500pos x ComplexF64)\n",
            pkE, nk, ptE, nf, 2 * nf * 50 * 500 * 16 / 1e9, 2 * nk * 50 * 500 * 16 / 1e9)
    end

    println("\n=== F. 训练权重可复现性（不同随机种子 -> 不同 W）===")
    T2, P2, _, _ = make_dataset(pk, pt, rel_l, nρ, npos; seed = 20260809)
    w2 = fit_shipped(tnodes, pnodes, T2, P2, nInterp; flag = flag)
    d1 = maximum(abs.(w * Tte .- Pte)) / maximum(abs.(Pte))
    d2 = maximum(abs.(w2 * Tte .- Pte)) / maximum(abs.(Pte))
    dw = maximum(abs.(w * Tte .- w2 * Tte)) / maximum(abs.(Pte))
    @printf("seed A 测试集最大相对误差=%.3e, seed B 测试集最大相对误差=%.3e, 两套 W 输出差异=%.3e\n",
        d1, d2, dw)

    println("\nDone.")
end

main()
