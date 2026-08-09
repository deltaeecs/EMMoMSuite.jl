# Fibonacci 网格能否用更少点数达到同 εi（EFIE 类，笛卡尔局部稀疏 L_loc=3）
# 对比 Lebedev (266->974) vs Fibonacci 更小规模
#
# 用法: julia --project=. scripts/lebedev_fibonacci_test.jl

using EMMoMSuite
using LinearAlgebra, SparseArrays, Random, Printf, Statistics
import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG
import EMMoMSuite.FastAlgorithms.Lebedev.SHInterp as SH
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncation_kernel

rel_l_for(pk) = begin
    f(x) = truncation_kernel(x) - (pk + 1) / 2
    lo, hi = 1e-4, 20.0
    for _ in 1:80
        mid = (lo + hi) / 2
        f(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    (lo + hi) / 2
end

function fibonacci_grid(n::Int)
    nodes = zeros(Float64, 3, n)
    ga = π * (3 - sqrt(5))
    for i in 1:n
        z = 1 - 2 * (i - 0.5) / n
        r = sqrt(1 - z * z)
        φ = ga * i
        nodes[1, i] = r * cos(φ)
        nodes[2, i] = r * sin(φ)
        nodes[3, i] = z
    end
    return nodes
end

paper_metric(W, Tte, Pte) =
    vec(maximum(abs.(W * Tte .- Pte); dims = 1) ./ maximum(abs.(Pte); dims = 1))

function data_on(nodes, pk, rel_l; ns = 3000, seed = 20260808)
    n = size(nodes, 2)
    Random.seed!(seed)
    poles = LSP.nodes2Poles(nodes)
    arm = min(0.12, 0.5 * rel_l)
    rscale = max(rel_l / 2 - arm, 1e-4)
    D = zeros(ComplexF64, n, 2, ns)
    for s in 1:ns
        rvec = DG.random_rvec() .* rscale
        geom = DG.random_source_geometry(rvec; arm_max = arm, off_max = 0.125 * arm)
        DG.evaluate_poles!(poles, view(D, :, :, s), geom)
    end
    return reshape(D, n * 2, :)
end

function main()
    pk, pt = 27, 53
    rel_l = rel_l_for(pk)
    println("pk=$pk -> pt=$pt, rel_l=$(round(rel_l, digits = 4))")

    # 参考：Lebedev
    n_lb_c, n_lb_f = 266, 974
    Tlb = data_on(LSP.get_t_nodes((pk - 1) ÷ 2), pk, rel_l)
    Plb = data_on(LSP.get_t_nodes((pt - 1) ÷ 2), pk, rel_l)
    flag = trunc(Int, 0.8 * size(Tlb, 2))
    Tte, Pte = Tlb[:, (flag+1):end], Plb[:, (flag+1):end]

    Wl = SH.interp_weights_cart_local(pk, pt; Lloc = 3, cap_rad = 0.6)
    @printf("Lebedev %d->%d  笛卡尔局部 L_loc=3: εi=%.3e\n", n_lb_c, n_lb_f,
        mean(paper_metric(Wl, Tte, Pte)))

    # Fibonacci 各规模
    for (nc, nf) in ((266, 974), (200, 740), (150, 560), (120, 450))
        tc = fibonacci_grid(nc)
        tf = fibonacci_grid(nf)
        poles_c = LSP.nodes2Poles(tc)
        poles_f = LSP.nodes2Poles(tf)
        Lloc = 3
        Yc = SH.realSHmatrix(tc, Lloc)
        Yf = SH.realSHmatrix(tf, Lloc)
        cap = 0.6
        W = SH._cart_core(tc, tf, poles_c, poles_f, Yc, Yf, Lloc; cap_rad = cap)
        # 数据要在 Fibonacci 节点上重新采样：细层用 Pte 的位置不同 -> 重新生成
        Tf = data_on(tc, pk, rel_l)
        Pf = data_on(tf, pk, rel_l)
        fflag = trunc(Int, 0.8 * size(Tf, 2))
        Tte_f, Pte_f = Tf[:, (fflag+1):end], Pf[:, (fflag+1):end]
        εi = paper_metric(W, Tte_f, Pte_f)
        @printf("Fibonacci %4d->%4d 笛卡尔局部 L_loc=3: εi=%.3e  nnz/行=%.1f\n",
            nc, nf, mean(εi), nnz(W) / (2nf))
    end
end

main()
