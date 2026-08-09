# 球面采样点效率对比：Lebedev vs Fibonacci 格点 vs GL×uniform 张量积
# 在相同点数 n 下比较：
#   - 标量 SH / 自旋加权 SH 合成矩阵条件数（决定插值鲁棒性）
#   - 最小角距（决定局部插值支撑成本）
#   - 限带标量场插值误差（稠密 SH 精确插值的实际容量）
#
# 用法: julia --project=. scripts/lebedev_grid_study.jl

using EMMoMSuite
using LinearAlgebra, Printf, Statistics, Random
import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.SHInterp as SH
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: octreeXWNCal

"""Fibonacci 格点：任意 n，准均匀"""
function fibonacci_grid(n::Int)
    nodes = zeros(Float64, 3, n)
    ga = π * (3 - sqrt(5))   # 黄金角
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

"""占位修正：直接构造笛卡尔节点"""
function gl_tensor_nodes(L::Int)
    Xcosθs, _ = octreeXWNCal(1.0, -1.0, L, :glq)
    Xθs = acos.(Xcosθs)
    Xϕs, _ = octreeXWNCal(0.0, 2π, L, :uni)
    pts = [begin
        sinθ, cosθ = sincos(θ)
        sinϕ, cosϕ = sincos(ϕ)
        [sinθ * cosϕ, sinθ * sinϕ, cosθ]
    end for ϕ in Xϕs for θ in Xθs]
    return reduce(hcat, pts)
end

function min_separation(nodes; nsample = 300)
    n = size(nodes, 2)
    qidx = sort(unique(rand(1:n, min(nsample, n))))
    mind = Inf
    for i in qidx
        for j in 1:n
            i == j && continue
            d = norm(nodes[:, i] .- nodes[:, j])
            mind = min(mind, 2asin(min(d / 2, 1.0)))
        end
    end
    return mind
end

function main()
    Random.seed!(42)
    for L in (6, 13)
        n_lb = LSP.p2nDict[2L + 1]          # Lebedev 点数
        n_fib = n_lb
        n_gl = 2 * (L + 1)^2
        println("\n===== L=$(L): Lebedev n=$(n_lb), Fibonacci n=$(n_fib), GL n=$(n_gl) =====")

        grids = [
            ("Lebedev", LSP.getlbSortedData(2L + 1)[1], n_lb),
            ("Fibonacci", fibonacci_grid(n_lb), n_fib),
            ("GL张量积", gl_tensor_nodes(L), n_gl),
        ]

        for (name, nodes, n) in grids
            Y = SH.realSHmatrix(nodes, L)
            κ = cond(Y)
            δ = min_separation(nodes)
            @printf("  %-10s n=%4d: κ(Y_%d)=%.2e, 最小角距=%.3f rad\n", name, n, L, κ, δ)
        end

        # 稠密 SH 精确插值容量：用 n 点网格对 degree<=L 随机标量场的重构误差
        # （取两个点数相近的网格比较：Lebedev vs Fibonacci，同一 L）
        nodes_lb = LSP.getlbSortedData(2L + 1)[1]
        nodes_fib = fibonacci_grid(n_lb)
        Ylb = SH.realSHmatrix(nodes_lb, L)
        Yfib = SH.realSHmatrix(nodes_fib, L)
        m = (L + 1)^2
        Random.seed!(7)
        for (name, Yg) in (("Lebedev", Ylb), ("Fibonacci", Yfib))
            # 用该网格样本重构随机限带场（自洽容量）
            worst = 0.0
            for _ in 1:50
                c = randn(m)
                f = Yg * c
                c_rec = pinv(Yg; rtol = 1e-12) * f
                worst = max(worst, norm(Yg * c_rec .- f) / norm(f))
            end
            @printf("    %-10s 自洽重构误差=%.2e\n", name, worst)
        end
        # 交叉：Lebedev 采样 -> Fibonacci 重构（细）与反向
        Random.seed!(9)
        worst_cross = 0.0
        for _ in 1:50
            c = randn(m)
            f_lb = Ylb * c
            f_fib = Yfib * c
            # Lebedev 样本 -> 用 SH 系数 -> Fibonacci 值
            c_rec = pinv(Ylb; rtol = 1e-12) * f_lb
            worst_cross = max(worst_cross, norm(Yfib * c_rec .- f_fib) / norm(f_fib))
        end
        @printf("    Lebedev→Fibonacci 插值误差=%.2e\n", worst_cross)
    end
end

main()
