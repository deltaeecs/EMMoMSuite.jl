# 分析 Lebedev 球面积分点结构：
#   - 八面体群（24 个真旋转）轨道分解
#   - 纬度环（|z| 分层）与 φ 分布规律（是否等间距）
#   - 网格旋转不变性验证、最小角距、正交性精确度
#
# 用法: julia --project=. scripts/lebedev_structure_analysis.jl

using EMMoMSuite
using LinearAlgebra, Printf, Statistics
import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP

"""八面体群的真旋转（24 个符号置换矩阵，行列式 +1）"""
function octahedral_rotations()
    mats = Matrix{Float64}[]
    for p1 in 1:3, p2 in 1:3, p3 in 1:3
        (p1 == p2 || p1 == p3 || p2 == p3) && continue
        for s1 in (-1.0, 1.0), s2 in (-1.0, 1.0), s3 in (-1.0, 1.0)
            M = zeros(3, 3)
            M[1, p1] = s1
            M[2, p2] = s2
            M[3, p3] = s3
            abs(det(M) - 1.0) < 1e-12 && push!(mats, M)
        end
    end
    return mats
end

"""每个节点的群规范形（对 24 个旋转取排序后的像），用于 O(n) 轨道分组"""
function group_orbit_representative(nodes)
    R = octahedral_rotations()
    n = size(nodes, 2)
    reps = Vector{Vector{Tuple{Float64,Float64,Float64}}}(undef, n)
    for i in 1:n
        imgs = [Rm * nodes[:, i] for Rm in R]
        key = sort([(round(x[1], digits = 6), round(x[2], digits = 6), round(x[3], digits = 6)) for x in imgs])
        reps[i] = key
    end
    return reps, R
end

"""按规范形分组，得到轨道（每轨道记录大小与首个节点）"""
function orbits_by_rep(reps)
    dict = Dict{Vector{Tuple{Float64,Float64,Float64}},Vector{Int}}()
    for (i, r) in enumerate(reps)
        push!(get!(dict, r, Int[]), i)
    end
    return collect(values(dict))
end

"""纬度环分析：按 |z| 分组，检查每环 φ 是否等间距"""
function ring_analysis(nodes)
    zs = abs.(nodes[3, :])
    zvals = sort(unique(round.(zs, digits = 8)))
    rows = String[]
    for zv in zvals
        idx = findall(abs.(zs .- zv) .< 1e-7)
        φs = sort(mod.(atan.(nodes[2, idx], nodes[1, idx]), 2π))
        m = length(φs)
        diffs = diff([φs; φs[1] + 2π])
        uniform = maximum(diffs) - minimum(diffs) < 1e-8
        push!(rows, @sprintf("  |z|=%.6f  n=%3d  φ等间距=%s  2π/n=%.6f  minΔφ=%.6f maxΔφ=%.6f",
            zv, m, uniform, 2π / m, minimum(diffs), maximum(diffs)))
    end
    return rows
end

"""最小角距（采样最近邻，kNN 暴力取 300 个查询点）"""
function min_angle(nodes; nsample = 300)
    n = size(nodes, 2)
    qidx = sort(unique(rand(1:n, min(nsample, n))))
    mind = Inf
    for i in qidx
        for j in 1:n
            i == j && continue
            d = norm(nodes[:, i] .- nodes[:, j])
            mind = min(mind, 2asin(min(d / 2, 1)))
        end
    end
    return mind
end

function main()
    println("=== Lebedev 点结构分析（nodesSorted 文件）===")
    R = octahedral_rotations()
    println("八面体群真旋转数 = ", length(R), "（应为 24）")

    for p in (3, 13, 27, 41, 65, 131)
        nodes, wts = LSP.getlbSortedData(p)
        n = size(nodes, 2)
        reps, _ = group_orbit_representative(nodes)
        orbs = orbits_by_rep(reps)
        orb_sizes = sort([length(o) for o in orbs])
        # 旋转不变性检查
        isinv = all(any(norm(Rm * nodes[:, i] .- nodes[:, j]) < 1e-8 for j in 1:n) for i in 1:n, Rm in R)
        # 环
        rings = ring_analysis(nodes)
        println("\n--- p=$p (n=$n) ---")
        println("轨道数=", length(orbs), " 轨道大小分布=", orb_sizes)
        println("24 旋转群不变性=", isinv, "  权重和/4π=", round(sum(wts) / 4π, digits = 12))
        println("纬度环数=", length(rings))
        foreach(println, rings[1:min(8, end)])
        length(rings) > 8 && println("  ... 其余 ", length(rings) - 8, " 环略")
    end

    println("\n=== 正交性精确度（Lebedev 权重下实球谐内积 ≈ δ）===")
    # 用简单复球谐 Y_l^0 的模方验证到阶 d：Σ w |Y_l^0|^2 = 1
    for p in (13, 27, 41)
        nodes, wts = LSP.getlbSortedData(p)
        x = nodes[3, :]
        println("p=$p: 对 l=0..min((p-1)÷2, 12) 检查 Σ w |P_l^0(x)|^2·(2l+1)/(4π) 与 1 的偏差")
        for l in 0:min(12, (p - 1) ÷ 2)
            # P_l^0 via recursion
            P0 = ones(length(x))
            P1 = x
            if l == 0
                Pl = P0
            elseif l == 1
                Pl = P1
            else
                Pprev, Pcur = P0, P1
                for k in 2:l
                    Pnext = ((2k - 1) .* x .* Pcur .- (k - 1) .* Pprev) ./ k
                    Pprev, Pcur = Pcur, Pnext
                end
                Pl = Pcur
            end
            val = sum(wts .* Pl .^ 2) * (2l + 1) / (4π)
            if abs(val - 1) > 1e-8
                @printf("  l=%2d  Σw|Y|^2=%.3e  <- 超精确度范围\n", l, val - 1)
            end
        end
        println("  （l ≤ (p-1)÷2 内全部通过则无输出）")
    end
end

main()
