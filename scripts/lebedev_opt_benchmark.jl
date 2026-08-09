# 复现原始插值权重管线 + 结构优化对比基准：
#   A. 原始方法：IDW(kNN) 初始化 + 逐行 pinv（仅实部约束）——忠实复现 runpinvCal
#   B. 球谐精确插值：W = Y_fine * pinv(Y_coarse)（限带函数的机器精度上界）
#   C. 局部约束最小二乘：支撑取角距帽内粗点，强制对 degree<=L_loc 球谐精确（行和=1）
#   D. 八面体群轨道压缩：代表性方向算权重 + 旋转复用（验证与朴素 C 完全一致 + 计时）
#
# 用法: julia --project=. scripts/lebedev_opt_benchmark.jl

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncation_kernel
using LinearAlgebra, SparseArrays, Random, Printf, Statistics
using SpecialFunctions: gamma

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.pinv2interpW as PW
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG

"""实球谐基（Condon-Shortley），节点 (3,n) -> Y (n, (Lmax+1)^2)"""
function realSHmatrix(nodes::AbstractMatrix, Lmax::Int)
    n = size(nodes, 2)
    Y = zeros(Float64, n, (Lmax + 1)^2)
    for i in 1:n
        x = Float64(nodes[3, i])
        s = sqrt(max(0.0, 1.0 - x * x))
        cphi = nodes[1, i] / max(s, 1e-300)
        sphi = nodes[2, i] / max(s, 1e-300)
        P = zeros(Lmax + 1, Lmax + 1)
        P[1, 1] = 1.0
        for l in 1:Lmax
            P[l+1, l+1] = -(2l - 1) * s * P[l, l]
            for m in 0:(l-1)
                pprev = l >= 2 ? P[l-1, m+1] : 0.0
                P[l+1, m+1] = ((2l - 1) * x * P[l, m+1] - (l - 1 + m) * pprev) / (l - m)
            end
        end
        col = 0
        for l in 0:Lmax, m in -l:l
            col += 1
            am = abs(m)
            nrm = sqrt((2l + 1) / (4π) * gamma(l - am + 1) / gamma(l + am + 1))
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

"""按现行 runpinvCal 布局生成 (pk, pt) 数据集（可缩小规模）"""
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

"""原始管线拟合（IDW + 逐行 pinv，仅实部约束）"""
function fit_original(tnodes, pnodes, T, P, nInterp; flag)
    xx2D = vcat(real(T[:, 1:flag]), imag(T[:, 1:flag]))
    yy2D = vcat(real(P[:, 1:flag]), imag(P[:, 1:flag]))
    w = PW.interpWeightsInitial(tnodes, pnodes; nInterp = nInterp)
    PW.pinv2W!(w, nInterp, xx2D, yy2D)
    return w
end

"""八面体群真旋转（24 个）"""
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

"""每个细点 -> (轨道代表下标, 使 M*r̂_rep = r̂_node 的旋转 M)"""
function orbit_reps(fine_nodes)
    R = octahedral_rotations()
    nf = size(fine_nodes, 2)
    reps = Vector{Int}(undef, nf)
    rots = Vector{Matrix{Float64}}(undef, nf)
    rep2nodes = Dict{Int,Vector{Int}}()
    assigned = falses(nf)
    for i in 1:nf
        assigned[i] && continue
        rep2nodes[i] = Int[]
        for j in i:nf
            assigned[j] && continue
            for Rm in R
                if norm(Rm * fine_nodes[:, i] .- fine_nodes[:, j]) < 1e-8
                    push!(rep2nodes[i], j)
                    reps[j] = i
                    rots[j] = Rm
                    assigned[j] = true
                    break
                end
            end
        end
    end
    return reps, rots, rep2nodes
end

"""局部约束最小二乘：每细点取角距帽内粗点为支撑，强制 degree<=L_loc 球谐精确"""
function local_constrained_weights(tnodes, pnodes, Lloc; cap_rad, grow = true)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    Yc = realSHmatrix(tnodes, Lloc)
    Yf = realSHmatrix(pnodes, Lloc)
    m = (Lloc + 1)^2
    W = spzeros(nf, nt)
    θ = cap_rad
    for i in 1:nf
        ang = [2asin(min(norm(tnodes[:, j] .- pnodes[:, i]) / 2, 1.0)) for j in 1:nt]
        S = findall(ang .< θ)
        if grow
            while length(S) < m
                θ *= 1.2
                S = findall(ang .< θ)
                θ >= π && break
            end
        end
        A = Yc[S, :]                    # |S| x m
        b = Yf[i, :]                    # m
        w = b' * pinv(A; rtol = 1e-12)  # 最小范数解: w A = b'（SVD，条件数稳健）
        W[i, S] = w[:]
    end
    return W
end

"""轨道压缩版：只对每轨道代表求解，旋转填充其余行"""
function orbit_compressed_weights(tnodes, pnodes, Lloc; cap_rad, grow = false)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    Yc = realSHmatrix(tnodes, Lloc)
    Yf = realSHmatrix(pnodes, Lloc)
    m = (Lloc + 1)^2
    reps, rots, rep2nodes = orbit_reps(pnodes)
    # 粗点坐标 -> 下标
    cdict = Dict{Tuple{Float64,Float64,Float64},Int}()
    for j in 1:nt
        cdict[(round(tnodes[1, j], digits = 8), round(tnodes[2, j], digits = 8), round(tnodes[3, j], digits = 8))] = j
    end
    W = spzeros(nf, nt)
    for (irep, nodes) in rep2nodes
        ang = [2asin(min(norm(tnodes[:, j] .- pnodes[:, irep]) / 2, 1.0)) for j in 1:nt]
        S = findall(ang .< cap_rad)
        while grow && length(S) < m
            cap_rad *= 1.2
            S = findall(ang .< cap_rad)
            cap_rad >= π && break
        end
        A = Yc[S, :]
        b = Yf[irep, :]
        w = b' * pinv(A; rtol = 1e-12)
        for j in nodes
            M = rots[j]
            cols = [cdict[(round(M[1,1]*tnodes[1,s] + M[1,2]*tnodes[2,s] + M[1,3]*tnodes[3,s], digits = 8),
                           round(M[2,1]*tnodes[1,s] + M[2,2]*tnodes[2,s] + M[2,3]*tnodes[3,s], digits = 8),
                           round(M[3,1]*tnodes[1,s] + M[3,2]*tnodes[2,s] + M[3,3]*tnodes[3,s], digits = 8))] for s in S]
            W[j, cols] = w[:]
        end
    end
    return W
end

function bandlimited_test(W, tnodes, pnodes, Lb; trials = 100, seed = 123)
    Yc = realSHmatrix(tnodes, Lb)
    Yf = realSHmatrix(pnodes, Lb)
    Random.seed!(seed)
    worst = 0.0
    for _ in 1:trials
        c = randn((Lb + 1)^2)
        fc, ff = Yc * c, Yf * c
        worst = max(worst, maximum(abs.(W * fc .- ff)) / maximum(abs.(ff)))
    end
    return worst
end

function main()
    cases = [
        (13, 27),
        (27, 53),
        (65, 131),
    ]
    for (pk, pt) in cases
        Lb = (pk - 1) ÷ 2
        nInterp = pk < 20 ? 9 : 8
        println("\n############ pk=$pk -> pt=$pt （Lb=$Lb）############")
        tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
        pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
        nt, nf = size(tnodes, 2), size(pnodes, 2)
        @printf("点数: 粗层 %d, 细层 %d\n", nt, nf)
        sweeps = if Lb <= 6
            [(min(Lb, 3), 0.6), (min(Lb, 4), 0.8), (Lb, 1.0), (Lb, 1.2)]
        elseif Lb <= 13
            [(3, 0.6), (4, 0.8), (6, 1.0), (13, 1.2)]
        else
            [(3, 0.35), (6, 0.5)]
        end

        # B: 球谐精确（机器精度上界）
        t0 = time()
        Yc = realSHmatrix(tnodes, Lb)
        Yf = realSHmatrix(pnodes, Lb)
        Wexact = Yf * pinv(Yc; rtol = 1e-13)
        @printf("B 球谐精确: 误差=%.3e, 非零=%d (稠密), 构造 %.1fs\n",
            bandlimited_test(Wexact, tnodes, pnodes, Lb), count(!iszero, Wexact), time() - t0)

        # C: 局部约束 LS，扫支撑角距
        println("C 局部约束 LS（degree<=L_loc 精确）:")
        for (Lloc, cap) in sweeps
            m = (Lloc + 1)^2
            m >= nt && continue
            t0 = time()
            Wloc = local_constrained_weights(tnodes, pnodes, Lloc; cap_rad = cap)
            err = bandlimited_test(Wloc, tnodes, pnodes, Lb)
            rowsum = maximum(abs.(Wloc * ones(nt) .- 1))
            @printf("  L_loc=%2d cap=%.2f 误差(deg<=%d)=%.3e nnz=%d (%.1f/行) 行和偏差=%.1e 构造 %.1fs\n",
                Lloc, cap, Lb, err, nnz(Wloc), nnz(Wloc) / nf, rowsum, time() - t0)
        end

        if nf <= 1000
            # D: 轨道压缩验证（固定支撑帽、无增长，保证 naive 与 orbit 支撑一致）
            LlocD = min(Lb, 4)
            capD = LlocD <= 3 ? 1.0 : 1.4
            t0 = time()
            Wnaive = local_constrained_weights(tnodes, pnodes, LlocD; cap_rad = capD, grow = false)
            tnaive = time() - t0
            t0 = time()
            Worb = orbit_compressed_weights(tnodes, pnodes, LlocD; cap_rad = capD)
            torb = time() - t0
            @printf("D 轨道压缩: max|W_naive-W_orbit|=%.2e, 构造 %.2fs vs naive %.2fs (加速 %.1fx), 轨道数=%d\n",
                maximum(abs.(Wnaive .- Worb)), torb, tnaive, tnaive / max(torb, 1e-9),
                length(unique(orbit_reps(pnodes)[1])))

            # A: 原始管线（数据集 + IDW + 逐行 pinv）
            rel_l = EMMoMSuite.Utilities.find_zero_bisection(x -> truncation_kernel(x) - (pk + 1) / 2, 0)
            nρ, npos = 40, 120
            T, P, _, _ = make_dataset(pk, pt, rel_l, nρ, npos; seed = 20260808)
            flag = trunc(Int, 0.8 * size(T, 2))
            Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]
            t0 = time()
            Worig = fit_original(tnodes, pnodes, T, P, nInterp; flag = flag)
            tfit = time() - t0
            err_bl = bandlimited_test(Worig[1:nf, 1:nt], tnodes, pnodes, Lb)
            err_efie = mean(PW.acc(Worig, Tte, Pte))
            Wexactv = [Wexact zeros(nf, nt); zeros(nf, nt) Wexact]
            err_sh_efie = mean(PW.acc(Wexactv, Tte, Pte))
            @printf("A 原始管线: EFIE留出=%.3e, 限带函数=%.3e, nnz=%d (%.1f/行), 拟合 %.1fs\n",
                err_efie, err_bl, nnz(Worig), nnz(Worig) / (2nf), tfit)
            @printf("  球谐精确在 EFIE 留出集: %.3e\n", err_sh_efie)
        else
            println("A/D 跳过（细层点数 > 1000：原始管线需 ~GB 级数据集/逐行 pinv，不可行——这正是要优化的点）")
        end
    end
end

main()
