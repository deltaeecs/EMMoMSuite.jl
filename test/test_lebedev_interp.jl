# Lebedev 插值权重新功能测试（轻量、无 MLFMA，供覆盖率套件使用）

using Test
using EMMoMSuite
using EMMoMSuite.FastAlgorithms.Lebedev
using EMMoMSuite.FastAlgorithms.Lebedev.SHInterp
import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncation_kernel
using LinearAlgebra, SparseArrays, Random, Statistics

@testset "Lebedev 插值权重" begin

    @testset "高阶节点（Fibonacci，p>131 无 Lebedev 数据集）" begin
        n133 = LSP.high_order_nodes(133)
        @test size(n133, 1) == 3
        @test all(abs.(sum(abs2, n133; dims = 1) .- 1) .< 1e-12)
        @test size(n133, 2) == round(Int, 4 / 3 * 66^2)   # n ≈ (4/3)τ²
        @test LSP.get_t_nodes(66) == n133                 # p=133 > 131 -> Fibonacci
        @test size(LSP.get_t_nodes(65), 2) == 5810        # p=131 仍是 Lebedev
    end

    @testset "实球谐正交性" begin
        nd, wts = LSP.getlbSortedData(13)
        Y = realSHmatrix(nd, 6)
        G = Y' * (wts .* Y)
        @test maximum(abs.(G - Matrix{Float64}(I, 49, 49))) < 1e-10
    end

    @testset "球谐精确插值：确定性 + 限带函数机器精度" begin
        W1 = interp_weights_exact(13, 27)
        W2 = interp_weights_exact(13, 27)
        @test W1 == W2
        tnodes = LSP.get_t_nodes(6)
        pnodes = LSP.get_t_nodes(13)
        Yc = realSHmatrix(tnodes, 6)
        Yf = realSHmatrix(pnodes, 6)
        Random.seed!(1)
        worst = 0.0
        for _ in 1:10
            c = randn(49)
            est = W1 * (Yc * c)
            worst = max(worst, maximum(abs.(est .- Yf * c)) / maximum(abs.(Yf * c)))
        end
        @test worst < 1e-10
    end

    @testset "局部约束与轨道压缩逐位一致" begin
        # 固定大支撑（cap=1.6 覆盖约一半球面），保证每点 |S| >= (Lloc+1)^2=16
        Wl = interp_weights_local(13, 27; Lloc = 3, cap_rad = 1.6, grow = false)
        Wo = interp_weights_local_orbit(13, 27; Lloc = 3, cap_rad = 1.6)
        @test maximum(abs.(Wl .- Wo)) < 1e-6
        @test maximum(abs.(Wl * ones(74) .- 1)) < 1e-8   # 行和 = 1（常数保持）
    end

    @testset "笛卡尔标量 SH 矢量插值（EFIE 类机器精度）" begin
        # 修正几何小数据集（源在盒内、粗/细层共享几何）
        f(x) = truncation_kernel(x) - (13 + 1) / 2
        lo, hi = 1e-4, 20.0
        for _ in 1:60
            mid = (lo + hi) / 2
            f(mid) > 0 ? (hi = mid) : (lo = mid)
        end
        rel_l = (lo + hi) / 2
        arm = min(0.12, 0.5 * rel_l)
        rscale = max(rel_l / 2 - arm, 1e-4)
        tnodes = LSP.get_t_nodes(6)
        pnodes = LSP.get_t_nodes(13)
        rC = LSP.nodes2Poles(tnodes)
        rF = LSP.nodes2Poles(pnodes)
        Random.seed!(20260808)
        T = zeros(ComplexF64, 74, 2, 8, 10)
        P = zeros(ComplexF64, 266, 2, 8, 10)
        for ir in 1:10, iρ in 1:8
            rvec = DG.random_rvec() .* rscale
            geom = DG.random_source_geometry(rvec; arm_max = arm, off_max = 0.125 * arm)
            DG.evaluate_poles!(rC, view(T, :, :, iρ, ir), geom)
            DG.evaluate_poles!(rF, view(P, :, :, iρ, ir), geom)
        end
        T = reshape(T, 148, :)
        P = reshape(P, 532, :)
        flag = trunc(Int, 0.8 * size(T, 2))
        Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]
        Wc = interp_weights_cart(13, 27)
        εi = maximum(abs.(Wc * Tte .- Pte); dims = 1) ./ maximum(abs.(Pte); dims = 1)
        @test mean(εi) < 1e-6
    end

    @testset "混合权重：确定性 + 数据一致性" begin
        Wh1 = interp_weights_hybrid(13, 27; Lloc = 3, support_scale = 1.5, nρ = 6, npos = 8)
        Wh2 = interp_weights_hybrid(13, 27; Lloc = 3, support_scale = 1.5, nρ = 6, npos = 8)
        @test Wh1 == Wh2                                   # 确定性（内部 RNG，不依赖全局状态）
        @test size(Wh1) == (532, 148)
        @test nnz(Wh1) / 532 <= 60                         # 稀疏（约 48/行）
        # 对"笛卡尔分量 degree<=3 的切向场"保持精确（混合权重的精确性约束对象）
        tnodes = LSP.get_t_nodes(6)
        pnodes = LSP.get_t_nodes(13)
        pC = LSP.nodes2Poles(tnodes)
        pF = LSP.nodes2Poles(pnodes)
        rC = hcat([p.r̂ for p in pC]...)
        rF = hcat([p.r̂ for p in pF]...)
        Random.seed!(2)
        for _ in 1:5
            J = randn(3)
            Fcart_c = J .- rC .* (vec(sum(rC .* J; dims = 1))')   # 切向，分量 degree<=2
            Fcart_f = J .- rF .* (vec(sum(rF .* J; dims = 1))')
            Fθc = [dot(pC[i].θhat, Fcart_c[:, i]) for i in 1:74]
            Fϕc = [dot(pC[i].ϕhat, Fcart_c[:, i]) for i in 1:74]
            Fθf = [dot(pF[i].θhat, Fcart_f[:, i]) for i in 1:266]
            Fϕf = [dot(pF[i].ϕhat, Fcart_f[:, i]) for i in 1:266]
            coarse = vcat(Fθc, Fϕc)
            exact = vcat(Fθf, Fϕf)
            est = Wh1 * coarse
            @test maximum(abs.(est .- exact)) / maximum(abs.(exact)) < 1e-6
        end
    end

    @testset "自旋加权球谐正交性" begin
        nd, wts = LSP.getlbSortedData(13)
        for s in (1, -1)
            S = spin_weighted_harmonics(nd, 6; s = s)
            G = S' * (wts .* S)
            @test maximum(abs.(G - Matrix{ComplexF64}(I, 48, 48))) < 1e-10
        end
    end

    @testset "修复后训练管线（hcat 复约束 + 共享几何）" begin
        # 小规模复现论文指标：k=8 时 εi 应远小于坏权重时代的 O(1)
        f(x) = truncation_kernel(x) - (13 + 1) / 2
        lo, hi = 1e-4, 20.0
        for _ in 1:60
            mid = (lo + hi) / 2
            f(mid) > 0 ? (hi = mid) : (lo = mid)
        end
        rel_l = (lo + hi) / 2
        arm = min(0.12, 0.5 * rel_l)
        rscale = max(rel_l / 2 - arm, 1e-4)
        tnodes = LSP.get_t_nodes(6)
        pnodes = LSP.get_t_nodes(13)
        nt, nf = size(tnodes, 2), size(pnodes, 2)
        rC = LSP.nodes2Poles(tnodes)
        rF = LSP.nodes2Poles(pnodes)
        Random.seed!(3)
        T = zeros(ComplexF64, nt, 2, 8, 10)
        P = zeros(ComplexF64, nf, 2, 8, 10)
        for ir in 1:10, iρ in 1:8
            rvec = DG.random_rvec() .* rscale
            geom = DG.random_source_geometry(rvec; arm_max = arm, off_max = 0.125 * arm)
            DG.evaluate_poles!(rC, view(T, :, :, iρ, ir), geom)
            DG.evaluate_poles!(rF, view(P, :, :, iρ, ir), geom)
        end
        T = reshape(T, nt * 2, :)
        P = reshape(P, nf * 2, :)
        flag = trunc(Int, 0.8 * size(T, 2))
        Aall = hcat(real(T[:, 1:flag]), imag(T[:, 1:flag]))
        Ball = hcat(real(P[:, 1:flag]), imag(P[:, 1:flag]))
        idxs = [partialsortperm([norm(tnodes[:, j] .- pnodes[:, i]) for j in 1:nt], 1:8) for i in 1:nf]
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
        Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]
        εi = maximum(abs.(W * Tte .- Pte); dims = 1) ./ maximum(abs.(Pte); dims = 1)
        @test mean(εi) < 1e-1
    end
end
