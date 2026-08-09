# 4λ 球真实 MLFMA 规模行为确认（新默认 :sh_auto：小对稠密精确 + 顶层对混合权重）
#   - 各层极点数 / 插值矩阵内存与 nnz
#   - 构造时间 / 单次 MVM / εq（vs 直接矩阵）
#   - GL 两段 Lagrange 对照
#
# 用法: julia --project=. scripts/lebedev_4lambda_benchmark.jl

using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.MLFMA
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.SHInterp as SH
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncation_kernel

function main()
    freq = 300e6
    radius = 4.0
    mesh = generate_sphere_mesh(radius, 24, 48)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    efie = EFIE(freq)
    @printf("PEC 球 r=4λ, N=%d\n", N)

    t0 = time()
    Z_direct = assemble_impedance_matrix(efie, basis)
    @printf("满阵 Z 装配: %.1fs (%.1f MB)\n", time() - t0, 16 * N * N / 2 / 1e6)

    Random.seed!(20260808)
    I = randn(ComplexF64, N)
    I ./= norm(I)

    eq_error(op) = begin
        y = op * I
        y_near = op.Z_near * I
        Zfr = (Z_direct - op.Z_near) * I
        norm(y - y_near - Zfr) / maximum(abs.(Zfr))
    end
    function time_mvm(op, nrep = 3)
        y = similar(I)
        mul!(y, op, I)
        t0 = time()
        for _ in 1:nrep
            mul!(y, op, I)
        end
        (time() - t0) / nrep
    end

    for (name, method) in (("GL 两段Lagrange", Val(:Lagrange2Step)), ("Lebedev :sh_auto(混合)", Val(:LbTrained1Step)))
        t0 = time()
        op = MLFMAOperator(efie, basis, 0.3, method, 4)
        tbuild = time() - t0
        @printf("\n--- %s: 构造 %.1fs, levels=%d ---\n", name, tbuild, op.octree.nLevels)
        @printf("  εq=%.3e  单次MVM=%.4fs\n", eq_error(op), time_mvm(op))
        tot_nnz = 0
        tot_poles = 0
        for id in sort(collect(keys(op.octree.levels)))
            lv = op.octree.levels[id]
            npol = length(lv.poles.Wθϕs)
            tot_poles += npol
            if isdefined(lv, :interpWθϕ) && lv.interpWθϕ !== nothing
                interp = lv.interpWθϕ
                if hasfield(typeof(interp), :θϕCSC)
                    nnzcount = nnz(interp.θϕCSC)
                    tot_nnz += nnzcount
                    @printf("  level %2d: L=%3d 极点数=%5d 插值 %s nnz=%d (%.2f MB)\n",
                        id, lv.L, npol, string(size(interp.θϕCSC, 1), "x", size(interp.θϕCSC, 2)),
                        nnzcount, nnzcount * 8 / 1e6)
                else
                    nnzcount = nnz(interp.θCSC) + nnz(interp.ϕCSC)
                    tot_nnz += nnzcount
                    @printf("  level %2d: L=%3d 极点数=%5d 插值 Lagrange nnz=%d (%.2f MB)\n",
                        id, lv.L, npol, nnzcount, nnzcount * 8 / 1e6)
                end
            else
                @printf("  level %2d: L=%3d 极点数=%5d（无插值矩阵）\n", id, lv.L, npol)
            end
        end
        @printf("  合计: 极点数=%d, 插值矩阵 %.2f MB\n", tot_poles, tot_nnz * 8 / 1e6)
    end

    println("\n--- 顶层对 (65->101) 权重层对比 ---")
    T, P, tnodes, pnodes = lebedev_data(65, 101)
    flag = trunc(Int, 0.8 * size(T, 2))
    Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]
    metric(W) = vec(maximum(abs.(W * Tte .- Pte); dims = 1) ./ maximum(abs.(Pte); dims = 1))
    t0 = time()
    Wh = SH.interp_weights_hybrid(65, 101; Lloc = 3, support_scale = 1.5, nρ = 15, npos = 30)
    @printf("混合 L_loc=3: εi=%.3e  nnz/行=%.1f  构造 %.1fs\n", mean(metric(Wh)), nnz(Wh) / (2 * 3470), time() - t0)
    Wd = SH.vectorize(SH.interp_weights_exact(65, 101))
    @printf("稠密精确:     εi=%.3e  矩阵 %.1f MB\n", mean(metric(Wd)), sizeof(Wd) / 1e6)
end

"""顶层对 EFIE 数据（修正几何，小样本）"""
function lebedev_data(pk, pt; nρ = 15, npos = 30, seed = 20260808)
    f(x) = truncation_kernel(x) - (pk + 1) / 2
    lo, hi = 1e-4, 20.0
    for _ in 1:80
        mid = (lo + hi) / 2
        f(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    rel_l = (lo + hi) / 2
    arm = min(0.12, 0.5 * rel_l)
    rscale = max(rel_l / 2 - arm, 1e-4)
    tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
    pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    Random.seed!(seed)
    rC = LSP.nodes2Poles(tnodes)
    rF = LSP.nodes2Poles(pnodes)
    T = zeros(ComplexF64, nt, 2, nρ, npos)
    P = zeros(ComplexF64, nf, 2, nρ, npos)
    for ir in 1:npos, iρ in 1:nρ
        rvec = DG.random_rvec() .* rscale
        geom = DG.random_source_geometry(rvec; arm_max = arm, off_max = 0.125 * arm)
        DG.evaluate_poles!(rC, view(T, :, :, iρ, ir), geom)
        DG.evaluate_poles!(rF, view(P, :, :, iρ, ir), geom)
    end
    return reshape(T, nt * 2, :), reshape(P, nf * 2, :), tnodes, pnodes
end

main()
