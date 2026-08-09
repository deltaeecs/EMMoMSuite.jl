# 真实 MLFMA 基准：GL 两段 Lagrange vs Lebedev 一步插值（修复后球谐权重）
# 对应论文 Fig.3 / Table I-II：
#   εq = |Z_MLFMA_far I − Z_far I| / max|Z_far I|   （Eq.8，Z_far 由满阵减近场得到）
#   插值矩阵内存/nnz、极点数（辐射/接收模式存储）、单次 MVM 时间
#
# 用法: julia --project=. scripts/lebedev_mlfma_benchmark.jl

using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.MLFMA
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP

"""εq（论文 Eq.8）：相对远场 MVM 误差，I 为单位随机电流向量"""
function eq_error(op, Z_direct, I)
    y = op * I
    y_near = op.Z_near * I
    Z_far_ref = (Z_direct - op.Z_near) * I
    return norm(y - y_near - Z_far_ref) / maximum(abs.(Z_far_ref))
end

"""插值矩阵内存 + 极点数统计（对应论文 Table I 的 Quadrature & Interpolation 列）"""
function interp_stats(op)
    total_nnz = 0
    total_bytes = 0.0
    total_poles = 0
    rows = String[]
    for id in sort(collect(keys(op.octree.levels)))
        lv = op.octree.levels[id]
        npol = length(lv.poles.Wθϕs)
        total_poles += npol
        if isdefined(lv, :interpWθϕ) && lv.interpWθϕ !== nothing
            interp = lv.interpWθϕ
            if hasfield(typeof(interp), :θϕCSC)
                nnzcount = nnz(interp.θϕCSC)
                sz = string(size(interp.θϕCSC, 1), "x", size(interp.θϕCSC, 2))
                kind = "Lebedev一步"
            else
                nnzcount = nnz(interp.θCSC) + nnz(interp.ϕCSC)
                sz = string(size(interp.θCSC, 1), "x", size(interp.θCSC, 2), " (θ+ϕ)")
                kind = "Lagrange两段"
            end
            total_nnz += nnzcount
            total_bytes += nnzcount * 8
            push!(rows, @sprintf("  level %d: L=%3d 极点数=%5d 插值[%s] %s nnz=%d (%.2f MB)",
                id, lv.L, npol, kind, sz, nnzcount, nnzcount * 8 / 1e6))
        else
            push!(rows, @sprintf("  level %d: L=%3d 极点数=%5d（无插值矩阵）", id, lv.L, npol))
        end
    end
    return total_nnz, total_bytes, total_poles, rows
end

function run_case(radius, ntheta, nphi, leaf; trials = 3)
    freq = 300e6
    λ = 1.0
    println("\n========== PEC 球 r=$(radius)λ, 网格 ($ntheta,$nphi), leaf=$(leaf)λ ==========")
    mesh = generate_sphere_mesh(radius, ntheta, nphi)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    efie = EFIE(freq)
    @printf("未知数 N=%d\n", N)

    # 满阵参考（远场 = 满阵 − 近场）
    t0 = time()
    Z_direct = assemble_impedance_matrix(efie, basis)
    @printf("满阵 Z 装配: %.1fs\n", time() - t0)

    Random.seed!(20260808)
    I = randn(ComplexF64, N)
    I ./= norm(I)

    for (name, method) in (("GL 两段Lagrange", Val(:Lagrange2Step)), ("Lebedev一步(:sh_auto)", Val(:LbTrained1Step)))
        t0 = time()
        op = MLFMAOperator(efie, basis, leaf, method, 4)
        tbuild = time() - t0
        nlevels = op.octree.nLevels

        # MVM 计时（多次取平均）
        y = similar(I)
        mul!(y, op, I)  # warmup
        tmv0 = time()
        for _ in 1:trials
            mul!(y, op, I)
        end
        tmvm = (time() - tmv0) / trials

        εq = eq_error(op, Z_direct, I)
        tot_nnz, tot_bytes, tot_poles, rows = interp_stats(op)

        @printf("\n--- %s: 构造 %.1fs, levels=%d, 单次MVM %.4fs ---\n", name, tbuild, nlevels, tmvm)
        @printf("  εq(远场MVM相对误差) = %.3e\n", εq)
        @printf("  总插值矩阵: nnz=%d, %.2f MB；极点数合计=%d\n", tot_nnz, tot_bytes / 1e6, tot_poles)
        foreach(println, rows)
    end
end

run_case(1.0, 12, 24, 0.3)
run_case(2.0, 16, 32, 0.3)
