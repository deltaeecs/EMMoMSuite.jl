# 同一真实 MLFMA 用例下，Lebedev 插值矩阵的精度-稀疏度-时间权衡：
#   a) :sh_auto（精确稠密）        b) :sh_local（L_loc=6 稀疏）
#   c) 修复后训练式（k=8/18，共享几何+复约束，论文方法修复版）
# 指标：εq（远场 MVM 相对误差）、插值矩阵内存、单次 MVM 时间、极点数。
#
# 用法: julia --project=. scripts/lebedev_mlfma_tradeoff.jl

using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.MLFMA
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: levelIntegralInfoCal, truncation_kernel
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG
import EMMoMSuite.FastAlgorithms.Lebedev.pinv2interpW as PW
import EMMoMSuite.FastAlgorithms.Lebedev.SHInterp as SH

"""truncation_kernel 单调递增，专用稳健二分求 rel_l"""
function rel_l_for(pk)
    f(x) = truncation_kernel(x) - (pk + 1) / 2
    lo, hi = 1e-4, 20.0
    for _ in 1:80
        mid = (lo + hi) / 2
        f(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    (lo + hi) / 2
end

eq_error(op, Z_direct, I) = begin
    y = op * I
    y_near = op.Z_near * I
    Z_far_ref = (Z_direct - op.Z_near) * I
    # 远场误差 = (MLFMA 总 - 近场) - (直接总 - 近场)
    norm(y - y_near - Z_far_ref) / maximum(abs.(Z_far_ref))
end

"""修复后训练权重（k 个最近点，hcat 全复约束）"""
function trained_weights(pk, pt, k; nρ = 40, npos = 120, seed = 20260808)
    Lb = (pk - 1) ÷ 2
    rel_l = rel_l_for(pk)
    tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
    pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
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
    T = zeros(ComplexF64, nt, 2, nρ, npos)
    P = zeros(ComplexF64, nf, 2, nρ, npos)
    for ir in axes(rbmrps, 2), iρ in axes(ρhats, 2)
        geom = DG.random_source_geometry(rbmrps[:, ir]; arm_max = 0.5 * rel_l, off_max = 0.125 * rel_l)
        DG.evaluate_poles!(rC, view(T, :, :, iρ, ir), geom)
        DG.evaluate_poles!(rF, view(P, :, :, iρ, ir), geom)
    end
    T = reshape(T, nt * 2, :)
    P = reshape(P, nf * 2, :)
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

"""原版坏权重复现：固定臂长 0.12λ + 两层各自随机几何 + vcat(real,imag) 只实部约束"""
function broken_weights(pk, pt, k; nρ = 40, npos = 120, seed = 20260808)
    Lb = (pk - 1) ÷ 2
    rel_l = rel_l_for(pk)
    tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
    pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
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
    T = zeros(ComplexF64, nt, 2, nρ, npos)
    P = zeros(ComplexF64, nf, 2, nρ, npos)
    for ir in axes(rbmrps, 2), iρ in axes(ρhats, 2)
        # 两层分别调用 -> 各自随机几何（原 bug）；固定臂长 0.12λ（不按层缩放）
        DG.generate_dataset_on_poles(rC, view(T, :, :, iρ, ir); rvec = rbmrps[:, ir], arm_max = 0.12, off_max = 0.03)
        DG.generate_dataset_on_poles(rF, view(P, :, :, iρ, ir); rvec = rbmrps[:, ir], arm_max = 0.12, off_max = 0.03)
    end
    T = reshape(T, nt * 2, :)
    P = reshape(P, nf * 2, :)
    flag = trunc(Int, 0.8 * size(T, 2))
    xx2D = vcat(real(T[:, 1:flag]), imag(T[:, 1:flag]))   # 原 vcat 布局
    yy2D = vcat(real(P[:, 1:flag]), imag(P[:, 1:flag]))
    w = PW.interpWeightsInitial(tnodes, pnodes; nInterp = k)
    PW.pinv2W!(w, k, xx2D, yy2D)
    return w
end

function swap_interp!(op, pkpt_pairs, new_weights)
    # new_weights: Dict{(pk,pt) => W}，W 为 2nf x 2nt 块矩阵
    for id in keys(op.octree.levels)
        lv = op.octree.levels[id]
        if !isdefined(lv, :interpWθϕ) || !hasfield(typeof(lv.interpWθϕ), :θϕCSC)
            continue
        end
        interp = lv.interpWθϕ
        nf = size(interp.θϕCSC, 1) ÷ 2
        nt = size(interp.θϕCSC, 2) ÷ 2
        pk = LSP.n2pDict[nt]
        pt = LSP.n2pDict[nf]
        if haskey(new_weights, (pk, pt))
            W = new_weights[(pk, pt)]
            @assert size(W) == size(interp.θϕCSC)
            interp.θϕCSC = convert(SparseMatrixCSC{Float64,Int}, W)
            interp.θϕCSCT = sparse(transpose(interp.θϕCSC))
        end
    end
    return op
end

function run_tradeoff(radius, ntheta, nphi, leaf; near_range = 4, nρ = 40, npos = 120)
    freq = 300e6
    mesh = generate_sphere_mesh(radius, ntheta, nphi)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    efie = EFIE(freq)
    @printf("\n########## PEC 球 r=%.1fλ, N=%d ##########\n", radius, N)

    Z_direct = assemble_impedance_matrix(efie, basis)
    Random.seed!(20260808)
    I = randn(ComplexF64, N); I ./= norm(I)

    op = MLFMAOperator(efie, basis, leaf, Val(:LbTrained1Step), near_range)
    pairs = Tuple{Int,Int}[]
    nfar3 = 0
    for id in keys(op.octree.levels)
        lv = op.octree.levels[id]
        nfar3 += sum(length(c.farneighbors) for c in lv.cubes)
        if isdefined(lv, :interpWθϕ) && hasfield(typeof(lv.interpWθϕ), :θϕCSC)
            nf = size(lv.interpWθϕ.θϕCSC, 1) ÷ 2
            nt = size(lv.interpWθϕ.θϕCSC, 2) ÷ 2
            push!(pairs, (LSP.n2pDict[nt], LSP.n2pDict[nf]))
        end
    end
    @printf("插值对: %s；非叶层远邻居总数=%d\n", pairs, nfar3)

    y = similar(I)
    nnz_tot(op) = sum(nnz(lv.interpWθϕ.θϕCSC) for lv in values(op.octree.levels) if
        isdefined(lv, :interpWθϕ) && hasfield(typeof(lv.interpWθϕ), :θϕCSC))
    time_mvm(op, I) = begin
        mul!(y, op, I)
        t0 = time()
        for _ in 1:3
            mul!(y, op, I)
        end
        (time() - t0) / 3
    end

    # a) 稠密精确（原样）
    @printf("a) 稠密精确: εq=%.3e  插值内存=%.1f MB  单次MVM=%.4fs\n",
        eq_error(op, Z_direct, I), nnz_tot(op) * 8 / 1e6, time_mvm(op, I))

    # b) :sh_local 稀疏
    local_weights = Dict{Tuple{Int,Int},Any}()
    for (pk, pt) in pairs
        local_weights[(pk, pt)] = SH.vectorize(SH.interp_weights_local(pk, pt; Lloc = min(6, (pk - 1) ÷ 2), cap_rad = 1.0))
    end
    swap_interp!(op, pairs, local_weights)
    @printf("b) :sh_local L_loc=6: εq=%.3e  插值内存=%.1f MB  单次MVM=%.4fs\n",
        eq_error(op, Z_direct, I), nnz_tot(op) * 8 / 1e6, time_mvm(op, I))

    # c) 修复后训练式 k=8 / k=18
    for k in (8, 18)
        tr_weights = Dict{Tuple{Int,Int},Any}()
        for (pk, pt) in pairs
            tr_weights[(pk, pt)] = trained_weights(pk, pt, k; nρ = nρ, npos = npos)
        end
        swap_interp!(op, pairs, tr_weights)
        @printf("c) 训练式 k=%2d: εq=%.3e  插值内存=%.1f MB  单次MVM=%.4fs\n",
            k, eq_error(op, Z_direct, I), nnz_tot(op) * 8 / 1e6, time_mvm(op, I))
    end

    # d) 对照组：原版坏权重（k=9）
    bad_weights = Dict{Tuple{Int,Int},Any}()
    for (pk, pt) in pairs
        bad_weights[(pk, pt)] = broken_weights(pk, pt, 9; nρ = nρ, npos = npos)
    end
    swap_interp!(op, pairs, bad_weights)
    @printf("d) 原版坏权重(k=9): εq=%.3e  插值内存=%.1f MB  单次MVM=%.4fs\n",
        eq_error(op, Z_direct, I), nnz_tot(op) * 8 / 1e6, time_mvm(op, I))

    # e) 对照组：零矩阵
    zero_weights = Dict{Tuple{Int,Int},Any}()
    for (pk, pt) in pairs
        nf = LSP.p2nDict[pt]
        nt = LSP.p2nDict[pk]
        zero_weights[(pk, pt)] = spzeros(2nf, 2nt)
    end
    swap_interp!(op, pairs, zero_weights)
    @printf("e) 零插值矩阵: εq=%.3e\n", eq_error(op, Z_direct, I))
end

function main()
    run_tradeoff(1.0, 12, 24, 0.3)
    run_tradeoff(2.0, 16, 32, 0.3; nρ = 25, npos = 60)
end

main()
