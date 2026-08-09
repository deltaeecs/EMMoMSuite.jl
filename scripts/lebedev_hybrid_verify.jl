# 混合权重（自适应支撑 + 笛卡尔标量 SH 精确性约束 + 数据拟合）验证：
#   1. EFIE 类 εi vs 训练式 k=8 与笛卡尔局部（L_loc / support_scale 扫描）
#   2. 确定性（同种子两次构造一致）
#   3. 真实 MLFMA（r=2λ）εq + 单次 MVM 时间
#   4. 高阶 65->131 构造时间
#
# 用法: julia --project=. scripts/lebedev_hybrid_verify.jl

using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.MLFMA
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

paper_metric(W, Tte, Pte) =
    vec(maximum(abs.(W * Tte .- Pte); dims = 1) ./ maximum(abs.(Pte); dims = 1))

function make_corrected(pk, pt; nρ = 40, npos = 120, seed = 20260808)
    rel_l = rel_l_for(pk)
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

function trained_weights(pk, pt, k; nρ = 40, npos = 120, seed = 20260808)
    tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
    pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    T, P, _, _ = make_corrected(pk, pt; nρ = nρ, npos = npos, seed = seed)
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
    for (pk, pt) in ((27, 53),)
        T, P, tnodes, pnodes = make_corrected(pk, pt)
        flag = trunc(Int, 0.8 * size(T, 2))
        Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]
        nt, nf = size(tnodes, 2), size(pnodes, 2)
        println("############ pk=$pk -> pt=$pt ############")

        Wt = trained_weights(pk, pt, 8)
        @printf("训练式 k=8:         εi=%.3e  nnz/行=%.1f\n", mean(paper_metric(Wt, Tte, Pte)), nnz(Wt) / (2nf))

        Wc = SH.interp_weights_cart_local(pk, pt; Lloc = 3, cap_rad = 0.6)
        @printf("笛卡尔局部 L_loc=3: εi=%.3e  nnz/行=%.1f\n", mean(paper_metric(Wc, Tte, Pte)), nnz(Wc) / (2nf))

        for Lloc in (2, 3, 4), ss in (1.5, 2.0)
            t0 = time()
            Wh = SH.interp_weights_hybrid(pk, pt; Lloc = Lloc, support_scale = ss, nρ = 25, npos = 60)
            εi = paper_metric(Wh, Tte, Pte)
            @printf("混合 L_loc=%d scale=%.1f: εi=%.3e  nnz/行=%.1f  构造 %.1fs\n",
                Lloc, ss, mean(εi), nnz(Wh) / (2nf), time() - t0)
        end

        # 确定性
        Wh1 = SH.interp_weights_hybrid(pk, pt; Lloc = 3, support_scale = 2.0, nρ = 25, npos = 60)
        Wh2 = SH.interp_weights_hybrid(pk, pt; Lloc = 3, support_scale = 2.0, nρ = 25, npos = 60)
        @printf("确定性(两次构造一致): %s\n", Wh1 == Wh2)
    end

    println("\n############ 高阶 65->131 混合权重 ############")
    pk2, pt2 = 65, 131
    t0 = time()
    Wh = SH.interp_weights_hybrid(pk2, pt2; Lloc = 3, support_scale = 2.0, nρ = 15, npos = 30)
    @printf("构造 %.1fs, nnz/行=%.1f\n", time() - t0, nnz(Wh) / (2 * 5810))

    println("\n############ 真实 MLFMA（PEC 球 r=2λ, N=1440）############")
    freq = 300e6
    mesh = generate_sphere_mesh(2.0, 16, 32)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    efie = EFIE(freq)
    Z_direct = assemble_impedance_matrix(efie, basis)
    Random.seed!(20260808)
    I = randn(ComplexF64, N)
    I ./= norm(I)

    eq_error(op) = begin
        y = op * I
        y_near = op.Z_near * I
        Zfr = (Z_direct - op.Z_near) * I
        norm(y - y_near - Zfr) / maximum(abs.(Zfr))
    end
    function time_mvm(op, nrep = 5)
        y = similar(I)
        mul!(y, op, I)
        t0 = time()
        for _ in 1:nrep
            mul!(y, op, I)
        end
        (time() - t0) / nrep
    end

    op = MLFMAOperator(efie, basis, 0.3, Val(:LbTrained1Step), 4)
    pairs = Tuple{Int,Int}[]
    for id in keys(op.octree.levels)
        lv = op.octree.levels[id]
        if isdefined(lv, :interpWθϕ) && hasfield(typeof(lv.interpWθϕ), :θϕCSC)
            nf = size(lv.interpWθϕ.θϕCSC, 1) ÷ 2
            nt = size(lv.interpWθϕ.θϕCSC, 2) ÷ 2
            push!(pairs, (LSP.n2pDict[nt], LSP.n2pDict[nf]))
        end
    end
    function swap!(Wdict)
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
            haskey(Wdict, (pk, pt)) || continue
            interp.θϕCSC = convert(SparseMatrixCSC{Float64,Int}, Wdict[(pk, pt)])
            interp.θϕCSCT = sparse(transpose(interp.θϕCSC))
        end
    end
    @printf("插值对: %s\n", pairs)
    Wt = Dict{Tuple{Int,Int},Any}()
    Wh = Dict{Tuple{Int,Int},Any}()
    for (pk, pt) in pairs
        Wt[(pk, pt)] = trained_weights(pk, pt, 8; nρ = 30, npos = 80)
        Wh[(pk, pt)] = SH.interp_weights_hybrid(pk, pt; Lloc = 3, support_scale = 2.0, nρ = 25, npos = 60)
    end
    swap!(Wt)
    @printf("训练式 k=8:  εq=%.3e  单次MVM=%.4fs\n", eq_error(op), time_mvm(op))
    swap!(Wh)
    @printf("混合 L_loc=3: εq=%.3e  单次MVM=%.4fs\n", eq_error(op), time_mvm(op))
end

main()
