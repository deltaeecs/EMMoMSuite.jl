# 性能对比：不同插值权重矩阵在真实 MLFMA 中的速度与精度权衡
#   - 插值矩阵微基准（27->53 与 65->131 对的单次 matvec 时间 / FLOPs）
#   - 真实 MLFMA（PEC 球 r=2λ, N=1440）单次 MVM 总时间 + εq
# 对比：GL 两段 Lagrange / 训练式 k=8 / 笛卡尔局部 L_loc=3,6 / 笛卡尔稠密
#
# 用法: julia --project=. scripts/lebedev_speed_compare.jl

using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.MLFMA
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: levelIntegralInfoCal, truncation_kernel
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG
import EMMoMSuite.FastAlgorithms.Lebedev.SHInterp as SH

rel_l_for(pk) = begin
    f(x) = truncation_kernel(x) - (pk + 1) / 2
    lo, hi = 1e-4, 20.0
    for _ in 1:80
        mid = (lo + hi) / 2
        f(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    (lo + hi) / 2
end

function trained_weights(pk, pt, k; nρ = 30, npos = 80, seed = 20260808)
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

function bench_matvec(name, W, ntrials = 50)
    n2 = size(W, 2)
    x = randn(ComplexF64, n2)
    y = similar(x, size(W, 1))
    mul!(y, W, x)
    t0 = time()
    for _ in 1:ntrials
        mul!(y, W, x)
    end
    t = (time() - t0) / ntrials
    nz = W isa SparseMatrixCSC ? nnz(W) : length(W)
    flops = 2 * nz
    @printf("  %-28s nnz=%9d (%.1f/行)  单次matvec %.3e s  有效GFLOP/s=%.1f\n",
        name, nz, nz / size(W, 1), t, flops / max(t, 1e-12) / 1e9)
end

function main()
    println("========== 插值矩阵微基准 ==========")
    for (pk, pt) in ((27, 53), (65, 131))
        println("--- pair $pk -> $pt ---")
        tnodes = LSP.get_t_nodes((pk - 1) ÷ 2)
        pnodes = LSP.get_t_nodes((pt - 1) ÷ 2)
        nt, nf = size(tnodes, 2), size(pnodes, 2)
        # 训练式（65->131 跳过：数据不可行）
        if pt <= 53
            Wt = trained_weights(pk, pt, 8)
            bench_matvec("训练式 k=8", Wt)
        end
        for (Lloc, cap) in ((3, 0.6), (6, 1.0))
            Wc = SH.interp_weights_cart_local_orbit(pk, pt; Lloc = Lloc, cap_rad = cap)
            bench_matvec("笛卡尔局部 L_loc=$Lloc", Wc)
        end
        Wd = SH.interp_weights_cart(pk, pt)
        bench_matvec("笛卡尔稠密", Wd)
    end

    println("\n========== 真实 MLFMA（PEC 球 r=2λ, N=1440, leaf=0.3, near_range=4）==========")
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
    @printf("插值对: %s\n", pairs)

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
            W = Wdict[(pk, pt)]
            interp.θϕCSC = convert(SparseMatrixCSC{Float64,Int}, W)
            interp.θϕCSCT = sparse(transpose(interp.θϕCSC))
        end
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

    # a) 稠密精确（原样 :sh_auto）
    @printf("a) 稠密精确:          εq=%.3e  单次MVM=%.4fs\n", eq_error(op), time_mvm(op))

    # b) 笛卡尔局部 L_loc=3 / 6
    for (Lloc, cap) in ((3, 0.6), (6, 1.0))
        Wdict = Dict{Tuple{Int,Int},Any}()
        for (pk, pt) in pairs
            Wdict[(pk, pt)] = SH.interp_weights_cart_local_orbit(pk, pt; Lloc = Lloc, cap_rad = cap)
        end
        swap!(Wdict)
        @printf("b) 笛卡尔局部 L_loc=%d: εq=%.3e  单次MVM=%.4fs\n", Lloc, eq_error(op), time_mvm(op))
    end

    # c) 训练式 k=8
    Wdict = Dict{Tuple{Int,Int},Any}()
    for (pk, pt) in pairs
        Wdict[(pk, pt)] = trained_weights(pk, pt, 8)
    end
    swap!(Wdict)
    @printf("c) 训练式 k=8:        εq=%.3e  单次MVM=%.4fs\n", eq_error(op), time_mvm(op))

    # d) GL 两段 Lagrange（对照）
    op_gl = MLFMAOperator(efie, basis, 0.3, Val(:Lagrange2Step), 4)
    @printf("d) GL 两段 Lagrange:  εq=%.3e  单次MVM=%.4fs\n", eq_error(op_gl), time_mvm(op_gl))
end

main()
