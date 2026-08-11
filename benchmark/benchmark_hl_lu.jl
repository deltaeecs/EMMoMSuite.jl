# benchmark/benchmark_hl_lu.jl - H2 H-LU 系统化用例基准
#
# 各用例（EFIE / CFIE / PMCHW / 低频）：H-LU vs BlockLU vs 稠密 LU
# - 分解时间、多 RHS 求解时间、压缩算子残差（‖Zc·X−B‖/‖B‖）
# - H 树节点数、低秩块数、存储估计 vs 稠密 N²
# 结果写 CSV 到 benchmark/results/（本地手动运行）。

using EMMoMSuite
using LinearAlgebra
using Printf
using SparseArrays
using DelimitedFiles
using EMMoMSuite.FastAlgorithms.ACAOperatorModule: ACAOperator
using EMMoMSuite.FastAlgorithms.MLACAOperatorModule: MLACAOperator
using EMMoMSuite.FastAlgorithms.BlockLUModule: block_lu, block_lu_solve
using EMMoMSuite.FastAlgorithms.HMatrixModule: hmatrix_from_mlaca, h_lu!, h_lu_solve, materialize, hmatrix_count

const HEADER = ["case", "method", "N", "t_factor_s", "t_solve_s", "residual",
                "n_nodes", "n_lowrank", "stored_entries", "compression_ratio"]

function hmatrix_stats(H)
    n_nodes = hmatrix_count(H)
    n_low = Ref(0)
    stored = Ref(0)
    function walk(n)
        if n.kind == :lowrank
            n_low[] += 1
            stored[] += size(n.U, 2) * (length(n.rows) + length(n.cols))
        elseif n.kind == :dense
            stored[] += length(n.rows) * length(n.cols)
        else
            for c in n.children
                walk(c)
            end
        end
    end
    walk(H)
    return n_nodes, n_low[], stored[]
end

function bench_case(; tag, op, N, nrhs = 4)
    rows = Vector{Vector{Any}}()
    X_true = randn(ComplexF64, N, nrhs)
    B = reduce(hcat, [op * X_true[:, j] for j in 1:nrhs])

    # 稠密 LU
    t_lu = @elapsed Xd = Matrix(op.Z_near) \ B  # 占位（Zc 稠密参照）
    # 用压缩算子的稠密参照
    H0 = hmatrix_from_mlaca(op)
    Zc = materialize(H0)
    t_dense = @elapsed Xd = Zc \ B
    push!(rows, Any[tag, "dense_lu", N, t_dense, 0.0, 0.0, 0, 0, N * N, 0.0])

    # BlockLU（叶层块 LU，精确）
    t_f = @elapsed F = block_lu(op)
    t_s = @elapsed Xb = block_lu_solve(F, B)
    res = norm(B - reduce(hcat, [op * Xb[:, j] for j in 1:nrhs])) / norm(B)
    stored_b = sum(length, F.blocks)^2  # 近似：对角稠密
    push!(rows, Any[tag, "blocklu", N, t_f, t_s, res, 0, 0, stored_b, 1 - stored_b / (N * N)])

    # H-LU（分层，精确）
    for (rname, rc) in (("hlu", false), ("hlu_rc", true))
        t_h = @elapsed begin
            H = hmatrix_from_mlaca(op)
            h_lu!(H; tol = 1e-4, recompress = rc)
        end
        t_hs = @elapsed Xh = h_lu_solve(H, B)
        res_h = norm(B - reduce(hcat, [op * Xh[:, j] for j in 1:nrhs])) / norm(B)
        n_nodes, n_low, stored_h = hmatrix_stats(H)
        push!(rows, Any[tag, rname, N, t_h, t_hs, res_h, n_nodes, n_low, stored_h, 1 - stored_h / (N * N)])
    end
    return rows
end

function main(; out = joinpath(@__DIR__, "results", "hl_lu_benchmark_2026-08-11.csv"))
    println("=============================================================")
    println("  H-LU systematic use-case benchmark")
    println("=============================================================")
    rows = Vector{Vector{Any}}()
    lambda = 299792458.0 / 300e6

    mesh = generate_sphere_mesh(0.5, 6, 10)
    basis = RWGBasis(mesh)
    N = num_basis(basis)

    println("\n[case] efie_mlaca (300 MHz, multi-level tree)")
    efie = EFIE(300e6)
    op = MLACAOperator(efie, basis, 0.125 * lambda; tol = 1e-4, near_range = 1)
    append!(rows, bench_case(; tag = "efie_mlaca", op, N))

    println("\n[case] cfie_aca (300 MHz, non-symmetric)")
    cfie = CFIE(300e6)
    op2 = ACAOperator(cfie, basis, 0.25 * lambda; tol = 1e-4, near_range = 1, symmetric = false)
    append!(rows, bench_case(; tag = "cfie_aca", op = op2, N))

    println("\n[case] pmchw_aca (300 MHz, 2N system)")
    pmchw = PMCHW(300e6, 4.0)
    op3 = ACAOperator(pmchw, basis, 0.25 * lambda; tol = 1e-4, near_range = 1)
    append!(rows, bench_case(; tag = "pmchw_aca", op = op3, N = 2N))

    println("\n[case] lowfreq_aca (30 MHz, leaf 0.03 lambda)")
    lambda_lf = 299792458.0 / 30e6
    lf = EFIE(30e6)
    op4 = ACAOperator(lf, basis, 0.03 * lambda_lf; tol = 1e-4, near_range = 1)
    append!(rows, bench_case(; tag = "lowfreq_aca", op = op4, N))

    println("\n[case] efie_aca_792 (300 MHz, larger blocks)")
    mesh2 = generate_sphere_mesh(0.5, 12, 24)
    basis2 = RWGBasis(mesh2)
    N2 = num_basis(basis2)
    op5 = ACAOperator(efie, basis2, 0.25 * lambda; tol = 1e-4, near_range = 1)
    append!(rows, bench_case(; tag = "efie_aca_792", op = op5, N = N2))

    for r in rows
        println(@sprintf("%-12s %-9s N=%-4d t_f=%.2fs t_s=%.3fs res=%.1e nodes=%d lr=%d stored=%d ratio=%.1f%%",
                         r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], 100 * r[10]))
    end

    mkpath(dirname(out))
    writedlm(out, vcat(permutedims(HEADER), reduce(vcat, [permutedims(r) for r in rows])), ',')
    println("\nCSV written: ", out)
end

main()
