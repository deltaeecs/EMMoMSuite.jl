# benchmark/run_full_fast_solvers_benchmark.jl — M5 综合基准
#
# 稠密 vs MLFMA vs ACA vs MLACA：
# - 记录装配/建立/求解时间、MatVec 相对误差、GMRES 相对残差、解相对误差、
#   近场 nnz、低秩存储、压缩率。
# - 用例 1：N=792 球体（叶层 0.25λ，near_range=1，各方法同一近/远场划分）。
# - 用例 2：N=252 球体（叶层 0.125λ，4 层八叉树）→ ACA vs MLACA 多层压缩对比。
# - 预条件扫描：Identity / BlockJacobi / ILU / SPAI 在 ACA 算子上的迭代数与时间。
# 结果写 CSV 到 benchmark/results/（该目录已被 gitignore，属本地生成物）。

using EMMoMSuite
using LinearAlgebra
using Printf
using SparseArrays
using DelimitedFiles
using IterativeSolvers

const HEADER = ["case", "method", "N", "leaf_lambda", "t_setup_s", "t_solve_s",
                "mv_rel_err", "gmres_rel_res", "sol_rel_err", "nnz_near",
                "n_blocks", "stored_lr", "total_stored", "compression_ratio", "iterations"]

function storage_stats(op, N)
    nnz_near = nnz(op.Z_near)
    if hasproperty(op, :blocks)
        stored_lr = sum(size(b.U, 2) * (length(b.rows) + length(b.cols)) for b in op.blocks)
        nblocks = length(op.blocks)
    else
        stored_lr = 0
        nblocks = 0
    end
    total = nnz_near + stored_lr
    return nnz_near, nblocks, stored_lr, total, 1 - total / (N * N)
end

function solve_with_log(op, V; Pl = IdentityPreconditioner())
    I, hist = gmres(op, V; Pl = Pl, abstol = 1e-6, reltol = 1e-8,
                    maxiter = 400, restart = 50, log = true)
    res = norm(V - op * I) / norm(V)
    return I, res, hist.iters
end

function bench_case(; freq = 300e6, radius = 0.5, n_theta, n_phi, leaf_factor,
                    tol = 1e-4, near_range = 1, tag)
    λ = 299792458.0 / freq
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    efie = EFIE(freq)
    leaf = leaf_factor * λ

    t_direct = @elapsed Z = assemble_impedance_matrix(efie, basis)
    src = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(efie, src, basis)
    I_dir = Z \ V
    x = randn(ComplexF64, N)

    rows = Vector{Vector{Any}}()
    push!(rows, Any[tag, "dense", N, leaf_factor, t_direct, 0.0, 0.0, 0.0, 0.0,
                    N * N, 0, 0, N * N, 0.0, 0])

    builds = [
        ("mlfma", () -> MLFMAOperator(efie, basis, leaf, Val(:Lagrange2Step), near_range)),
        ("aca", () -> ACAOperator(efie, basis, leaf; tol = tol, near_range = near_range)),
        ("mlaca", () -> MLACAOperator(efie, basis, leaf; tol = tol, near_range = near_range)),
    ]
    for (name, build) in builds
        t_setup = @elapsed op = build()
        mv_err = norm(op * x - Z * x) / norm(Z * x)
        P = ILUPreconditioner(op)
        t_solve = @elapsed (I, res, iters) = solve_with_log(op, V; Pl = P)
        sol_err = norm(I - I_dir) / norm(I_dir)
        nnz_near, nblocks, stored_lr, total, ratio = storage_stats(op, N)
        push!(rows, Any[tag, name, N, leaf_factor, t_setup, t_solve, mv_err, res,
                        sol_err, nnz_near, nblocks, stored_lr, total, ratio, iters])
    end
    return rows
end

function bench_preconditioners(; freq = 300e6, radius = 0.5, n_theta = 8, n_phi = 12,
                               leaf_factor = 0.25, tol = 1e-4, near_range = 1)
    λ = 299792458.0 / freq
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    efie = EFIE(freq)
    leaf = leaf_factor * λ
    op = ACAOperator(efie, basis, leaf; tol = tol, near_range = near_range)
    src = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(efie, src, basis)

    rows = Vector{Vector{Any}}()
    for (pname, Pl) in [
        ("identity", IdentityPreconditioner()),
        ("blockjacobi", BlockJacobiPreconditioner(op)),
        ("ilu", ILUPreconditioner(op)),
        ("spai", SPAIPreconditioner(op)),
    ]
        t_solve = @elapsed (I, res, iters) = solve_with_log(op, V; Pl = Pl)
        push!(rows, Any["precond", pname, num_basis(basis), leaf_factor, 0.0, t_solve,
                        0.0, res, 0.0, 0, 0, 0, 0, 0.0, iters])
    end
    return rows
end

function main(; out = joinpath(@__DIR__, "results", "fast_solvers_benchmark_2026-08-11.csv"))
    println("=============================================================")
    println("  Fast Solvers Benchmark: dense / MLFMA / ACA / MLACA (EFIE)")
    println("=============================================================")
    rows = Vector{Vector{Any}}()
    append!(rows, bench_case(; n_theta = 12, n_phi = 24, leaf_factor = 0.25, tag = "case1"))
    append!(rows, bench_case(; n_theta = 8, n_phi = 12, leaf_factor = 0.125, tag = "case2"))
    append!(rows, bench_preconditioners())

    for r in rows
        println(@sprintf("%-10s %-11s N=%-5d leaf=%.3f setup=%.3fs solve=%.3fs mv=%.1e res=%.1e sol=%.1e nnz=%d blk=%d lr=%d ratio=%.1f%% it=%d",
                         r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9],
                         r[10], r[11], r[12], 100 * r[14], r[15]))
    end

    mkpath(dirname(out))
    writedlm(out, vcat(permutedims(HEADER), reduce(vcat, [permutedims(r) for r in rows])), ',')
    println("CSV written: ", out)
end

main()
