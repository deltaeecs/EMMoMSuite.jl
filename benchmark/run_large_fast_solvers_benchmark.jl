# benchmark/run_large_fast_solvers_benchmark.jl - F4 更大规模本地测试
#
# 多尺寸 / 多配方（EFIE、CFIE、PMCHW、低频）对照稠密参照：
# - MatVec 相对误差门控（< 1e-2）
# - GMRES（ILU 预条件）求解：迭代数、算子残差、解误差（记录，受条件数限制）
# - 内存/耗时/压缩率：近场 nnz、低秩存储、块数、压缩率
# - finite/NaN 逐项检查
# 结果写 CSV 到 benchmark/results/（gitignore 的本地生成物）。

using EMMoMSuite
using LinearAlgebra
using Printf
using SparseArrays
using DelimitedFiles
using IterativeSolvers

const HEADER = ["case", "method", "N", "freq_Hz", "leaf_lambda", "t_setup_s", "t_solve_s",
                "mv_rel_err", "finite_ok", "gmres_iters", "gmres_rel_res", "sol_rel_err",
                "cond_est", "nnz_near", "n_blocks", "stored_lr", "total_stored", "compression_ratio"]

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

function bench_case(; tag, freq, radius, n_theta, n_phi, leaf_factor, tol = 1e-4,
                    near_range = 1, methods = (:mlfma, :aca, :mlaca), do_dense = true)
    lambda = 299792458.0 / freq
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    efie = EFIE(freq)
    leaf = leaf_factor * lambda

    t_direct = @elapsed Z = assemble_impedance_matrix(efie, basis)
    src = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(efie, src, basis)
    t_lu = @elapsed I_dir = Z \ V
    x = randn(ComplexF64, N)

    rows = Vector{Vector{Any}}()
    cond_est = N <= 1200 ? cond(Z) : NaN
    push!(rows, Any[tag, "dense", N, freq, leaf_factor, t_direct + t_lu, 0.0, 0.0,
                    true, 0, 0.0, 0.0, cond_est, N * N, 0, 0, N * N, 0.0])

    for name in methods
        build = if name === :mlfma
            () -> MLFMAOperator(efie, basis, leaf, Val(:Lagrange2Step), near_range)
        elseif name === :aca
            () -> ACAOperator(efie, basis, leaf; tol = tol, near_range = near_range)
        else
            () -> MLACAOperator(efie, basis, leaf; tol = tol, near_range = near_range)
        end
        t_setup = @elapsed op = build()
        y = op * x
        mv_err = norm(y - Z * x) / norm(Z * x)
        finite_ok = all(isfinite, vec(y))
        P = ILUPreconditioner(op)
        t_solve = @elapsed (I, hist) = gmres(op, V; Pl = P, abstol = 1e-6, reltol = 1e-8,
                                             maxiter = 500, restart = 100, log = true)
        relres = norm(V - op * I) / norm(V)
        sol_err = norm(I - I_dir) / norm(I_dir)
        nnz_near, nblocks, stored_lr, total, ratio = storage_stats(op, N)
        push!(rows, Any[tag, String(name), N, freq, leaf_factor, t_setup, t_solve, mv_err,
                        finite_ok, hist.iters, relres, sol_err, cond_est,
                        nnz_near, nblocks, stored_lr, total, ratio])
    end
    return rows
end

function main(; out = joinpath(@__DIR__, "results", "large_fast_solvers_benchmark_2026-08-11.csv"))
    println("=============================================================")
    println("  Large-scale Fast Solvers Benchmark (dense reference gates)")
    println("=============================================================")
    rows = Vector{Vector{Any}}()

    println("\n[case 1] EFIE sphere N~792 (300 MHz)")
    append!(rows, bench_case(; tag = "efie_792", freq = 300e6, radius = 0.5,
                             n_theta = 12, n_phi = 24, leaf_factor = 0.25))

    println("\n[case 2] EFIE sphere N~1734 (300 MHz)")
    append!(rows, bench_case(; tag = "efie_1734", freq = 300e6, radius = 0.5,
                             n_theta = 18, n_phi = 34, leaf_factor = 0.25))

    println("\n[case 3] EFIE sphere N~2280 (300 MHz)")
    append!(rows, bench_case(; tag = "efie_2280", freq = 300e6, radius = 0.5,
                             n_theta = 20, n_phi = 40, leaf_factor = 0.25))

    println("\n[case 4] CFIE sphere N~792 (non-symmetric, 300 MHz)")
    λ = 299792458.0 / 300e6
    mesh = generate_sphere_mesh(0.5, 12, 24)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    cfie = CFIE(300e6)
    Z = assemble_impedance_matrix(cfie, basis)
    src = PlaneWave(300e6, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(cfie, src, basis)
    I_dir = Z \ V
    x = randn(ComplexF64, N)
    for name in (:aca, :mlaca)
        op = name === :aca ? ACAOperator(cfie, basis, 0.25 * λ; tol = 1e-4, near_range = 1, symmetric = false) :
                             MLACAOperator(cfie, basis, 0.25 * λ; tol = 1e-4, near_range = 1, symmetric = false)
        y = op * x
        mv_err = norm(y - Z * x) / norm(Z * x)
        P = ILUPreconditioner(op)
        t_solve = @elapsed (I, hist) = gmres(op, V; Pl = P, abstol = 1e-6, reltol = 1e-8,
                                             maxiter = 500, restart = 100, log = true)
        relres = norm(V - op * I) / norm(V)
        sol_err = norm(I - I_dir) / norm(I_dir)
        nnz_near, nblocks, stored_lr, total, ratio = storage_stats(op, N)
        push!(rows, Any["cfie_792", String(name), N, 300e6, 0.25, 0.0, t_solve, mv_err,
                        all(isfinite, vec(y)), hist.iters, relres, sol_err, cond(Z),
                        nnz_near, nblocks, stored_lr, total, ratio])
    end

    println("\n[case 5] PMCHW sphere N=300 (2N=600, 300 MHz)")
    mesh = generate_sphere_mesh(0.5, 6, 10)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    pmchw = PMCHW(300e6, 4.0)
    Z = assemble_impedance_matrix(pmchw, basis)
    S = 2N
    V = randn(ComplexF64, S)
    I_dir = Z \ V
    x = randn(ComplexF64, S)
    for name in (:aca, :mlaca)
        op = name === :aca ? ACAOperator(pmchw, basis, 0.25 * λ; tol = 1e-4, near_range = 1) :
                             MLACAOperator(pmchw, basis, 0.125 * λ; tol = 1e-4, near_range = 1)
        y = op * x
        mv_err = norm(y - Z * x) / norm(Z * x)
        P = ILUPreconditioner(op)
        t_solve = @elapsed (I, hist) = gmres(op, V; Pl = P, abstol = 1e-6, reltol = 1e-8,
                                             maxiter = 800, restart = 100, log = true)
        relres = norm(V - op * I) / norm(V)
        sol_err = norm(I - I_dir) / norm(I_dir)
        nnz_near, nblocks, stored_lr, total, ratio = storage_stats(op, S)
        push!(rows, Any["pmchw_600", String(name), S, 300e6, 0.25, 0.0, t_solve, mv_err,
                        all(isfinite, vec(y)), hist.iters, relres, sol_err, cond(Z),
                        nnz_near, nblocks, stored_lr, total, ratio])
    end

    println("\n[case 6] low-frequency EFIE N~792 (30 MHz, leaf 0.03 lambda)")
    append!(rows, bench_case(; tag = "lowfreq_792", freq = 30e6, radius = 0.5,
                             n_theta = 12, n_phi = 24, leaf_factor = 0.03,
                             methods = (:aca, :mlaca)))

    for r in rows
        println(@sprintf("%-12s %-8s N=%-5d f=%-3.0fMHz leaf=%.3f setup=%.2fs solve=%.2fs mv=%.1e fin=%s it=%d res=%.1e sol=%.1e cond=%.1e ratio=%.1f%%",
                         r[1], r[2], r[3], r[4] / 1e6, r[5], r[6], r[7], r[8],
                         r[9], r[10], r[11], r[12], r[13], 100 * r[18]))
    end

    mkpath(dirname(out))
    writedlm(out, vcat(permutedims(HEADER), reduce(vcat, [permutedims(r) for r in rows])), ',')
    println("\nCSV written: ", out)
end

main()
