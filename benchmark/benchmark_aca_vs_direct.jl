# benchmark/benchmark_aca_vs_direct.jl — M2/Task5 ACA vs Direct 基准
#
# 输出：N、稠密装配/求解时间、ACA 建立时间、MatVec 相对误差、GMRES 残差、
# 解相对误差、近场 nnz、远场块数、低秩存储、压缩率；写 CSV 到 benchmark/results/。

using EMMoMSuite
using LinearAlgebra
using Printf
using SparseArrays
using DelimitedFiles
using IterativeSolvers

function run_benchmark(;
    freq::Real = 300e6,
    radius::Real = 0.5,
    n_theta::Int = 12,
    n_phi::Int = 24,
    leaf_factor::Real = 0.25,
    tol::Real = 1e-4,
    near_range::Int = 1,
    out::String = joinpath(@__DIR__, "results", "aca_benchmark_2026-08-11.csv"),
)
    λ = 299792458.0 / freq
    println("======================================================")
    println("   Benchmark: ACA vs Direct Solver (EFIE, RWG sphere)")
    println("======================================================")
    println(@sprintf("freq=%.2f MHz, radius=%.4f m, leaf=%.3f λ, tol=%.0e, near_range=%d",
                     freq / 1e6, radius, leaf_factor, tol, near_range))

    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println(@sprintf("Unknowns: %d", N))

    efie = EFIE(freq)

    # 1. Direct
    t_direct = @elapsed Z = assemble_impedance_matrix(efie, basis)
    println(@sprintf("[1/3] Direct assembly: %.4f s", t_direct))
    src = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(efie, src, basis)
    I_dir = Z \ V

    # 2. ACA
    leaf = leaf_factor * λ
    t_aca_setup = @elapsed op = ACAOperator(efie, basis, leaf; tol = tol, near_range = near_range)
    println(@sprintf("[2/3] ACA setup: %.4f s", t_aca_setup))

    # MatVec accuracy
    x = randn(ComplexF64, N)
    y = op * x
    mv_err = norm(y - Z * x) / norm(Z * x)
    println(@sprintf("      MatVec rel err vs direct: %.3e", mv_err))

    # Solve with ILU(Z_near) preconditioner
    P = ILUPreconditioner(op.Z_near, τ = 0.01)
    t_aca_solve = @elapsed begin
        I_aca = gmres(op, V; Pl = P, abstol = 1e-6, reltol = 1e-8, maxiter = 400, restart = 50)
    end
    println(@sprintf("      ACA GMRES solve: %.4f s", t_aca_solve))
    res = norm(V - op * I_aca) / norm(V)
    sol_err = norm(I_aca - I_dir) / norm(I_dir)
    println(@sprintf("      Relative residual: %.3e", res))
    println(@sprintf("      Solution rel err vs direct: %.3e", sol_err))

    # 3. Storage / compression
    nnz_near = nnz(op.Z_near)
    nblocks = length(op.blocks)
    stored_lr = 0
    for blk in op.blocks
        k = size(blk.U, 2)
        stored_lr += k * (length(blk.rows) + length(blk.cols))
    end
    total_stored = nnz_near + stored_lr
    full = N * N
    ratio = 1 - total_stored / full
    println(@sprintf("[3/3] near nnz=%d (%d blocks), stored low-rank=%d, total=%d / %d",
                     nnz_near, nblocks, stored_lr, total_stored, full))
    println(@sprintf("      Compression ratio vs dense: %.2f%%", 100 * ratio))

    # 4. CSV
    mkpath(dirname(out))
    header = ["N", "freq_Hz", "t_direct_assembly_s", "t_aca_setup_s", "t_aca_solve_s",
              "mv_rel_err", "gmres_rel_res", "sol_rel_err", "nnz_near", "n_blocks",
              "stored_lr", "total_stored", "compression_ratio"]
    row = [N, freq, t_direct, t_aca_setup, t_aca_solve, mv_err, res, sol_err,
           nnz_near, nblocks, stored_lr, total_stored, ratio]
    if isfile(out)
        data = readdlm(out, ',')
        vcat(data, permutedims(row)) |> x -> writedlm(out, x, ',')
    else
        writedlm(out, vcat(permutedims(header), permutedims(row)), ',')
    end
    println("      CSV written: ", out)

    return (N = N, t_direct = t_direct, t_aca_setup = t_aca_setup, t_aca_solve = t_aca_solve,
            mv_err = mv_err, res = res, sol_err = sol_err, nnz_near = nnz_near,
            nblocks = nblocks, stored_lr = stored_lr, total_stored = total_stored, ratio = ratio)
end

run_benchmark()
