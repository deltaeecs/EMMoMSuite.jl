# P1 性能门基准：FFTSpectral vs Lagrange2Step
# 输出：每层 φ 步耗时（稀疏 ϕCSC vs 批量 FFT）与 matvec 时间/内存分列
using EMMoMSuite
using LinearAlgebra, Random
using Printf
import EMMoMSuite.FastAlgorithms.MLFMA.Interpolation as Interp

function run_case(radius, nθ, nφ, leaf_ratio; nrep = 8)
    freq = 300e6
    λ = 299792458.0 / freq
    mesh = generate_sphere_mesh(radius, nθ, nφ)
    basis = RWGBasis(mesh)
    efie = EFIE(freq)
    leaf = leaf_ratio * λ
    op_lag = MLFMAOperator(efie, basis, leaf, Val(:Lagrange2Step))
    op_fft = MLFMAOperator(efie, basis, leaf, Val(:FFTSpectral))
    N = num_basis(basis)
    println("case: r=", radius, "λ mesh=(", nθ, ",", nφ, ") leaf=", leaf_ratio, "λ N=", N,
            " nLevels=", op_fft.octree.nLevels)
    println("level | τ | nθ | M1→M2 | nCubes | φ步 sparse(ms) | φ步 fft-batch(ms) | 比值")
    for lv in 3:op_fft.octree.nLevels
        info_l = op_lag.octree.levels[lv].interpWθϕ
        info_f = op_fft.octree.levels[lv].interpWθϕ
        nc = op_fft.octree.levels[lv].nCubes
        L = op_fft.octree.levels[lv].L
        agg = randn(ComplexF64, info_f.nθ * info_f.M1, 2, nc)
        out = zeros(ComplexF64, info_f.nθ * info_f.M2, 2, nc)
        # 稀疏路径总耗时（每盒×2 极化）
        y = info_l.ϕCSC * agg[:, 1, 1]
        t_sp = minimum(@elapsed(info_l.ϕCSC * agg[:, 1, 1]) for _ in 1:200)
        t_sp_total = t_sp * nc * 2 * 1000
        # 批量 FFT 路径
        Interp.fft_interp_phi_batch!(out, agg, info_f)
        t_ff = minimum(@elapsed(Interp.fft_interp_phi_batch!(out, agg, info_f)) for _ in 1:30)
        println(@sprintf("%d | %d | %d | %d→%d | %d | %.3f | %.3f | %.2f",
            lv, L, info_f.nθ, info_f.M1, info_f.M2, nc, t_sp_total, t_ff * 1000,
            t_sp_total / (t_ff * 1000)))
    end
    Random.seed!(9)
    x = randn(ComplexF64, N)
    op_lag * x; op_fft * x
    a_lag = @allocated op_lag * x
    a_fft = @allocated op_fft * x
    for _ in 1:3
        op_lag * x; op_fft * x
    end
    t_lag = minimum(@elapsed(op_lag * x) for _ in 1:nrep)
    t_fft = minimum(@elapsed(op_fft * x) for _ in 1:nrep)
    y_lag = op_lag * x
    y_fft = op_fft * x
    println(@sprintf("matvec: lag=%.1f ms  fft=%.1f ms  speedup=%.3f  alloc lag=%.0f MB  fft=%.0f MB  rel diff=%.2e",
        t_lag * 1000, t_fft * 1000, t_lag / t_fft, a_lag / 1e6, a_fft / 1e6,
        norm(y_fft - y_lag) / norm(y_lag)))
end

run_case(1.5, 16, 32, 0.25)
