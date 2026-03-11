"""
test_pmchw_gate_s_planewave_trajectory_medium.jl — Phase 15 PMCHW medium plane-wave dense weak/strong GMRES 轨迹

目的：
  - 在同一 medium dielectric sphere 夹具上，用 PlaneWave 物理激励比较 PMCHW dense weak/strong 的 GMRES 轨迹；
  - 把“conditioning 放大”从单点 Gate S 推进到 checkpoint 级正式门禁；
  - 为后续 deeper restart / Arnoldi 诊断提供一个更贴近实际散射工况的 PMCHW 主线基线。

说明：
  - 该测试运行成本较高，因此不并入默认 runtests；
  - 该测试只比较 PMCHW 自身 weak/strong 离散在同一物理 RHS 下相对 LU 的行为，不引入 N-Muller 对照。
"""

using Test
using EMSuite
using LinearAlgebra
using IterativeSolvers

function make_pmchw_planewave_gate_s_fixture(; freq = 120e6, eps_r = 4.0, mu_r = 1.0, radius = 0.1, n_theta = 6, n_phi = 10)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r, mu_r)
    shell = PMCHWBlockOperator(pmchw, basis)
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    rhs = excitation_vector(pmchw, source, basis)
    return pmchw, basis, shell, rhs
end

function solve_weak_against_lu(shell, rhs; reltol = 1e-5, maxiter = 250)
    A = weak_form(shell)
    x_lu = A \ rhs
    x_gmres, hist = gmres(A, rhs; reltol = reltol, maxiter = maxiter, log = true)
    return (
        rel_err = norm(x_gmres - x_lu) / (norm(x_lu) + 1e-30),
        rel_res = norm(rhs - A * x_gmres) / (norm(rhs) + 1e-30),
        niters = length(hist.data[:resnorm]),
    )
end

function solve_strong_against_lu(shell, rhs; reltol = 1e-5, maxiter = 250)
    A = strong_form(shell)
    rhs_sf = strong_form_rhs(shell, rhs)
    coeffs_lu = A \ rhs_sf
    coeffs_gmres, hist = gmres(A, rhs_sf; reltol = reltol, maxiter = maxiter, log = true)
    x_lu = recover_trial_coefficients(shell, coeffs_lu)
    x_gmres = recover_trial_coefficients(shell, coeffs_gmres)
    return (
        rel_err = norm(x_gmres - x_lu) / (norm(x_lu) + 1e-30),
        rel_res = norm(rhs - weak_form(shell) * x_gmres) / (norm(rhs) + 1e-30),
        niters = length(hist.data[:resnorm]),
    )
end

@testset "PMCHW medium plane-wave dense weak/strong GMRES trajectory" begin
    pmchw, basis, shell, rhs = make_pmchw_planewave_gate_s_fixture()
    N = num_basis(basis)
    checkpoints = [50, 100, 150, 200, 250]
    weak_rows = Dict{Int,NamedTuple}()
    strong_rows = Dict{Int,NamedTuple}()

    for checkpoint in checkpoints
        weak_rows[checkpoint] = solve_weak_against_lu(shell, rhs; reltol = 1e-5, maxiter = checkpoint)
        strong_rows[checkpoint] = solve_strong_against_lu(shell, rhs; reltol = 1e-5, maxiter = checkpoint)
    end

    @info "PMCHW medium plane-wave dense weak/strong GMRES trajectory" N checkpoints weak_rel=[weak_rows[k].rel_err for k in checkpoints] strong_rel=[strong_rows[k].rel_err for k in checkpoints] weak_res=[weak_rows[k].rel_res for k in checkpoints] strong_res=[strong_rows[k].rel_res for k in checkpoints] weak_iters=[weak_rows[k].niters for k in checkpoints] strong_iters=[strong_rows[k].niters for k in checkpoints]

    @test N == 150
    @test norm(rhs) > 0

    for checkpoint in checkpoints
        @test strong_rows[checkpoint].rel_err < weak_rows[checkpoint].rel_err
        @test strong_rows[checkpoint].rel_err / (weak_rows[checkpoint].rel_err + 1e-30) > 0.99
        @test abs(strong_rows[checkpoint].rel_res / (weak_rows[checkpoint].rel_res + 1e-30) - 1.0) < 0.1
        @test weak_rows[checkpoint].niters == checkpoint
        @test strong_rows[checkpoint].niters <= checkpoint
    end

    @test strong_rows[250].rel_err > 0.9
    @test strong_rows[250].rel_res > 1e-3
    @test weak_rows[250].rel_err > 0.9
    @test weak_rows[250].rel_res > 1e-3
end