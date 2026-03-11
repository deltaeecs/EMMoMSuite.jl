"""
test_pmchw_nmuller_planewave_restart_sweep_medium.jl — Phase 15 plane-wave dense restart 敏感性对照

目的：
  - 在同一 medium dielectric sphere 夹具上，显式扫描 dense GMRES 的 `restart` 参数；
  - 判断此前 PMCHW 的 plane-wave 轨迹失真有多少来自默认 `restart=20`；
  - 在隔离 restart 影响后，再检查 PMCHW 相对 N-Muller 是否仍保留 formulation gap。

说明：
  - 该测试运行成本较高，因此不并入默认 runtests；
  - 它只比较各 formulation 在各自 physical RHS 下相对各自 LU 的行为，
    不把未校准的端口语义混入门禁。
"""

using Test
using EMSuite
using LinearAlgebra
using IterativeSolvers

function make_restart_sweep_fixture(; freq = 120e6, eps_r = 4.0, mu_r = 1.0, radius = 0.1, n_theta = 6, n_phi = 10)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r, mu_r)
    nmuller = NMuller(freq, eps_r, mu_r)
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    return pmchw, nmuller, basis, source
end

function solve_gmres_against_lu(A::AbstractMatrix, rhs::AbstractVector; restart::Int, reltol = 1e-5, maxiter = 250)
    x_lu = A \ rhs
    x_gmres, hist = gmres(A, rhs; restart = restart, reltol = reltol, maxiter = maxiter, log = true)
    return (
        rel_err = norm(x_gmres - x_lu) / (norm(x_lu) + 1e-30),
        rel_res = norm(rhs - A * x_gmres) / (norm(rhs) + 1e-30),
        niters = length(hist.data[:resnorm]),
    )
end

function sweep_restarts(A::AbstractMatrix, rhs::AbstractVector, restarts)
    rows = Dict{Int, NamedTuple}()
    for restart in restarts
        rows[restart] = solve_gmres_against_lu(A, rhs; restart = restart)
    end
    return rows
end

@testset "PMCHW vs N-Muller medium plane-wave dense restart sweep" begin
    pmchw, nmuller, basis, source = make_restart_sweep_fixture()
    N = num_basis(basis)
    Z_pmchw = assemble_impedance_matrix(pmchw, basis)
    Z_nmuller = assemble_impedance_matrix(nmuller, basis)
    rhs_pmchw = excitation_vector(pmchw, source, basis)
    rhs_nmuller = excitation_vector(nmuller, source, basis)
    restarts = [20, 50, 100, 150, 250]

    pmchw_rows = sweep_restarts(Z_pmchw, rhs_pmchw, restarts)
    nmuller_rows = sweep_restarts(Z_nmuller, rhs_nmuller, restarts)

    @info "PMCHW vs N-Muller medium plane-wave dense restart sweep" N restarts pmchw_rel=[pmchw_rows[k].rel_err for k in restarts] pmchw_res=[pmchw_rows[k].rel_res for k in restarts] pmchw_iters=[pmchw_rows[k].niters for k in restarts] nmuller_rel=[nmuller_rows[k].rel_err for k in restarts] nmuller_res=[nmuller_rows[k].rel_res for k in restarts] nmuller_iters=[nmuller_rows[k].niters for k in restarts]

    @test N == 150
    @test norm(rhs_pmchw) > 0
    @test norm(rhs_nmuller) > 0

    @test pmchw_rows[20].niters == 250
    @test pmchw_rows[20].rel_err > 0.9
    @test pmchw_rows[20].rel_res > 1e-3

    @test pmchw_rows[150].niters == 250
    @test pmchw_rows[150].rel_err < 0.35
    @test pmchw_rows[150].rel_res < 1e-4

    @test pmchw_rows[250].niters < 200
    @test pmchw_rows[250].rel_err < 0.1
    @test pmchw_rows[250].rel_res < 1e-5
    @test pmchw_rows[250].rel_err / (pmchw_rows[20].rel_err + 1e-30) < 0.1
    @test pmchw_rows[250].rel_res / (pmchw_rows[20].rel_res + 1e-30) < 0.01

    @test nmuller_rows[20].niters < 100
    @test nmuller_rows[250].rel_err < 0.02
    @test nmuller_rows[250].rel_res < 1e-5
    @test isapprox(nmuller_rows[100].rel_err, nmuller_rows[250].rel_err; rtol = 1e-10, atol = 1e-12)
    @test isapprox(nmuller_rows[100].rel_res, nmuller_rows[250].rel_res; rtol = 1e-10, atol = 1e-12)
    @test nmuller_rows[100].niters == nmuller_rows[250].niters

    @test pmchw_rows[250].rel_err > 4 * nmuller_rows[250].rel_err
    @test pmchw_rows[250].rel_res < 2 * nmuller_rows[250].rel_res
end