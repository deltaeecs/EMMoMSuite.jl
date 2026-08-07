"""
test_nmuller_planewave_gmres_trajectory_medium.jl — Phase 15 PMCHW vs N-Muller 中尺度 plane-wave GMRES 轨迹对照

目的：
  - 把已有的 random-RHS dense trajectory 对照推进到物理激励 RHS；
  - 在同一 medium dielectric sphere 夹具上，使用 plane-wave excitation 比较 PMCHW 与 N-Muller
    相对各自 LU 的 GMRES 轨迹；
  - 验证“NMuller 更像 formulation / conditioning 改善”不依赖随机 probe。

说明：
  - 该测试运行成本较高，因此不并入默认 runtests；
  - 它只比较各 formulation 在各自 physical RHS 下相对各自 LU 的收敛行为，
    不把未校准的 N-Muller 端口语义混入验收门。
"""

using Test
using EMMoMSuite
using LinearAlgebra
using IterativeSolvers

function make_nmuller_planewave_fixture(; freq = 120e6, eps_r = 4.0, mu_r = 1.0, radius = 0.1, n_theta = 6, n_phi = 10)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r, mu_r)
    nmuller = NMuller(freq, eps_r, mu_r)
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    return pmchw, nmuller, basis, source
end

function solve_gmres_against_lu(A::AbstractMatrix, rhs::AbstractVector; reltol = 1e-5, maxiter = 250)
    x_lu = A \ rhs
    x_gmres, hist = gmres(A, rhs; reltol = reltol, maxiter = maxiter, log = true)
    final_rel_res = norm(rhs - A * x_gmres) / (norm(rhs) + 1e-30)
    return (
        rel_err = norm(x_gmres - x_lu) / (norm(x_lu) + 1e-30),
        rel_res = final_rel_res,
        niters = length(hist.data[:resnorm]),
    )
end

@testset "PMCHW vs N-Muller medium plane-wave GMRES trajectory" begin
    pmchw, nmuller, basis, source = make_nmuller_planewave_fixture()
    N = num_basis(basis)
    Z_pmchw = assemble_impedance_matrix(pmchw, basis)
    Z_nmuller = assemble_impedance_matrix(nmuller, basis)
    rhs_pmchw = excitation_vector(pmchw, source, basis)
    rhs_nmuller = excitation_vector(nmuller, source, basis)

    checkpoints = [50, 100, 150, 200, 250]
    pmchw_rows = Dict{Int,NamedTuple}()
    nmuller_rows = Dict{Int,NamedTuple}()

    for checkpoint in checkpoints
        pmchw_rows[checkpoint] = solve_gmres_against_lu(Z_pmchw, rhs_pmchw; reltol = 1e-5, maxiter = checkpoint)
        nmuller_rows[checkpoint] = solve_gmres_against_lu(Z_nmuller, rhs_nmuller; reltol = 1e-5, maxiter = checkpoint)
    end

    @info "PMCHW vs N-Muller medium plane-wave GMRES trajectory" N checkpoints pmchw_rel=[pmchw_rows[k].rel_err for k in checkpoints] nmuller_rel=[nmuller_rows[k].rel_err for k in checkpoints] pmchw_res=[pmchw_rows[k].rel_res for k in checkpoints] nmuller_res=[nmuller_rows[k].rel_res for k in checkpoints] pmchw_iters=[pmchw_rows[k].niters for k in checkpoints] nmuller_iters=[nmuller_rows[k].niters for k in checkpoints]

    @test N == 150
    @test norm(rhs_pmchw) > 0
    @test norm(rhs_nmuller) > 0

    for checkpoint in checkpoints
        @test nmuller_rows[checkpoint].rel_err < pmchw_rows[checkpoint].rel_err
        @test nmuller_rows[checkpoint].rel_res < pmchw_rows[checkpoint].rel_res
        @test pmchw_rows[checkpoint].niters == checkpoint
    end

    @test nmuller_rows[50].niters == 50
    @test nmuller_rows[100].niters < 100
    @test nmuller_rows[250].niters < pmchw_rows[250].niters
    @test nmuller_rows[250].rel_err < 0.02
    @test nmuller_rows[250].rel_res < 1e-4
    @test pmchw_rows[250].rel_err > 0.9
    @test pmchw_rows[250].rel_res > 1e-3
end