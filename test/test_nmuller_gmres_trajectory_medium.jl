"""
test_nmuller_gmres_trajectory_medium.jl — Phase 15 PMCHW vs N-Muller 中尺度 GMRES 轨迹对照

目的：
  - 在同一 medium dense dielectric sphere 夹具上，比较 PMCHW 与 N-Muller 的 GMRES 收敛轨迹；
  - 把“NMuller 更像 formulation 改善而不是 backend 改善”的判断升级为正式专门回归；
  - 固定多个迭代预算点，直接比较两种 formulation 相对各自 LU 的误差与残差。

说明：
  - 该测试运行成本较高，因此不并入默认 runtests；
  - 它的定位是 Phase 15 的 formulation / conditioning 专门门禁，而不是日常快速回归。
"""

using Test
using EMMoMSuite
using LinearAlgebra
using IterativeSolvers
using Random

function make_nmuller_trajectory_fixture(; freq = 120e6, eps_r = 4.0, mu_r = 1.0, radius = 0.1, n_theta = 6, n_phi = 10)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r, mu_r)
    nmuller = NMuller(freq, eps_r, mu_r)
    return pmchw, nmuller, basis
end

function solve_gmres_against_lu(A::AbstractMatrix, rhs::AbstractVector; reltol = 1e-5, maxiter = 250)
    x_lu = A \ rhs
    x_gmres, hist = gmres(A, rhs; reltol = reltol, maxiter = maxiter, log = true)
    resnorms = hist.data[:resnorm]
    return (
        rel_err = norm(x_gmres - x_lu) / (norm(x_lu) + 1e-30),
        resnorm = resnorms[end],
        niters = length(resnorms),
    )
end

@testset "PMCHW vs N-Muller medium dense GMRES trajectory" begin
    pmchw, nmuller, basis = make_nmuller_trajectory_fixture()
    N = num_basis(basis)
    Z_pmchw = assemble_impedance_matrix(pmchw, basis)
    Z_nmuller = assemble_impedance_matrix(nmuller, basis)

    Random.seed!(42)
    rhs = ComplexF64.(randn(2N), randn(2N))
    rhs ./= norm(rhs)

    checkpoints = [50, 100, 150, 200, 250]
    pmchw_rows = Dict{Int,NamedTuple}()
    nmuller_rows = Dict{Int,NamedTuple}()

    for checkpoint in checkpoints
        pmchw_rows[checkpoint] = solve_gmres_against_lu(Z_pmchw, rhs; reltol = 1e-5, maxiter = checkpoint)
        nmuller_rows[checkpoint] = solve_gmres_against_lu(Z_nmuller, rhs; reltol = 1e-5, maxiter = checkpoint)
    end

    @info "PMCHW vs N-Muller medium dense GMRES trajectory" N checkpoints pmchw_rel=[pmchw_rows[k].rel_err for k in checkpoints] nmuller_rel=[nmuller_rows[k].rel_err for k in checkpoints] pmchw_res=[pmchw_rows[k].resnorm for k in checkpoints] nmuller_res=[nmuller_rows[k].resnorm for k in checkpoints]

    @test N == 150
    for checkpoint in checkpoints
        @test nmuller_rows[checkpoint].rel_err < pmchw_rows[checkpoint].rel_err
        @test nmuller_rows[checkpoint].resnorm < pmchw_rows[checkpoint].resnorm
        @test nmuller_rows[checkpoint].niters == checkpoint
        @test pmchw_rows[checkpoint].niters == checkpoint
    end

    @test nmuller_rows[250].rel_err < 0.1
    @test pmchw_rows[250].rel_err > 0.9
    @test nmuller_rows[250].resnorm < 0.1
    @test pmchw_rows[250].resnorm > 0.5
end