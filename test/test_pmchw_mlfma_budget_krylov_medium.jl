"""
test_pmchw_mlfma_budget_krylov_medium.jl — Phase 15 PMCHW MLFMA 长 Krylov budget 专门门禁

目的：
  - 在 `N=540` medium 夹具上，把 budget 子流正式接到长 Krylov 子流；
  - 固定 `restart=300, maxiter=600, reltol=1e-6` 的 strong-form 求解；
  - 比较代表性 budget 中的 `default` 与 `loose_near`，锁定 backend fidelity 在长 Krylov 下会进入主导层。

说明：
  - 该测试运行成本很高，因此不并入默认 `runtests.jl`；
  - 它的定位是 Phase 15 的专门长 Krylov 验收门，而不是日常快速回归。
"""

using Test
using EMSuite
using IterativeSolvers
using LinearAlgebra
using SparseArrays
using EMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator

function make_pmchw_budget_krylov_medium_fixture()
    mesh = generate_sphere_mesh(0.5, 10, 20)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(300e6, 4.0)
    feed = DeltaGapSource(pmchw.freq, [1], 1.0 + 0im)
    Z_dense = assemble_impedance_matrix(pmchw, basis)
    rhs = excitation_vector(pmchw, feed, basis)
    Zin_lu = input_impedance(pmchw, feed, Z_dense \ rhs, basis)
    dense_shell = PMCHWBlockOperator(pmchw, basis)
    return pmchw, basis, feed, Z_dense, rhs, Zin_lu, dense_shell
end

function solve_strong(shell, rhs; restart, maxiter, reltol)
    coeffs_sf, hist = gmres(
        strong_form(shell),
        strong_form_rhs(shell, rhs);
        restart = restart,
        maxiter = maxiter,
        reltol = reltol,
        log = true,
    )
    return recover_trial_coefficients(shell, coeffs_sf), hist.data[:resnorm][end], length(hist.data[:resnorm])
end

function evaluate_long_budget_case(pmchw, basis, feed, Z_dense, rhs, Zin_lu, dense_result, name, budget)
    op = PMCHWMLFMAOperator(pmchw, basis, 0.10; budget = budget)
    shell = PMCHWBlockOperator(pmchw, basis, MatrixFreePMCHWBackend(op); block_source = Z_dense)
    I_mlfma, res_mlfma, niters_mlfma = solve_strong(shell, rhs; restart = 300, maxiter = 600, reltol = 1e-6)
    Zin_mlfma = input_impedance(pmchw, feed, I_mlfma, basis)

    return (
        name = name,
        leaf_size_eff = op.leaf_size_eff,
        near_range = op.near_range,
        near_density = nnz(op.Z_near) / (length(rhs)^2),
        Zin_mlfma = Zin_mlfma,
        Zin_re_err_pct = abs(real(Zin_mlfma) - real(Zin_lu)) / (max(abs(real(Zin_lu)), 1.0) + 1e-30) * 100,
        Zin_im_err_ohm = abs(imag(Zin_mlfma) - imag(Zin_lu)),
        Zin_gap_re_ohm = abs(real(Zin_mlfma) - real(dense_result.Zin)),
        Zin_gap_im_ohm = abs(imag(Zin_mlfma) - imag(dense_result.Zin)),
        res_mlfma = res_mlfma,
        niters_mlfma = niters_mlfma,
    )
end

@testset "15.BG2 PMCHW MLFMA medium long-Krylov budget gate" begin
    pmchw, basis, feed, Z_dense, rhs, Zin_lu, dense_shell = make_pmchw_budget_krylov_medium_fixture()

    I_dense, res_dense, niters_dense = solve_strong(dense_shell, rhs; restart = 300, maxiter = 600, reltol = 1e-6)
    Zin_dense = input_impedance(pmchw, feed, I_dense, basis)
    dense_result = (
        Zin = Zin_dense,
        Zin_re_err_pct = abs(real(Zin_dense) - real(Zin_lu)) / (max(abs(real(Zin_lu)), 1.0) + 1e-30) * 100,
        Zin_im_err_ohm = abs(imag(Zin_dense) - imag(Zin_lu)),
        res_dense = res_dense,
        niters_dense = niters_dense,
    )

    budgets = [
        (
            "default",
            PMCHWMLFMAErrorBudget(Float64),
        ),
        (
            "loose_near",
            PMCHWMLFMAErrorBudget(Float64;
                near_range_scale = 4.0,
                min_near_range = 4,
                max_near_range = 16,
            ),
        ),
    ]

    results = [
        evaluate_long_budget_case(pmchw, basis, feed, Z_dense, rhs, Zin_lu, dense_result, name, budget)
        for (name, budget) in budgets
    ]
    by_name = Dict(result.name => result for result in results)
    default = by_name["default"]
    loose = by_name["loose_near"]

    @info "PMCHW medium long-Krylov budget gate" N=num_basis(basis) Zin_lu Zin_dense dense_re_err_pct=dense_result.Zin_re_err_pct dense_im_err_ohm=dense_result.Zin_im_err_ohm dense_res=dense_result.res_dense default_near_density=default.near_density default_gap_re_ohm=default.Zin_gap_re_ohm default_gap_im_ohm=default.Zin_gap_im_ohm default_re_err_pct=default.Zin_re_err_pct default_im_err_ohm=default.Zin_im_err_ohm default_res=default.res_mlfma loose_near_density=loose.near_density loose_gap_re_ohm=loose.Zin_gap_re_ohm loose_gap_im_ohm=loose.Zin_gap_im_ohm loose_re_err_pct=loose.Zin_re_err_pct loose_im_err_ohm=loose.Zin_im_err_ohm loose_res=loose.res_mlfma

    @test num_basis(basis) == 540

    @test dense_result.Zin_re_err_pct < 8.0
    @test dense_result.Zin_im_err_ohm < 5.0
    @test dense_result.res_dense < 1e-4

    @test default.near_density > 0.85
    @test loose.near_density < 0.35
    @test default.near_density > loose.near_density

    @test default.Zin_gap_re_ohm < 0.1
    @test default.Zin_gap_im_ohm < 0.05
    @test loose.Zin_gap_re_ohm < 1.0
    @test loose.Zin_gap_im_ohm < 0.1
    @test loose.Zin_gap_re_ohm > 5 * default.Zin_gap_re_ohm

    @test abs(default.Zin_re_err_pct - dense_result.Zin_re_err_pct) < 0.1
    @test abs(loose.Zin_re_err_pct - dense_result.Zin_re_err_pct) < 0.3
    @test default.res_mlfma < 1e-4
    @test loose.res_mlfma < 1e-4
end