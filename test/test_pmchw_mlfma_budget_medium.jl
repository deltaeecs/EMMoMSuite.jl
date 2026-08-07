"""
test_pmchw_mlfma_budget_medium.jl — Phase 15 PMCHW MLFMA budget 中尺度专门门禁

目的：
  - 在真正能拉开 near/far 划分的 `N=540` 夹具上，固定三组代表性 budget；
  - 锁定 budget 对 `near_density` 与 `matvec fidelity` 的影响；
  - 同时验证在固定 `GMRES(100)` strong-form 路径下，预算调节本身不会显著改变 `Z_in`。

说明：
  - 该测试运行成本较高，因此不并入默认 `runtests.jl`；
  - 它的定位是 Phase 15 的专门 budget 门，而不是日常快速回归。
"""

using Test
using EMMoMSuite
using IterativeSolvers
using LinearAlgebra
using Random
using SparseArrays
using EMMoMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator

function make_pmchw_budget_medium_fixture()
    mesh = generate_sphere_mesh(0.5, 10, 20)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(300e6, 4.0)
    feed = DeltaGapSource(pmchw.freq, [1], 1.0 + 0im)
    Z_dense = assemble_impedance_matrix(pmchw, basis)
    rhs = excitation_vector(pmchw, feed, basis)
    I_ref = Z_dense \ rhs
    Zin_ref = input_impedance(pmchw, feed, I_ref, basis)
    return pmchw, basis, feed, Z_dense, rhs, Zin_ref
end

function relcorr(a, b)
    rel = norm(a - b) / (norm(b) + 1e-30)
    corr = abs(dot(a, b)) / ((norm(a) * norm(b)) + 1e-30)
    return rel, corr
end

function evaluate_budget_case(pmchw, basis, feed, Z_dense, rhs, Zin_ref, name, budget)
    op = PMCHWMLFMAOperator(pmchw, basis, 0.10; budget = budget)

    Random.seed!(42)
    x = randn(ComplexF64, 2 * num_basis(basis))
    x ./= norm(x)
    y_mlfma = zeros(ComplexF64, length(x))
    mul!(y_mlfma, op, x)
    y_dense = Z_dense * x
    rel_matvec, corr_matvec = relcorr(y_mlfma, y_dense)

    shell = PMCHWBlockOperator(
        pmchw,
        basis,
        MatrixFreePMCHWBackend(op);
        block_source = Z_dense,
    )
    coeffs_sf, hist = gmres(
        strong_form(shell),
        strong_form_rhs(shell, rhs);
        reltol = 1e-4,
        maxiter = 100,
        log = true,
    )
    I_budget = recover_trial_coefficients(shell, coeffs_sf)
    Zin_budget = input_impedance(pmchw, feed, I_budget, basis)

    return (
        name = name,
        near_range = op.near_range,
        leaf_size_eff = op.leaf_size_eff,
        nnz_near = nnz(op.Z_near),
        near_density = nnz(op.Z_near) / (length(rhs)^2),
        rel_matvec = rel_matvec,
        corr_matvec = corr_matvec,
        Zin = Zin_budget,
        Zin_re_err_pct = abs(real(Zin_budget) - real(Zin_ref)) / (max(abs(real(Zin_ref)), 1.0) + 1e-30) * 100,
        Zin_im_err_ohm = abs(imag(Zin_budget) - imag(Zin_ref)),
        resnorm = hist.data[:resnorm][end],
    )
end

@testset "15.BG1 PMCHW MLFMA medium budget gate" begin
    pmchw, basis, feed, Z_dense, rhs, Zin_ref = make_pmchw_budget_medium_fixture()

    cases = [
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
        (
            "fixed_leaf_0p04_nr9",
            PMCHWMLFMAErrorBudget(Float64;
                fixed_leaf_size_eff = 0.04,
                fixed_near_range = 9,
            ),
        ),
    ]

    results = [evaluate_budget_case(pmchw, basis, feed, Z_dense, rhs, Zin_ref, name, budget) for (name, budget) in cases]
    by_name = Dict(result.name => result for result in results)

    default = by_name["default"]
    loose = by_name["loose_near"]
    fixed = by_name["fixed_leaf_0p04_nr9"]

    re_values = [real(result.Zin) for result in results]
    im_values = [imag(result.Zin) for result in results]
    resnorms = [result.resnorm for result in results]

    @info "PMCHW medium budget gate" N=num_basis(basis) Zin_ref default_near_density=default.near_density default_rel_matvec=default.rel_matvec loose_near_density=loose.near_density loose_rel_matvec=loose.rel_matvec fixed_near_density=fixed.near_density fixed_rel_matvec=fixed.rel_matvec re_spread=(maximum(re_values) - minimum(re_values)) im_spread=(maximum(im_values) - minimum(im_values)) resnorm_spread=(maximum(resnorms) - minimum(resnorms))

    @test num_basis(basis) == 540

    @test default.near_density > 0.85
    @test loose.near_density < 0.35
    @test fixed.near_density < 0.27
    @test default.near_density > loose.near_density > fixed.near_density

    @test default.rel_matvec < 5e-4
    @test loose.rel_matvec < 1e-3
    @test fixed.rel_matvec < 1e-3
    @test default.corr_matvec > 0.99999
    @test loose.corr_matvec > 0.99999
    @test fixed.corr_matvec > 0.99999
    @test default.rel_matvec < loose.rel_matvec
    @test default.rel_matvec < fixed.rel_matvec

    @test maximum(re_values) - minimum(re_values) < 0.2
    @test maximum(im_values) - minimum(im_values) < 0.2
    @test maximum(resnorms) - minimum(resnorms) < 2e-5
end