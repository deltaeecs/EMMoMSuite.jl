# test/test_mlfma_params.jl - M4 MLFMA parameter configurability (TDD)
#
# Goal: nInterp (interpolation points) and precision_digits (truncation d0)
# are no longer hardcoded; configurable via MLFMAOperator/build_octree kwargs
# with defaults preserving current behavior (6 / 9.0).

using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using LinearAlgebra

@testset "MLFMA configurable parameters" begin
    mesh = generate_sphere_mesh(0.5, 8, 12)
    basis = RWGBasis(mesh)
    efie = EFIE(300e6)
    λ = 299792458.0 / 300e6
    leaf = 0.25 * λ

    op_default = MLFMAOperator(efie, basis, leaf)
    op_interp = MLFMAOperator(efie, basis, leaf; nInterp = 8)
    op_prec = MLFMAOperator(efie, basis, leaf; precision_digits = 3.0)

    # nInterp takes effect: inter-level interpolation matrix changes
    leaf_lvl_d = op_default.octree.levels[op_default.octree.nLevels]
    leaf_lvl_i = op_interp.octree.levels[op_interp.octree.nLevels]
    @test leaf_lvl_d.interpWθϕ != leaf_lvl_i.interpWθϕ

    # precision_digits takes effect: smaller truncation -> fewer top-level poles
    nPoles_d = length(op_default.octree.levels[2].poles.r̂sθsϕs)
    nPoles_p = length(op_prec.octree.levels[2].poles.r̂sθsϕs)
    @test nPoles_p < nPoles_d

    # Accuracy gate: MatVec still matches dense reference
    Z = assemble_impedance_matrix(efie, basis)
    x = randn(ComplexF64, num_basis(basis))
    y = op_interp * x
    @test norm(y - Z * x) / norm(Z * x) < 5e-2
    y2 = op_prec * x
    @test norm(y2 - Z * x) / norm(Z * x) < 5e-2
end
