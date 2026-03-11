"""
test_pmchw_mlfma_operator_medium.jl — Phase 15 中尺度正确性回归

目的：避免小矩阵 (N=54) 对 PMCHWMLFMAOperator 正确性的证明力度不足。

本文件固定一个仍可做 Direct 对照、但明显更大的球面夹具：
  r = 0.5m, lat_divs = 10, lon_divs = 20  -> N = 540, 2N = 1080

验证目标分两层：
1. 算子正确性：MLFMA matvec 与 Direct matvec 对齐。
2. 迭代一致性：在相同 GMRES 参数下，MLFMA operator 与 dense matrix 的 GMRES
   结果应一致，从而把更大规模上的偏差归因到“收敛不足/预条件不足”，而非 MLFMA matvec 本身。
"""

using Test
using EMSuite
using LinearAlgebra
using Random
using SparseArrays
using IterativeSolvers
using EMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator

function make_medium_pmchw_fixture()
    mesh = generate_sphere_mesh(0.5, 10, 20)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(300e6, 4.0)
    return pmchw, basis
end

function relcorr(a, b)
    rel = norm(a - b) / (norm(b) + 1e-30)
    corr = abs(dot(a, b)) / ((norm(a) * norm(b)) + 1e-30)
    return rel, corr
end

@testset "15.M1 Medium-scale Gate C matvec" begin
    pmchw, basis = make_medium_pmchw_fixture()
    N = num_basis(basis)

    op = PMCHWMLFMAOperator(pmchw, basis, 0.10)
    Z = assemble_impedance_matrix(pmchw, basis)

    Random.seed!(42)
    x = randn(ComplexF64, 2N)
    x ./= norm(x)

    y_mlfma = zeros(ComplexF64, 2N)
    mul!(y_mlfma, op, x)
    y_direct = Z * x

    rel, corr = relcorr(y_mlfma, y_direct)
    @info "Medium-scale Gate C" N rel corr nnz_near=nnz(op.Z_near) nlevels0=op.octree0.nLevels nlevels1=op.octree1.nLevels

    @test N == 540
    @test rel < 2e-3
    @test corr > 0.9999
end

@testset "15.M2 Medium-scale GMRES parity" begin
    pmchw, basis = make_medium_pmchw_fixture()
    N = num_basis(basis)

    op = PMCHWMLFMAOperator(pmchw, basis, 0.10)
    Z = assemble_impedance_matrix(pmchw, basis)
    feed = DeltaGapSource(pmchw.freq, [1], 1.0 + 0im)
    V = excitation_vector(pmchw, feed, basis)

    I_dense, hist_dense = gmres(Z, V; reltol = 1e-4, maxiter = 200, log = true)
    I_mlfma, hist_mlfma = gmres(op, V; reltol = 1e-4, maxiter = 200, log = true)

    Zin_dense = input_impedance(pmchw, feed, I_dense, basis)
    Zin_mlfma = input_impedance(pmchw, feed, I_mlfma, basis)

    rel_I, corr_I = relcorr(I_mlfma, I_dense)
    rel_re = abs(real(Zin_mlfma) - real(Zin_dense)) / (abs(real(Zin_dense)) + 1e-30)
    rel_im = abs(imag(Zin_mlfma) - imag(Zin_dense)) / (abs(imag(Zin_dense)) + 1e-30)
    res_dense = hist_dense.data[:resnorm][end]
    res_mlfma = hist_mlfma.data[:resnorm][end]
    res_gap = abs(res_mlfma - res_dense) / (abs(res_dense) + 1e-30)

    @info "Medium-scale GMRES parity" N rel_I corr_I rel_re rel_im res_dense res_mlfma res_gap

    @test N == 540
    @test rel_I < 1e-2
    @test corr_I > 0.999
    @test rel_re < 1e-3
    @test rel_im < 1e-3
    @test res_gap < 1e-2
end