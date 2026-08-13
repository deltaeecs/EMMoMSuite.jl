"""
test_pmchw_mlfma_operator.jl — Phase 15 步骤 15.7

TDD RED 测试：PMCHWMLFMAOperator 结构、构造函数与 mul!

测试覆盖：
  15.8  assemble_near_field_pmchw — 2N×2N 稀疏矩阵，4 块正确
  15.9  aggregate_leaf_pmchw!    — 聚合叶结点（J-Pass / M-Pass）
  15.10 disaggregate_leaf_pmchw_j! / _m! — 四块接收核函数
  15.11 PMCHWMLFMAOperator struct + 构造函数 + mul! (4 遍远场)

验收门限（B2）：|Z_in_MLFMA - Z_in_Direct| / |Z_in_Direct| < 5%

夹具参数：r=0.5m, lat_divs=4, lon_divs=6, leaf_size=0.1m
  → N=54, nnz_near=4216/11664 (36.1%)，真正的 MLFMA 近/远场分离
"""

using Test
using EMMoMSuite
using LinearAlgebra
using SparseArrays
using Random
using IterativeSolvers
using EMMoMSuite.FastAlgorithms.MLFMA: build_octree, aggregate!, aggregate_leaf!, disaggregate_leaf!
using EMMoMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: aggregate_leaf_pmchw!, disaggregate_leaf_pmchw_j!, disaggregate_leaf_pmchw_m!
using EMMoMSuite.IntegralEquations.PMCHWModule: assemble_K_pmchw_offdiag
using EMMoMSuite.FastAlgorithms.MLFMA.Aggregation: aggregate_upward!
using EMMoMSuite.FastAlgorithms.MLFMA.Disaggregation: disaggregate_downward!
using EMMoMSuite.FastAlgorithms.MLFMA.Translation: translate!

# ─────────────────────────────────────────────────────────────────────────────
# 公用测试夹具
# ─────────────────────────────────────────────────────────────────────────────

function make_pmchw_test_fixture(; freq = 300e6, eps_r = 4.0,
                                    lat_divs = 4, lon_divs = 6)
    # r=0.5m, leaf_size=0.1m → octree 中叶片数 ~20+，nnz_near≈36% < 100%
    # 确保远场确实有贡献（MLFMA 非平凡）
    mesh  = generate_sphere_mesh(0.5, lat_divs, lon_divs)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r)
    return pmchw, basis
end

function clear_agg!(oct)
    for (_, lv) in oct.levels
        if isdefined(lv, :aggS)
            fill!(lv.aggS, zero(eltype(lv.aggS)))
        end
        if isdefined(lv, :disaggG)
            fill!(lv.disaggG, zero(eltype(lv.disaggG)))
        end
    end
end

function up_translate_down!(oct)
    for levelID in (oct.nLevels - 1):-1:2
        aggregate_upward!(oct.levels[levelID], oct.levels[levelID + 1])
    end
    for levelID in 2:oct.nLevels
        translate!(oct.levels[levelID])
    end
    for levelID in 2:(oct.nLevels - 1)
        disaggregate_downward!(oct.levels[levelID], oct.levels[levelID + 1])
    end
end

function compare_vectors(a, b)
    rel = norm(a - b) / (norm(b) + 1e-30)
    ratio = norm(a) / (norm(b) + 1e-30)
    corr = abs(dot(a, b)) / ((norm(a) * norm(b)) + 1e-30)
    return rel, ratio, corr
end

struct SharedEFIEProbe <: AbstractIntegralOperator
    freq::Float64
    k::Float64
    eta::Float64
    factor::ComplexF64
end

function near_pairs_ej(octree, sorted_ids, N)
    leaf_level = octree.levels[octree.nLevels]
    pairs = Set{Tuple{Int,Int}}()

    for i_cube in 1:leaf_level.nCubes
        cube = leaf_level.cubes[i_cube]
        isempty(cube.bfInterval) && continue
        my_bfs = [sorted_ids[s] for s in cube.bfInterval]

        for neigh_idx in cube.neighbors
            neigh_cube = leaf_level.cubes[neigh_idx]
            isempty(neigh_cube.bfInterval) && continue
            neigh_bfs = [sorted_ids[s] for s in neigh_cube.bfInterval]

            for i in my_bfs, j in neigh_bfs
                push!(pairs, (i, j))
            end
        end
    end

    return pairs
end

function leaf_far_pairs_from_octree(octree, sorted_ids)
    leaf_level = octree.levels[octree.nLevels]
    pairs = Set{Tuple{Int,Int}}()

    for i_cube in 1:leaf_level.nCubes
        cube = leaf_level.cubes[i_cube]
        isempty(cube.bfInterval) && continue
        test_bfs = [sorted_ids[s] for s in cube.bfInterval]

        for j_cube in cube.farneighbors
            src_cube = leaf_level.cubes[j_cube]
            isempty(src_cube.bfInterval) && continue
            src_bfs = [sorted_ids[s] for s in src_cube.bfInterval]
            for i in test_bfs, j in src_bfs
                push!(pairs, (i, j))
            end
        end
    end

    return pairs
end

function direct_m_far_reference(pmchw, basis, octree, sorted_ids, kmode::Symbol)
    k_real = kmode === :k0 ? pmchw.k0 : abs(real(pmchw.k1))
    eta_real = kmode === :k0 ? pmchw.eta0 : abs(real(pmchw.eta1))
    factor_hm = kmode === :k0 ? (im * pmchw.k0 / (pmchw.eta0 * 16 * pi)) : (im * pmchw.k1 / (pmchw.eta1 * 16 * pi))

    Z_em = assemble_K_pmchw_offdiag(basis, k_real)
    Z_hm = assemble_impedance_matrix(efie_from_keta(k_real, eta_real, factor_hm), basis)

    for (i, j) in near_pairs_ej(octree, sorted_ids, num_basis(basis))
        Z_em[i, j] = 0
        Z_hm[i, j] = 0
    end

    return Z_em, Z_hm
end

function shared_k1_core_metrics(pmchw, basis; leaf_size = 0.10, near_range = 4)
    N = num_basis(basis)
    k1r = real(pmchw.k1)
    eta1r = abs(real(pmchw.eta1))
    factor1 = im * pmchw.k1 * pmchw.eta1 / (16 * pi)
    probe, octree, sorted_ids = build_shared_k1_probe_context(pmchw, basis; leaf_size = leaf_size, near_range = near_range)

    Random.seed!(42)
    x = randn(ComplexF64, N)
    x ./= norm(x)

    aggregate!(octree, basis, probe, x, sorted_ids)
    for levelID in 2:octree.nLevels
        translate!(octree.levels[levelID])
    end
    for levelID in 2:(octree.nLevels - 1)
        disaggregate_downward!(octree.levels[levelID], octree.levels[levelID + 1])
    end

    y_far = zeros(ComplexF64, N)
    disaggregate_leaf!(octree.levels[octree.nLevels], basis, probe, y_far, sorted_ids)
    y_far .*= (4 * probe.factor)

    Z_ref = assemble_impedance_matrix(efie_from_keta(k1r, eta1r, factor1), basis)
    for (i, j) in near_pairs_ej(octree, sorted_ids, N)
        Z_ref[i, j] = 0
    end
    y_ref = Z_ref * x
    rel, ratio, corr = compare_vectors(y_far, y_ref)

    return rel, ratio, corr, octree.nLevels
end

function build_shared_k1_probe_context(pmchw, basis; leaf_size = 0.10, near_range = 4)
    k1r = real(pmchw.k1)
    eta1r = abs(real(pmchw.eta1))
    freq1 = 299792458.0 * k1r / (2 * pi)
    factor1 = im * pmchw.k1 * pmchw.eta1 / (16 * pi)
    probe = SharedEFIEProbe(freq1, k1r, eta1r, factor1)

    centers = reduce(hcat, [bf.center for bf in basis.functions])
    octree, sorted_ids = build_octree(centers, leaf_size; λ = 299792458.0 / freq1, near_range = near_range)
    return probe, octree, sorted_ids
end

function run_shared_k1_to_leaf(pmchw, basis; leaf_size = 0.10, near_range = 4, seed = 42)
    probe, octree, sorted_ids = build_shared_k1_probe_context(pmchw, basis; leaf_size = leaf_size, near_range = near_range)
    Random.seed!(seed)
    x = randn(ComplexF64, num_basis(basis))
    x ./= norm(x)

    aggregate!(octree, basis, probe, x, sorted_ids)
    for levelID in 2:octree.nLevels
        translate!(octree.levels[levelID])
    end
    for levelID in 2:(octree.nLevels - 1)
        disaggregate_downward!(octree.levels[levelID], octree.levels[levelID + 1])
    end

    return probe, octree, sorted_ids, x
end

function direct_far_ej_reference(pmchw, basis, op, which::Symbol)
    N = num_basis(basis)
    if which === :k0
        k = pmchw.k0
        eta = pmchw.eta0
        factor = im * pmchw.k0 * pmchw.eta0 / (16 * pi)
    elseif which === :k1
        k = real(pmchw.k1)
        eta = abs(real(pmchw.eta1))
        factor = im * pmchw.k1 * pmchw.eta1 / (16 * pi)
    else
        error("unsupported kernel selector: $which")
    end

    Z = assemble_impedance_matrix(efie_from_keta(k, eta, factor), basis)
    ZnEJ = op.Z_near[1:N, 1:N]
    I, J, _ = findnz(ZnEJ)
    for t in eachindex(I)
        Z[I[t], J[t]] = 0
    end
    return Z
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD0: Gate D0 — Dense PMCHW 基线先验复验（Dense -> Mie）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD0 Gate D0 Dense PMCHW 基线复验" begin

    freq = 300e6
    radius = 0.5
    eps_r = 4.0
    mu_r = 1.0

    mesh = generate_sphere_mesh(radius, 8, 12)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    pmchw = PMCHW(freq, eps_r, mu_r)

    @test N > 0

    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(pmchw, source, basis)
    Z = assemble_impedance_matrix(pmchw, basis)
    I = Z \ V

    theta_deg = collect(0.0:10.0:180.0)
    theta_rad = deg2rad.(theta_deg)
    phi_obs = [0.0, pi / 2]

    rcs_comp, _, _ = radarCrossSection(theta_rad, phi_obs, I, basis, pmchw.k0, pmchw.eta0)
    rcs_e_mie, rcs_h_mie, _ = calculate_mie_rcs_dielectric_sphere(radius, freq, theta_rad, eps_r, mu_r)

    rcs_e_dense_dB = 10 .* log10.(max.(vec(rcs_comp[1, :, 1]), 1e-100))
    rcs_h_dense_dB = 10 .* log10.(max.(vec(rcs_comp[2, :, 2]), 1e-100))
    rcs_e_mie_dB = 10 .* log10.(max.(rcs_e_mie, 1e-100))
    rcs_h_mie_dB = 10 .* log10.(max.(rcs_h_mie, 1e-100))

    acc_e = compute_rcs_accuracy(rcs_e_dense_dB, rcs_e_mie_dB, theta_deg, "PMCHW Dense E-plane"; threshold = 1.5)
    acc_h = compute_rcs_accuracy(rcs_h_dense_dB, rcs_h_mie_dB, theta_deg, "PMCHW Dense H-plane"; threshold = 2.0)

    @info "Gate D0 Dense PMCHW vs Mie" e_rmse=acc_e.rmse_dB e_bs=acc_e.backscatter_err_dB h_rmse=acc_h.rmse_dB h_bs=acc_h.backscatter_err_dB

    @test acc_e.pass
    @test acc_e.backscatter_err_dB < 0.2
    @test acc_h.pass
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD1: Gate D1 — 近场逐元素对齐（Z_near vs Dense PMCHW）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD1 Gate D1 近场逐元素对齐" begin

    pmchw, basis = make_pmchw_test_fixture()
    op_mlfma = PMCHWMLFMAOperator(pmchw, basis, 0.10)

    Zn = op_mlfma.Z_near
    Z_full = assemble_impedance_matrix(pmchw, basis)
    I, J, V = findnz(Zn)

    dense_vals = ComplexF64[Z_full[i, j] for (i, j) in zip(I, J)]
    rel_errs = abs.(V .- dense_vals) ./ (abs.(dense_vals) .+ 1e-30)
    max_err = maximum(rel_errs)
    mean_err = sum(rel_errs) / length(rel_errs)
    worst_idx = argmax(rel_errs)

    @info "Gate D1 near-field element check" max_err mean_err worst_i=I[worst_idx] worst_j=J[worst_idx]

    @test max_err < 1e-3
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD2: Gate D2 — EJ k1 最坏列扫描（far path only）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD2 Gate D2 EJ k1 最坏列扫描" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)
    op = PMCHWMLFMAOperator(pmchw, basis, 0.10)
    ZEJ1_far = direct_far_ej_reference(pmchw, basis, op, :k1)

    worst_col = 0
    worst_rel = -1.0
    worst_corr = 0.0
    worst_row = 0

    for col in 1:N
        x2N = zeros(ComplexF64, 2N)
        x2N[col] = 1.0 + 0im

        yk1 = zeros(ComplexF64, 2N)
        clear_agg!(op.octree1)
        aggregate_leaf_pmchw!(op.octree1, basis, x2N, op.sorted_ids1, 1:N, pmchw.k1)
        up_translate_down!(op.octree1)
        disaggregate_leaf_pmchw_j!(op.octree1, basis, pmchw, yk1, op.sorted_ids1, :k1)

        y = yk1[1:N]
        yref = ZEJ1_far[:, col]
        rel, _, corr = compare_vectors(y, yref)

        if rel > worst_rel
            worst_rel = rel
            worst_col = col
            worst_corr = corr
            elem_rel = abs.(y .- yref) ./ (abs.(yref) .+ 1e-30)
            worst_row = argmax(elem_rel)
        end
    end

    @info "Gate D2 EJ k1 worst column" worst_col worst_row worst_rel worst_corr

    @test worst_rel < 0.15
    @test worst_corr > 0.95
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD2S: Gate D2 Shared-Core 分界测试（generic EFIE k1 far path）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD2S Gate D2 shared core EFIE k1 对照" begin

    pmchw, basis = make_pmchw_test_fixture()
    rel, ratio, corr, nLevels = shared_k1_core_metrics(pmchw, basis)

    @info "Gate D2 shared core EFIE k1" rel ratio corr nLevels

    @test_broken rel < 0.15
    @test_broken corr > 0.95
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD2L: Gate D2 层级分界（4-level 红, 3-level 绿）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD2L Gate D2 level boundary" begin

    pmchw, basis = make_pmchw_test_fixture()
    rel4, _, corr4, nlevels4 = shared_k1_core_metrics(pmchw, basis; leaf_size = 0.10)
    rel3, _, corr3, nlevels3 = shared_k1_core_metrics(pmchw, basis; leaf_size = 0.20)

    @info "Gate D2 level boundary" rel4 corr4 nlevels4 rel3 corr3 nlevels3

    @test nlevels4 == 4
    @test nlevels3 == 3
    @test_broken rel4 < 0.15
    @test_broken corr4 > 0.95
    @test rel3 < 0.15
    @test corr3 > 0.95
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD2P: Gate D2 中间层对齐（4-level level-3 vs 3-level leaf）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD2P Gate D2 level-3 parity" begin

    pmchw, basis = make_pmchw_test_fixture()
    _, oct4, _, _ = run_shared_k1_to_leaf(pmchw, basis; leaf_size = 0.10)
    _, oct3, _, _ = run_shared_k1_to_leaf(pmchw, basis; leaf_size = 0.20)

    field4 = oct4.levels[3].disaggG
    field3 = oct3.levels[3].disaggG

    rel = norm(field4 - field3) / (norm(field3) + 1e-30)
    corr = abs(sum(conj.(field4) .* field3)) / ((norm(field4) * norm(field3)) + 1e-30)

    @info "Gate D2 level-3 parity" rel corr cubes4=oct4.levels[3].nCubes cubes3=oct3.levels[3].nCubes

    @test rel < 1e-3
    @test corr > 0.999
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD2A: Gate D2 Aggregation 包装对照（PMCHW J/k1 vs shared EFIE）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD2A Gate D2 aggregation parity" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)
    probe, octree, sorted_ids = build_shared_k1_probe_context(pmchw, basis)

    Random.seed!(42)
    x = randn(ComplexF64, N)
    x ./= norm(x)

    aggregate!(octree, basis, probe, x, sorted_ids)
    agg_generic = copy(octree.levels[octree.nLevels].aggS)

    x2N = zeros(ComplexF64, 2N)
    x2N[1:N] .= x
    aggregate_leaf_pmchw!(octree, basis, x2N, sorted_ids, 1:N, pmchw.k1)
    agg_pmchw = copy(octree.levels[octree.nLevels].aggS)

    rel = norm(agg_pmchw - agg_generic) / (norm(agg_generic) + 1e-30)
    @info "Gate D2 aggregation parity" rel

    @test rel < 1e-12
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD2R: Gate D2 Receive 包装对照（PMCHW J/k1 vs shared EFIE）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD2R Gate D2 receive parity" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)
    probe, octree, sorted_ids = build_shared_k1_probe_context(pmchw, basis)

    Random.seed!(42)
    x = randn(ComplexF64, N)
    x ./= norm(x)

    aggregate!(octree, basis, probe, x, sorted_ids)
    for levelID in 2:octree.nLevels
        translate!(octree.levels[levelID])
    end
    for levelID in 2:(octree.nLevels - 1)
        disaggregate_downward!(octree.levels[levelID], octree.levels[levelID + 1])
    end

    y_generic = zeros(ComplexF64, N)
    disaggregate_leaf!(octree.levels[octree.nLevels], basis, probe, y_generic, sorted_ids)
    y_generic .*= (4 * probe.factor)

    y_pmchw = zeros(ComplexF64, 2N)
    disaggregate_leaf_pmchw_j!(octree, basis, pmchw, y_pmchw, sorted_ids, :k1)

    rel, ratio, corr = compare_vectors(y_pmchw[1:N], y_generic)
    @info "Gate D2 receive parity" rel ratio corr

    @test rel < 1e-12
    @test corr > 1 - 1e-12
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD2RM: Gate D2 M-pass receive parity（当前预期 RED）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD2RM Gate D2 M-pass receive parity" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)
    op = PMCHWMLFMAOperator(pmchw, basis, 0.10)

    Random.seed!(43)
    x_m = zeros(ComplexF64, 2N)
    x_m[N+1:2N] .= randn(ComplexF64, N)
    x_m ./= norm(x_m)
    x_coeff = x_m[N+1:2N]

    for (kmode, octree, sorted_ids) in ((:k0, op.octree0, op.sorted_ids0), (:k1, op.octree1, op.sorted_ids1))
        k = kmode === :k0 ? pmchw.k0 : pmchw.k1
        Z_em, Z_hm = direct_m_far_reference(pmchw, basis, octree, sorted_ids, kmode)
        y_ref = vcat(Z_em * x_coeff, Z_hm * x_coeff)

        clear_agg!(octree)
        aggregate_leaf_pmchw!(octree, basis, x_m, sorted_ids, (N + 1):(2N), k)
        up_translate_down!(octree)

        y_pmchw = zeros(ComplexF64, 2N)
        disaggregate_leaf_pmchw_m!(octree, basis, pmchw, y_pmchw, sorted_ids, kmode)

        rel, ratio, corr = compare_vectors(y_pmchw, y_ref)
        rel_e, _, corr_e = compare_vectors(y_pmchw[1:N], y_ref[1:N])
        rel_h, _, corr_h = compare_vectors(y_pmchw[N+1:2N], y_ref[N+1:2N])
        @info "Gate D2 M-pass receive parity" kmode rel ratio corr rel_e corr_e rel_h corr_h

        @test rel < 0.15
        @test corr > 0.95
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD2T: Gate D2 leaf translation 主导性
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD2T Gate D2 leaf translation dominance" begin

    pmchw, basis = make_pmchw_test_fixture()
    probe, octree, sorted_ids = build_shared_k1_probe_context(pmchw, basis; leaf_size = 0.10)
    Random.seed!(42)
    x = randn(ComplexF64, num_basis(basis))
    x ./= norm(x)

    aggregate!(octree, basis, probe, x, sorted_ids)
    for levelID in 2:octree.nLevels
        translate!(octree.levels[levelID])
    end

    leaf_translation = copy(octree.levels[octree.nLevels].disaggG)
    for levelID in 2:(octree.nLevels - 1)
        disaggregate_downward!(octree.levels[levelID], octree.levels[levelID + 1])
    end

    leaf_total = octree.levels[octree.nLevels].disaggG
    leaf_parent = leaf_total - leaf_translation
    ratio = norm(leaf_parent) / (norm(leaf_translation) + 1e-30)

    @info "Gate D2 leaf translation dominance" ratio norm_leaf_translation=norm(leaf_translation) norm_leaf_parent=norm(leaf_parent)

    @test ratio < 0.05
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD2U: Gate D2 leaf translation vs 仅 leaf-farneighbor dense 对照
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD2U Gate D2 leaf translation isolated" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)
    k1r = real(pmchw.k1)
    eta1r = abs(real(pmchw.eta1))
    fac1 = im * pmchw.k1 * pmchw.eta1 / (16 * pi)

    probe, octree, sorted_ids = build_shared_k1_probe_context(pmchw, basis; leaf_size = 0.10)
    Random.seed!(42)
    x = randn(ComplexF64, N)
    x ./= norm(x)

    aggregate!(octree, basis, probe, x, sorted_ids)
    for levelID in 2:octree.nLevels
        translate!(octree.levels[levelID])
    end

    y_leaf = zeros(ComplexF64, N)
    disaggregate_leaf!(octree.levels[octree.nLevels], basis, probe, y_leaf, sorted_ids)
    y_leaf .*= (4 * probe.factor)

    Z = assemble_impedance_matrix(efie_from_keta(k1r, eta1r, fac1), basis)
    keep_pairs = leaf_far_pairs_from_octree(octree, sorted_ids)
    Z_leaf = zeros(ComplexF64, size(Z))
    for (i, j) in keep_pairs
        Z_leaf[i, j] = Z[i, j]
    end
    y_ref = Z_leaf * x

    rel, ratio, corr = compare_vectors(y_leaf, y_ref)
    @info "Gate D2 leaf translation isolated" rel ratio corr npairs=length(keep_pairs)

    @test_broken rel < 0.15
    @test_broken corr > 0.95
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GD2V: Gate D2 near-range 阈值分界
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GD2V Gate D2 near-range threshold" begin

    pmchw, basis = make_pmchw_test_fixture()
    rel4, _, corr4, _ = shared_k1_core_metrics(pmchw, basis; leaf_size = 0.10, near_range = 4)
    rel7, _, corr7, _ = shared_k1_core_metrics(pmchw, basis; leaf_size = 0.10, near_range = 7)

    @info "Gate D2 near-range threshold" rel4 corr4 rel7 corr7

    @test_broken rel4 < 0.15
    @test_broken corr4 > 0.95
    @test rel7 < 0.15
    @test corr7 > 0.95
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.11: struct 存在性 + size
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.11 PMCHWMLFMAOperator 结构与构造" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)

    # leaf_size = 0.1m（约 λ/10，300 MHz λ≈1 m）
    leaf_size = 0.10

    op_mlfma = PMCHWMLFMAOperator(pmchw, basis, leaf_size)

    @testset "类型正确" begin
        @test op_mlfma isa PMCHWMLFMAOperator
    end

    @testset "size 为 2N×2N" begin
        @test size(op_mlfma) == (2N, 2N)
        @test size(op_mlfma, 1) == 2N
        @test size(op_mlfma, 2) == 2N
    end

    @testset "eltype 为 ComplexF64" begin
        @test eltype(op_mlfma) == ComplexF64
    end

    @testset "octree0 / octree1 均已建立" begin
        @test isdefined(op_mlfma, :octree0)
        @test isdefined(op_mlfma, :octree1)
    end

    @testset "sorted_ids 长度为 N" begin
        @test length(op_mlfma.sorted_ids0) == N
        @test length(op_mlfma.sorted_ids1) == N
    end

    @testset "Z_near 尺寸为 2N×2N" begin
        @test size(op_mlfma.Z_near) == (2N, 2N)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GA: Gate A — 结构不变量（Theory -> Implementation）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GA Gate A 结构不变量" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)
    op_mlfma = PMCHWMLFMAOperator(pmchw, basis, 0.10)
    Zn = op_mlfma.Z_near

    @testset "A1 HJ + EM = 0（近场结构不变量）" begin
        EM_block = Zn[1:N, N+1:2N]
        HJ_block = Zn[N+1:2N, 1:N]
        @test norm(HJ_block + EM_block, Inf) < 1e-8 * (norm(EM_block, Inf) + 1e-30)
    end

    @testset "A2 近远场非平凡划分" begin
        @test nnz(Zn) < (2N)^2
    end

    @testset "A3 双树拓扑索引一致性" begin
        @test sort(op_mlfma.sorted_ids0) == collect(1:N)
        @test sort(op_mlfma.sorted_ids1) == collect(1:N)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.GB: Gate B — 分 Pass 对齐（可检验链路）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.GB Gate B EJ pass 对齐" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)
    op = PMCHWMLFMAOperator(pmchw, basis, 0.10)

    # Direct EJ k0/k1 reference with PMCHW near-mask removed.
    k0 = pmchw.k0
    eta0 = pmchw.eta0
    k1r = real(pmchw.k1)
    eta1r = abs(real(pmchw.eta1))
    fac0 = im * pmchw.k0 * pmchw.eta0 / (16 * pi)
    fac1 = im * pmchw.k1 * pmchw.eta1 / (16 * pi)
    ZEJ0 = assemble_impedance_matrix(efie_from_keta(k0, eta0, fac0), basis)
    ZEJ1 = assemble_impedance_matrix(efie_from_keta(k1r, eta1r, fac1), basis)

    ZnEJ = op.Z_near[1:N, 1:N]
    I, J, _ = findnz(ZnEJ)
    for t in eachindex(I)
        ZEJ0[I[t], J[t]] = 0
        ZEJ1[I[t], J[t]] = 0
    end

    Random.seed!(42)
    xJ = randn(ComplexF64, N)
    xJ ./= norm(xJ)
    x2N = zeros(ComplexF64, 2N)
    x2N[1:N] .= xJ

    yk0 = zeros(ComplexF64, 2N)
    clear_agg!(op.octree0)
    aggregate_leaf_pmchw!(op.octree0, basis, x2N, op.sorted_ids0, 1:N, pmchw.k0)
    up_translate_down!(op.octree0)
    disaggregate_leaf_pmchw_j!(op.octree0, basis, pmchw, yk0, op.sorted_ids0, :k0)

    yk1 = zeros(ComplexF64, 2N)
    clear_agg!(op.octree1)
    aggregate_leaf_pmchw!(op.octree1, basis, x2N, op.sorted_ids1, 1:N, pmchw.k1)
    up_translate_down!(op.octree1)
    disaggregate_leaf_pmchw_j!(op.octree1, basis, pmchw, yk1, op.sorted_ids1, :k1)

    yk0E = yk0[1:N]
    yk1E = yk1[1:N]
    yk0_true = ZEJ0 * xJ
    yk1_true = ZEJ1 * xJ

    rel0, ratio0, corr0 = compare_vectors(yk0E, yk0_true)
    rel1, ratio1, corr1 = compare_vectors(yk1E, yk1_true)

    @info "Gate B EJ k0" rel0 ratio0 corr0
    @info "Gate B EJ k1" rel1 ratio1 corr1

    @test rel0 < 0.15
    @test corr0 > 0.95

    @test rel1 < 0.15
    @test corr1 > 0.95
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.8: assemble_near_field_pmchw — 2N×2N 近场矩阵
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.8 assemble_near_field_pmchw 近场矩阵" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)
    leaf_size = 0.10

    op_mlfma = PMCHWMLFMAOperator(pmchw, basis, leaf_size)
    Zn = op_mlfma.Z_near

    @testset "稀疏矩阵类型" begin
        @test Zn isa SparseMatrixCSC
    end

    @testset "4 块均非全零" begin
        @test nnz(Zn[1:N,   1:N])   > 0   # EJ 块
        @test nnz(Zn[1:N,   N+1:2N]) > 0  # EM 块
        @test nnz(Zn[N+1:2N, 1:N])   > 0  # HJ 块
        @test nnz(Zn[N+1:2N, N+1:2N]) > 0 # HM 块
    end

    @testset "近场矩阵为真正稀疏（nnz < (2N)^2，远场有贡献）" begin
        # 核心非平凡性检查：若失败说明所有基函数在同一叶片内
        @test nnz(Zn) < (2N)^2
        @info "nnz_near=$(nnz(Zn)) < full=$(( 2N)^2) ($(round(100*nnz(Zn)/(2N)^2, digits=1))%)"
    end

    @testset "HJ = -EM（结构不变量）" begin
        # Z^HJ = -Z^EM（精确成立于近场部分）
        EJ_block = Zn[1:N,  1:N]
        EM_block = Zn[1:N,  N+1:2N]
        HJ_block = Zn[N+1:2N, 1:N]
        @test norm(HJ_block + EM_block, Inf) < 1e-8 * norm(EM_block, Inf)
    end

    @testset "近场矩阵与 Direct PMCHW 矩阵近邻部分相符（相对误差 < 1%）" begin
        Z_full = assemble_impedance_matrix(pmchw, basis)
        for i in 1:min(5, N)
            @test abs(Zn[i,i] - Z_full[i,i]) / (abs(Z_full[i,i]) + 1e-30) < 0.01
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.11: mul! 基本性质
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.11 mul! 基本性质" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)
    leaf_size = 0.10

    op_mlfma = PMCHWMLFMAOperator(pmchw, basis, leaf_size)

    @testset "零输入 → 零输出" begin
        x = zeros(ComplexF64, 2N)
        y = zeros(ComplexF64, 2N)
        mul!(y, op_mlfma, x)
        @test norm(y) < 1e-12
    end

    @testset "输出向量长度为 2N" begin
        x = randn(ComplexF64, 2N)
        y = zeros(ComplexF64, 2N)
        mul!(y, op_mlfma, x)
        @test length(y) == 2N
    end

    @testset "线性性：mul!(y, A, 2x) == 2*mul!(y, A, x)" begin
        x  = randn(ComplexF64, 2N)
        y1 = zeros(ComplexF64, 2N)
        y2 = zeros(ComplexF64, 2N)
        mul!(y1, op_mlfma, 2x)
        mul!(y2, op_mlfma, x)
        @test norm(y1 - 2y2) / (norm(y2) + 1e-30) < 1e-6
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.11: MLFMA vs Direct — 矩阵向量积近似一致（宽容差验证）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.11 MLFMA mul! vs Direct 矩阵向量积（相对误差 < 10%）" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)
    leaf_size = 0.10

    op_mlfma = PMCHWMLFMAOperator(pmchw, basis, leaf_size)

    # ── 非平凡性前置检查：近场矩阵 nnz 必须严格小于完整矩阵 (2N)² ──────────
    nnz_near = nnz(op_mlfma.Z_near)
    full_nnz = (2N)^2
    @test nnz_near < full_nnz   # 远场确实贡献，MLFMA 非平凡
    @info "MLFMA 非平凡性：nnz_near=$nnz_near < full_nnz=$full_nnz ($(round(100*nnz_near/full_nnz, digits=1))%)"

    Z_direct  = assemble_impedance_matrix(pmchw, basis)

    Random.seed!(42)
    x = randn(ComplexF64, 2N)
    x ./= norm(x)

    y_mlfma  = zeros(ComplexF64, 2N)
    mul!(y_mlfma, op_mlfma, x)

    y_direct = Z_direct * x

    rel_err = norm(y_mlfma - y_direct) / norm(y_direct)
    @info "MLFMA matvec 相对误差：$(round(rel_err*100, digits=2))%"
    @test rel_err < 0.10    # 10% 宽容差
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.11: B2 天线输入阻抗 MLFMA vs Direct（相对误差 < 5%）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.11 B2: PMCHW MLFMA Z_in vs Direct（误差 < 5%）" begin

    pmchw, basis = make_pmchw_test_fixture(freq = 300e6, eps_r = 4.0)
    N = num_basis(basis)
    leaf_size = 0.10

    op_mlfma = PMCHWMLFMAOperator(pmchw, basis, leaf_size)

    feed  = DeltaGapSource(pmchw.freq, [1], 1.0 + 0im)
    V_2N  = excitation_vector(pmchw, feed, basis)

    Z_direct = assemble_impedance_matrix(pmchw, basis)
    I_direct = Z_direct \ V_2N
    Z_in_direct = input_impedance(pmchw, feed, I_direct, basis)

    # restart=30（默认 min(30,N)）在本夹具下 GMRES 残差停滞（稠密矩阵亦然，
    # 实测 restart=2N 后 109 次收敛到机器精度、Zin 误差 ~0.2%），故用全空间 restart。
    I_mlfma, hist = gmres(
        op_mlfma, V_2N;
        restart = 2 * num_basis(basis), reltol = 1e-4, maxiter = 200, log = true,
    )
    Z_in_mlfma  = input_impedance(pmchw, feed, I_mlfma, basis)

    @info "Z_in_direct = $(round(Z_in_direct, digits=3))"
    @info "Z_in_mlfma  = $(round(Z_in_mlfma, digits=3))"

    re_err = abs(real(Z_in_mlfma) - real(Z_in_direct)) / (abs(real(Z_in_direct)) + 1e-30)
    @test re_err < 0.05
    @test real(Z_in_mlfma) > 0.0
end

@testset "PMCHWMLFMAOperator budget interface preserves defaults and exposes overrides" begin
    pmchw, basis = make_pmchw_test_fixture()

    op_default = PMCHWMLFMAOperator(pmchw, basis, 0.10)
    @test op_default.budget isa PMCHWMLFMAErrorBudget
    @test op_default.near_range == clamp(
        round(Int, 0.10 / op_default.leaf_size_eff * op_default.budget.near_range_scale),
        op_default.budget.min_near_range,
        op_default.budget.max_near_range,
    )
    @test op_default.budget.L_min == 0

    budget = PMCHWMLFMAErrorBudget(Float64;
        fixed_near_range = 9,
        fixed_leaf_size_eff = 0.04,
        L_min = 4,
    )
    op_budgeted = PMCHWMLFMAOperator(pmchw, basis, 0.10; budget = budget)

    @test op_budgeted.near_range == 9
    @test op_budgeted.leaf_size_eff ≈ 0.04
    @test op_budgeted.budget.L_min == 4
    @test op_budgeted.octree0.nLevels >= op_default.octree0.nLevels
    @test op_budgeted.octree1.nLevels >= op_default.octree1.nLevels
end
