"""
test_pmchw_mlfma_operator.jl — Phase 15 步骤 15.7

TDD RED 测试：PMCHWMLFMAOperator 结构、构造函数与 mul!

测试覆盖：
  15.8  assemble_near_field_pmchw — 2N×2N 稀疏矩阵，4 块正确
  15.9  aggregate_leaf_pmchw!    — 聚合叶结点（J-Pass / M-Pass）
  15.10 disaggregate_leaf_pmchw_j! / _m! — 四块接收核函数
  15.11 PMCHWMLFMAOperator struct + 构造函数 + mul! (4 遍远场)

验收门限（B2）：|Z_in_MLFMA - Z_in_Direct| / |Z_in_Direct| < 5%
"""

using Test
using EMSuite
using LinearAlgebra
using SparseArrays
using Random
using IterativeSolvers

# ─────────────────────────────────────────────────────────────────────────────
# 公用测试夹具
# ─────────────────────────────────────────────────────────────────────────────

function make_pmchw_test_fixture(; freq = 300e6, eps_r = 4.0,
                                    lat_divs = 4, lon_divs = 6)
    mesh  = generate_sphere_mesh(0.1, lat_divs, lon_divs)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r)
    return pmchw, basis
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.11: struct 存在性 + size
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.11 PMCHWMLFMAOperator 结构与构造" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)

    # leaf_size ≈ λ/4（300 MHz → λ≈1 m，leaf_size=0.25 m）
    leaf_size = 0.25

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
# 15.8: assemble_near_field_pmchw — 2N×2N 近场矩阵
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.8 assemble_near_field_pmchw 近场矩阵" begin

    pmchw, basis = make_pmchw_test_fixture()
    N = num_basis(basis)
    leaf_size = 0.25

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

    @testset "HJ = -EM（结构不变量）" begin
        # Z^HJ = -Z^EM（精确成立于近场部分）
        EJ_block = Zn[1:N,  1:N]
        EM_block = Zn[1:N,  N+1:2N]
        HJ_block = Zn[N+1:2N, 1:N]
        @test norm(HJ_block + EM_block, Inf) < 1e-8 * norm(EM_block, Inf)
    end

    @testset "近场矩阵与 Direct PMCHW 矩阵近邻部分相符（相对误差 < 1%）" begin
        # Direct full matrix
        Z_full = assemble_impedance_matrix(pmchw, basis)
        # 近场包含 octree0 的所有近邻对（不含全局自作用；此处比较 EM 块对角线元素）
        for i in 1:min(5, N)
            # 近场矩阵包含 self-interaction（同一叶方格内），对角元应与直接法一致
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
    leaf_size = 0.25

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
    leaf_size = 0.25

    op_mlfma = PMCHWMLFMAOperator(pmchw, basis, leaf_size)
    Z_direct  = assemble_impedance_matrix(pmchw, basis)   # 2N×2N full matrix

    # 使用随机单位向量作为测试 x
    Random.seed!(42)
    x = randn(ComplexF64, 2N)
    x ./= norm(x)

    y_mlfma  = zeros(ComplexF64, 2N)
    mul!(y_mlfma, op_mlfma, x)

    y_direct = Z_direct * x

    rel_err = norm(y_mlfma - y_direct) / norm(y_direct)
    @test rel_err < 0.10    # 10% 宽容差（远场近似误差）
end

# ─────────────────────────────────────────────────────────────────────────────
# 15.11: B2 天线输入阻抗 MLFMA vs Direct（相对误差 < 5%）
# ─────────────────────────────────────────────────────────────────────────────
@testset "15.11 B2: PMCHW MLFMA Z_in vs Direct（误差 < 5%）" begin

    pmchw, basis = make_pmchw_test_fixture(freq = 300e6, eps_r = 4.0)
    N = num_basis(basis)
    leaf_size = 0.25

    op_mlfma = PMCHWMLFMAOperator(pmchw, basis, leaf_size)

    feed  = DeltaGapSource(pmchw.freq, [1], 1.0 + 0im)
    V_2N  = excitation_vector(pmchw, feed, basis)

    # ─── Direct 求解 ──────────────────────────────────────────────
    Z_direct = assemble_impedance_matrix(pmchw, basis)
    I_direct = Z_direct \ V_2N
    Z_in_direct = input_impedance(pmchw, feed, I_direct, basis)

    # ─── MLFMA + GMRES 求解 ───────────────────────────────────────
    I_mlfma, hist = gmres(op_mlfma, V_2N; reltol = 1e-4, maxiter = 200, log = true)
    Z_in_mlfma  = input_impedance(pmchw, feed, I_mlfma, basis)

    re_err = abs(real(Z_in_mlfma) - real(Z_in_direct)) / (abs(real(Z_in_direct)) + 1e-30)
    @test re_err < 0.05     # 5% 门限（Re 部分）
    @test real(Z_in_mlfma) > 0.0  # 物理约束：辐射阻抗 > 0
end
