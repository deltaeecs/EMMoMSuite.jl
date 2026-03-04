"""
test_pmchw.jl — Phase 22 PMCHWT 均匀介质体 SIE 测试

测试覆盖：
  22.1  assemble_K_offdiag : 参数化 K 算子（无质量矩阵项）
  22.2  PMCHW 结构体 + 矩阵装配
  22.3  激励向量 excitation_pmchw
"""

using Test
using EMSuite
using LinearAlgebra
using StaticArrays

# ────────────────────────────────────────────────────────
# 工具：小型球面网格（用于单元测试，与大 Mie 验证分离）
# ────────────────────────────────────────────────────────
function make_small_sphere(; n_theta = 6, n_phi = 8, radius = 1.0)
    return generate_sphere_mesh(radius, n_theta, n_phi)
end

# ============================================================
# 22.1  assemble_K_offdiag — 参数化 K 算子
# ============================================================
@testset "22.1 assemble_K_offdiag" begin
    mesh  = make_small_sphere(n_theta = 4, n_phi = 6)
    basis = RWGBasis(mesh)
    N     = num_basis(basis)
    @test N > 0

    c0    = 299792458.0
    freq  = 100e6
    k0    = 2π * freq / c0

    K = assemble_K_offdiag(basis, k0)

    # 尺寸正确
    @test size(K) == (N, N)

    # 非零（K 算子有实质贡献）
    @test norm(K) > 0

    # 对角元不是严格零：每条基函数 m 的自交叉 t+_m ↔ t-_m 仍被计入（三角形不同）
    # 检查方式：仅验证 K 矩阵整体非零、复数，以及形状正确（上面已检查）
    # （mass-matrix 自项 ±0.5*I 被跳过，通过 PMCHW 结构不变量 Z^EM+Z^HJ=0 来验证）

    # K 矩阵应为复数值（虚部非零）
    @test norm(imag(K)) > 0
end

# ============================================================
# 22.2  PMCHW 构造与矩阵装配
# ============================================================
@testset "22.2 PMCHW Matrix Assembly" begin
    mesh  = make_small_sphere(n_theta = 4, n_phi = 6)
    basis = RWGBasis(mesh)
    N     = num_basis(basis)

    freq  = 100e6
    eps_r = 4.0  # 无损介质
    mu_r  = 1.0

    pmchw = PMCHW(freq, eps_r, mu_r)
    Z = assemble_impedance_matrix(pmchw, basis)

    # 尺寸应为 2N × 2N
    @test size(Z) == (2N, 2N)

    # 1. Z^EJ 子块（左上）应复对称 Z = Zᵀ（非 Hermitian，因为 G 是对称核）
    Z_EJ = Z[1:N, 1:N]
    @test norm(Z_EJ - transpose(Z_EJ)) / norm(Z_EJ) < 1e-8

    # 2. Z^HM 子块（右下）应复对称 Z = Zᵀ
    Z_HM = Z[N+1:2N, N+1:2N]
    @test norm(Z_HM - transpose(Z_HM)) / norm(Z_HM) < 1e-8

    # 3. 构造不变量：Z^EM = -Z^HJ（两块是同一 K 求和，符号相反）
    Z_EM = Z[1:N, N+1:2N]
    Z_HJ = Z[N+1:2N, 1:N]
    @test norm(Z_EM + Z_HJ) / (norm(Z_EM) + 1e-30) < 1e-10

    # 4. 各块非零
    @test norm(Z_EJ) > 0
    @test norm(Z_HM) > 0
    @test norm(Z_EM) > 0

    # 5. EJ 和 HM 量级比：factor_EJ ~ jkη, factor_HM ~ jk/η，比值 ~ η₀² ≈ 1.4e5
    #    有效性检查：两块量级之比应在 [η₀²/100, η₀²×100] 范围内
    eta0_sq = (sqrt(4π * 1e-7 / (1.0 / (299792458.0^2 * 4π * 1e-7))))^2
    ratio = norm(Z_EJ) / norm(Z_HM)
    @test ratio > eta0_sq / 100.0 && ratio < eta0_sq * 100.0
end

@testset "22.2 PMCHW PEC Limit (eps_r large)" begin
    # 当 eps_r 极大时，PMCHW 应退化为 PEC 行为
    # 此测试仅验证：eps_r 很大时不报错，矩阵正常构成
    mesh  = make_small_sphere(n_theta = 4, n_phi = 6)
    basis = RWGBasis(mesh)
    N     = num_basis(basis)

    freq  = 100e6
    pmchw_large = PMCHW(freq, 1e6, 1.0)  # "PEC-like" eps_r
    Z = assemble_impedance_matrix(pmchw_large, basis)
    @test size(Z) == (2N, 2N)
    @test !any(isnan, Z)
    @test !any(isinf, Z)
end

@testset "22.2 PMCHW lossy medium" begin
    mesh  = make_small_sphere(n_theta = 4, n_phi = 6)
    basis = RWGBasis(mesh)
    N     = num_basis(basis)

    freq  = 100e6
    # 有损介质 eps_r = 4 + 0.4j
    eps_r_complex = 4.0 + 0.4im
    pmchw_lossy = PMCHW(freq, eps_r_complex, 1.0)
    Z = assemble_impedance_matrix(pmchw_lossy, basis)
    @test size(Z) == (2N, 2N)
    @test !any(isnan, Z)
    # Z^HM 应该不再实对称（因为复数 eps_r），但矩阵本身有效
    @test norm(Z) > 0
end

# ============================================================
# 22.3  激励向量 excitation_pmchw
# ============================================================
@testset "22.3 excitation_pmchw" begin
    mesh  = make_small_sphere(n_theta = 4, n_phi = 6)
    basis = RWGBasis(mesh)
    N     = num_basis(basis)

    freq  = 100e6
    c0    = 299792458.0
    mu0   = 4π * 1e-7
    eps0  = 1.0 / (c0^2 * mu0)
    eta0  = sqrt(mu0 / eps0)

    # +z 方向传播（theta=pi/2, phi=0），z 极化平面波
    # k_hat = (1,0,0), E^inc 极化为 (0,0,1)（与传播方向垂直）
    # H^inc = k̂ × E / η₀ = (1,0,0)×(0,0,1)/η₀ = (0,1,0)/η₀ ≠ 0
    src = PlaneWave(freq, π/2, 0.0, [0.0, 0.0, 1.0])

    pmchw = PMCHW(freq, 4.0, 1.0)
    V = excitation_vector(pmchw, src, basis)

    # 长度为 2N
    @test length(V) == 2N

    # V_E（前 N 个）和 V_H（后 N 个）均非零
    V_E = V[1:N]
    V_H = V[N+1:2N]
    @test norm(V_E) > 0
    @test norm(V_H) > 0

    # V_H = ∫ f · H^inc dS，H^inc = (k̂ × E^inc)/η₀
    # 因此 |V_H|/|V_E| ≈ 1/η₀ ≈ 2.65e-3（H 场比 E 场小 η₀ 倍）
    # 注意：旧实现错用 MFIE 激励（η₀ × ∫ f·(n̂×H^inc)），使 ratio ≈ O(1)，物理不正确
    ratio = norm(V_H) / norm(V_E)
    @test 1e-4 < ratio < 0.05  # 正确物理：ratio ≈ 1/η₀ ≈ 2.65e-3
end

# ============================================================
# 22.2 + 22.3  端到端求解（不验证 RCS 精度，只保证流程畅通）
# ============================================================
@testset "22.2/22.3 End-to-end solve (no Mie check)" begin
    mesh  = make_small_sphere(n_theta = 4, n_phi = 6)
    basis = RWGBasis(mesh)
    N     = num_basis(basis)

    freq  = 100e6
    eps_r = 4.0
    pmchw = PMCHW(freq, eps_r, 1.0)
    Z = assemble_impedance_matrix(pmchw, basis)

    src = PlaneWave(freq, π/2, 0.0, [0.0, 0.0, 1.0])  # z 极化，与 k̂=[1,0,0] 垂直
    V   = excitation_vector(pmchw, src, basis)

    # LU 求解
    I_coeff = Z \ V

    # 系数非零
    @test norm(I_coeff) > 0
    @test !any(isnan, I_coeff)

    # 系数分解：I^J（前 N）和 I^M（后 N）
    I_J = I_coeff[1:N]
    I_M = I_coeff[N+1:2N]
    @test norm(I_J) > 0
    @test norm(I_M) > 0
end
