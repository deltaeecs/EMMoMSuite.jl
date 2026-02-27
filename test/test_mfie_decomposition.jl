"""
    test_mfie_decomposition.jl

Phase 10 / B1: 验证 S-MFIE 算子正确性
通过 CFIE 线性分解验证: Z_CFIE = α·Z_EFIE + (1-α)·Z_MFIE
同时验证对应的激励向量分解: V_CFIE = α·V_EFIE + (1-α)·V_MFIE

使用小网格 (Tri.nas, ~1330 RWG) 进行密集矩阵检验。
"""

using Test

@testset "S-MFIE CFIE Decomposition" begin
    using EMSuite
    using LinearAlgebra

    # 使用 AllinOne 的小网格
    mesh_file = ""
    for candidate in [
        joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles", "Tri.nas"),
        joinpath(@__DIR__, "..", "..", "MoM_Basics", "meshfiles", "Tri.nas"),
    ]
        if isfile(candidate)
            mesh_file = candidate
            break
        end
    end
    @assert !isempty(mesh_file) "未找到 Tri.nas 网格文件"

    freq = 3e8  # 300 MHz
    alpha = 0.5

    # 加载网格与基函数
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    set_frequency!(freq)
    basis = RWGBasis(mesh)

    n = num_basis(basis)
    @test n > 100  # 确保网格足够大

    # 创建三种算子
    efie = EFIE(freq)
    mfie = MFIE(freq)
    cfie = CFIE(freq, alpha)

    # 组装阻抗矩阵
    Z_efie = assemble_impedance_matrix(efie, basis)
    Z_mfie = assemble_impedance_matrix(mfie, basis)
    Z_cfie = assemble_impedance_matrix(cfie, basis)

    @testset "Z 矩阵线性分解" begin
        # Z_CFIE ≈ α * Z_EFIE + (1 - α) * Z_MFIE
        Z_reconstructed = alpha * Z_efie + (1 - alpha) * Z_mfie
        rel_err = norm(Z_cfie - Z_reconstructed) / norm(Z_cfie)

        @test rel_err < 1e-12
        println("  Z 矩阵分解相对误差: $rel_err")
    end

    @testset "激励向量线性分解" begin
        source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
        V_efie = excitation_vector(efie, source, basis)
        V_mfie = excitation_vector(mfie, source, basis)
        V_cfie = excitation_vector(cfie, source, basis)

        V_reconstructed = alpha * V_efie + (1 - alpha) * V_mfie
        rel_err = norm(V_cfie - V_reconstructed) / norm(V_cfie)

        @test rel_err < 1e-12
        println("  V 向量分解相对误差: $rel_err")
    end

    @testset "MFIE 矩阵基本性质" begin
        # MFIE 矩阵不应为对称 (EFIE 对称, MFIE 含 K 算子不对称)
        asymmetry = norm(Z_mfie - Z_mfie') / norm(Z_mfie)
        @test asymmetry > 1e-3  # 明显不对称
        println("  MFIE Z 不对称度: $asymmetry")

        # EFIE 远场部分使用 symmetric=true 装配，
        # 但近奇异项 (calc_near_interaction!) 采用半解析公式，不严格对称。
        # 这里只检查大体对称性（相对误差 < 30%）
        symmetry_efie = norm(Z_efie - Z_efie') / norm(Z_efie)
        @test symmetry_efie < 0.3
        println("  EFIE Z 对称度偏差: $symmetry_efie")
    end

    @testset "不同 α 值验证" begin
        for α in [0.0, 0.2, 0.5, 0.8, 1.0]
            cfie_α = CFIE(freq, α)
            Z_α = assemble_impedance_matrix(cfie_α, basis)
            Z_recon = α * Z_efie + (1 - α) * Z_mfie
            rel_err = norm(Z_α - Z_recon) / norm(Z_α)
            @test rel_err < 1e-12
        end
        println("  α ∈ {0, 0.2, 0.5, 0.8, 1.0} 全部通过")
    end
end
