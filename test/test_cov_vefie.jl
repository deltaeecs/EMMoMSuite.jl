# test_cov_vefie.jl — VEFIE 模块轻量覆盖（小体积网格 direct）
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.IntegralEquations.VEFIEModule:
    VEFIE, vefie_element_interaction, vefie_mass_matrix, vefie_mass_matrix_cached,
    precompute_vefie_basis, vefie_element_interaction_kernel
using LinearAlgebra

@testset "VEFIE 轻量覆盖" begin
    vol = generate_box_tet_mesh(0.5, 0.5, 0.5, 1, 1, 1)
    basis = SWGBasis(vol)
    perms = fill(ComplexF64(2.0), num_elements(vol))
    freq = 2e9
    vefie = VEFIE(freq, perms)
    N = num_basis(basis)

    # 元素级核与质量矩阵
    tets = get_tetrahedra_info(vol, basis, perms)
    M = vefie_mass_matrix(vefie, tets[1])
    @test size(M) == (4, 4)
    Z_el = vefie_element_interaction(vefie, tets[1], tets[1])
    Z_el isa Tuple && (Z_el = Z_el[1])
    @test size(Z_el) == (4, 4)
    cache = precompute_vefie_basis(vefie, tets)
    @test cache isa Vector
    Mc = vefie_mass_matrix_cached(vefie, tets[1], cache[1])
    @test size(Mc) == (4, 4)
    Zk = vefie_element_interaction_kernel(vefie, tets[1], tets[2], cache[1], cache[2])
    @test size(Zk) == (4, 4)

    # 直接装配
    Z = assemble_impedance_matrix(vefie, basis)
    @test size(Z) == (N, N)
    @test norm(Z) > 0

    # 激励（PWC 路径）
    pwc = PWCBasis(vol)
    src = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(vefie, src, pwc)
    @test length(V) == num_basis(pwc)
end
