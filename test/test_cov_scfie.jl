# test_cov_scfie.jl — SCFIE 模块轻量覆盖（小网格 direct，无 MLFMA）
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.IntegralEquations.SCFIEModule: SCFIE, scfie_coupling_interaction, scfie_sv_only_interaction
using LinearAlgebra

@testset "SCFIE 轻量覆盖" begin
    vol = generate_box_tet_mesh(0.5, 0.5, 0.5, 1, 1, 1)
    surf = extract_surface(vol)
    basis_surf = RWGBasis(surf)
    basis_vol = SWGBasis(vol)
    perms = fill(ComplexF64(2.0), num_elements(vol))
    freq = 2e9
    scfie = SCFIE(freq, perms; alpha = 0.5)
    n_surf = num_basis(basis_surf)
    n_vol = num_basis(basis_vol)
    n_total = n_surf + n_vol

    # 元素级核
    tris = get_triangles_info(surf, basis_surf)
    tetras = get_tetrahedra_info(vol, basis_vol, perms)
    Z_sv, Z_vs = scfie_coupling_interaction(scfie, tris[1], tetras[1])
    @test size(Z_sv) == (3, 4)
    @test size(Z_vs) == (4, 3)
    Z_opt = scfie_sv_only_interaction(scfie, tris[1], tetras[1])
    @test norm(Z_opt - Z_sv) / norm(Z_sv) < 1e-12

    # 直接装配（含全部子块）
    Z = assemble_impedance_matrix(scfie, basis_surf, basis_vol)
    @test size(Z) == (n_total, n_total)
    @test norm(Z) > 0
    Z_SS = Z[1:n_surf, 1:n_surf]
    Z_SV = Z[1:n_surf, n_surf+1:end]
    Z_VS = Z[n_surf+1:end, 1:n_surf]
    Z_VV = Z[n_surf+1:end, n_surf+1:end]
    @test norm(Z_SS) > 0
    @test norm(Z_SV) > 0
    @test norm(Z_VS) > 0
    @test norm(Z_VV) > 0

    # 激励与求解（SCFIE 支持 DeltaGapSource）
    src = DeltaGapSource(freq, [1], 1.0 + 0im)
    V = excitation_vector(scfie, src, basis_surf, basis_vol)
    @test length(V) == n_total
    I = Z \ V
    @test all(isfinite, I)
end
