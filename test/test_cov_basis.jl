# test_cov_basis.jl — BasisFunctions 覆盖率补测
using Test
using EMMoMSuite
using EMMoMSuite.BasisFunctions
using EMMoMSuite.Geometry
using LinearAlgebra, Random

@testset "RWG/SWG/PWC 基函数 evaluate 与 info" begin
    smesh = generate_sphere_mesh(0.5, 4, 8)
    rwg = RWGBasis(smesh)
    @test count_unknowns(rwg) == num_basis(rwg)
    tris = get_triangles_info(smesh, rwg)
    @test tris isa Vector
    @test length(tris) == num_elements(smesh)
    ti = get_triangle_info(smesh, rwg, 1)
    @test ti isa EMMoMSuite.Geometry.TriangleInfo
    # 混合网格（surface + volume）
    vol = generate_box_tet_mesh(1.0, 1.0, 1.0, 1, 1, 1)
    surf = extract_surface(vol)
    comp = CompositeMesh(surf, vol)
    @test num_vertices(comp) > 0
    @test EMMoMSuite.BasisFunctions.RWGBasis(comp) isa RWGBasis
    # 馈电边选择
    edges = select_gap_feed_edges(rwg)
    @test edges isa Vector
end

@testset "SWG/PWC 构造与未知数" begin
    vol = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)
    swg = SWGBasis(vol)
    @test num_basis(swg) > 0
    pwc = PWCBasis(vol)
    @test num_basis(pwc) > 0
    @test count_unknowns(swg) == num_basis(swg)
    # 六面体
    nodes = zeros(3, 8)
    for i in 1:8
        nodes[1, i] = (i - 1) % 2
        nodes[2, i] = ((i - 1) ÷ 2) % 2
        nodes[3, i] = (i - 1) ÷ 4
    end
    hexes = collect(reshape(1:8, 8, 1))
    hm = HexahedraMesh(1, nodes, hexes, [1])
    ph = PWCHexBasis(hm)
    @test num_basis(ph) > 0
end

@testset "基函数 info 构造（Tet/Hex）" begin
    vol = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)
    # 带 permittivity 的体积 info（VEFIE 路径）
    perm = fill(ComplexF64(2.0), num_elements(vol))
    info = get_tetrahedra_info(vol, SWGBasis(vol), perm)
    @test info isa Vector
    @test length(info) == num_elements(vol)
    @test info[1] isa EMMoMSuite.Geometry.TetrahedraInfo
end
