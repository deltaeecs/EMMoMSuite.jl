# test_cov_geometry.jl — Geometry 模块覆盖率补测（轻量冒烟）
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using LinearAlgebra, StaticArrays, Random

@testset "Mesh 类型构造与查询" begin
    nodes = [0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0]
    tris = [1 2; 2 4; 3 3]
    tm = TriangleMesh(2, nodes, tris, [1, 2])
    @test num_vertices(tm) == 4
    @test num_elements(tm) == 2
    @test dimension(tm) == 2
    @test vertices(tm) == nodes
    @test elements(tm) == tris
    # 四面体网格
    tn = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    tets = [1 2; 2 3; 3 4; 4 1]
    tet_mesh = TetrahedraMesh(2, tn, tets, [1, 1])
    @test num_elements(tet_mesh) == 2
    @test dimension(tet_mesh) == 3
    # 六面体网格
    hn = zeros(3, 8)
    for i in 1:8
        hn[1, i] = (i - 1) % 2
        hn[2, i] = ((i - 1) ÷ 2) % 2
        hn[3, i] = (i - 1) ÷ 4
    end
    hexes = [1 5; 2 6; 3 7; 4 8; 5 1; 6 2; 7 3; 8 4]
    hm = HexahedraMesh(1, hn, hexes, [1])
    @test num_elements(hm) == 2
    @test isfinite(tet_volume(tn[:, 1], tn[:, 2], tn[:, 3], tn[:, 4]))
    @test isfinite(hex_volume(hn[:, 1], hn[:, 2], hn[:, 3], hn[:, 4], hn[:, 5], hn[:, 6], hn[:, 7], hn[:, 8]))
end

@testset "Mesh 生成器" begin
    for m in (generate_rectangle_mesh(1.0, 1.0, 3, 3),
              generate_sphere_mesh(0.5, 4, 8),
              generate_cylinder_mesh(0.5, 1.0, 4, 4),
              generate_ellipsoid_mesh(1.0, 0.5, 0.5, 4, 8),
              generate_torus_mesh(1.0, 0.3, 8, 6))
        @test m isa TriangleMesh
        @test num_vertices(m) > 0
        @test num_elements(m) > 0
    end
    vol = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)
    @test vol isa TetrahedraMesh
    @test num_elements(vol) > 0
end

@testset "Mesh 变换/修复/质量/边界" begin
    mesh = generate_rectangle_mesh(1.0, 1.0, 4, 4)
    t = translate_mesh(mesh, [1.0, 2.0, 3.0])
    @test num_vertices(t) == num_vertices(mesh)
    s = scale_mesh(mesh, 2.0)
    @test s isa TriangleMesh
    r = rotate_mesh(mesh, [0, 0, 1], pi / 2)
    @test r isa TriangleMesh
    tr = transform_mesh(mesh, I(3), [0.0, 0.0, 0.0])
    @test tr isa TriangleMesh
    merged = merge_meshes([mesh, mesh])
    @test num_elements(merged) == 2 * num_elements(mesh)
    # 修复
    dup = TriangleMesh(num_elements(mesh), vertices(mesh), elements(mesh), fill(1, num_elements(mesh)))
    fixed = remove_duplicate_nodes(dup)
    @test num_vertices(fixed) <= num_vertices(dup)
    @test fix_element_orientation(mesh) isa TriangleMesh
    degen = detect_degenerates(mesh)
    @test degen isa Vector
    # 质量
    q = mesh_quality(mesh)
    @test q isa MeshQualityReport
    # 边界提取（体积网格）
    vol = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)
    surf = extract_surface(vol)
    @test surf isa TriangleMesh
    @test num_elements(surf) > 0
end

@testset "坐标变换与求积" begin
    v = sphere2cart(1.0, pi / 2, 0.0)
    @test norm(v) ≈ 1.0
    info = EMMoMSuite.Geometry.r̂θϕInfo(pi / 2, 0.0)
    @test norm(info.r̂) ≈ 1.0
    r = EMMoMSuite.Geometry.r̂θϕInfo{Float64}([1.0, 0.0, 0.0])
    @test r.r̂ ≈ [1.0, 0.0, 0.0]
    # 观测坐标旋转
    obs = reshape(
        [EMMoMSuite.Geometry.r̂θϕInfo(pi / 2, 0.0), EMMoMSuite.Geometry.r̂θϕInfo(pi / 3, pi / 4)],
        1, 2,
    )
    R = SMatrix{3,3,Float64}(I)
    l = globalObs2LocalObs(obs, R)
    @test size(l) == (1, 2)
    g = localObs2GlobalObs(l, R)
    @test size(g) == (1, 2)
    # 求积规则
    for gq in (GaussQuadratureInfo(:Triangle, 4),
               GaussQuadratureInfo(:Tetrahedron, 4),
               GaussQuadratureInfo(:Hexahedron, 8),
               GaussQuadratureInfo(:Quadrangle, 4))
        @test length(gq.weight) > 0
    end
    @test length(gaussQuadratureHexa1D(3)) == 2
end

@testset "TetrahedraInfo 与 STL/Hex I/O" begin
    # STL 读写
    smesh = generate_sphere_mesh(0.5, 4, 8)
    path = joinpath(tempdir(), "cov_$(getpid()).stl")
    write_stl_mesh(path, smesh)
    m2 = read_stl_mesh(path)
    @test num_elements(m2) > 0
    rm(path; force = true)
end

@testset "BRep 实体与 CSG" begin
    b = box_solid(1.0, 1.0, 1.0)
    @test b isa BRepSolid
    @test solid_volume(b) ≈ 1.0 atol = 1e-6
    @test solid_surface_area(b) ≈ 6.0 atol = 1e-6
    @test check_manifold(b)
    tri = convert_to_triangle_mesh(b)
    @test tri isa TriangleMesh
    @test num_elements(tri) > 0
    b2 = box_solid(0.5, 0.5, 0.5)
    inter = intersect_solids(b, b2)
    @test solid_volume(inter) ≈ 0.125 atol = 1e-4
    node = CSGNode(:intersect, b, b2)
    @test csg_volume(node) ≈ 0.125 atol = 1e-4
end
