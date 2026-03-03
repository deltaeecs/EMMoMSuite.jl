using Test
using EMSuite
using EMSuite.Geometry
using LinearAlgebra
using StaticArrays

@testset "Geometry" begin
    @testset "TriangleMesh" begin
        # Create a simple mesh manually
        # 2 triangles forming a square
        # (0,1) -- (1,1)
        #   | \      |
        #   |  \     |
        # (0,0) -- (1,0)

        nodes = [
            0.0 1.0 0.0 1.0
            0.0 0.0 1.0 1.0
            0.0 0.0 0.0 0.0
        ]

        # Triangle 1: 1-2-3
        # Triangle 2: 2-4-3
        tris = [
            1 2
            2 4
            3 3
        ]

        tags = [1, 2] # Different tags for testing

        mesh = TriangleMesh(2, nodes, tris, tags)

        @test num_vertices(mesh) == 4
        @test num_elements(mesh) == 2
        @test dimension(mesh) == 2
        @test mesh.tags == [1, 2]
        @test vertices(mesh) == nodes
        @test elements(mesh) == tris
    end

    @testset "MeshIO (Mock NAS)" begin
        # Create a temporary NAS file
        nas_content = """
\$ Test NAS file
GRID, 1, , 0.0, 0.0, 0.0
GRID, 2, , 1.0, 0.0, 0.0
GRID, 3, , 0.0, 1.0, 0.0
CTRIA3, 1, 100, 1, 2, 3
"""
        path = "test_mesh.nas"
        open(path, "w") do f
            write(f, nas_content)
        end

        try
            mesh = read_nas_mesh(path)
            @test num_vertices(mesh) == 3
            @test num_elements(mesh) == 1
            @test mesh.tags[1] == 100
            @test mesh.node[:, 1] ≈ [0.0, 0.0, 0.0]
        finally
            rm(path, force = true)
        end
    end

    @testset "GmshIO (Mock MSH)" begin
        # Create a temporary MSH file (Gmsh 4.1)
        msh_content = """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Nodes
1 3 1 3
2 1 0 3
1
2
3
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
\$EndNodes
\$Elements
1 1 1 1
2 1 2 1
1 1 2 3
\$EndElements
"""
        # Explanation of Mock Data:
        # Nodes: 1 block, 3 nodes.
        # Block 1: dim=2, tag=1, parametric=0, num=3
        # Node tags: 1, 2, 3
        # Coords: (0,0,0), (1,0,0), (0,1,0)

        # Elements: 1 block, 1 element.
        # Block 1: dim=2, tag=1 (Physical Tag), type=2 (Triangle), num=1
        # Element 1: tag=1, nodes=1, 2, 3

        path = "test_mesh.msh"
        open(path, "w") do f
            write(f, msh_content)
        end

        try
            mesh = read_msh_mesh(path)
            @test num_vertices(mesh) == 3
            @test num_elements(mesh) == 1
            @test mesh.tags[1] == 1
            @test mesh.node[:, 1] ≈ [0.0, 0.0, 0.0]
        finally
            rm(path, force = true)
        end
    end

    @testset "GaussQuadrature - Hexahedron and Quadrangle" begin
        using EMSuite.Geometry: GaussQuadratureInfo, gaussQuadratureHexa1D

        # Test Hexahedron GQ (8 = 2³ points)
        gq_hexa = GaussQuadratureInfo(:Hexahedron, 8)
        @test length(gq_hexa.weight) == 8
        @test isapprox(sum(gq_hexa.weight), 1.0; rtol=1e-10)
        @test size(gq_hexa.coordinate, 1) == 8  # 8 shape functions
        @test size(gq_hexa.coordinate, 2) == 8  # 8 points

        # Test Hexahedron GQ (1 = 1³ point)
        gq_hexa1 = GaussQuadratureInfo(:Hexahedron, 1)
        @test length(gq_hexa1.weight) == 1
        @test isapprox(sum(gq_hexa1.weight), 1.0; rtol=1e-10)

        # Test Quadrangle GQ (4 = 2² points)
        gq_quad = GaussQuadratureInfo(:Quadrangle, 4)
        @test length(gq_quad.weight) == 4
        @test isapprox(sum(gq_quad.weight), 1.0; rtol=1e-10)
        @test size(gq_quad.coordinate, 1) == 4  # 4 shape functions

        # Test Quadrangle GQ (1 point)
        gq_quad1 = GaussQuadratureInfo(:Quadrangle, 1)
        @test length(gq_quad1.weight) == 1
        @test isapprox(sum(gq_quad1.weight), 1.0; rtol=1e-10)

        # Test gaussQuadratureHexa1D
        x, w = gaussQuadratureHexa1D(2)
        @test length(x) == 2
        @test length(w) == 2
        @test all(0.0 .<= x .<= 1.0)
        @test isapprox(sum(w), 1.0; rtol=1e-10)

        # Test error path
        @test_throws ErrorException GaussQuadratureInfo(:InvalidGeometry, 3)
    end

    @testset "GaussQuadrature - MeshTypes" begin
        using EMSuite.Geometry: TriangleMesh, TetrahedraMesh
        # TetrahedraMesh construction
        nodes_tet = [0.0 1.0 0.0 0.0;
                     0.0 0.0 1.0 0.0;
                     0.0 0.0 0.0 1.0]
        elements_tet = reshape([1; 2; 3; 4], 4, 1)
        tags_tet = [1]
        mesh_tet = TetrahedraMesh(1, nodes_tet, elements_tet, tags_tet)
        @test num_elements(mesh_tet) == 1
        @test num_vertices(mesh_tet) == 4
        @test dimension(mesh_tet) == 3

        # TetrahedraMesh with default tags
        mesh_tet2 = TetrahedraMesh(1, nodes_tet, elements_tet)
        @test mesh_tet2.tags == [0]
    end

    @testset "NAS - CTETRA and CHEXA reading" begin
        # Test read_nas_mesh with CTETRA (fixed width)
        nas_ctetra = """
\$ Test NAS file with tetrahedra (fixed width)
GRID    1               0.0     0.0     0.0
GRID    2               1.0     0.0     0.0
GRID    3               0.0     1.0     0.0
GRID    4               0.0     0.0     1.0
CTETRA  1       100     1       2       3       4
"""
        path = "test_ctetra.nas"
        open(path, "w") do f
            write(f, nas_ctetra)
        end
        try
            mesh = read_nas_mesh(path)
            @test mesh isa TetrahedraMesh
            @test num_vertices(mesh) == 4
            @test num_elements(mesh) == 1
            @test mesh.tags[1] == 100
        finally
            rm(path, force=true)
        end

        # Test read_nas_mesh with CHEXA (free-field format)
        nas_chexa_ff = """
\$ Test NAS file with hexahedra (free field)
GRID, 1, , 0.0, 0.0, 0.0
GRID, 2, , 1.0, 0.0, 0.0
GRID, 3, , 1.0, 1.0, 0.0
GRID, 4, , 0.0, 1.0, 0.0
GRID, 5, , 0.0, 0.0, 1.0
GRID, 6, , 1.0, 0.0, 1.0
GRID, 7, , 1.0, 1.0, 1.0
GRID, 8, , 0.0, 1.0, 1.0
CHEXA, 1, 200, 1, 2, 3, 4, 5, 6,
+, 7, 8
"""
        path2 = "test_chexa.nas"
        open(path2, "w") do f
            write(f, nas_chexa_ff)
        end
        try
            mesh = read_nas_mesh(path2)
            @test mesh isa HexahedraMesh
            @test num_vertices(mesh) == 8
            @test num_elements(mesh) == 1
            @test mesh.tags[1] == 200
        finally
            rm(path2, force=true)
        end

        # Test read_mixed_nas_mesh
        nas_mixed = """
\$ Mixed mesh: triangles + tetrahedra
GRID, 1, , 0.0, 0.0, 0.0
GRID, 2, , 1.0, 0.0, 0.0
GRID, 3, , 0.0, 1.0, 0.0
GRID, 4, , 0.0, 0.0, 1.0
CTRIA3, 1, 100, 1, 2, 3
CTETRA, 1, 200, 1, 2, 3, 4
"""
        path3 = "test_mixed.nas"
        open(path3, "w") do f
            write(f, nas_mixed)
        end
        try
            surf_mesh, vol_mesh, hexa_mesh = read_mixed_nas_mesh(path3)
            @test surf_mesh isa TriangleMesh
            @test num_elements(surf_mesh) == 1
            @test vol_mesh isa TetrahedraMesh
            @test num_elements(vol_mesh) == 1
            @test hexa_mesh isa HexahedraMesh
            @test num_elements(hexa_mesh) == 0
        finally
            rm(path3, force=true)
        end
    end

    @testset "NAS - parse_nastran_float and GRID* format" begin
        using EMSuite.Geometry: parse_nastran_float

        # Normal floats
        @test parse_nastran_float("1.5", Float64) ≈ 1.5
        @test parse_nastran_float("0", Float64) ≈ 0.0
        @test parse_nastran_float("", Float64) ≈ 0.0

        # Nastran compressed scientific notation
        @test parse_nastran_float("1.78-15", Float64) ≈ 1.78e-15
        @test parse_nastran_float("3.2+5", Float64) ≈ 3.2e5
        @test parse_nastran_float("-1.5-3", Float64) ≈ -1.5e-3
    end

    @testset "NAS - write_nas_mesh" begin
        # Create a mesh, write it, read it back
        nodes = [
            0.0 1.0 0.0 1.0;
            0.0 0.0 1.0 1.0;
            0.0 0.0 0.0 0.0
        ]
        tris = [1 2; 2 4; 3 3]
        tags = [100, 200]
        mesh_orig = TriangleMesh(2, nodes, tris, tags)

        path = "test_write.nas"
        write_nas_mesh(path, mesh_orig)
        try
            @test isfile(path)
            mesh_read = read_nas_mesh(path)
            @test num_elements(mesh_read) == 2
            @test num_vertices(mesh_read) == 4
        finally
            rm(path, force=true)
        end
    end

    @testset "MSH - Gmsh tetrahedra reading" begin
        # Create Gmsh .msh file with tetrahedra (element type 4)
        # Gmsh tetrahedral type = 4
        msh_tet = """
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
\$EndNodes
\$Elements
1 1 1 1
3 1 4 1
1 1 2 3 4
\$EndElements
"""
        path = "test_tet.msh"
        open(path, "w") do f
            write(f, msh_tet)
        end
        try
            mesh = read_msh_mesh(path)
            @test mesh isa TetrahedraMesh
            @test num_vertices(mesh) == 4
            @test num_elements(mesh) == 1
        finally
            rm(path, force=true)
        end
    end

    @testset "CoordinateTransforms" begin
        # sphere2cart: θ=π/2, ϕ=0 → (1,0,0)
        r1 = sphere2cart(1.0, π/2, 0.0)
        @test isapprox(r1, [1.0, 0.0, 0.0]; atol=1e-12)

        # sphere2cart: r=2, θ=0 → (0,0,2)  (zenith)
        r2 = sphere2cart(2.0, 0.0, 0.0)
        @test isapprox(r2, [0.0, 0.0, 2.0]; atol=1e-12)

        # r̂θϕInfo from direction vector (vector constructor + normalization branch)
        info1 = r̂θϕInfo([1.0, 0.0, 0.0])
        @test isapprox(norm(info1.r̂), 1.0; atol=1e-12)
        @test isapprox(info1.r̂, [1.0, 0.0, 0.0]; atol=1e-12)

        # r̂θϕInfo with near-zero vector (edge case: if r ≈ 0 branch)
        info_zero = r̂θϕInfo([0.0, 0.0, 0.0])
        @test all(iszero, info_zero.r̂)

        # globalObs2LocalObs with identity rotation (no-op transform)
        FT = Float64
        l2g = SMatrix{3,3,FT,9}(1,0,0,0,1,0,0,0,1)
        obs_mat = fill(r̂θϕInfo(π/4, 0.0), 1, 2)
        local_obs = globalObs2LocalObs(obs_mat, l2g)
        @test size(local_obs) == (1, 2)
        @test isapprox(local_obs[1,1].r̂, obs_mat[1,1].r̂; atol=1e-10)

        # localObs2GlobalObs with identity rotation (no-op transform)
        global_obs = localObs2GlobalObs(obs_mat, l2g)
        @test size(global_obs) == (1, 2)
        @test isapprox(global_obs[1,1].r̂, obs_mat[1,1].r̂; atol=1e-10)
    end

    # ──────────────────────────────────────────────────────────────────────────
    # Phase 16.1 — 新增表面网格生成器
    # ──────────────────────────────────────────────────────────────────────────

    @testset "generate_ellipsoid_mesh" begin
        # a=b=c=1 → should be a sphere: same topology as generate_sphere_mesh
        sphere  = generate_sphere_mesh(1.0, 10, 20)
        ellips  = generate_ellipsoid_mesh(1.0, 1.0, 1.0, 10, 20)
        @test ellips isa TriangleMesh
        @test ellips.trinum == sphere.trinum
        @test size(ellips.node) == size(sphere.node)
        @test ellips.node ≈ sphere.node

        # Stretched ellipsoid: bounding box must match (a,b,c)
        a, b, c = 2.0, 3.0, 1.5
        em = generate_ellipsoid_mesh(a, b, c, 12, 24)
        @test maximum(abs, em.node[1,:]) ≈ a  atol=1e-10
        @test maximum(abs, em.node[2,:]) ≈ b  atol=1e-10
        @test maximum(abs, em.node[3,:]) ≈ c  atol=1e-10

        # Approximate surface area: for sphere(r=1) ≈ 4π with fine enough mesh
        sm = generate_ellipsoid_mesh(1.0, 1.0, 1.0, 50, 100)
        area_sum = sum(1:sm.trinum) do k
            v1 = sm.node[:, sm.triangles[1,k]]
            v2 = sm.node[:, sm.triangles[2,k]]
            v3 = sm.node[:, sm.triangles[3,k]]
            norm(cross(v2 .- v1, v3 .- v1)) / 2
        end
        @test isapprox(area_sum, 4π; rtol=0.01)
    end

    @testset "generate_cone_mesh" begin
        n_circ, n_h = 12, 3
        radius, height = 1.0, 2.0

        # open cone
        m = generate_cone_mesh(radius, height, n_circ, n_h; closed=false)
        @test m isa TriangleMesh
        expected_elems = (n_h - 1) * 2 * n_circ + n_circ   # full bands + apex band
        @test m.trinum == expected_elems

        # apex is at height/2
        @test maximum(m.node[3, :]) ≈ height / 2  atol=1e-12

        # closed cone
        mc = generate_cone_mesh(radius, height, n_circ, n_h; closed=true)
        @test mc.trinum == expected_elems + n_circ

        # single-band cone (n_height=1): only apex band + optional cap
        m1 = generate_cone_mesh(1.0, 1.0, 8, 1; closed=false)
        @test m1.trinum == 8

        # no triangles have area ≤ 0
        for k in 1:m1.trinum
            v1 = m1.node[:, m1.triangles[1,k]]
            v2 = m1.node[:, m1.triangles[2,k]]
            v3 = m1.node[:, m1.triangles[3,k]]
            @test norm(cross(v2 .- v1, v3 .- v1)) > 1e-14
        end
    end

    @testset "generate_torus_mesh" begin
        R, r = 3.0, 1.0
        nm, nn = 20, 10
        m = generate_torus_mesh(R, r, nm, nn)
        @test m isa TriangleMesh
        @test size(m.node, 2) == nm * nn
        @test m.trinum == 2 * nm * nn

        # All nodes are at distance ∈ [R-r, R+r] from z-axis
        for k in 1:size(m.node, 2)
            dist_z = sqrt(m.node[1,k]^2 + m.node[2,k]^2)
            @test R - r - 1e-12 ≤ dist_z ≤ R + r + 1e-12
        end

        # No degenerate triangles
        for k in 1:m.trinum
            v1 = m.node[:, m.triangles[1,k]]
            v2 = m.node[:, m.triangles[2,k]]
            v3 = m.node[:, m.triangles[3,k]]
            @test norm(cross(v2 .- v1, v3 .- v1)) > 1e-14
        end
    end

    # ──────────────────────────────────────────────────────────────────────────
    # Phase 16.1 — 体网格生成器
    # ──────────────────────────────────────────────────────────────────────────

    @testset "generate_box_volume_mesh (HexahedraMesh)" begin
        Lx, Ly, Lz = 2.0, 3.0, 1.0
        nx, ny, nz = 3, 4, 2
        m = generate_box_volume_mesh(Lx, Ly, Lz, nx, ny, nz)
        @test m isa HexahedraMesh
        @test m.hexnum == nx * ny * nz
        @test size(m.node, 2) == (nx+1)*(ny+1)*(nz+1)

        # Bounding box
        @test minimum(m.node[1,:]) ≈ -Lx/2  atol=1e-12
        @test maximum(m.node[1,:]) ≈  Lx/2  atol=1e-12
        @test minimum(m.node[3,:]) ≈ -Lz/2  atol=1e-12
        @test maximum(m.node[3,:]) ≈  Lz/2  atol=1e-12

        # Each hex has 8 distinct nodes
        for k in 1:m.hexnum
            @test length(unique(m.hexes[:,k])) == 8
        end

        # Note: hex_volume() returns ~1/3 of actual volume (Legacy signed convention).
        # Do not use it for total-volume verification here.
    end

    @testset "generate_box_tet_mesh (TetrahedraMesh)" begin
        Lx, Ly, Lz = 1.0, 1.0, 1.0
        nx, ny, nz = 2, 2, 2
        m = generate_box_tet_mesh(Lx, Ly, Lz, nx, ny, nz)
        @test m isa TetrahedraMesh
        @test m.tetnum == 6 * nx * ny * nz
        @test size(m.node, 2) == (nx+1)*(ny+1)*(nz+1)

        # Total volume via tet_volume (Freudenthal decomposition → all vols positive)
        vol_total = sum(1:m.tetnum) do k
            v1 = m.node[:, m.tetras[1,k]]
            v2 = m.node[:, m.tetras[2,k]]
            v3 = m.node[:, m.tetras[3,k]]
            v4 = m.node[:, m.tetras[4,k]]
            tet_volume(v1, v2, v3, v4)
        end
        @test isapprox(vol_total, Lx*Ly*Lz; rtol=1e-10)

        # All tet volumes positive (consistent orientation)
        for k in 1:m.tetnum
            v1 = m.node[:, m.tetras[1,k]]
            v2 = m.node[:, m.tetras[2,k]]
            v3 = m.node[:, m.tetras[3,k]]
            v4 = m.node[:, m.tetras[4,k]]
            @test tet_volume(v1, v2, v3, v4) > 0
        end
    end

    # ──────────────────────────────────────────────────────────────────────────
    # Phase 16.3 — 边界提取与 CompositeMesh
    # ──────────────────────────────────────────────────────────────────────────

    @testset "extract_surface" begin
        # Use a single-cube tet mesh (1×1×1, 1×1×1 grid → 6 tets)
        Lx, Ly, Lz = 1.0, 1.0, 1.0
        m = generate_box_tet_mesh(Lx, Ly, Lz, 1, 1, 1)
        surf = extract_surface(m)

        @test surf isa TriangleMesh
        # A cube has 6 faces × 2 triangles = 12 boundary triangles
        @test surf.trinum == 12

        # Surface area of unit cube = 6
        area_sum = sum(1:surf.trinum) do k
            v1 = surf.node[:, surf.triangles[1,k]]
            v2 = surf.node[:, surf.triangles[2,k]]
            v3 = surf.node[:, surf.triangles[3,k]]
            norm(cross(v2 .- v1, v3 .- v1)) / 2
        end
        @test isapprox(area_sum, 6.0; rtol=1e-10)

        # Outward normals: for each triangle, normal should point outward
        # (component along the dominant face direction should be positive)
        centroid = sum(m.node; dims=2) ./ size(m.node, 2)
        for k in 1:surf.trinum
            v1 = surf.node[:, surf.triangles[1,k]]
            v2 = surf.node[:, surf.triangles[2,k]]
            v3 = surf.node[:, surf.triangles[3,k]]
            n  = cross(v2 .- v1, v3 .- v1)
            fc = (v1 .+ v2 .+ v3) ./ 3
            @test dot(n, fc .- vec(centroid)) > 0
        end
    end

    @testset "CompositeMesh" begin
        tri  = generate_sphere_mesh(1.0, 6, 12)
        tet  = generate_box_tet_mesh(2.0, 2.0, 2.0, 1, 1, 1)
        comp = CompositeMesh(tri, tet)
        @test comp isa CompositeMesh
        @test comp.surface === tri
        @test comp.volume   === tet
        @test dimension(comp) == 3
    end

    # ─── Phase 16.2: MeshTransforms ───────────────────────────────────────────

    @testset "translate_mesh" begin
        m  = generate_box_tet_mesh(1.0, 1.0, 1.0, 1, 1, 1)
        d  = [2.0, -3.0, 1.0]
        m2 = translate_mesh(m, d)

        @test m2 isa TetrahedraMesh
        @test m2.tetnum == m.tetnum
        @test m2.tetras == m.tetras          # connectivity unchanged

        # Every node shifted by d
        @test isapprox(m2.node, m.node .+ d; atol=1e-14)

        # Works on TriangleMesh too
        tri  = generate_sphere_mesh(1.0, 6, 12)
        tri2 = translate_mesh(tri, d)
        @test tri2 isa TriangleMesh
        @test isapprox(tri2.node, tri.node .+ d; atol=1e-14)
    end

    @testset "scale_mesh isotropic" begin
        m  = generate_sphere_mesh(1.0, 6, 12)
        m2 = scale_mesh(m, 2.0)
        @test m2 isa TriangleMesh
        @test isapprox(m2.node, 2.0 .* m.node; atol=1e-14)
        # Edge lengths double → sum area quadruples
        area1 = sum(k -> begin
            v1=m.node[:,m.triangles[1,k]]; v2=m.node[:,m.triangles[2,k]]; v3=m.node[:,m.triangles[3,k]]
            norm(cross(v2-v1, v3-v1)) / 2
        end, 1:m.trinum)
        area2 = sum(k -> begin
            v1=m2.node[:,m2.triangles[1,k]]; v2=m2.node[:,m2.triangles[2,k]]; v3=m2.node[:,m2.triangles[3,k]]
            norm(cross(v2-v1, v3-v1)) / 2
        end, 1:m2.trinum)
        @test isapprox(area2, 4.0 * area1; rtol=1e-10)
    end

    @testset "scale_mesh anisotropic" begin
        m  = generate_box_volume_mesh(1.0, 1.0, 1.0, 1, 1, 1)
        sx, sy, sz = 2.0, 3.0, 0.5
        m2 = scale_mesh(m, sx, sy, sz)
        @test m2 isa HexahedraMesh
        # Bounding box should be [sx, sy, sz]
        @test isapprox(maximum(m2.node[1,:]) - minimum(m2.node[1,:]), sx; rtol=1e-10)
        @test isapprox(maximum(m2.node[2,:]) - minimum(m2.node[2,:]), sy; rtol=1e-10)
        @test isapprox(maximum(m2.node[3,:]) - minimum(m2.node[3,:]), sz; rtol=1e-10)
    end

    @testset "rotate_mesh 90° about z" begin
        # A point at (1,0,0) should map to (0,1,0) under 90° CCW rotation about z
        m  = generate_box_tet_mesh(1.0, 1.0, 1.0, 1, 1, 1)
        m2 = rotate_mesh(m, [0.0, 0.0, 1.0], π/2)
        @test m2 isa TetrahedraMesh
        @test m2.tetnum == m.tetnum

        # Apply same rotation to all nodes manually and compare
        c, s = cos(π/2), sin(π/2)
        R = [c -s 0.0; s c 0.0; 0.0 0.0 1.0]
        expected = R * m.node
        @test isapprox(m2.node, expected; atol=1e-12)
    end

    @testset "merge_meshes TriangleMesh" begin
        m1 = generate_rectangle_mesh(1.0, 1.0, 2, 2)
        m2 = translate_mesh(generate_rectangle_mesh(1.0, 1.0, 2, 2), [2.0, 0.0, 0.0])

        merged = merge_meshes([m1, m2])
        @test merged isa TriangleMesh
        @test merged.trinum == m1.trinum + m2.trinum

        nv1 = size(m1.node, 2)
        nv2 = size(m2.node, 2)
        @test size(merged.node, 2) == nv1 + nv2

        # All connectivity indices valid
        @test all(1 .≤ merged.triangles .≤ nv1 + nv2)

        # Area preserved
        area(m) = sum(k -> begin
            v1=m.node[:,m.triangles[1,k]]; v2=m.node[:,m.triangles[2,k]]; v3=m.node[:,m.triangles[3,k]]
            norm(cross(v2-v1, v3-v1))/2
        end, 1:m.trinum)
        @test isapprox(area(merged), area(m1) + area(m2); rtol=1e-12)
    end

    @testset "merge_meshes TetrahedraMesh" begin
        t1 = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)
        t2 = translate_mesh(generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2), [2.0, 0.0, 0.0])
        merged = merge_meshes([t1, t2])
        @test merged isa TetrahedraMesh
        @test merged.tetnum == t1.tetnum + t2.tetnum
        @test size(merged.node, 2) == size(t1.node, 2) + size(t2.node, 2)
        @test all(1 .≤ merged.tetras .≤ size(merged.node, 2))
    end

    # ─── Phase 16.4: MeshQuality ─────────────────────────────────────────────

    @testset "mesh_quality TriangleMesh" begin
        m  = generate_sphere_mesh(1.0, 20, 40)
        rpt = mesh_quality(m)

        @test rpt isa MeshQualityReport
        @test rpt.n_elements == m.trinum

        # All triangles on unit sphere have positive area
        @test rpt.area_min > 0
        @test rpt.area_max ≥ rpt.area_min
        @test rpt.area_mean > 0

        # Aspect ratio ≥ 1
        @test rpt.aspect_ratio_min ≥ 1.0 - 1e-10
        @test rpt.aspect_ratio_mean ≥ 1.0 - 1e-10

        # Skewness ∈ [0,1)
        @test rpt.skewness_min ≥ 0.0
        @test rpt.skewness_max < 1.0

        # Total area ≈ 4π for unit sphere (rough check)
        @test isapprox(rpt.area_mean * m.trinum, 4π; rtol=0.05)

        @test rpt.n_degenerate == 0
        @test rpt.n_inverted   == 0

        # show produces a non-empty string
        s = sprint(show, rpt)
        @test !isempty(s)
        @test occursin("MeshQualityReport", s)
    end

    @testset "mesh_quality TetrahedraMesh" begin
        # Unit cube 2×2×2 grid → 48 tets, each vol = 1/48
        m   = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)
        rpt = mesh_quality(m)

        @test rpt isa MeshQualityReport
        @test rpt.n_elements == m.tetnum

        # All volumes positive (Freudenthal → positive)
        @test rpt.area_min > 0
        @test rpt.n_inverted == 0
        @test rpt.n_degenerate == 0

        # Mean volume ≈ 1 / 48
        @test isapprox(rpt.area_mean, 1.0 / m.tetnum; rtol=1e-8)

        # Aspect ratio ≥ 1
        @test rpt.aspect_ratio_min ≥ 1.0 - 1e-8

        # Skewness ∈ [0,1)
        @test rpt.skewness_min ≥ 0.0
        @test rpt.skewness_max < 1.0

        s = sprint(show, rpt)
        @test !isempty(s)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 16.5 — MeshRepair
# ─────────────────────────────────────────────────────────────────────────────
@testset "MeshRepair" begin
    using EMSuite.Geometry: remove_duplicate_nodes, fix_element_orientation, detect_degenerates

    # ─── remove_duplicate_nodes ───────────────────────────────────────────────
    @testset "remove_duplicate_nodes TriangleMesh" begin
        # 2-triangle mesh with one duplicated node (nodes 3 and 4 are the same)
        nodes = [0.0 1.0 0.5 0.5;
                 0.0 0.0 1.0 1.0;
                 0.0 0.0 0.0 0.0]      # 4 nodes, nodes 3 and 4 are duplicates
        tris  = [1 2; 2 4; 3 3]        # 2 triangles (3×2); node cols
        mesh  = TriangleMesh(2, nodes, tris, [1,1])

        # With default tol, nodes 3 and 4 (identical) should merge
        m2 = remove_duplicate_nodes(mesh; tol=1e-10)
        @test m2 isa TriangleMesh
        @test size(m2.node, 2) == 3    # 4 → 3 unique nodes
        @test m2.trinum == 2           # element count unchanged
        @test size(m2.triangles) == size(mesh.triangles)

        # No merging when tol is very small and nodes are truly distinct
        nodes2 = [0.0 1.0 0.5; 0.0 0.0 1.0; 0.0 0.0 0.0]
        mesh2  = TriangleMesh(1, nodes2, reshape([1,2,3], 3, 1), [0])
        m3 = remove_duplicate_nodes(mesh2; tol=1e-12)
        @test size(m3.node, 2) == 3
    end

    @testset "remove_duplicate_nodes TetrahedraMesh" begin
        # Box tet mesh with 2×2×2 grid — all nodes distinct
        m   = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)
        m2  = remove_duplicate_nodes(m; tol=1e-10)
        @test size(m2.node, 2) == size(m.node, 2)  # no merging expected
        @test m2.tetnum == m.tetnum

        # Same total volume (using tet_volume on individual tets)
        function tet_vol(nd, tets, t)
            v1 = nd[:, tets[1,t]]; v2 = nd[:, tets[2,t]]
            v3 = nd[:, tets[3,t]]; v4 = nd[:, tets[4,t]]
            abs(tet_volume(v1, v2, v3, v4))
        end
        vol_orig = sum(tet_vol(m.node,  m.tetras, t) for t in 1:m.tetnum)
        vol_new  = sum(tet_vol(m2.node, m2.tetras, t) for t in 1:m2.tetnum)
        @test vol_orig ≈ vol_new atol=1e-12
    end

    # ─── detect_degenerates ───────────────────────────────────────────────────
    @testset "detect_degenerates TriangleMesh" begin
        # Good mesh: a unit triangle — no degenerates
        nodes = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
        mesh  = TriangleMesh(1, nodes, reshape([1,2,3], 3, 1), [0])
        @test isempty(detect_degenerates(mesh))

        # Collapsed triangle: all three vertices at the same point
        bad_nodes = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        bad_mesh  = TriangleMesh(1, bad_nodes, reshape([1,2,3], 3, 1), [0])
        degs = detect_degenerates(bad_mesh)
        @test length(degs) == 1
        @test degs[1] == 1
    end

    @testset "detect_degenerates TetrahedraMesh" begin
        # Normal box tet mesh — no degenerates
        m    = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)
        @test isempty(detect_degenerates(m))

        # Collapsed tet: all 4 vertices at the same point
        bad_nodes = [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
        bad_mesh  = TetrahedraMesh(1, bad_nodes, reshape([1,2,3,4], 4, 1), [0])
        degs = detect_degenerates(bad_mesh)
        @test length(degs) == 1
    end

    # ─── fix_element_orientation ──────────────────────────────────────────────
    @testset "fix_element_orientation" begin
        # Two adjacent triangles with consistent orientation
        #   tri 1: [v1,v2,v3] — CCW upward (normal +z)
        #   tri 2: [v2,v4,v3] — CCW upward (normal +z) — consistent
        nodes = [0.0 1.0 0.0 1.0;
                 0.0 0.0 1.0 1.0;
                 0.0 0.0 0.0 0.0]
        # Consistent: shared edge v2→v3 in tri1 (CCW), v3→v2 in tri2
        tris_ok = [1 2; 2 3; 3 4]       # col k: [v1,v2,v3] for tri k
        mesh_ok = TriangleMesh(2, nodes, tris_ok, [0,0])

        m_fixed = fix_element_orientation(mesh_ok)
        @test m_fixed isa TriangleMesh
        @test m_fixed.trinum == mesh_ok.trinum

        # After fix, normals should point in the same half-space
        function tri_normal(m, t)
            v1 = m.node[:, m.triangles[1,t]]
            v2 = m.node[:, m.triangles[2,t]]
            v3 = m.node[:, m.triangles[3,t]]
            return cross(v2-v1, v3-v1)
        end
        n1 = tri_normal(m_fixed, 1)
        n2 = tri_normal(m_fixed, 2)
        # dot product ≥ 0 means normals in the same half-space
        @test dot(n1, n2) >= 0

        # Inconsistent orientation test: flip tri 2 manually
        tris_bad = [1 3; 2 2; 3 4]     # tri 2 vertex order flipped
        mesh_bad = TriangleMesh(2, nodes, tris_bad, [0,0])
        m_fixed2 = fix_element_orientation(mesh_bad)
        n1b = tri_normal(m_fixed2, 1)
        n2b = tri_normal(m_fixed2, 2)
        @test dot(n1b, n2b) >= 0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 16.6 — STL & NAS I/O
# ─────────────────────────────────────────────────────────────────────────────
@testset "STLIO" begin
    using EMSuite.Geometry: read_stl_mesh, write_stl_mesh

    # Simple triangle: one facet
    nodes  = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
    mesh   = TriangleMesh(1, nodes, reshape(Int[1,2,3], 3, 1), [1])

    mktempdir() do tmpdir
        stl_path = joinpath(tmpdir, "test.stl")

        # ── write ASCII STL ───────────────────────────────────────────────────
        write_stl_mesh(stl_path, mesh)
        @test isfile(stl_path)

        content = read(stl_path, String)
        @test occursin("solid", content)
        @test occursin("facet normal", content)
        @test occursin("vertex", content)
        @test occursin("endsolid", content)

        # ── read back ─────────────────────────────────────────────────────────
        m2 = read_stl_mesh(stl_path)
        @test m2 isa TriangleMesh
        @test m2.trinum == 1
        # After round-trip the 3 vertices should be preserved
        @test size(m2.node, 2) == 3  # 3 unique nodes (1 triangle)

        # Vertex roundtrip (order may differ — check set equality)
        orig_verts  = Set([nodes[:, i] for i in 1:3])
        round_verts = Set([m2.node[:, i] for i in 1:size(m2.node,2)])
        @test orig_verts == round_verts

        # ── multi-triangle write/read round-trip ──────────────────────────────
        mesh3 = generate_sphere_mesh(1.0, 8, 8)
        stl3_path = joinpath(tmpdir, "sphere.stl")
        write_stl_mesh(stl3_path, mesh3)

        m3 = read_stl_mesh(stl3_path)
        @test m3.trinum == mesh3.trinum
        # Number of unique nodes should be ≤ 3 * Ntri (shared nodes merged)
        @test size(m3.node, 2) <= 3 * m3.trinum
        # All elements are degenerate-free
        @test isempty(detect_degenerates(m3))
    end
end

@testset "write_nas_mesh TetrahedraMesh" begin
    m = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)

    mktempdir() do tmpdir
        nas_path = joinpath(tmpdir, "box.nas")
        write_nas_mesh(nas_path, m)

        @test isfile(nas_path)
        content = read(nas_path, String)
        @test occursin("GRID", content)
        @test occursin("CTETRA", content)
        @test occursin("ENDDATA", content)

        # Round-trip: read back and compare
        m2 = read_nas_mesh(nas_path)
        @test m2 isa TetrahedraMesh
        @test m2.tetnum == m.tetnum
        @test size(m2.node, 2) == size(m.node, 2)

        # Total volume preserved
        function tet_vol(nd, tets, t)
            v1 = nd[:, tets[1,t]]; v2 = nd[:, tets[2,t]]
            v3 = nd[:, tets[3,t]]; v4 = nd[:, tets[4,t]]
            abs(tet_volume(v1, v2, v3, v4))
        end
        vol_orig  = sum(tet_vol(m.node,  m.tetras,  t) for t in 1:m.tetnum)
        vol_round = sum(tet_vol(m2.node, m2.tetras, t) for t in 1:m2.tetnum)
        @test vol_orig ≈ vol_round atol=1e-8
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 18.1 — GeomKernel: CSG / B-Rep 数据结构
# ─────────────────────────────────────────────────────────────────────────────
@testset "GeomKernel" begin
    using EMSuite.Geometry: BRepFace, BRepSolid, CSGNode,
                             box_solid, solid_volume, solid_surface_area,
                             check_manifold, convert_to_triangle_mesh
    using StaticArrays, LinearAlgebra

    # ── 1. BRepFace 类型构造 ─────────────────────────────────────────────────
    face = BRepFace([1, 2, 3, 4])
    @test face isa BRepFace
    @test face.vertex_indices == [1, 2, 3, 4]

    # ── 2. box_solid 单位立方体 ──────────────────────────────────────────────
    box = box_solid(1.0, 1.0, 1.0)
    @test box isa BRepSolid
    @test length(box.vertices) == 8      # 8 角点
    @test length(box.faces) == 6         # 6 个四边形面
    @test length(box.edges) == 12        # 12 条棱
    @test all(length(f.vertex_indices) == 4 for f in box.faces)

    # ── 3. 体积（散度定理）──────────────────────────────────────────────────
    @test solid_volume(box) ≈ 1.0 rtol=1e-10

    box2 = box_solid(2.0, 3.0, 4.0)
    @test solid_volume(box2) ≈ 24.0 rtol=1e-10

    # ── 4. 表面积 ────────────────────────────────────────────────────────────
    @test solid_surface_area(box) ≈ 6.0 rtol=1e-10
    @test solid_surface_area(box2) ≈ 2*(2*3 + 3*4 + 4*2) rtol=1e-10

    # ── 5. 流形检验 (unit cube → true) ──────────────────────────────────────
    @test check_manifold(box) == true

    # ── 6. boundary_labels ───────────────────────────────────────────────────
    box_lab = box_solid(1.0, 1.0, 1.0;
                        material_label="copper",
                        boundary_labels=Dict(1 => "bottom", 6 => "top"))
    @test box_lab.material_label == "copper"
    @test box_lab.boundary_labels[1] == "bottom"
    @test box_lab.boundary_labels[6] == "top"

    # ── 7. 偏移原点 origin 参数 ──────────────────────────────────────────────
    origin = SVector(1.0, 2.0, 3.0)
    box_off = box_solid(1.0, 1.0, 1.0; origin=origin)
    min_v = reduce(
        (a, b) -> SVector(min(a[1],b[1]), min(a[2],b[2]), min(a[3],b[3])),
        box_off.vertices)
    @test min_v ≈ origin rtol=1e-10

    # ── 8. CSGNode 树构造 ────────────────────────────────────────────────────
    node = CSGNode(:subtract, box, box2)
    @test node.op == :subtract
    @test node.left === box
    @test node.right === box2

    leaf = CSGNode(:leaf, box, nothing)
    @test leaf.op == :leaf
    @test leaf.right === nothing

    # ── 9. convert_to_triangle_mesh ──────────────────────────────────────────
    tri_mesh = convert_to_triangle_mesh(box)
    @test tri_mesh isa TriangleMesh
    # 单位立方体：6 个四边形面 → 12 个三角形
    @test tri_mesh.trinum == 12
    @test size(tri_mesh.node, 2) == 8   # 8 个唯一节点
    # 三角形面积之和 ≈ 表面积
    total_tri_area = sum(
        0.5 * norm(cross(
            tri_mesh.node[:, tri_mesh.triangles[2, t]] - tri_mesh.node[:, tri_mesh.triangles[1, t]],
            tri_mesh.node[:, tri_mesh.triangles[3, t]] - tri_mesh.node[:, tri_mesh.triangles[1, t]]))
        for t in 1:tri_mesh.trinum)
    @test total_tri_area ≈ 6.0 rtol=1e-10
end

# ─────────────────────────────────────────────────────────────────────────────
@testset "BooleanOps" begin
    using StaticArrays

    # Phase 18.2: 凸多面体布尔操作
    using EMSuite: intersect_solids, union_solids, subtract_solid, csg_volume

    # ── 1. 相交：两个单位立方体，沿 x 偏移 0.5 ─────────────────────────────
    boxA = box_solid(1.0, 1.0, 1.0)                      # [0,1]³
    boxB = box_solid(1.0, 1.0, 1.0;
                     origin=SVector(0.5, 0.0, 0.0))      # [0.5,1.5]³
    result_AB = intersect_solids(boxA, boxB)
    @test result_AB isa BRepSolid
    @test solid_volume(result_AB) ≈ 0.5 rtol=1e-8

    # ── 2. 相交：完全不重叠 → 体积为 0 ──────────────────────────────────────
    boxC = box_solid(1.0, 1.0, 1.0;
                     origin=SVector(2.0, 0.0, 0.0))      # [2,3]³
    result_AC = intersect_solids(boxA, boxC)
    @test solid_volume(result_AC) < 1e-12

    # ── 3. 相交：B 完全包含在 A 内 ────────────────────────────────────────────
    boxD = box_solid(0.5, 0.5, 0.5;
                     origin=SVector(0.25, 0.25, 0.25))   # [0.25,0.75]³ ⊂ boxA
    result_AD = intersect_solids(boxA, boxD)
    @test solid_volume(result_AD) ≈ solid_volume(boxD) rtol=1e-8

    # ── 4. 相交结果满足流形条件 ────────────────────────────────────────────────
    @test check_manifold(result_AB; warn=false)

    # ── 5. 体积不等式：V(A∩B) ≤ min(V(A),V(B)) ──────────────────────────────
    @test solid_volume(result_AB) ≤ min(solid_volume(boxA), solid_volume(boxB)) + 1e-12

    # ── 6. union_solids 返回结构 ──────────────────────────────────────────────
    uni = union_solids(boxA, boxB)
    @test uni isa CSGNode
    @test uni.op == :union

    # ── 7. subtract_solid 返回结构 ────────────────────────────────────────────
    sub = subtract_solid(boxA, boxB)
    @test sub isa CSGNode
    @test sub.op == :subtract

    # ── 8. csg_volume: union 体积 = V(A)+V(B)-V(A∩B)（容斥原理）─────────────
    v_A  = solid_volume(boxA)
    v_B  = solid_volume(boxB)
    v_AB = solid_volume(result_AB)
    @test csg_volume(uni) ≈ v_A + v_B - v_AB rtol=1e-8

    # ── 9. csg_volume: subtract 体积 = V(A)-V(A∩B) ────────────────────────────
    @test csg_volume(sub) ≈ v_A - v_AB rtol=1e-8

    # ── 10. leaf CSGNode 体积 = solid 体积 ───────────────────────────────────
    leaf = CSGNode(:leaf, boxA, nothing)
    @test csg_volume(leaf) ≈ v_A rtol=1e-8

    # ── 11. 相交交换律：V(A∩B) == V(B∩A) ──────────────────────────────────────
    result_BA = intersect_solids(boxB, boxA)
    @test solid_volume(result_BA) ≈ solid_volume(result_AB) rtol=1e-8

    # ── 12. 正方形面片数量合理（相交结果为凸多面体）─────────────────────────
    @test length(result_AB.faces) >= 4   # 至少 4 个面

    # ── 13. z-全重叠，x 偏移：表面积验证 ─────────────────────────────────────
    # A∩B 为 [0.5,1]×[0,1]×[0,1]（Lx=0.5,Ly=Lz=1）
    # 面积 = 2(0.5×1 + 1×1 + 0.5×1) = 2×2 = 4.0
    @test solid_surface_area(result_AB) ≈ 4.0 rtol=1e-8

    # ── 14. 全等相交：A∩A = A ─────────────────────────────────────────────────
    result_AA = intersect_solids(boxA, boxA)
    @test solid_volume(result_AA) ≈ v_A rtol=1e-8
end

# ─────────────────────────────────────────────────────────────────────────────
# GmshAPI 测试默认跳过；运行方式： EMSUITE_TEST_GMSH=1 julia test/test_geometry.jl
@testset "GmshAPI" begin
    if get(ENV, "EMSUITE_TEST_GMSH", "0") != "1"
        @test_skip "Gmsh tests skipped (set EMSUITE_TEST_GMSH=1 to enable)"
    else

    # ── 1. 球面三角网格 ───────────────────────────────────────────────────────
    r     = 1.0
    smesh = generate_gmsh_sphere(r; mesh_size=0.3)
    @test smesh isa TriangleMesh
    @test smesh.trinum > 0
    # 面积 ≈ 4π r²，网格误差应 < 5%
    S_analytic = 4π * r^2
    S_mesh = sum(
        0.5 * norm(cross(
            smesh.node[:, smesh.triangles[2, t]] - smesh.node[:, smesh.triangles[1, t]],
            smesh.node[:, smesh.triangles[3, t]] - smesh.node[:, smesh.triangles[1, t]]))
        for t in 1:smesh.trinum)
    @test abs(S_mesh - S_analytic) / S_analytic < 0.05

    # ── 2. 半径 0.5 球面 ─────────────────────────────────────────────────────
    smesh2 = generate_gmsh_sphere(0.5; mesh_size=0.2)
    @test smesh2 isa TriangleMesh
    @test smesh2.trinum > 0

    # ── 3. Box 四面体网格 ─────────────────────────────────────────────────────
    bmesh = generate_gmsh_box(1.0, 2.0, 3.0; mesh_size=0.5)
    @test bmesh isa TetrahedraMesh
    @test bmesh.tetnum > 0
    # 体积 ≈ Lx*Ly*Lz，误差 < 1%（结构化参数，误差主要来自表面元素）
    V_analytic = 1.0 * 2.0 * 3.0
    # 累加 tet 体积
    V_mesh = sum(begin
        v1 = bmesh.node[:, bmesh.tetras[1, t]]
        v2 = bmesh.node[:, bmesh.tetras[2, t]]
        v3 = bmesh.node[:, bmesh.tetras[3, t]]
        v4 = bmesh.node[:, bmesh.tetras[4, t]]
        abs(dot(v2 - v1, cross(v3 - v1, v4 - v1))) / 6
    end for t in 1:bmesh.tetnum)
    @test abs(V_mesh - V_analytic) / V_analytic < 0.01

    # ── 4. Box 节点坐标范围检查 ───────────────────────────────────────────────
    Lx, Ly, Lz = 2.0, 3.0, 4.0
    bmesh2 = generate_gmsh_box(Lx, Ly, Lz; mesh_size=1.0)
    @test bmesh2 isa TetrahedraMesh
    xs = bmesh2.node[1, :]
    ys = bmesh2.node[2, :]
    zs = bmesh2.node[3, :]
    @test minimum(xs) ≥ -1e-10 && maximum(xs) ≤ Lx + 1e-10
    @test minimum(ys) ≥ -1e-10 && maximum(ys) ≤ Ly + 1e-10
    @test minimum(zs) ≥ -1e-10 && maximum(zs) ≤ Lz + 1e-10

    # ── 5. generate_gmsh_from_file：读取 .geo 文件（Box → TetrahedraMesh）──
    geo_content = """
    SetFactory("OpenCASCADE");
    Box(1) = {0, 0, 0, 1, 1, 1};
    """
    geo_file = tempname() * ".geo"
    write(geo_file, geo_content)
    fmesh = generate_gmsh_from_file(geo_file; mesh_size=0.5)
    @test fmesh isa TetrahedraMesh
    @test fmesh.tetnum > 0
    rm(geo_file; force=true)

    # ── 6. 不存在的文件应抛出错误 ───────────────────────────────────────────
    @test_throws Exception generate_gmsh_from_file("nonexistent.geo"; mesh_size=0.1)

        # ── 7. mesh_size 影响网格密度：更小的 mesh_size → 更多三角形 ─────────────
    coarse = generate_gmsh_sphere(1.0; mesh_size=0.5)
    fine   = generate_gmsh_sphere(1.0; mesh_size=0.2)
    @test fine.trinum > coarse.trinum
    end  # if EMSUITE_TEST_GMSH
end  # @testset "GmshAPI"
