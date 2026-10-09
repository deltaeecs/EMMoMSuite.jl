# test_volume_material_binding.jl
# Volume (tetrahedral) geometry→material mapping per docs/dev/volume_material_binding.md:
# - GmshIO: Physical Volume tag remap + region name map (read_msh_volume)
# - GmshAPI: generate_gmsh_volume with Physical Volume tags
# - MeshMaterialBind: bind_regions + validate_bindings gating

using Test
using EMMoMSuite

const _FIXDIR = joinpath(@__DIR__, "fixtures")
isdir(_FIXDIR) || mkpath(_FIXDIR)

# ── synthetic .msh 4.1 with two Physical Volumes and two tet blocks ─────────
const _MSH_FULL = raw"""
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
2
3 1 "outer"
3 2 "core"
$EndPhysicalNames
$Entities
0 0 0 2
1 0 0 0 0 0 0 1 1 0
2 0 0 0 0 0 0 1 2 0
$EndEntities
$Nodes
2 6 1 6
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
3 2 0 2
5
6
0 0 2
1 1 2
$EndNodes
$Elements
2 2 1 2
3 1 4 1
1 1 2 3 4
3 2 4 1
2 5 6 2 1
$EndElements
"""

# variant: volume entity 2 belongs to NO Physical Volume → tet tag must be 0
const _MSH_UNGROUPED = raw"""
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
1
3 1 "outer"
$EndPhysicalNames
$Entities
0 0 0 2
1 0 0 0 0 0 0 1 1 0
2 0 0 0 0 0 0 0 0
$EndEntities
$Nodes
2 6 1 6
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
3 2 0 2
5
6
0 0 2
1 1 2
$EndNodes
$Elements
2 2 1 2
3 1 4 1
1 1 2 3 4
3 2 4 1
2 5 6 2 1
$EndElements
"""

function _write_fixture(name::String, content::String)
    path = joinpath(_FIXDIR, name)
    open(path, "w") do io
        print(io, content)
    end
    return path
end

@testset "volume material binding" begin
    @testset "read_msh_volume: Physical Volume remap + regions" begin
        path = _write_fixture("two_regions.msh", _MSH_FULL)
        mesh, regions = read_msh_volume(path)
        @test mesh isa TetrahedraMesh
        @test mesh.tetnum == 2
        @test regions == Dict("outer" => 1, "core" => 2)
        @test mesh.tags[1] == 1   # entity 1 → Physical Volume 1
        @test mesh.tags[2] == 2   # entity 2 → Physical Volume 2
    end

    @testset "read_msh_volume: ungrouped volume entity → tag 0" begin
        path = _write_fixture("ungrouped_region.msh", _MSH_UNGROUPED)
        mesh, regions = read_msh_volume(path)
        @test regions == Dict("outer" => 1)
        @test mesh.tags[1] == 1
        @test mesh.tags[2] == 0   # not in any Physical Volume
    end

    @testset "bind_regions + validate_bindings (full coverage)" begin
        path = _write_fixture("two_regions.msh", _MSH_FULL)
        mesh, regions = read_msh_volume(path)
        eps_out = 4.0 + 0.0im
        eps_core = 2.0 + 0.0im
        bm = bind_regions(mesh, regions, Dict{String,ComplexF64}(
            "outer" => eps_out, "core" => eps_core))
        @test bm isa BoundMesh
        @test validate_bindings(bm)
        @test element_material(bm, 1) == eps_out
        @test element_material(bm, 2) == eps_core
    end

    @testset "validate_bindings rejects unbound tags" begin
        path = _write_fixture("ungrouped_region.msh", _MSH_UNGROUPED)
        mesh, regions = read_msh_volume(path)
        bm = bind_regions(mesh, regions, Dict{String,ComplexF64}("outer" => 4.0 + 0.0im))
        @test !validate_bindings(bm)   # tag 0 has no binding
        @test_throws KeyError element_material(bm, 2)
    end

    @testset "bind_regions: unknown region errors" begin
        path = _write_fixture("two_regions.msh", _MSH_FULL)
        mesh, regions = read_msh_volume(path)
        @test_throws ErrorException bind_regions(
            mesh, regions, Dict{String,ComplexF64}("nope" => 1.0 + 0.0im))
    end

    @testset "generate_gmsh_volume: sphere Physical Volume" begin
        geo = normpath(joinpath(@__DIR__, "..", "benchmark", "cases", "geo",
                                "sphere_r0p15.geo"))
        mesh, regions = generate_gmsh_volume(geo; mesh_size = 0.15)
        @test mesh isa TetrahedraMesh
        @test mesh.tetnum > 0
        @test regions == Dict("body" => regions["body"])   # single named region
        ptag = regions["body"]
        @test ptag >= 1
        @test all(t -> t == ptag, mesh.tags)
        bm = bind_regions(mesh, regions, Dict{String,ComplexF64}("body" => 4.0 + 0.0im))
        @test validate_bindings(bm)
    end

    @testset "legacy read_msh_mesh keeps entity tags" begin
        path = _write_fixture("two_regions.msh", _MSH_FULL)
        mesh = read_msh_mesh(path)
        @test mesh isa TetrahedraMesh
        @test mesh.tags == [1, 2]   # entity tags (same here, but no remap pass)
    end

    @testset "generic build_volume_system (material protocol)" begin
        # dielectric bindings → VEFIE on the volume mesh, permittivities per tet
        path = _write_fixture("two_regions.msh", _MSH_FULL)
        mesh, regions = read_msh_volume(path)
        bm = bind_regions(mesh, regions, Dict{String,ComplexF64}(
            "outer" => 4.0 + 0.0im, "core" => 2.0 + 0.0im))
        sys = build_volume_system(bm, 3e8)
        @test sys.kind === :vefie
        @test sys.op isa VEFIE
        @test sys.mesh === bm.mesh
        @test sys.permittivities == [4.0 + 0.0im, 2.0 + 0.0im]
        @test !is_all_pec(bm)
        @test volume_permittivities(bm) == sys.permittivities

        # PEC-like bindings (custom material type via the protocol) → auto
        # fallback to surface EFIE on the extracted boundary
        struct _TestPEC end
        EMMoMSuite.is_pec_material(::_TestPEC) = true
        bm_pec = bind_regions(mesh, regions, Dict{String,_TestPEC}(
            "outer" => _TestPEC(), "core" => _TestPEC()))
        @test is_all_pec(bm_pec)
        sys_pec = build_volume_system(bm_pec, 3e8)
        @test sys_pec.kind === :efie_fallback
        @test sys_pec.op isa EFIE
        @test sys_pec.mesh isa TriangleMesh
        @test isempty(sys_pec.permittivities)
        @test_throws ErrorException build_volume_system(bm_pec, 3e8; ie = "VEFIE")

        # far-field reconstruction from RCS
        rcs = [0.0 10.0; 20.0 30.0]
        ff = farfield_from_rcs(rcs)
        @test size(ff) == (2, 2, 2)
        @test ff[1, 2, 2] ≈ sqrt(10^((30.0 + 30.0) / 10.0))
    end
end
