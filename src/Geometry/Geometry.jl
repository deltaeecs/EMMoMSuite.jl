module Geometry

using ..CoreModule
using StaticArrays
using LinearAlgebra

include("MeshTypes.jl")
include("TriangleInfoConstructor.jl")
include("TetrahedraInfo.jl")
include("MeshIO.jl")
include("MeshGen.jl")
include("MeshBoundary.jl")
include("GmshIO.jl")
include("CoordinateTransforms.jl")
include("GaussQuadrature.jl")

export TriangleMesh, TetrahedraMesh, HexahedraMesh, TriangleInfo, TetrahedraInfo
export HexahedraInfo, Quads4Hexa
export CompositeMesh
export HEXA_FACE_VERTEX_IDS, HEXA_OPP_FACE, HEXA_OPP_FACE_VERTEX_IDS
export hex_volume, tet_volume, area
export get_free_vns, gq3d_to_face2d_idx, construct_gq3d_index_map, set_delta_kappa!
export read_nas_mesh,
    read_mixed_nas_mesh,
    write_nas_mesh,
    read_msh_mesh,
    generate_rectangle_mesh,
    generate_cylinder_mesh,
    generate_sphere_mesh,
    generate_ellipsoid_mesh,
    generate_cone_mesh,
    generate_torus_mesh,
    generate_box_volume_mesh,
    generate_box_tet_mesh,
    extract_surface
export globalObs2LocalObs, localObs2GlobalObs
export GaussQuadratureInfo, GaussQuadratureInfoStruct, get_global_quad_points
export gaussQuadratureHexa, gaussQuadratureQuad, gaussQuadratureHexa1D
export get_tetrahedra_info
export r̂θϕInfo, sphere2cart

end
