module Geometry

using ..CoreModule
using StaticArrays
using LinearAlgebra

include("MeshTypes.jl")
include("TriangleInfoConstructor.jl")
include("TetrahedraInfo.jl")
include("MeshIO.jl")
include("MeshGen.jl")
include("GmshIO.jl")
include("CoordinateTransforms.jl")
include("GaussQuadrature.jl")

export TriangleMesh, TetrahedraMesh, HexahedraMesh, TriangleInfo, TetrahedraInfo
export read_nas_mesh, read_mixed_nas_mesh, write_nas_mesh, read_msh_mesh, generate_rectangle_mesh, generate_cylinder_mesh, generate_sphere_mesh
export globalObs2LocalObs, localObs2GlobalObs
export GaussQuadratureInfo, GaussQuadratureInfoStruct, get_global_quad_points
export get_tetrahedra_info
export r̂θϕInfo, sphere2cart

end
