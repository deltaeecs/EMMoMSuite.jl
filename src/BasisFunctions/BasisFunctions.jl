module BasisFunctions

using ..CoreModule
using ..Geometry
using StaticArrays
using LinearAlgebra

include("RWG.jl")
include("SWG.jl")
include("PWC.jl")
include("RBF.jl")
include("BasisUtilities.jl")

export RWGBasis, RWG
export SWGBasis, SWG, evaluate_swg
export PWCBasis, PWCHexBasis, PWC
export RBFBasis, RBF
export get_triangle_info, get_triangles_info, get_tetrahedra_info, get_hexahedra_info

# Re-export utilities
# using .BasisUtilities # Removed because BasisUtilities is now included directly
# export count_unknowns, get_triangle_info # Already exported by include? No, include just defines them.
# I need to export them from BasisFunctions module.

export num_basis, support, evaluate

end
