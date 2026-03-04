module MLFMA

include("Interpolation.jl")
include("Level.jl")
include("Octree.jl")
include("Precomputations.jl")
include("Translation.jl")
include("Disaggregation.jl")
include("OctreeBuilder.jl")
include("Aggregation.jl")
include("MLFMAOperator.jl")
include("PMCHWMLFMAOperator.jl")

using .Interpolation
using .Level
using .Octree
using .Precomputations
using .Translation
using .Disaggregation
using .OctreeBuilder
using .Aggregation
using .MLFMAOperatorModule
using .PMCHWMLFMAOperatorModule

export AbstractPolesInfo, AbstractInterpInfo
export CubeInfo, AbstractLevel, LevelInfo
export OctreeInfo
export build_octree
export aggregate!
export compute_translation_factors!, translate!
export disaggregate_downward!, disaggregate_leaf!
export MLFMAOperator, mul!, get_leaf_intervals
export MLFMAOperatorMPI
export PMCHWMLFMAOperator, assemble_near_field_pmchw

end
