"""
    MLFMA

多层快速多极算法（Multilevel Fast Multipole Algorithm）实现模块：
包含八叉树构建（`Octree` / `OctreeBuilder`）、聚合-转移-散播
（`Aggregation` / `Translation` / `Disaggregation`）、插值（`Interpolation`）、
MLFMA 算子（`MLFMAOperator` / `MLFMAOperatorMPI`）与
PMCHW-MLFMA 算子（`PMCHWMLFMAOperator`）。
"""
module MLFMA

include("Interpolation.jl")
include("Level.jl")
include("SphericalHarmonics.jl")
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
using .SphericalHarmonics
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
export PMCHWMLFMAErrorBudget, PMCHWMLFMAOperator, PMCHWMLFMAOperatorMPI, assemble_near_field_pmchw

end
