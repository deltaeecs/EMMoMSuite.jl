module PostProcessing

using LinearAlgebra
using StaticArrays
using Statistics
using ProgressMeter
using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..Utilities

include("RadiationIntegral.jl")
include("CurrentOnGeos.jl")
include("RCS.jl")
include("FarField.jl")
include("NearField.jl")
include("FieldCut.jl")
include("AntennaMetrics.jl")
include("Absorption.jl")
# Phase 21: 高级后处理与快速算法
include("SolverResult.jl")
include("NearFieldAdvanced.jl")
include("FarFieldPattern.jl")
include("RCSBatch.jl")
include("MLFMAFastPost.jl")

using .RadiationIntegral
using .CurrentOnGeos
using .RCS
using .FarField
using .NearField
using .FieldCut
using .AntennaMetrics
using .Absorption
using .SolverResultModule
using .NearFieldAdvancedModule
using .FarFieldPatternModule
using .RCSBatchModule
using .MLFMAFastPostModule

export radarCrossSection, farField, geoElectricJCal, geoVolumeCurrentCal, calculate_near_field
export r̂θϕInfo, raditionalIntegralNθϕCal
export field_cut_line, field_cut_plane
export antenna_directivity, input_impedance, beam_metrics
export absorbed_power, sar
# Phase 21 exports
export SolverResult
export NearFieldGrid, NearFieldLine
export FarFieldPattern
export gain, gain_db, hpbw, side_lobe_level
export axial_ratio, co_cross_decompose, xpd
export RCSResult, rcs_frequency_response, rcs_angular_pattern
export MLFMACache, invalidate_cache!, validate_cache, solve_multi_rhs

end
