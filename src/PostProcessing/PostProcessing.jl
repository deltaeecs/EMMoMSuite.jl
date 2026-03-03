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

using .RadiationIntegral
using .CurrentOnGeos
using .RCS
using .FarField
using .NearField
using .FieldCut

export radarCrossSection, farField, geoElectricJCal, geoVolumeCurrentCal, calculate_near_field
export r̂θϕInfo, raditionalIntegralNθϕCal
export field_cut_line, field_cut_plane

end
