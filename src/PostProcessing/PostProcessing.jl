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

using .RadiationIntegral
using .CurrentOnGeos
using .RCS
using .FarField
using .NearField

export radarCrossSection, farField, geoElectricJCal, calculate_near_field
export r̂θϕInfo, raditionalIntegralNθϕCal

end
