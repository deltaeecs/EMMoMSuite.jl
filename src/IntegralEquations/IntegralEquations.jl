module IntegralEquations

using ..CoreModule
using ..Geometry
using ..BasisFunctions

include("Kernels.jl")
using .Kernels

include("Impedance.jl")
using .Impedance

include("EFIE.jl")
using .EFIEModule

include("MFIE.jl")
using .MFIEModule

include("CFIE.jl")
using .CFIEModule

include("VEFIE.jl")
using .VEFIEModule

include("SCFIE.jl")
using .SCFIEModule

include("Excitation.jl")
using .Excitation

export EFIE, MFIE, CFIE, VEFIE, SCFIE, assemble_impedance_matrix
export green_function_free_space

end
