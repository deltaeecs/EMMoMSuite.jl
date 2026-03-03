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

include("PMCHW.jl")
using .PMCHWModule

include("Excitation.jl")
using .Excitation

export EFIE, MFIE, CFIE, VEFIE, SCFIE, PMCHW, assemble_impedance_matrix
export green_function_free_space
export assemble_K_offdiag, efie_from_keta, excitation_vector

end
