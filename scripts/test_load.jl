using Pkg
Pkg.activate(joinpath(@__DIR__, "../../MoM_Kernels"))
Pkg.instantiate()
using MoM_Kernels
println("Success loading MoM_Kernels")
