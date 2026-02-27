using Pkg
Pkg.activate(".")
push!(LOAD_PATH, @__DIR__)
include("scripts/verification/verify_VEFIE_mie.jl")
verify_vefie_mie()
