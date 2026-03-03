"""
    module MaterialsModule

Phase 19 — Electromagnetic material library and dispersion models.

Subfiles:
  MaterialLibrary.jl  — 19.3: MaterialModel hierarchy, library CRUD, built-ins
  DispersionModels.jl — 19.4: Debye / Drude / Lorentz dispersive models
"""
module MaterialsModule

using StaticArrays
using LinearAlgebra
using JLD2

include("MaterialLibrary.jl")
include("DispersionModels.jl")

end # module MaterialsModule
