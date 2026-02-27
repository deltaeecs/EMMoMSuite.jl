module Materials

using ..Constants

export AbstractMaterial
export PEC, PMC, Dielectric, ImpedanceSurface
export permittivity, permeability, impedance

"""
    AbstractMaterial

Abstract base type for all materials.
"""
abstract type AbstractMaterial end

"""
    PEC

Perfect Electric Conductor.
"""
struct PEC <: AbstractMaterial end

"""
    PMC

Perfect Magnetic Conductor.
"""
struct PMC <: AbstractMaterial end

"""
    Dielectric{T<:Real}

Homogeneous dielectric material.
"""
struct Dielectric{T<:Real} <: AbstractMaterial
    εr::Complex{T}  # Relative permittivity
    μr::Complex{T}  # Relative permeability
end

Dielectric(εr::Real, μr::Real) = Dielectric(complex(εr), complex(μr))
Dielectric(εr::Number) = Dielectric(complex(εr), complex(1.0))

"""
    ImpedanceSurface{T<:Real}

Surface with specified impedance.
"""
struct ImpedanceSurface{T<:Real} <: AbstractMaterial
    Zs::Complex{T}  # Surface impedance
end

# Interface functions
permittivity(::PEC) = Inf
permeability(::PEC) = Constants.mu0

permittivity(::PMC) = Constants.eps0
permeability(::PMC) = Inf

permittivity(m::Dielectric) = m.εr * Constants.eps0
permeability(m::Dielectric) = m.μr * Constants.mu0
impedance(m::Dielectric) = sqrt(permeability(m) / permittivity(m))

end
