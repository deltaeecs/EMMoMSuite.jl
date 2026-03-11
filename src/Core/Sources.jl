module Sources

using ..CoreModule
using ..Constants
using LinearAlgebra

import ..CoreModule: incident_field

export PlaneWave, DeltaGapSource, incident_field

struct PlaneWave <: AbstractSource
    frequency::Float64
    theta::Float64
    phi::Float64
    polarization::Vector{Float64}
    k::Float64

    function PlaneWave(frequency, theta, phi, polarization)
        k = 2 * pi * frequency / Constants.c0
        st, ct = sincos(theta)
        sp, cp = sincos(phi)
        k_dir = [st * cp, st * sp, ct]

        pol_vec = Vector{Float64}(polarization)
        pol_transverse = pol_vec .- dot(pol_vec, k_dir) .* k_dir
        norm(pol_transverse) > sqrt(eps(Float64)) ||
            throw(ArgumentError("PlaneWave polarization must have a non-zero component transverse to the propagation direction"))

        new(frequency, theta, phi, normalize(pol_transverse), k)
    end
end

struct DeltaGapSource <: AbstractSource
    frequency::Float64
    edge_indices::Vector{Int}
    voltage::ComplexF64
end

function incident_field(source::PlaneWave, r::AbstractVector)
    # k vector direction (propagation direction)
    # k = k0 * (sin theta cos phi, sin theta sin phi, cos theta)

    st, ct = sincos(source.theta)
    sp, cp = sincos(source.phi)
    k_dir = [st * cp, st * sp, ct]

    # Phase factor exp(-j * k * r . k_dir)
    phase = exp(-im * source.k * dot(r, k_dir))
    return source.polarization * phase
end

end
