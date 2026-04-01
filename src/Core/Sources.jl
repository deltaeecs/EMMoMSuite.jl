module Sources

using ..CoreModule
using ..Constants
using LinearAlgebra

import ..CoreModule: incident_field

export PlaneWave, DeltaGapSource, incident_field

"""
    PlaneWave <: AbstractSource

Plane wave electromagnetic excitation source.

# Mathematical Formulation

A plane wave is a uniform electromagnetic wave with constant phase in planes
perpendicular to the propagation direction. The electric field is given by:

```math
\\mathbf{E}(\\mathbf{r}) = \\mathbf{E}_0 e^{-j\\mathbf{k} \\cdot \\mathbf{r}}
```

where:
- `\\mathbf{E}_0` is the polarization vector (transverse to k)
- `\\mathbf{k} = k_0 \\hat{\\mathbf{k}}` is the wave vector (magnitude ``k_0 = 2\\pi/\\lambda``)
- Direction: ``\\hat{\\mathbf{k}} = (\\sin\\theta\\cos\\phi, \\sin\\theta\\sin\\phi, \\cos\\theta)``

# Fields

- `frequency::Float64`: Operating frequency [Hz]
- `theta::Float64`: Polar angle of propagation direction [radians]
  - θ = 0: propagation along +z axis
  - θ = π/2: propagation in xy-plane
  - θ = π: propagation along -z axis
- `phi::Float64`: Azimuthal angle of propagation direction [radians]
  - φ = 0: propagation in xz-plane
  - φ = π/2: propagation in yz-plane
- `polarization::Vector{Float64}`: Electric field polarization (3-vector, automatically normalized and made transverse to k)
- `k::Float64`: Wavenumber magnitude ``k_0 = 2\\pi f / c_0`` [rad/m] (computed internally)

# Constructor

    PlaneWave(frequency, theta, phi, polarization)

Constructs a plane wave with specified parameters. The polarization vector is automatically:
1. Projected onto the plane transverse to the propagation direction
2. Normalized to unit magnitude

Throws `ArgumentError` if polarization has zero transverse component.

# Examples

```julia
# z-directed plane wave with x-polarization
pw = PlaneWave(1e9, 0.0, 0.0, [1.0, 0.0, 0.0])

# Oblique incidence (θ=45°, φ=30°) with y-polarization
pw = PlaneWave(300e6, π/4, π/6, [0.0, 1.0, 0.0])

# Evaluate incident field at point r = [0.1, 0.2, 0.3]
E_inc = incident_field(pw, [0.1, 0.2, 0.3])
```

# See Also

- [`incident_field`](@ref): Evaluate plane wave field at a point
- [`DeltaGapSource`](@ref): Delta-gap voltage source for ports
"""
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

"""
    DeltaGapSource <: AbstractSource

Delta-gap voltage source for lumped port excitation.

# Mathematical Formulation

A delta-gap source models a localized voltage excitation across an edge
or gap in the mesh. It is typically used for:
- Microstrip antenna feeds
- Coaxial cable ports
- Lumped element excitations

The voltage ``V`` is applied uniformly across the specified edges, creating
a highly localized electric field in the gap region.

# Fields

- `frequency::Float64`: Operating frequency [Hz]
- `edge_indices::Vector{Int}`: Indices of edges (RWG basis functions) where voltage is applied
- `voltage::ComplexF64`: Applied voltage [V] (complex to support phase)

# Constructor

    DeltaGapSource(frequency, edge_indices, voltage)

# Examples

```julia
# Single-edge delta gap with 1V excitation
dg = DeltaGapSource(1e9, [42], 1.0 + 0.0im)

# Multi-edge differential port (edges 10 and 11 with opposite polarity)
dg_pos = DeltaGapSource(2.4e9, [10], 0.5 + 0.0im)
dg_neg = DeltaGapSource(2.4e9, [11], -0.5 + 0.0im)
```

# See Also

- [`PlaneWave`](@ref): Plane wave excitation
- [`LumpedPort`](@ref): Structured lumped port definition (in Ports module)
"""
struct DeltaGapSource <: AbstractSource
    frequency::Float64
    edge_indices::Vector{Int}
    voltage::ComplexF64
end

"""
    incident_field(source::AbstractSource, r::AbstractVector) -> Vector{ComplexF64}

Evaluate the incident electric field at a point in space.

# Arguments

- `source::AbstractSource`: Excitation source (e.g., PlaneWave, DeltaGapSource)
- `r::AbstractVector`: Observation point coordinates [x, y, z] in meters

# Returns

- `Vector{ComplexF64}`: Complex electric field vector ``\\mathbf{E}_{inc}(\\mathbf{r})`` [V/m]

# Details

For a plane wave, the field is:
```math
\\mathbf{E}_{inc}(\\mathbf{r}) = \\mathbf{E}_0 e^{-j\\mathbf{k} \\cdot \\mathbf{r}}
```

For a delta-gap source, the field is typically non-zero only in the immediate
vicinity of the gap (implementation-specific).

# Examples

```julia
# Plane wave incident field
pw = PlaneWave(1e9, π/4, 0.0, [1.0, 0.0, 0.0])
r = [0.1, 0.2, 0.3]
E = incident_field(pw, r)

# E is a 3-element complex vector
@assert length(E) == 3
@assert eltype(E) == ComplexF64
```

# See Also

- [`PlaneWave`](@ref): Plane wave source definition
- [`DeltaGapSource`](@ref): Delta-gap source definition
"""
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
