module FarField

using StaticArrays
using LinearAlgebra
using Base.Threads

using ...Geometry
using ...BasisFunctions
using ...Utilities.Parameters
using ..CurrentOnGeos
using ..RadiationIntegral: r̂θϕInfo, raditionalIntegralNθϕCal, ∠Info
using ...CoreModule: num_elements

export farField

"""
    farField(θs_obs, ϕs_obs, ICoeff, basis, source) -> Array{Complex, 3}

Calculate far-field radiation pattern from computed surface currents.

# Arguments

- `θs_obs::Vector{FT}`: Observation theta angles [radians], range [0, π]
  - θ = 0: +z direction (zenith)
  - θ = π/2: xy-plane (horizon)
  - θ = π: -z direction (nadir)
- `ϕs_obs::Vector{FT}`: Observation phi angles [radians], range [0, 2π]
  - φ = 0: +x direction
  - φ = π/2: +y direction
- `ICoeff::Vector{Complex}`: Current coefficients (RWG basis expansion)
  - Obtained from MoM solve: `ICoeff = Z \\ V`
- `basis::RWGBasis`: RWG basis functions on surface mesh
- `source`: Excitation source (used for normalization, if applicable)

# Returns

- `Array{Complex{FT}, 3}`: Far-field electric field components [V/m]
  - Dimensions: `[2, length(θs_obs), length(ϕs_obs)]`
  - `farE[1, i, j]`: E_θ component at (θs_obs[i], ϕs_obs[j])
  - `farE[2, i, j]`: E_φ component at (θs_obs[i], ϕs_obs[j])

# Mathematical Background

The far-field electric field is computed via radiation integral:

```math
\\mathbf{E}(r, \\theta, \\phi) = -j\\frac{k_0 \\eta_0}{4\\pi} e^{-jk_0 r} 
\\int_S \\mathbf{J}(\\mathbf{r}') e^{jk_0 \\hat{r} \\cdot \\mathbf{r}'} dS'
```

where only the transverse components (E_θ, E_φ) are retained in far field.

# Examples

```julia
# After solving EFIE
Z = assemble_impedance_matrix(efie)
V = compute_excitation_vector(efie, plane_wave)
I = Z \\ V

# Compute far-field pattern (θ: 0-180°, φ=0° cut)
θ_obs = range(0, π, length=181)
φ_obs = [0.0]
E_far = farField(θ_obs, φ_obs, I, basis, plane_wave)

# Extract E_θ component (co-pol for z-polarized excitation)
E_theta = E_far[1, :, 1]
E_phi = E_far[2, :, 1]

# Convert to magnitude (dB)
E_mag = sqrt.(abs2.(E_theta) .+ abs2.(E_phi))
E_dB = 20 .* log10.(E_mag ./ maximum(E_mag))
```

# Performance

- Uses multi-threading via `@threads` for parallel evaluation
- Complexity: O(N_obs × N_tri) where N_obs = length(θs_obs) × length(ϕs_obs)
- For large meshes (>10k triangles), consider MLFMA-accelerated far-field

# See Also

- [`compute_rcs`](@ref): Radar cross section (magnitude of far-field)
- [`raditionalIntegralNθϕCal`](@ref): Internal radiation integral computation
- [`RWGBasis`](@ref): Surface basis functions

# Notes

- Far-field approximation valid for r >> λ and r >> D (D = antenna size)
- Returns spherical components (E_θ, E_φ), NOT Cartesian (E_x, E_y, E_z)
- Use `r̂θϕInfo` utilities to convert to Cartesian if needed
"""
function farField(
    θs_obs::Vector{FT},
    ϕs_obs::Vector{FT},
    ICoeff::Vector{CT},
    basis::RWGBasis{IT,FT},
    source,
) where {IT<:Integer,FT<:Real,CT<:Complex{FT}}

    # Calculate currents on geometry
    Jtris = geoElectricJCal(ICoeff, basis)

    k0 = get_k0()
    eta0 = get_eta0()
    jk0 = im * k0
    div4pi = 1.0 / (4 * pi)

    Nθ_obs = length(θs_obs)
    Nϕ_obs = length(ϕs_obs)
    nobs = Nθ_obs * Nϕ_obs

    # Precompute angle info
    θsobsInfo = [∠Info(θ) for θ in θs_obs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs_obs]

    # Create grid of observation directions
    r̂θsϕs = [r̂θϕInfo(θ, ϕ) for θ in θsobsInfo, ϕ in ϕsobsInfo]
    r̂θsϕs_flat = vec(r̂θsϕs)

    farEθsϕs = zeros(Complex{FT}, 2, Nθ_obs, Nϕ_obs)
    farE_flat = reshape(farEθsϕs, 2, nobs)

    # Construct TriangleInfo for radiation integral
    mesh = basis.mesh
    ntri = num_elements(mesh)
    trianglesInfo = Vector{TriangleInfo{IT,FT}}(undef, ntri)
    for t = 1:ntri
        trianglesInfo[t] = TriangleInfo(mesh, t)
    end

    @threads for ii = 1:nobs
        r_info = r̂θsϕs_flat[ii]

        # Calculate Radiation Integral
        Nθϕ = raditionalIntegralNθϕCal(r_info, trianglesInfo, Jtris)
        Nθϕ = raditionalIntegralNθϕCal(r_info, trianglesInfo, Jtris)

        # Calculate Far Field
        # E_far = (-j k0 eta0 / 4pi) * N
        # Legacy code uses div4pi.
        factor = -jk0 * eta0 * div4pi
        farEθϕ = factor .* Nθϕ

        # Add source field if applicable
        # farEθϕ .+= sourceFarEfield(source, r_info)

        farE_flat[:, ii] = farEθϕ
    end

    return farEθsϕs
end

end
