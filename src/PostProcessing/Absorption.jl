"""
    Absorption.jl

Specific Absorption Rate (SAR) and absorbed power calculation for lossy
dielectric bodies simulated with volume integral equations (VEFIE).

Exported:
  absorbed_power  — total and per-element absorbed power (W, W/m³)
  sar             — specific absorption rate (W/kg)
"""
module Absorption

using LinearAlgebra

using ...BasisFunctions
using ...CoreModule: num_elements
using ...Utilities.Parameters
using ..CurrentOnGeos: geoVolumeCurrentCal

export absorbed_power, sar

# Physical constants (kept local to avoid module-level globals)
const _c0_abs  = 299_792_458.0
const _mu0_abs = 4π * 1e-7
const _eps0_abs = 1.0 / (_c0_abs^2 * _mu0_abs)

# ─────────────────────────────────────────────────────────────────────────────
# absorbed_power
# ─────────────────────────────────────────────────────────────────────────────

"""
    absorbed_power(basis, I_coeffs, permittivities) → (; P_total, P_density)

Compute the time-averaged absorbed power inside a lossy dielectric body.

## Physics

The time-averaged volumetric absorbed power density at element `t` is

    P_density[t] = (ω/2) ε₀ ε_r''[t] |E[t]|²

where ``ε_r''[t] = -Im(ε_r[t]) ≥ 0`` (engineering convention, lossy medium has
negative imaginary part).

From the VEFIE expansion the polarisation current density at element centroid is

    **J_pol**[t]  =  geoVolumeCurrentCal(I_coeffs, basis, permittivities)[t]

and the total E field follows from

    **E**[t]  =  **J_pol**[t] / ( jω ε₀ (ε_r[t] − 1) )

Substituting:

    P_density[t] = −Im(ε_r[t]) · |**J_pol**[t]|² / ( 2 ω ε₀ |ε_r[t] − 1|² )

For vacuum elements where `|ε_r − 1|` is negligibly small the density is set to
zero.

## Arguments
- `basis`          : `SWGBasis` or `PWCBasis` — volume basis for the VEFIE.
- `I_coeffs`       : MoM solution vector (expansion coefficients of the VEFIE).
- `permittivities` : Complex relative permittivity per element (length = num_elements).

## Returns
`NamedTuple` with
- `P_total    :: Float64` — total power absorbed [W]
- `P_density  :: Vector{Float64}` — per-element absorbed power density [W/m³]
"""
function absorbed_power(
    basis::Union{SWGBasis{IT,FT}, PWCBasis{IT,FT}},
    I_coeffs::Vector{Complex{FT}},
    permittivities::Vector{Complex{FT}},
) where {IT,FT}
    omega = get_omega()
    nt    = num_elements(basis.mesh)

    # Polarisation current density at each element centroid
    J_pol = geoVolumeCurrentCal(I_coeffs, basis, permittivities)

    # Element volumes
    V = _element_volumes(basis)

    P_density = Vector{FT}(undef, nt)

    for t in 1:nt
        eps_r       = permittivities[t]
        abs2_epsr_m1 = abs2(eps_r - one(Complex{FT}))

        if abs2_epsr_m1 < 1e-30
            # Vacuum or near-vacuum: no loss
            P_density[t] = zero(FT)
            continue
        end

        eps_pp = FT(-imag(eps_r))          # ε_r'' = -Im(ε_r) ≥ 0 for lossy

        # |J_pol[t]|² = Σ |J_i|²
        J2 = sum(abs2, J_pol[t])     # real, non-negative

        P_density[t] = eps_pp * J2 / (2 * FT(omega) * FT(_eps0_abs) * abs2_epsr_m1)
    end

    P_total = dot(P_density, V)

    return (; P_total, P_density)
end

# ─────────────────────────────────────────────────────────────────────────────
# sar
# ─────────────────────────────────────────────────────────────────────────────

"""
    sar(basis, I_coeffs, permittivities, mass_densities) → (; SAR_total, SAR_per_element)

Compute the Specific Absorption Rate (SAR) in W/kg.

    SAR[t] = P_density[t] / ρ[t]

where `ρ[t]` is the mass density [kg/m³] of element `t`.

    SAR_total = P_total / total_mass

## Arguments
- `basis`           : `SWGBasis` or `PWCBasis`.
- `I_coeffs`        : MoM solution vector.
- `permittivities`  : Complex relative permittivity per element.
- `mass_densities`  : Mass density [kg/m³] per element (length = num_elements).

## Returns
`NamedTuple` with
- `SAR_total        :: Float64` — whole-body SAR [W/kg]
- `SAR_per_element  :: Vector{Float64}` — per-element SAR [W/kg]
"""
function sar(
    basis::Union{SWGBasis{IT,FT}, PWCBasis{IT,FT}},
    I_coeffs::Vector{Complex{FT}},
    permittivities::Vector{Complex{FT}},
    mass_densities::Vector{FT},
) where {IT,FT}
    (; P_total, P_density) = absorbed_power(basis, I_coeffs, permittivities)

    V   = _element_volumes(basis)
    nt  = length(P_density)

    SAR_per_element = Vector{FT}(undef, nt)
    for t in 1:nt
        ρ = mass_densities[t]
        SAR_per_element[t] = ρ > 0 ? P_density[t] / ρ : zero(FT)
    end

    total_mass = dot(mass_densities, V)
    SAR_total  = total_mass > 0 ? P_total / total_mass : zero(FT)

    return (; SAR_total, SAR_per_element)
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _element_volumes(basis) → Vector{FT}

Return the volume of each mesh element for SWG or PWC bases.
For PWCBasis the stored `volume` field is used directly.
For SWGBasis the volume is computed from node coordinates.
"""
function _element_volumes(basis::PWCBasis{IT,FT}) where {IT,FT}
    return FT[basis.functions[t].volume for t in 1:length(basis.functions)]
end

function _element_volumes(basis::SWGBasis{IT,FT}) where {IT,FT}
    mesh   = basis.mesh
    nodes  = mesh.node
    tetras = mesh.tetras
    nt     = num_elements(mesh)
    V      = Vector{FT}(undef, nt)
    for t in 1:nt
        v1 = nodes[:, tetras[1, t]]
        v2 = nodes[:, tetras[2, t]]
        v3 = nodes[:, tetras[3, t]]
        v4 = nodes[:, tetras[4, t]]
        V[t] = FT(abs(dot(v2 - v1, cross(v3 - v1, v4 - v1))) / 6)
    end
    return V
end

end  # module Absorption
