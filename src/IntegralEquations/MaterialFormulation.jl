# MaterialFormulation.jl
#
# Generic (case-independent) material → integral-equation formulation layer.
#
# Design rule: the library owns the algorithmic decisions (which operator, which
# basis, which post-processing dispatch); callers (benchmark cases, scripts,
# tests) only describe geometry + materials and call these functions.  Nothing
# in this file knows about a specific case format.
#
# Material protocol — a caller's material model type implements:
#   is_pec_material(m)          -> Bool   (default: false)
#   relative_permittivity(m)    -> Complex (default: complex(m.eps_r))
#   relative_permeability(m)    -> Real    (default: m.mu_r or 1.0)
# `BoundMesh` bindings may hold ANY material type implementing the protocol.

"""
    is_pec_material(m) -> Bool

Material-protocol hook: is `m` a perfect electric conductor?
Default implementation returns `false`; PEC-like material models override it.
"""
is_pec_material(m) = false

"""
    relative_permittivity(m) -> Complex

Material-protocol hook: complex relative permittivity of `m`.
A bare number is its own permittivity; structured material models default to
`complex(m.eps_r)` and override for loss or field-dependent ε.
"""
relative_permittivity(m::Number) = ComplexF64(m)
relative_permittivity(m) = complex(m.eps_r)

"""
    relative_permeability(m) -> Real

Material-protocol hook: relative permeability of `m`.
Default implementation returns `m.mu_r` when present, otherwise `1.0`.
"""
function relative_permeability(m)
    return hasproperty(m, :mu_r) ? m.mu_r : 1.0
end

"""
    volume_permittivities(bm::BoundMesh) -> Vector{ComplexF64}

Per-tetrahedron complex relative permittivity vector read from the material
bindings of `bm` (tetrahedron `t` uses tag `t`).
"""
function volume_permittivities(bm::BoundMesh)
    return ComplexF64[relative_permittivity(element_material(bm, t))
                      for t in 1:bm.mesh.tetnum]
end

"""
    is_all_pec(bm::BoundMesh) -> Bool

`true` when every bound region material is a PEC (material protocol).
"""
is_all_pec(bm::BoundMesh) = all(is_pec_material, values(bm.material_map))

"""
    build_volume_system(bm::BoundMesh, freq; ie = "auto")
        -> (op, mesh, kind::Symbol, permittivities)

Generic volume (tetrahedral) formulation dispatch — the single branch point
from a bound tetrahedral mesh to a solver operator.  `kind` is one of:
- `:vefie`         — dielectric volume: `VEFIE` operator on `bm.mesh`
  (solve with `SWGBasis(mesh)`; pass `permittivities` to
  `excitation_vector` and `radarCrossSection`).
- `:efie_fallback` — all-PEC volume: degenerate to the surface `EFIE` on
  `extract_surface(bm.mesh)` (solve with `RWGBasis(mesh)`).

`ie = "auto"` selects `:efie_fallback` iff every bound region is PEC;
`ie = "VEFIE"` / `ie = "EFIE"` force one branch (error when inconsistent with
the bindings).  A mix of PEC and dielectric regions is rejected: model PEC
boundaries as surface interfaces instead.
"""
function build_volume_system(bm::BoundMesh, freq::Real; ie::String = "auto")
    sym = Symbol(ie)
    sym in (:auto, :VEFIE, :EFIE) ||
        error("unsupported ie \"$ie\" for volume system")
    isempty(bm.material_map) && error("no material bound to mesh")

    allpec = is_all_pec(bm)
    if sym === :EFIE || (sym === :auto && allpec)
        surface = extract_surface(bm.mesh)
        return (op = EFIE(freq), mesh = surface, kind = :efie_fallback,
                permittivities = ComplexF64[])
    end

    anypec = any(is_pec_material, values(bm.material_map))
    anypec && error("mixed PEC + dielectric volume regions are not supported; " *
                    "model PEC regions as surfaces (interface) instead")
    permittivities = volume_permittivities(bm)
    return (op = VEFIE(freq, permittivities), mesh = bm.mesh, kind = :vefie,
            permittivities = permittivities)
end

"""
    derive_formulation(has_dielectric::Bool, ie::String) -> Symbol

Core surface-formulation rule shared by all callers:
- `ie = "auto"` → `:PMCHW` when any side is dielectric, `:EFIE` otherwise;
- `ie = "PMCHW" | "CFIE" | "EFIE"` selects the named formulation.
"""
function derive_formulation(has_dielectric::Bool, ie::String)
    sym = Symbol(ie)
    if sym === :auto
        return has_dielectric ? :PMCHW : :EFIE
    end
    sym in (:PMCHW, :CFIE, :EFIE) || error("unknown formulation :$sym")
    return sym
end

"""
    farfield_from_rcs(rcs_dB) -> Array{Float64,3}

Derive a far-field magnitude array from bistatic RCS values (dBsm): volume
currents have no `farField` overload, so the magnitude is reconstructed as
`sqrt(10^((RCS_dB + 30)/10))` in the `[2, nθ, nϕ]` layout used by plotting
(slot 1 holds the magnitude; slot 2 stays zero).
"""
function farfield_from_rcs(rcs_dB)
    rcs_lin = 10.0 .^ ((rcs_dB .+ 30.0) ./ 10.0)   # dBsm → m²
    mag = sqrt.(rcs_lin)
    nθ, nϕ = size(rcs_dB)
    FF = zeros(Float64, 2, nθ, nϕ)
    FF[1, :, :] .= mag
    return FF
end
