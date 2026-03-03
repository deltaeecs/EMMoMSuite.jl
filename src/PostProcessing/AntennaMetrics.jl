module AntennaMetrics

using StaticArrays
using LinearAlgebra

using ...BasisFunctions
using ...CoreModule
using ...Utilities.Parameters
using ..FarField: farField

export antenna_directivity, input_impedance, beam_metrics

# ─────────────────────────────────────────────────────────────────────────────
# antenna_directivity
# ─────────────────────────────────────────────────────────────────────────────

"""
    antenna_directivity(θs, ϕs, ICoeff, basis; P_input=nothing, source=nothing)

Compute the antenna directivity `D(θ,ϕ)` over the specified far-field grid.

## Formulation

The radiation intensity per unit solid angle:
```
U(θ,ϕ) = (|E_θ|² + |E_ϕ|²) / (2 η₀)
```
where `E_θ/E_ϕ` are the far-field angular functions (units V, excluding
the `e^{-jkr}/r` propagation factor).

Total radiated power (numerical integration with trapezoidal rule):
```
P_rad = ∬ U(θ,ϕ) sin(θ) dθ dϕ
```

Directivity:
```
D(θ,ϕ) = 4π U(θ,ϕ) / P_rad
```

If `P_input` (the power delivered to the antenna input, in watts) is
provided, the gain and radiation efficiency are also returned:
```
η_eff = P_rad / P_input
G(θ,ϕ) = η_eff · D(θ,ϕ)
```

## Arguments
- `θs`        : Vector of elevation angles (radians), length Nθ (must be ≥ 2 for integration).
- `ϕs`        : Vector of azimuth angles (radians), length Nϕ (must be ≥ 2 for integration).
- `ICoeff`    : Current-coefficient vector (MoM solution).
- `basis`     : RWG (or similar) basis object.
- `P_input`   : Optional real input power (W). When given, gain and η_eff are returned.
- `source`    : Source object forwarded to `farField`. Pass `nothing` for scatter problems.

## Returns
Named tuple with fields:
- `:D`      — directivity matrix (Nθ × Nϕ, dimensionless).
- `:P_rad`  — total radiated power (W).
- `:G`      — gain matrix (Nθ × Nϕ) — *only* when `P_input` is supplied.
- `:η_eff`  — radiation efficiency — *only* when `P_input` is supplied.
"""
function antenna_directivity(
    θs::Vector{FT},
    ϕs::Vector{FT},
    ICoeff::Vector{<:Complex},
    basis;
    P_input::Union{Nothing,Real} = nothing,
    source = nothing,
) where {FT<:Real}

    Nθ = length(θs)
    Nϕ = length(ϕs)
    Nθ >= 2 || throw(ArgumentError("θs must have ≥ 2 elements for integration"))
    Nϕ >= 2 || throw(ArgumentError("ϕs must have ≥ 2 elements for integration"))

    η₀ = get_eta0()

    # ── Far-field pattern (2 × Nθ × Nϕ) ──────────────────────────────────────
    farE = farField(θs, ϕs, ICoeff, basis, source)  # shape (2, Nθ, Nϕ)

    # ── Radiation intensity U(θ,ϕ) [W/sr] ────────────────────────────────────
    U = (abs2.(farE[1, :, :]) .+ abs2.(farE[2, :, :])) ./ (2 * η₀)
    # U is Nθ × Nϕ

    # ── Numerical integration: P_rad = ∫∫ U sinθ dθ dϕ (trapezoidal) ─────────
    P_rad = _trapz2d(U, θs, ϕs)

    # ── Directivity ───────────────────────────────────────────────────────────
    D = P_rad > 0 ? 4π .* U ./ P_rad : zeros(FT, Nθ, Nϕ)

    if P_input !== nothing
        P_input_r = Float64(P_input)::Float64
        η_eff = P_rad / P_input_r
        G = η_eff .* D
        return (; D, G, P_rad, η_eff)
    else
        return (; D, P_rad)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# input_impedance
# ─────────────────────────────────────────────────────────────────────────────

"""
    input_impedance(source::DeltaGapSource, ICoeff, basis::RWGBasis) → Complex

Compute the input impedance of a delta-gap-excited antenna.

## Formulation

For a delta-gap source applied to a set of excited edges, the total
input current flowing through the feed is:
```
I_in = Σ ICoeff[n],   n ∈ source.edge_indices
```

The input impedance is then:
```
Z_in = V_applied / I_in
```

where `V_applied = source.voltage`.

## Arguments
- `source`  : `DeltaGapSource` carrying the applied voltage and edge indices.
- `ICoeff`  : Current-coefficient vector (MoM solution).
- `basis`   : `RWGBasis` (used for bound-checking `num_basis`).

## Returns
`Z_in` as `ComplexF64`.
"""
function input_impedance(
    source::DeltaGapSource,
    ICoeff::Vector{<:Complex},
    basis::RWGBasis,
)
    N = num_basis(basis)
    I_in = zero(ComplexF64)
    for idx in source.edge_indices
        if 1 <= idx <= N
            I_in += ICoeff[idx]
        else
            @warn "DeltaGapSource edge index $idx out of bounds (1:$N); skipped"
        end
    end
    iszero(I_in) && error("input_impedance: total feed current is zero — check DeltaGapSource.edge_indices")
    return ComplexF64(source.voltage) / ComplexF64(I_in)
end

# ─────────────────────────────────────────────────────────────────────────────
# beam_metrics
# ─────────────────────────────────────────────────────────────────────────────

"""
    beam_metrics(θs, pattern_dB; peak_threshold_dB=-3.0)

Compute scalar beam metrics from a 1-D antenna pattern (dB scale).

## Returns
Named tuple:
- `:peak_angle`  — angle (same units as `θs`) of the pattern maximum.
- `:HPBW`        — half-power beam width (angular span where pattern ≥ peak + `peak_threshold_dB`).
- `:SLL_dB`      — side-lobe level relative to the peak (negative number; 0 if no side lobe found).

## Arguments
- `θs`               : 1-D vector of angles (radians or degrees — consistent with `pattern_dB`).
- `pattern_dB`        : 1-D vector of gain / RCS values in dB.
- `peak_threshold_dB` : Threshold for half-power; default `-3.0`.
"""
function beam_metrics(
    θs::AbstractVector{<:Real},
    pattern_dB::AbstractVector{<:Real};
    peak_threshold_dB::Real = -3.0,
)
    length(θs) == length(pattern_dB) || throw(DimensionMismatch(
        "θs ($(length(θs))) and pattern_dB ($(length(pattern_dB))) must have the same length"))
    length(θs) >= 3 || throw(ArgumentError("beam_metrics requires ≥ 3 points"))

    FT = float(eltype(θs))

    # ── Peak ──────────────────────────────────────────────────────────────────
    peak_val, peak_idx = findmax(pattern_dB)
    peak_angle = FT(θs[peak_idx])

    # ── HPBW — angular span of the main-lobe above `peak_val + threshold` ────
    threshold = peak_val + peak_threshold_dB      # e.g. peak - 3 dB
    in_main = pattern_dB .>= threshold

    hp_angles = θs[in_main]
    HPBW = isempty(hp_angles) ? FT(0) : FT(last(hp_angles) - first(hp_angles))

    # ── SLL — maximum of side-lobe LOCAL MAXIMA relative to peak ─────────────
    # A "side lobe peak" is a point outside the main lobe that is ≥ both neighbors.
    sll_abs = FT(-Inf)          # absolute peak of side lobes (dB)
    n = length(pattern_dB)

    # Interior points
    for i in 2:n-1
        in_main[i] && continue
        if pattern_dB[i] >= pattern_dB[i-1] && pattern_dB[i] >= pattern_dB[i+1]
            sll_abs = max(sll_abs, FT(pattern_dB[i]))
        end
    end
    # Left endpoint
    if n >= 1 && !in_main[1]
        is_local = (n == 1) || (pattern_dB[1] >= pattern_dB[2])
        is_local && (sll_abs = max(sll_abs, FT(pattern_dB[1])))
    end
    # Right endpoint
    if n >= 2 && !in_main[n]
        pattern_dB[n] >= pattern_dB[n-1] && (sll_abs = max(sll_abs, FT(pattern_dB[n])))
    end

    SLL_dB = isinf(sll_abs) ? FT(0) : FT(sll_abs - peak_val)

    return (; peak_angle, HPBW, SLL_dB)
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helper: 2-D trapezoidal integration
#   ∫∫ f(θ,ϕ) sin(θ) dθ dϕ  over a uniform-ish grid
# ─────────────────────────────────────────────────────────────────────────────
function _trapz2d(
    f::Matrix{<:Real},
    θs::AbstractVector{<:Real},
    ϕs::AbstractVector{<:Real},
)
    Nθ, Nϕ = size(f)
    (Nθ == length(θs) && Nϕ == length(ϕs)) ||
        throw(DimensionMismatch("f size $(size(f)) ≠ ($(length(θs)), $(length(ϕs)))"))

    # Trapezoidal weights in θ and ϕ
    wθ = _trapz_weights(θs)
    wϕ = _trapz_weights(ϕs)

    result = 0.0
    for iϕ in 1:Nϕ, iθ in 1:Nθ
        result += f[iθ, iϕ] * sin(θs[iθ]) * wθ[iθ] * wϕ[iϕ]
    end
    return result
end

function _trapz_weights(xs::AbstractVector{<:Real})
    n = length(xs)
    w = similar(xs, Float64)
    if n == 1
        w[1] = 0.0
        return w
    end
    w[1]   = (xs[2] - xs[1]) / 2
    w[end] = (xs[end] - xs[end-1]) / 2
    for i in 2:n-1
        w[i] = (xs[i+1] - xs[i-1]) / 2
    end
    return w
end

end  # module AntennaMetrics
