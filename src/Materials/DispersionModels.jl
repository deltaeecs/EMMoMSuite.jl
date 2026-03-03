"""
    DispersionModels.jl

Phase 19.4 — Frequency-dispersive material models for electromagnetics.

Supported models:
  DebyeModel    — relaxation / polar dielectric  (Debye relaxation)
  DrudeModel    — free-electron metal            (Drude plasma)
  LorentzModel  — bound-electron resonance       (Lorentz oscillator)

All models implement `eval_permittivity(model, ω::Float64) → ComplexF64`
and a vectorised overload `eval_permittivity(model, ωs::Vector{Float64})`.

Sign convention: e^{-jωt} (standard CEM).  Under this convention a lossy
medium has `imag(εr) > 0`, i.e. the imaginary part of the complex permittivity
ε = ε' - jε'' is negative when written in the ε'' (positive) form.
"""

# ─────────────────────────────────────────────────────────────────────────────
# Pole structs
# ─────────────────────────────────────────────────────────────────────────────

"""
    DebyePole(Δε, τ)

A single Debye relaxation pole with static permittivity contribution `Δε`
(dimensionless) and relaxation time `τ` (seconds).
"""
struct DebyePole
    Δε :: Float64
    τ  :: Float64
end

"""
    LorentzPole(Δε, ω0, δ)

A single Lorentz oscillator pole.

# Fields
- `Δε` : static permittivity contribution (dimensionless)
- `ω0` : undamped resonance frequency (rad/s)
- `δ`  : damping rate (rad/s)
"""
struct LorentzPole
    Δε :: Float64
    ω0 :: Float64
    δ  :: Float64
end

# ─────────────────────────────────────────────────────────────────────────────
# Model structs
# ─────────────────────────────────────────────────────────────────────────────

"""
    DebyeModel(ε_inf, poles)

Multi-pole Debye relaxation model for polar dielectrics:
```
ε(ω) = ε_inf + Σ_k Δε_k / (1 + jωτ_k)
```

# Arguments
- `ε_inf` : high-frequency (optical) permittivity limit
- `poles` : `Vector{DebyePole}` — one pole per relaxation process
"""
struct DebyeModel <: MaterialModel
    ε_inf :: Float64
    poles :: Vector{DebyePole}
end

"""
    DrudeModel(ε_inf, ωp, γ)

Drude free-electron model for metals / plasmas:
```
ε(ω) = ε_inf - ωp² / (ω² + jωγ)
```

# Arguments
- `ε_inf` : background (inter-band) permittivity
- `ωp`    : plasma frequency (rad/s)
- `γ`     : free-electron damping / collision rate (rad/s)
"""
struct DrudeModel <: MaterialModel
    ε_inf :: Float64
    ωp    :: Float64
    γ     :: Float64
end

"""
    LorentzModel(ε_inf, poles)

Multi-pole Lorentz oscillator model for bound-electron resonances:
```
ε(ω) = ε_inf + Σ_k Δε_k ω0_k² / (ω0_k² - ω² + 2jω δ_k)
```

# Arguments
- `ε_inf` : high-frequency background permittivity
- `poles` : `Vector{LorentzPole}`
"""
struct LorentzModel <: MaterialModel
    ε_inf :: Float64
    poles :: Vector{LorentzPole}
end

# ─────────────────────────────────────────────────────────────────────────────
# Scalar evaluation  (ω in rad/s)
# ─────────────────────────────────────────────────────────────────────────────

"""
    eval_permittivity(model::DebyeModel, ω::Float64) → ComplexF64

Evaluate the complex relative permittivity at angular frequency `ω` (rad/s).
"""
function eval_permittivity(model::DebyeModel, ω::Float64)
    ε = ComplexF64(model.ε_inf)
    @inbounds for p in model.poles
        ε += p.Δε / (1.0 + im * ω * p.τ)
    end
    return ε
end

"""
    eval_permittivity(model::DrudeModel, ω::Float64) → ComplexF64

Evaluate the complex relative permittivity at angular frequency `ω` (rad/s).
"""
function eval_permittivity(model::DrudeModel, ω::Float64)
    return ComplexF64(model.ε_inf) - model.ωp^2 / (ω^2 + im * ω * model.γ)
end

"""
    eval_permittivity(model::LorentzModel, ω::Float64) → ComplexF64

Evaluate the complex relative permittivity at angular frequency `ω` (rad/s).
"""
function eval_permittivity(model::LorentzModel, ω::Float64)
    ε  = ComplexF64(model.ε_inf)
    ω2 = ω^2
    @inbounds for p in model.poles
        ω02 = p.ω0^2
        ε  += p.Δε * ω02 / (ω02 - ω2 + 2.0im * ω * p.δ)
    end
    return ε
end

# ─────────────────────────────────────────────────────────────────────────────
# Vectorised evaluation (fast broadcast)
# ─────────────────────────────────────────────────────────────────────────────

"""
    eval_permittivity(model::MaterialModel, ωs::Vector{Float64}) → Vector{ComplexF64}

Vectorised evaluation of complex relative permittivity.
Equivalent to `[eval_permittivity(model, ω) for ω in ωs]` but avoids
per-element function-call overhead for large frequency sweeps.
"""
function eval_permittivity(model::MaterialModel, ωs::Vector{Float64})
    return [eval_permittivity(model, ω) for ω in ωs]
end

export DebyePole, LorentzPole
export DebyeModel, DrudeModel, LorentzModel
export eval_permittivity
