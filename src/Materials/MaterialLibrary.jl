"""
    MaterialLibrary.jl

Phase 19.3 — Electromagnetic material library with named entries, frequency
validity ranges, version metadata, and JLD2 serialisation.

Exported types:
  MaterialModel      — abstract base for all Phase-19 material models
  Isotropic          — homogeneous isotropic medium (εr, μr as ComplexF64)
  Anisotropic        — bi-anisotropic tensor medium (3×3 SMatrix)
  MaterialEntry      — single library entry (name, model, freq range, source)
  MaterialLibrary    — mutable Dict-backed catalogue

Exported functions:
  add_material!      — insert / overwrite an entry
  get_material       — look up model by name + frequency (range-checked)
  save_library       — serialize to JLD2
  load_library       — deserialize from JLD2
  load_builtin_library — return pre-populated standard library
"""

using StaticArrays
using LinearAlgebra
using JLD2

# ─────────────────────────────────────────────────────────────────────────────
# Abstract base
# ─────────────────────────────────────────────────────────────────────────────

"""
    MaterialModel

Abstract base type for electromagnetic material parameter models used in the
Phase-19 material library.  All subtypes carry enough data to evaluate ε(ω)
and μ(ω) at arbitrary frequency.
"""
abstract type MaterialModel end

# ─────────────────────────────────────────────────────────────────────────────
# Concrete frequency-independent types
# ─────────────────────────────────────────────────────────────────────────────

"""
    Isotropic(εr, μr=1)

Frequency-independent homogeneous isotropic medium described by complex
relative permittivity `εr` and relative permeability `μr` (both `ComplexF64`).

# Special values
- `Isotropic(complex(Inf), complex(1.0))` — Perfect Electric Conductor (PEC)
- `Isotropic(1.0, 1.0)` — free space / vacuum
"""
struct Isotropic <: MaterialModel
    εr :: ComplexF64
    μr :: ComplexF64
end

Isotropic(εr::Number, μr::Number = 1.0) = Isotropic(ComplexF64(εr), ComplexF64(μr))

"""
    Anisotropic(ε_tensor, μ_tensor)

Frequency-independent bi-anisotropic medium described by 3×3 complex relative
tensors.  Accepts any `AbstractMatrix` and promotes to `SMatrix{3,3,ComplexF64,9}`.

!!! note
    Requires `size(ε_tensor) == (3,3)` and `size(μ_tensor) == (3,3)`.
"""
struct Anisotropic <: MaterialModel
    ε_tensor :: SMatrix{3,3,ComplexF64,9}
    μ_tensor :: SMatrix{3,3,ComplexF64,9}
end

function Anisotropic(ε::AbstractMatrix, μ::AbstractMatrix)
    @assert size(ε) == (3, 3) "Anisotropic: ε_tensor must be 3×3, got $(size(ε))"
    @assert size(μ) == (3, 3) "Anisotropic: μ_tensor must be 3×3, got $(size(μ))"
    return Anisotropic(
        SMatrix{3,3,ComplexF64}(ComplexF64.(ε)),
        SMatrix{3,3,ComplexF64}(ComplexF64.(μ)),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Library entry and catalogue
# ─────────────────────────────────────────────────────────────────────────────

"""
    MaterialEntry(name, model, frequency_range, version, source)

A single entry in a `MaterialLibrary`.

# Fields
- `name`            : unique string identifier (e.g. `"fr4_standard"`)
- `model`           : the `MaterialModel` (frequency-independent)
- `frequency_range` : `(f_min, f_max)` in Hz over which the model is valid
- `version`         : `VersionNumber` tracking data provenance
- `source`          : human-readable provenance string
"""
struct MaterialEntry
    name            :: String
    model           :: MaterialModel
    frequency_range :: Tuple{Float64,Float64}   # Hz [f_min, f_max]
    version         :: VersionNumber
    source          :: String
end

"""
    MaterialLibrary()

Mutable catalogue of `MaterialEntry` objects, keyed by material name.

Use `add_material!`, `get_material`, `save_library`, and `load_library`
to manage contents.
"""
mutable struct MaterialLibrary
    entries :: Dict{String,MaterialEntry}
end
MaterialLibrary() = MaterialLibrary(Dict{String,MaterialEntry}())

# ─────────────────────────────────────────────────────────────────────────────
# CRUD
# ─────────────────────────────────────────────────────────────────────────────

"""
    add_material!(lib::MaterialLibrary, entry::MaterialEntry) → lib

Insert or overwrite `entry` in `lib`.  Returns `lib` for chaining.
"""
function add_material!(lib::MaterialLibrary, entry::MaterialEntry)
    lib.entries[entry.name] = entry
    return lib
end

"""
    get_material(lib::MaterialLibrary, name::String, freq::Real) → MaterialModel

Return the `MaterialModel` for `name` at frequency `freq` (Hz).

Throws `KeyError` if the material is not found; throws `DomainError` if `freq`
is outside the declared `frequency_range`.

!!! note
    Extrapolation outside the declared range is intentionally disallowed.
    Adjust `frequency_range` if your application requires a wider range.
"""
function get_material(lib::MaterialLibrary, name::String, freq::Real)
    haskey(lib.entries, name) ||
        throw(KeyError("Material $(repr(name)) not found in library"))
    entry = lib.entries[name]
    f = Float64(freq)
    flo, fhi = entry.frequency_range
    if f < flo || f > fhi
        throw(DomainError(
            f,
            "Frequency $f Hz is outside the declared range [$flo, $fhi] Hz " *
            "for material \"$(name)\". Extrapolation is not supported.",
        ))
    end
    return entry.model
end

# ─────────────────────────────────────────────────────────────────────────────
# Serialisation (JLD2)
# ─────────────────────────────────────────────────────────────────────────────

"""
    save_library(lib::MaterialLibrary, path::String; format::Symbol = :jld2)

Serialize `lib` to disk.  Only `:jld2` format is currently supported.
"""
function save_library(lib::MaterialLibrary, path::String; format::Symbol = :jld2)
    format === :jld2 ||
        error("save_library: unsupported format $(format). Only :jld2 is supported.")
    jldsave(path; library = lib)
    return nothing
end

"""
    load_library(path::String) → MaterialLibrary

Deserialize a `MaterialLibrary` from a JLD2 file written by `save_library`.
"""
function load_library(path::String)
    jldopen(path, "r") do fh
        return fh["library"]
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Built-in materials
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_builtin_library() → MaterialLibrary

Return a fresh `MaterialLibrary` pre-populated with five standard RF/microwave
materials:

| Name              | εr (relative)       | tanδ   |
|-------------------|---------------------|--------|
| `pec`             | ∞ (conductor)       | 0      |
| `vacuum`          | 1.0                 | 0      |
| `fr4_standard`    | 4.4                 | 0.020  |
| `rogers_ro4003c`  | 3.55                | 0.0027 |
| `silicon`         | 11.9                | 0.005  |

The loss tangent `tanδ` is encoded as `εr_complex = εr_real * (1 - j*tanδ)`
following the `e^{-jωt}` CEM time convention.
"""
function load_builtin_library()
    lib = MaterialLibrary()

    # PEC — perfect electric conductor (sentinel: real(εr) = Inf)
    add_material!(lib, MaterialEntry(
        "pec",
        Isotropic(complex(Inf, 0.0), complex(1.0, 0.0)),
        (0.0, Inf),
        v"1.0.0",
        "model",
    ))

    # Vacuum / free space
    add_material!(lib, MaterialEntry(
        "vacuum",
        Isotropic(1.0, 1.0),
        (0.0, Inf),
        v"1.0.0",
        "model",
    ))

    # FR-4 standard laminate: εr = 4.4, tanδ = 0.020
    add_material!(lib, MaterialEntry(
        "fr4_standard",
        Isotropic(4.4 * (1.0 - 0.020im), 1.0),
        (1e6, 20e9),
        v"1.0.0",
        "reference:IPC-4101",
    ))

    # Rogers RO4003C high-frequency laminate: εr = 3.55, tanδ = 0.0027
    add_material!(lib, MaterialEntry(
        "rogers_ro4003c",
        Isotropic(3.55 * (1.0 - 0.0027im), 1.0),
        (1e6, 40e9),
        v"1.0.0",
        "reference:Rogers_RO4003C_datasheet",
    ))

    # Silicon (lightly doped): εr = 11.9, tanδ = 0.005
    add_material!(lib, MaterialEntry(
        "silicon",
        Isotropic(11.9 * (1.0 - 0.005im), 1.0),
        (1e6, 100e9),
        v"1.0.0",
        "reference:Pozar2011",
    ))

    return lib
end

export MaterialModel, Isotropic, Anisotropic
export MaterialEntry, MaterialLibrary
export add_material!, get_material
export save_library, load_library, load_builtin_library
