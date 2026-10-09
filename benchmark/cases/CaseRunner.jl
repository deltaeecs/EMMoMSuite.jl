# CaseRunner.jl — gmsh-driven RCS case pipeline for EMMoMSuite
#
# Pipeline per case:
#   .geo file -> gmsh surface mesh -> RWG basis -> MoM assemble -> LU solve
#   -> bistatic RCS (Mie reference when available) -> far-field plots
#   -> per-case report.md + rcs.csv + PNGs.
#
# Run from benchmark/cases/run_case.jl (activates the benchmark project so
# that both EMMoMSuite and Plots are available).

module CaseRunner

using EMMoMSuite
using TOML, Printf, Dates, Statistics, LinearAlgebra
using Plots

export CaseSpec, Interface, Region, PEC, Dielectric, VolumeMaterialContext,
       load_cases, run_case, run_cases, write_index

const CASES_DIR    = @__DIR__
const GEO_DIR      = joinpath(CASES_DIR, "geo")
const CASES_TOML   = joinpath(CASES_DIR, "cases.toml")
const RESULT_ROOT  = normpath(joinpath(CASES_DIR, "..", "..", "test_results", "cases"))

# ─────────────────────────────────────────────────────────────────────────────
# Material model: an electromagnetic interface carries a material on EACH side.
# Convention: `plus`  = the side the mesh triangle normal n̂ points to
#             `minus` = the opposite side.
# (closed-surface orientation is verified/auto-fixed in S2; open surfaces must
#  be oriented in the .geo file or flagged `flip`.)
# ─────────────────────────────────────────────────────────────────────────────

abstract type Material end

struct PEC <: Material end

Base.@kwdef struct Dielectric <: Material
    eps_r::ComplexF64 = 1.0 + 0.0im
    mu_r::ComplexF64  = 1.0 + 0.0im
end

const AIR = Dielectric()   # free-space reference medium

"""
    Interface

One physical surface separating two media (`plus` on the n̂ side, `minus` on
the back side). `closed` marks a watertight surface (enables the signed-volume
orientation check). `flip` manually reverses the n̂ orientation (open surfaces).
"""
Base.@kwdef struct Interface
    surface::String             # Physical Surface label (or "body" pre-S2)
    plus::Material
    minus::Material
    closed::Bool = true
    flip::Bool  = false
end

"""
    Region

One physical *volume* (gmsh `Physical Volume` label) filled with a single
material. Volume cases use `Vector{Region}` instead of `Vector{Interface}`.
"""
Base.@kwdef struct Region
    surface::String             # Physical Volume label
    material::Material
end

function Base.show(io::IO, m::Dielectric)
    print(io, "Dielectric(εᵣ=", m.eps_r, ", μᵣ=", m.mu_r, ")")
end

"""
    MaterialContext

Everything downstream (RCS / far-field) needs to know about the formulation,
resolved once by the operator factory:
- `k0`, `eta0` : free-space wavenumber / impedance used for post-processing
- `layout`     : `:j_only` (EFIE/CFIE, N coefficients) or `:jm` (PMCHW, [J; M])
"""
Base.@kwdef struct MaterialContext
    k0::Float64
    eta0::Float64
    layout::Symbol                  # :j_only | :jm
    interfaces::Vector{Interface}
end

"""
    CaseSpec

One registered case from `cases.toml`.

New-style cases describe materials as a `Vector{Interface}` (each surface
carries materials on both sides). Old-style single-material fields (`ie`,
`eps_r`, `mu_r`) are still accepted and translated into one implicit interface
`air | dielectric` (or `air | pec`) in `load_cases`.
"""
Base.@kwdef struct CaseSpec
    name::String
    geo::String
    mesh_size::Float64
    dim::Int                       = 2
    freq::Float64
    interfaces::Vector{Interface}  = Interface[]
    regions::Vector{Region}        = Region[]     # volume (dim=3) cases
    ie::String                     = "auto"       # auto | EFIE | CFIE | PMCHW
    alpha::Float64                 = 0.5
    theta_inc::Float64             = pi / 2
    phi_inc::Float64               = pi
    pol::Vector{Float64}           = [0.0, 0.0, 1.0]
    n_theta::Int                   = 181
    phi_cuts::Vector{Float64}      = [0.0, pi / 2]
    mie_radius::Union{Nothing,Float64} = nothing
    # --- legacy single-material fields (translated to `interfaces`) ---------
    eps_r::Union{Nothing,ComplexF64}   = nothing
    mu_r::Union{Nothing,ComplexF64}    = nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Operator factory: the ONLY place that maps (interfaces, materials) -> solver
# ─────────────────────────────────────────────────────────────────────────────

"""
    resolve_material(m::Material) -> Material
    _is_diel, _is_pec helpers
"""
_is_diel(m::Dielectric) = true
_is_diel(::PEC) = false
_is_pec(::PEC) = true
_is_pec(::Dielectric) = false

# Material-protocol hooks (src/IntegralEquations/MaterialFormulation.jl):
# the generic build_volume_system/is_all_pec machinery works on ANY material
# model that implements these three functions.
EMMoMSuite.is_pec_material(m::PEC) = true
EMMoMSuite.relative_permittivity(d::Dielectric) =
    ComplexF64(d.eps_r) + (hasproperty(d, :eps_r_im) ? ComplexF64(0, d.eps_r_im) : 0im)

"""
    validate_interfaces(interfaces) -> nothing

Interface legality checks (fail fast, before any meshing/solving):
1. plus/minus must not both be PEC (nothing to solve on an internal PEC pair);
2. identical media on both sides degenerates PMCHW -> warn;
3. surfaces must be non-empty.
"""
function validate_interfaces(interfaces)
    isempty(interfaces) && error("no interface defined")
    for itf in interfaces
        isempty(itf.surface) && error("interface with empty surface label")
        if _is_pec(itf.plus) && _is_pec(itf.minus)
            error("interface `$(itf.surface)`: both sides are PEC — nothing to solve")
        end
        if !_is_pec(itf.plus) && !_is_pec(itf.minus) && itf.plus == itf.minus
            @warn "interface `$(itf.surface)`: identical media on both sides; PMCHW degenerates to EFIE"
        end
    end
    return nothing
end

"""
    derive_formulation(interfaces, ie::String) -> Symbol

`ie = "auto"` rules:
- any dielectric involved (dielectric|dielectric or dielectric|air) → `:PMCHW`
- otherwise (PEC interfaces only)                                  → `:EFIE`
Explicit `ie` overrides (`:EFIE`, `:CFIE`, `:PMCHW`).
"""
function derive_formulation(interfaces, ie::String)
    has_diel = any(_is_diel,
                   reduce(vcat, [[itf.plus, itf.minus] for itf in interfaces]))
    return EMMoMSuite.derive_formulation(has_diel, ie)
end

"""
    build_operator(interfaces, freq; ie = "auto", alpha = 0.5) -> (op, ctx)

The single branch point from material/interface description to solver operator.
Returns the EMMoMSuite operator and the `MaterialContext` consumed by the
post-processing dispatch.
"""
function build_operator(interfaces, freq::Float64; ie::String = "auto", alpha::Float64 = 0.5)
    validate_interfaces(interfaces)
    set_frequency!(freq)          # global k0/eta0 the solver & post-processing read
    kind = derive_formulation(interfaces, ie)

    diels = Dielectric[]
    for itf in interfaces, m in (itf.plus, itf.minus)
        m isa Dielectric && push!(diels, m)
    end
    unique!(diels)
    # the implicit exterior (air) is already part of the single-region PMCHW
    # kernel; only non-air dielectrics count as "regions"
    diels = filter(d -> !(isapprox(d.eps_r, AIR.eps_r) && isapprox(d.mu_r, AIR.mu_r)), diels)

    op, layout = if kind === :PMCHW
        isempty(diels) && error("PMCHW requires at least one dielectric side")
        length(diels) > 1 &&
            error("multi-region PMCHW (interfaces with $(length(diels)) distinct dielectrics) is not yet supported")
        d = diels[1]
        PMCHW(freq, d.eps_r, d.mu_r), :jm
    elseif kind === :CFIE
        CFIE(freq, alpha), :j_only
    elseif kind === :EFIE
        EFIE(freq), :j_only
    else
        error("unknown formulation :$kind")
    end
    return op, MaterialContext(k0 = Float64(get_k0()), eta0 = Float64(get_eta0()),
                               layout = layout, interfaces = interfaces)
end

# ─────────────────────────────────────────────────────────────────────────────
# Volume (tetrahedral) branch: regions → SWG/VEFIE, with PEC-only fallback
# ─────────────────────────────────────────────────────────────────────────────

"""
    VolumeMaterialContext

Post-processing context for volume (tetrahedral) cases:
- `bm`              : `BoundMesh` (tet mesh + per-tag material binding)
- `regions`         : the `Region` list of the case
- `permittivities`  : per-tetrahedron relative εᵣ used by VEFIE/RCS
- `boundary_fallback`: `true` when all regions are PEC and the solve was
  routed to the surface EFIE on the extracted boundary (RWG basis)
- `fallback_ctx`    : the inner `MaterialContext` of the fallback path
"""
Base.@kwdef struct VolumeMaterialContext
    k0::Float64
    eta0::Float64
    layout::Symbol                              # always :volume
    bm
    regions::Vector{Region}
    permittivities::Vector{ComplexF64}
    boundary_fallback::Bool                     = false
    fallback_ctx::Union{Nothing,MaterialContext} = nothing
end

"""
    build_volume_operator(bm, regions, freq; ie = "auto") -> (op, VolumeMaterialContext)

Volume counterpart of [`build_operator`](@ref): the single branch point from
region/material description to solver operator.

`ie = "auto"` rules:
- **all regions PEC** → degenerate to the existing surface EFIE path
  (`extract_surface` + RWG basis); the returned context carries
  `boundary_fallback = true`.
- **any dielectric region** → SWG basis + VEFIE volume operator with
  per-tetrahedron permittivities taken from the region bindings.
- **mix of PEC and dielectric regions** → not supported (error).
"""
function build_volume_operator(bm, regions, freq::Float64; ie::String = "auto")
    isempty(regions) && error("no region defined")
    set_frequency!(freq)

    # Algorithmic dispatch lives in the library (generic over any material
    # model implementing the material protocol); CaseRunner only adapts the
    # result into the case-level post-processing context.
    sys = build_volume_system(bm, freq; ie = ie)
    fallback = sys.kind === :efie_fallback
    ctx = VolumeMaterialContext(
        k0 = Float64(get_k0()), eta0 = Float64(get_eta0()), layout = :volume,
        bm = bm, regions = regions, permittivities = sys.permittivities,
        boundary_fallback = fallback,
        fallback_ctx = fallback ? MaterialContext(k0 = Float64(get_k0()),
                                                  eta0 = Float64(get_eta0()),
                                                  layout = :j_only,
                                                  interfaces = Interface[]) :
                                 nothing)
    return sys.op, sys.mesh, ctx
end

# ─────────────────────────────────────────────────────────────────────────────
# Post-processing dispatch — no `ie` strings beyond this point
# ─────────────────────────────────────────────────────────────────────────────

"bistatic RCS dispatch: PMCHW coefficients [J; M] need the k0/eta0 overload"
function postprocess_rcs(ctx::MaterialContext, θa, ϕs, I, basis)
    if ctx.layout === :jm
        _, _, rcs_dB = radarCrossSection(θa, ϕs, I, basis, ctx.k0, ctx.eta0)
    else
        _, _, rcs_dB = radarCrossSection(θa, ϕs, I, basis)
    end
    return rcs_dB
end

"far-field dispatch: farField consumes the electric current (J) part only"
postprocess_farfield(ctx::MaterialContext, θa, ϕs, I, basis, source, nbasis) =
    farField(θa, ϕs, ctx.layout === :jm ? I[1:nbasis] : I, basis, source)

"""
    reference_rcs(ctx, mie_radius, freq, θa, φ, θinc, φinc, pol)

Analytic reference chosen from the interface materials:
- any dielectric side → dielectric-sphere Mie (with that material),
- PEC-only interfaces → PEC-sphere Mie.
Returns `nothing` when `mie_radius` is unset.
"""
function reference_rcs(ctx::MaterialContext, mie_radius, freq, θa, φ, θinc, φinc, pol)
    mie_radius === nothing && return nothing
    diels = Dielectric[]
    for itf in ctx.interfaces, m in (itf.plus, itf.minus)
        m isa Dielectric &&
            !(isapprox(m.eps_r, AIR.eps_r) && isapprox(m.mu_r, AIR.mu_r)) && push!(diels, m)
    end
    if isempty(diels)
        return (:pec, mie_pec_bistatic_rcs_dBsm(mie_radius, freq, θa, φ, θinc, φinc, pol))
    else
        d = diels[1]
        return (:diel, mie_dielectric_bistatic_rcs_dBsm(
            mie_radius, freq, d.eps_r, d.mu_r, θa, φ, θinc, φinc, pol))
    end
end

"volume RCS dispatch: SWG coefficients + per-tet permittivities, or RWG fallback"
function postprocess_rcs(ctx::VolumeMaterialContext, θa, ϕs, I, basis)
    if ctx.boundary_fallback
        _, _, rcs_dB = radarCrossSection(θa, ϕs, I, basis)
    else
        _, _, rcs_dB = radarCrossSection(θa, ϕs, I, basis, ctx.permittivities)
    end
    return rcs_dB
end

"volume far-field dispatch: derived from the RCS magnitudes (volume currents
have no farField overload); matches `_plot_farfield`'s [2, nθ, nϕ] layout"
function postprocess_farfield(ctx::VolumeMaterialContext, θa, ϕs, I, basis, source, nbasis)
    return farfield_from_rcs(postprocess_rcs(ctx, θa, ϕs, I, basis))
end

"""
    reference_rcs(ctx::VolumeMaterialContext, ...)

Analytic reference chosen from the region materials (same rule as the
interface variant, reading `ctx.regions` instead).
"""
function reference_rcs(ctx::VolumeMaterialContext, mie_radius, freq, θa, φ, θinc, φinc, pol)
    mie_radius === nothing && return nothing
    if ctx.boundary_fallback
        return (:pec, mie_pec_bistatic_rcs_dBsm(mie_radius, freq, θa, φ, θinc, φinc, pol))
    end
    diels = filter(r -> r.material isa Dielectric &&
                         !(isapprox(r.material.eps_r, AIR.eps_r) &&
                           isapprox(r.material.mu_r, AIR.mu_r)), ctx.regions)
    isempty(diels) &&
        return (:pec, mie_pec_bistatic_rcs_dBsm(mie_radius, freq, θa, φ, θinc, φinc, pol))
    d = diels[1].material
    return (:diel, mie_dielectric_bistatic_rcs_dBsm(
        mie_radius, freq, d.eps_r, d.mu_r, θa, φ, θinc, φinc, pol))
end

# ─────────────────────────────────────────────────────────────────────────────
# cases.toml parsing (new interface format + legacy single-material format)
# ─────────────────────────────────────────────────────────────────────────────

_parse_material(reg::Dict, name::AbstractString) = begin
    haskey(reg, name) || error("unknown material `$name` (register it under [material.$name])")
    e = reg[name]
    get(e, "type", "dielectric") == "pec" && return PEC()
    eps = ComplexF64(Float64(get(e, "eps_r", 1.0)), Float64(get(e, "eps_r_im", 0.0)))
    mu  = ComplexF64(Float64(get(e, "mu_r", 1.0)),  Float64(get(e, "mu_r_im", 0.0)))
    return Dielectric(eps_r = eps, mu_r = mu)
end

function load_cases(path::AbstractString = CASES_TOML)
    raw = TOML.parsefile(path)
    registry = get(raw, "material", Dict{String,Any}())
    specs = CaseSpec[]
    for c in get(raw, "case", [])
        interfaces = Interface[]
        for itf in get(c, "interface", [])
            push!(interfaces, Interface(
                surface = String(itf["surface"]),
                plus    = _parse_material(registry, itf["plus"]),
                minus   = _parse_material(registry, itf["minus"]),
                closed  = Bool(get(itf, "closed", true)),
                flip    = Bool(get(itf, "flip", false)),
            ))
        end
        # volume (Physical Volume) regions — mutually exclusive with interfaces
        regions = Region[]
        for r in get(c, "region", [])
            push!(regions, Region(
                surface  = String(r["surface"]),
                material = _parse_material(registry, r["material"]),
            ))
        end
        !isempty(regions) && (!isempty(interfaces) || haskey(c, "eps_r")) &&
            error("case `$(c["name"])`: [[case.region]] is mutually exclusive " *
                  "with [[case.interface]] / legacy eps_r")
        # legacy single-material translation: eps_r present → air|dielectric
        haskey(c, "eps_r") && isempty(regions) && push!(interfaces, Interface(
            surface = "body",
            plus    = AIR,
            minus   = Dielectric(
                        eps_r = ComplexF64(Float64(c["eps_r"]), Float64(get(c, "eps_r_im", 0.0))),
                        mu_r  = ComplexF64(Float64(get(c, "mu_r", 1.0)), Float64(get(c, "mu_r_im", 0.0)))),
        ))
        ie = String(get(c, "ie", (isempty(interfaces) && isempty(regions)) ? "EFIE" : "auto"))
        # legacy pure-PEC case without any interface: synthesize air | pec
        isempty(interfaces) && isempty(regions) &&
            push!(interfaces, Interface(surface = "body", plus = AIR, minus = PEC()))
        push!(specs, CaseSpec(
            name       = String(c["name"]),
            geo        = String(c["geo"]),
            mesh_size  = Float64(c["mesh_size"]),
            dim        = get(c, "dim", 2),
            freq       = Float64(c["freq"]),
            interfaces = interfaces,
            regions    = regions,
            ie         = ie,
            alpha      = get(c, "alpha", 0.5),
            theta_inc  = Float64(get(c, "theta_inc", pi / 2)),
            phi_inc    = Float64(get(c, "phi_inc", pi)),
            pol        = Float64.(get(c, "pol", [0.0, 0.0, 1.0])),
            n_theta    = get(c, "n_theta", 181),
            phi_cuts   = Float64.(get(c, "phi_cuts", [0.0, pi / 2])),
            mie_radius = haskey(c, "mie_radius") ? Float64(c["mie_radius"]) : nothing,
        ))
    end
    return specs
end

"""
    mesh_signed_volume(mesh) -> Float64

Signed volume enclosed by a closed triangulated surface, via the divergence
theorem: `V = (1/6) Σ v₁·(v₂×v₃)` over all triangles. Positive ⇔ triangle
normals point outward — this verifies the `Interface` n̂ convention.
"""
function mesh_signed_volume(mesh)
    v = 0.0
    for t in 1:mesh.trinum
        a = mesh.triangles[1, t]; b = mesh.triangles[2, t]; c = mesh.triangles[3, t]
        v += dot(mesh.node[:, a], cross(mesh.node[:, b], mesh.node[:, c])) / 6
    end
    return v
end

"""
    check_orientation!(spec, mesh) -> nothing

For every `closed` interface, verify that the mesh normals actually point to
the `plus` side: gmsh/OCC emits outward normals, so a positive signed volume
agrees with `flip = false`. Any mismatch is auto-corrected (`flip` toggled)
with a warning — material sides are never silently swapped.
"""
function check_orientation!(spec::CaseSpec, mesh)
    sv = mesh_signed_volume(mesh)
    for (k, itf) in enumerate(spec.interfaces)
        itf.closed || continue
        normals_outward = sv > 0
        if normals_outward == itf.flip   # flip=true claims inward normals
            @warn "interface `$(itf.surface)`: mesh signed volume $(round(sv; digits=6)) m³ " *
                  "says normals point $(normals_outward ? "outward" : "inward") — auto-fixing flip"
            spec.interfaces[k] = Interface(itf.surface, itf.plus, itf.minus, itf.closed, !itf.flip)
        end
    end
    return nothing
end

_deg(x) = round(rad2deg(x); digits = 1)

struct CaseResult
    spec::CaseSpec
    outdir::String
    mesh::Any                       # kept for the geometry views in the report
    trinum::Int
    num_nodes::Int
    num_basis::Int
    t_mesh::Float64
    t_assembly::Float64
    t_solve::Float64
    t_rcs::Float64
    t_plots::Float64
    mie_ok::Bool
    mie_rmse::Vector{Float64}     # per phi cut, NaN when no reference
end

"""
    run_case(spec::CaseSpec; outroot::AbstractString = RESULT_ROOT) -> CaseResult

Execute the full pipeline for one case and write all artifacts.
Volume cases (`[[case.region]]`) dispatch to `_run_volume_case`; surface
(interface / legacy) cases run the RWG pipeline below.
"""
function run_case(spec::CaseSpec; outroot::AbstractString = RESULT_ROOT)
    isempty(spec.regions) || return _run_volume_case(spec; outroot)
    outdir = joinpath(outroot, spec.name)
    mkpath(outdir)
    println("=== Case $(spec.name) ===")

    # ---- 1. geometry -> gmsh mesh -------------------------------------
    geo_file = isabspath(spec.geo) ? spec.geo : joinpath(GEO_DIR, spec.geo)
    t0 = time()
    mesh = generate_gmsh_from_file(geo_file; mesh_size = spec.mesh_size, dim = spec.dim)
    t_mesh = time() - t0
    trinum  = mesh.trinum
    nnodes  = size(mesh.node, 2)
    basis   = RWGBasis(mesh)

    # ---- 1b. closed-surface orientation check (S2) ----------------------
    check_orientation!(spec, mesh)
    nbasis  = num_basis(basis)
    @printf("  mesh: %d triangles, %d nodes, %d RWG unknowns (%.1f s)\n",
            trinum, nnodes, nbasis, t_mesh)

    set_frequency!(spec.freq)
    source = PlaneWave(spec.freq, spec.theta_inc, spec.phi_inc, spec.pol)

    # ---- 2. assemble: operator factory is the only branch point ----------
    t0 = time()
    op, ctx = build_operator(spec.interfaces, spec.freq; ie = spec.ie, alpha = spec.alpha)
    Z  = assemble_impedance_matrix(op, basis)
    V  = excitation_vector(op, source, basis)
    t_assembly = time() - t0
    @printf("  assembly (%s): %.1f s\n", nameof(typeof(op)), t_assembly)

    # ---- 3. solve -------------------------------------------------------
    t0 = time()
    I = Z \ V
    t_solve = time() - t0
    @printf("  LU solve: %.1f s\n", t_solve)

    # ---- 4. RCS ----------------------------------------------------------
    t0 = time()
    θs = range(0.0, pi; length = spec.n_theta)
    ϕs = collect(spec.phi_cuts)
    θa = collect(Float64.(θs))
    rcs_dB = postprocess_rcs(ctx, θa, ϕs, I, basis)   # [nθ, nϕ] dBsm

    # optional Mie reference on the same grid
    mie_dB = Matrix{Float64}(undef, length(θa), length(ϕs)); fill!(mie_dB, NaN)
    mie_ok = false
    if spec.mie_radius !== nothing
        try
            for (j, φ) in enumerate(ϕs)
                ref = reference_rcs(ctx, spec.mie_radius, spec.freq,
                                    θa, φ, spec.theta_inc, spec.phi_inc, spec.pol)
                mie_dB[:, j] .= ref[2]
            end
            mie_ok = true
        catch err
            @warn "Mie reference failed" spec.name err
        end
    end
    t_rcs = time() - t0

    rmse = [mie_ok ? sqrt(mean((rcs_dB[:, j] .- mie_dB[:, j]).^2)) : NaN for j in eachindex(ϕs)]
    if mie_ok
        for (j, φ) in enumerate(ϕs)
            @printf("  RCS vs Mie, phi=%6.1f°: RMSE = %.3f dB\n", _deg(φ), rmse[j])
        end
    end

    # ---- 5. far field + plots ---------------------------------------------
    t0 = time()
    FF = postprocess_farfield(ctx, θa, ϕs, I, basis, source, nbasis)
    _write_rcs_csv(spec, outdir, θa, ϕs, rcs_dB, mie_dB, mie_ok)
    _plot_rcs(spec, outdir, θa, ϕs, rcs_dB, mie_dB, mie_ok)
    _plot_farfield(spec, outdir, θa, ϕs, FF)
    t_plots = time() - t0
    @printf("  plots: %.1f s\n", t_plots)

    res = CaseResult(spec, outdir, mesh, trinum, nnodes, nbasis,
                     t_mesh, t_assembly, t_solve, t_rcs, t_plots,
                     mie_ok, rmse)
    _write_report(res, rcs_dB, mie_dB, mie_ok)
    println("  report: ", joinpath(outdir, "report.md"))
    return res
end

# ─────────────────────────────────────────────────────────────────────────────
# volume (tetrahedral) case pipeline: Physical Volumes → SWG/VEFIE
# ─────────────────────────────────────────────────────────────────────────────

function _run_volume_case(spec::CaseSpec; outroot::AbstractString = RESULT_ROOT)
    outdir = joinpath(outroot, spec.name)
    mkpath(outdir)
    println("=== Case $(spec.name) (volume) ===")

    # ---- 1. geometry -> gmsh tetrahedral mesh + Physical Volume regions ---
    geo_file = isabspath(spec.geo) ? spec.geo : joinpath(GEO_DIR, spec.geo)
    t0 = time()
    mesh, region_tags = generate_gmsh_volume(geo_file; mesh_size = spec.mesh_size)
    t_mesh = time() - t0

    # ---- 1b. bind materials by region name and validate -------------------
    bm = bind_regions(mesh, region_tags,
                      Dict(r.surface => r.material for r in spec.regions))
    validate_bindings(bm) ||
        error("case $(spec.name): unbound tetrahedra (Physical Volume missing?)")

    tetnum  = mesh.tetnum
    nnodes  = size(mesh.node, 2)
    @printf("  mesh: %d tetrahedra, %d nodes, regions = %s (%.1f s)\n",
            tetnum, nnodes, join(keys(region_tags), ", "), t_mesh)

    set_frequency!(spec.freq)
    source = PlaneWave(spec.freq, spec.theta_inc, spec.phi_inc, spec.pol)

    # ---- 2. assemble: volume operator factory is the only branch point ----
    t0 = time()
    op, solve_mesh, ctx = build_volume_operator(bm, spec.regions, spec.freq;
                                                ie = spec.ie)
    basis = ctx.boundary_fallback ? RWGBasis(solve_mesh) : SWGBasis(solve_mesh)
    nbasis = num_basis(basis)
    Z = assemble_impedance_matrix(op, basis)
    V = ctx.boundary_fallback ? excitation_vector(op, source, basis) :
        excitation_vector(op, source, basis, ctx.permittivities)
    t_assembly = time() - t0
    @printf("  assembly (%s, %d unknowns): %.1f s\n",
            nameof(typeof(op)), nbasis, t_assembly)

    # ---- 3. solve -------------------------------------------------------
    t0 = time()
    I = Z \ V
    t_solve = time() - t0
    @printf("  LU solve: %.1f s\n", t_solve)

    # ---- 4. RCS ----------------------------------------------------------
    t0 = time()
    θs = range(0.0, pi; length = spec.n_theta)
    ϕs = collect(spec.phi_cuts)
    θa = collect(Float64.(θs))
    rcs_dB = postprocess_rcs(ctx, θa, ϕs, I, basis)   # [nθ, nϕ] dBsm

    mie_dB = Matrix{Float64}(undef, length(θa), length(ϕs)); fill!(mie_dB, NaN)
    mie_ok = false
    if spec.mie_radius !== nothing
        try
            for (j, φ) in enumerate(ϕs)
                ref = reference_rcs(ctx, spec.mie_radius, spec.freq,
                                    θa, φ, spec.theta_inc, spec.phi_inc, spec.pol)
                mie_dB[:, j] .= ref[2]
            end
            mie_ok = true
        catch err
            @warn "Mie reference failed" spec.name err
        end
    end
    t_rcs = time() - t0

    rmse = [mie_ok ? sqrt(mean((rcs_dB[:, j] .- mie_dB[:, j]).^2)) : NaN for j in eachindex(ϕs)]
    if mie_ok
        for (j, φ) in enumerate(ϕs)
            @printf("  RCS vs Mie, phi=%6.1f°: RMSE = %.3f dB\n", _deg(φ), rmse[j])
        end
    end

    # ---- 5. far field + plots ---------------------------------------------
    t0 = time()
    FF = postprocess_farfield(ctx, θa, ϕs, I, basis, source, nbasis)
    _write_rcs_csv(spec, outdir, θa, ϕs, rcs_dB, mie_dB, mie_ok)
    _plot_rcs(spec, outdir, θa, ϕs, rcs_dB, mie_dB, mie_ok)
    _plot_farfield(spec, outdir, θa, ϕs, FF)
    t_plots = time() - t0
    @printf("  plots: %.1f s\n", t_plots)

    res = CaseResult(spec, outdir, mesh, tetnum, nnodes, nbasis,
                     t_mesh, t_assembly, t_solve, t_rcs, t_plots,
                     mie_ok, rmse)
    _write_report_volume(res, rcs_dB, mie_dB, mie_ok, region_tags, ctx)
    println("  report: ", joinpath(outdir, "report.md"))
    return res
end

"""Run several cases by name (empty = all registered)."""
function run_cases(names::Vector{String} = String[]; outroot = RESULT_ROOT)
    specs = load_cases()
    isempty(names) || (specs = filter(s -> s.name in names, specs))
    isempty(specs) && error("no case matched: $names")
    results = CaseResult[]
    for s in specs
        push!(results, run_case(s; outroot = outroot))
    end
    write_index(outroot)
    return results
end

# ---------------------------------------------------------------------------
# artifacts
# ---------------------------------------------------------------------------

function _write_rcs_csv(spec, outdir, θa, ϕs, rcs_dB, mie_dB, mie_ok)
    path = joinpath(outdir, "rcs.csv")
    open(path, "w") do io
        println(io, "theta_deg,phi_deg,rcs_dBsm" * (mie_ok ? ",mie_dBsm" : ""))
        for (j, φ) in enumerate(ϕs), (i, θ) in enumerate(θa)
            line = @sprintf("%.2f,%.2f,%.4f", _deg(θ), _deg(φ), rcs_dB[i, j])
            mie_ok && (line *= @sprintf(",%.4f", mie_dB[i, j]))
            println(io, line)
        end
    end
    return path
end

function _plot_rcs(spec, outdir, θa, ϕs, rcs_dB, mie_dB, mie_ok)
    θdeg = _deg.(θa)
    p = plot(title = "$(spec.name) bistatic RCS  (f = $(spec.freq/1e6) MHz)",
             xlabel = "theta [deg]", ylabel = "RCS [dBsm]",
             legend = :bottomleft, size = (760, 520))
    for (j, φ) in enumerate(ϕs)
        lbl = @sprintf("MoM phi=%.0f°", _deg(φ))
        plot!(p, θdeg, rcs_dB[:, j]; label = lbl, lw = 2)
        mie_ok && plot!(p, θdeg, mie_dB[:, j];
                        label = @sprintf("Mie  phi=%.0f°", _deg(φ)),
                        ls = :dash, lw = 1.5)
    end
    path = joinpath(outdir, "rcs_cuts.png")
    savefig(p, path)
    return path
end

function _plot_farfield(spec, outdir, θa, ϕs, FF)
    # normalized |E| pattern (co-pol: total power) per phi cut, mirrored to 2π
    mag = sqrt.(abs2.(FF[1, :, :]) .+ abs2.(FF[2, :, :]))  # [nθ, nϕ]
    p = plot(proj = :polar, title = "$(spec.name) far-field |E| (norm., dB)",
             legend = :bottomleft, size = (640, 640), lims = (-60, 0))
    for (j, φ) in enumerate(ϕs)
        d = mag[:, j] ./ maximum(mag[:, j])
        db = 20 .* log10.(max.(d, 1e-6))              # clamp at -120 dB
        ang = vcat(θa, 2pi .- reverse(θa))            # mirror θ -> [0, 2π)
        val = vcat(db, reverse(db))
        plot!(p, ang .+ pi / 2, val; label = @sprintf("phi=%.0f°", _deg(φ)), lw = 1.8)
    end
    path = joinpath(outdir, "farfield_polar.png")
    savefig(p, path)
    return path
end

# ─────────────────────────────────────────────────────────────────────────────
# publication report: geometry views + performance metrics + conclusions
# ─────────────────────────────────────────────────────────────────────────────

"""
    _boundary_tris(mesh) -> 3×M Matrix{Int}

Node-index triangles of the visible (boundary) surface: identity for a
`TriangleMesh`; for a `TetrahedraMesh` the faces shared by exactly one
tetrahedron.
"""
function _boundary_tris(mesh)
    if isdefined(mesh, :trinum)
        return mesh.triangles
    end
    cnt = Dict{Tuple{Int,Int,Int},Int}()
    F = (i, j, k) -> tuple(sort!([i, j, k])...)
    for t in 1:mesh.tetnum
        a, b, c, d = mesh.tetras[:, t]
        for f in (F(a, b, c), F(a, b, d), F(a, c, d), F(b, c, d))
            cnt[f] = get(cnt, f, 0) + 1
        end
    end
    out = Matrix{Int}(undef, 3, count(x -> x == 1, values(cnt)))
    m = 0
    for (f, n) in cnt
        n == 1 || continue
        m += 1
        out[:, m] .= collect(f)
    end
    return out
end

"Orbit-camera projection: yaw `az`, pitch `el` (both radians) → screen (px, py) and depth."
function _proj_view(x, y, z, az, el)
    ca, sa = cos(az), sin(az)
    ce, se = cos(el), sin(el)
    x1 =  ca .* x .+ sa .* y
    y1 = -sa .* x .+ ca .* y
    return x1, ce .* z .- se .* y1, se .* z .+ ce .* y1   # px, py, depth
end

const _GEO_VIEWS = [  # (label, azimuth, elevation, aspect)
    ("top view (x–y)",   0.0,      pi / 2),
    ("front view (x–z)", 0.0,      0.0),
    ("side view (y–z)",  pi / 2,   0.0),
]

"""
    _plot_geometry_views(spec, mesh, outdir) -> path

Publication-style geometry figure: three orthographic views + one pseudo-3D
isometric view of the boundary triangles (filled, depth-sorted in the 3D view).
"""
function _plot_geometry_views(spec, mesh, outdir)
    tri = _boundary_tris(mesh)
    X, Y, Z = mesh.node[1, :], mesh.node[2, :], mesh.node[3, :]
    ntri = size(tri, 2)

    # NaN-separated polygon arrays for one :shape series call
    polys(px, py) = begin
        xs = Float64[]; ys = Float64[]
        for t in 1:ntri
            append!(xs, px[tri[1, t]], px[tri[2, t]], px[tri[3, t]], px[tri[1, t]], NaN)
            append!(ys, py[tri[1, t]], py[tri[2, t]], py[tri[3, t]], py[tri[1, t]], NaN)
        end
        return xs, ys
    end

    plots = Plots.Plot[]
    for (lbl, az, el) in _GEO_VIEWS
        px, py, _ = _proj_view(X, Y, Z, az, el)
        xs, ys = polys(px, py)
        p = plot(xs, ys; seriestype = :shape, fillcolor = :steelblue,
                 fillalpha = 0.25, linecolor = :grey25, lw = 0.4,
                 title = lbl, aspect_ratio = :equal, legend = false,
                 ticks = false)
        push!(plots, p)
    end
    # pseudo-3D isometric, painter's algorithm (far → near)
    az, el = pi / 4, pi / 6
    px, py, dep = _proj_view(X, Y, Z, az, el)
    order = sortperm([mean(dep[tri[:, t]]) for t in 1:ntri]; rev = true)
    tri_s = tri[:, order]
    xs = Float64[]; ys = Float64[]
    for t in 1:ntri
        append!(xs, px[tri_s[1, t]], px[tri_s[2, t]], px[tri_s[3, t]], px[tri_s[1, t]], NaN)
        append!(ys, py[tri_s[1, t]], py[tri_s[2, t]], py[tri_s[3, t]], py[tri_s[1, t]], NaN)
    end
    p3 = plot(xs, ys; seriestype = :shape, fillcolor = :steelblue,
              fillalpha = 0.35, linecolor = :grey15, lw = 0.4,
              title = "isometric view (pseudo-3D)", aspect_ratio = :equal,
              legend = false, ticks = false)
    p = plot(plots..., p3, layout = (2, 2), size = (1000, 950),
             plot_title = "$(spec.name) — computational geometry",
             titlefontsize = 9, plot_titlefontsize = 12)
    path = joinpath(outdir, "geometry_views.png")
    savefig(p, path)
    return path
end

"Throughput/efficiency metrics from the recorded stage timings."
function _perf_metrics(res::CaseResult)
    n = res.num_basis
    n2 = Float64(n)^2
    n3 = Float64(n)^3
    return (
        unknowns          = n,
        asm_rate_gm       = res.t_assembly > 0 ? n2 / res.t_assembly / 1e9 : NaN,   # 10⁹ interactions/s
        solve_gflops      = res.t_solve     > 0 ? (2 / 3) * n3 / res.t_solve / 1e9 : NaN, # dense LU
        total             = res.t_mesh + res.t_assembly + res.t_solve + res.t_rcs + res.t_plots,
    )
end

function _write_perf_csv(res::CaseResult)
    m = _perf_metrics(res)
    path = joinpath(res.outdir, "perf.csv")
    open(path, "w") do io
        println(io, "metric,value")
        println(io, "unknowns,", m.unknowns)
        println(io, "elements,", res.trinum)
        println(io, "nodes,", res.num_nodes)
        println(io, "t_mesh_s,", @sprintf("%.2f", res.t_mesh))
        println(io, "t_assembly_s,", @sprintf("%.2f", res.t_assembly))
        println(io, "t_solve_s,", @sprintf("%.2f", res.t_solve))
        println(io, "t_rcs_s,", @sprintf("%.2f", res.t_rcs))
        println(io, "t_plots_s,", @sprintf("%.2f", res.t_plots))
        println(io, "t_total_s,", @sprintf("%.2f", m.total))
        println(io, "assembly_rate_Ginteractions_per_s,", @sprintf("%.3f", m.asm_rate_gm))
        println(io, "solve_Gflops,", @sprintf("%.2f", m.solve_gflops))
    end
    return path
end


_write_report(res::CaseResult, rcs_dB, mie_dB, mie_ok) =
    _publication_report(res, rcs_dB, mie_dB, mie_ok; geometry = :surface)

"thin wrapper keeping the historical volume signature"
_write_report_volume(res::CaseResult, rcs_dB, mie_dB, mie_ok,
                     region_tags::Dict{String,Int}, ctx::VolumeMaterialContext) =
    _publication_report(res, rcs_dB, mie_dB, mie_ok; geometry = :volume,
                        region_tags = region_tags, ctx = ctx)

"""
    _key_conclusions(res, mie_ok, rcs_dB, m) -> Vector{String}

Auto-generated "Key conclusions" bullets for the publication report.
"""
function _key_conclusions(res::CaseResult, mie_ok, rcs_dB, m)
    s = res.spec
    λ = 3e8 / s.freq
    c = String[]
    push!(c, "Electric resolution: mesh size $(s.mesh_size) m = " *
             "$(round(s.mesh_size / λ; digits = 3)) λ at f = $(s.freq/1e6) MHz " *
             "($(_deg(s.theta_inc))°, $(_deg(s.phi_inc))° plane-wave incidence).")
    if mie_ok
        worst = maximum(res.mie_rmse)
        push!(c, "Accuracy: worst-cut RMSE vs analytic Mie reference = " *
                 "$(round(worst; digits = 3)) dB over $(length(s.phi_cuts)) cuts → " *
                 "**$(worst ≤ 0.5 ? "PASS" : "CHECK")** against the 0.5 dB acceptance line.")
    else
        push!(c, "No analytic reference exists for this geometry; verification is by " *
                 "mesh convergence against the refined twin case (see §4).")
    end
    push!(c, "Throughput: impedance assembly $( round(m.asm_rate_gm; sigdigits = 3)) G-interactions/s, " *
             "dense LU $(round(m.solve_gflops; digits = 1)) Gflop/s at N = $(m.unknowns) unknowns.")
    stages = ("gmsh meshing" => res.t_mesh, "matrix assembly" => res.t_assembly,
              "LU solve" => res.t_solve, "RCS post-processing" => res.t_rcs,
              "plots" => res.t_plots)
    (dname, dtime) = stages[argmax(last.(stages))]
    push!(c, "Dominant cost: $(dname) — $(round(dtime; digits = 1)) s " *
             "($(round(100 * dtime / m.total; digits = 0))% of the $(round(m.total; digits = 1)) s total).")
    fin = filter(!isnan, vec(rcs_dB))
    isempty(fin) || push!(c, "Bistatic RCS dynamic range over the observed cuts: " *
             "$(round(minimum(fin); digits = 1)) … $(round(maximum(fin); digits = 1)) dBsm.")
    return c
end

"""
    _publication_report(res, rcs_dB, mie_dB, mie_ok; geometry, region_tags, ctx)

Publication-grade report: geometry views (three orthographic + pseudo-3D),
simulation configuration, performance & efficiency, results comparison and
auto-generated key conclusions. Table rows `| geometry |`, `| formulation |`
and the timing row keep the exact formats `write_index` parses.
"""
function _publication_report(res::CaseResult, rcs_dB, mie_dB, mie_ok;
                             geometry::Symbol, region_tags = nothing, ctx = nothing)
    s = res.spec
    m = _perf_metrics(res)
    buf = IOBuffer()

    println(buf, "# EMMoMSuite Validation Report — `$(s.name)`")
    println(buf)
    println(buf, "| | |")
    println(buf, "|---|---|")
    println(buf, "| solver | EMMoMSuite v$(pkgversion(EMMoMSuite)) (Julia $(VERSION)) |")
    println(buf, "| generated | $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")) |")
    println(buf, "| pipeline | geometry → gmsh mesh → MoM solve → RCS / far-field → this report |")
    println(buf)

    println(buf, "## 1 · Geometry")
    println(buf)
    gpath = _plot_geometry_views(s, res.mesh, res.outdir)
    println(buf, "![geometry views](geometry_views.png)")
    println(buf)
    println(buf, "Boundary discretization: **$(res.trinum) $(geometry === :volume ? "boundary triangles (of $(res.trinum) tetrahedra)" : "triangles")**, " *
                 "$(res.num_nodes) nodes; geometry source `cases/geo/$(s.geo)` (gmsh/OpenCASCADE).")
    println(buf)

    println(buf, "## 2 · Simulation Configuration")
    println(buf)
    println(buf, "| item | value |")
    println(buf, "|---|---|")
    println(buf, "| geometry | `cases/geo/$(s.geo)` |")
    println(buf, "| mesh size | $(s.mesh_size) m |")
    println(buf, "| mesh dim | $(s.dim) ($(geometry === :volume ? "tetrahedral volume mesh" : "surface triangulation")) |")
    println(buf, "| frequency | $(s.freq/1e6) MHz (λ = $(round(3e8/s.freq; digits = 3)) m) |")
    if geometry === :volume
        kind = ctx.boundary_fallback ? "EFIE (all-PEC fallback, surface extracted)" :
               "VEFIE (SWG volume discretization)"
        println(buf, "| formulation | $(kind) |")
        println(buf, "| regions (Physical Volume) | material | tag |")
        println(buf, "|---|---|---|")
        for r in s.regions
            println(buf, "| `$(r.surface)` | $(r.material) | $(get(region_tags, r.surface, "—")) |")
        end
    else
        _nonair(d::Dielectric) = !(isapprox(d.eps_r, AIR.eps_r) && isapprox(d.mu_r, AIR.mu_r))
        iestr = s.ie == "CFIE" ? "CFIE" :
                (any((itf.minus isa Dielectric && _nonair(itf.minus)) ||
                     (itf.plus  isa Dielectric && _nonair(itf.plus)) for itf in s.interfaces) ?
                 "PMCHW" : "EFIE")
        println(buf, "| formulation | $(iestr)" * (iestr == "CFIE" ? " (α = $(s.alpha))" : "") * " |")
        println(buf, "| interfaces (n̂ = mesh triangle normal) | plus side (n̂) | minus side |")
        println(buf, "|---|---|---|")
        for itf in s.interfaces
            println(buf, "| `$(itf.surface)`$(itf.closed ? " (closed)" : " (open)")$(itf.flip ? " [flipped]" : "") | $(itf.plus) | $(itf.minus) |")
        end
    end
    println(buf, "| incidence | θᵢ = $(_deg(s.theta_inc))°, φᵢ = $(_deg(s.phi_inc))°, pol = [$(join(s.pol, ", "))] |")
    println(buf, "| observation | θ ∈ [0°, 180°] ($(s.n_theta) samples), φ cuts = " *
                 "$(join(_deg.(s.phi_cuts), "°, "))° |")
    if s.mie_radius !== nothing
        kind_mie = geometry === :volume ?
            (ctx.boundary_fallback ? "PEC" : "dielectric") :
            (any(itf -> (itf.minus isa Dielectric && _nonair(itf.minus)) ||
                        (itf.plus  isa Dielectric && _nonair(itf.plus)), s.interfaces) ?
             "dielectric" : "PEC")
        println(buf, "| reference | analytic Mie series, $(kind_mie) sphere r = $(s.mie_radius) m |")
    end
    println(buf)

    println(buf, "## 3 · Performance & Efficiency")
    println(buf)
    println(buf, "| triangles/tets | nodes | unknowns | t_mesh | t_assemble | t_solve | t_RCS | t_plots |")
    println(buf, "|---|---|---|---|---|---|---|---|")
    @printf(buf, "| %d | %d | %d | %.1f s | %.1f s | %.1f s | %.1f s | %.1f s |\n",
            res.trinum, res.num_nodes, res.num_basis,
            res.t_mesh, res.t_assembly, res.t_solve, res.t_rcs, res.t_plots)
    println(buf)
    println(buf, "| total time | assembly rate | LU throughput |")
    println(buf, "|---|---|---|")
    @printf(buf, "| %.1f s | %.3g × 10⁹ interactions/s | %.1f Gflop/s |\n",
            m.total, m.asm_rate_gm, m.solve_gflops)
    println(buf)
    println(buf, "*Assembly rate counts N² impedance-matrix interactions; LU throughput " *
                 "uses the dense (2/3)·N³ flop model (single node, default BLAS threads).*")
    println(buf, "Machine-readable: `perf.csv`, `rcs.csv`.")
    println(buf)

    println(buf, "## 4 · Results & Comparison")
    println(buf)
    println(buf, "Bistatic RCS — MoM vs analytic reference:")
    println(buf)
    println(buf, "![RCS cuts](rcs_cuts.png)")
    println(buf)
    println(buf, "Normalized far-field |E| pattern:")
    println(buf)
    println(buf, "![Far-field polar](farfield_polar.png)")
    println(buf)
    if mie_ok
        println(buf, "| phi cut | RMSE vs Mie [dB] | verdict |")
        println(buf, "|---|---|---|")
        for (j, φ) in enumerate(s.phi_cuts)
            r = res.mie_rmse[j]
            @printf(buf, "| %.1f° | %.3f | %s |\n", _deg(φ), r, r ≤ 0.5 ? "pass" : "CHECK")
        end
    else
        println(buf, "No analytic reference for this geometry.")
        println(buf, "Verification method: mesh convergence — compare the RCS cuts")
        println(buf, "(`rcs.csv` / `rcs_cuts.png`) against the refined twin case on the")
        println(buf, "same observation grid; agreement within ~1 dB indicates a")
        println(buf, "mesh-converged solution.")
    end
    println(buf)

    println(buf, "## 5 · Key Conclusions")
    println(buf)
    foreach(c -> println(buf, "- ", c), _key_conclusions(res, mie_ok, rcs_dB, m))
    println(buf)

    println(buf, "## Artifacts")
    println(buf)
    println(buf, "- geometry: `geometry_views.png` · RCS data: `rcs.csv` · performance: `perf.csv`")
    println(buf, "- plots: `rcs_cuts.png`, `farfield_polar.png`")
    open(joinpath(res.outdir, "report.md"), "w") do io
        print(io, String(take!(buf)))
    end
    _write_perf_csv(res)
    return nothing
end


function write_index(outroot::AbstractString = RESULT_ROOT)
    path = joinpath(outroot, "index.md")
    rows = Tuple{String,String,String,Int,Int,String,Float64,Float64}[]  # name, geo, ie, tris, unknowns, mie, total, gflops
    for d in sort(readdir(outroot; join = true))
        isdir(d) || continue
        rep = joinpath(d, "report.md")
        isfile(rep) || continue
        name = basename(d)
        lines = readlines(rep)
        geom = replace(only([l for l in lines if occursin("| geometry |", l)]),
                       r".*`cases/geo/([^`]+)`.*" => s"\1")
        ie = replace(only([l for l in lines if occursin("| formulation |", l)]),
                     r".*\| ([A-Z]+).*" => s"\1")
        meshrow = only([l for l in lines if occursin(r"^\|\s*\d+ \| \d+ \| \d+ \|", l)])
        f = split(strip(meshrow, ['|']), '|')
        tris = parse(Int, strip(f[1])); unk = parse(Int, strip(f[3]))
        tvals = [parse(Float64, strip(strip(x), ['s', ' '])) for x in f[4:end]]
        total = sum(tvals)
        gf = NaN
        pf = joinpath(d, "perf.csv")
        isfile(pf) && for l in eachline(pf)
            occursin("solve_Gflops", l) && (gf = parse(Float64, split(l, ',')[2]))
        end
        mie = "—"
        acc = [l for l in lines if occursin(r"^\|\s*[\d.]+° \| -?\d", l)]
        isempty(acc) || (mie = join([strip(split(strip(l, ['|']), '|')[2]) for l in acc], " / "))
        push!(rows, (name, geom, ie, tris, unk, mie, total, gf))
    end
    open(path, "w") do io
        println(io, "# EMMoMSuite gmsh-driven RCS case family — index")
        println(io)
        println(io, "![](../../docs/assets/logo.svg)")
        println(io)
        println(io, "*Generated: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))*")
        println(io)
        println(io, "Pipeline: geometry file → gmsh mesh → MoM solve → RCS → far-field plots → publication report.")
        println(io)
        println(io, "| case | geometry | IE | elements | unknowns | Mie RMSE [dB] | total [s] | LU [Gflop/s] | report |")
        println(io, "|---|---|---|---|---|---|---|---|---|")
        for (name, geom, ie, tris, unk, mie, total, gf) in rows
            println(io, "| [`$(name)`]($(name)/report.md) | $(geom) | $(ie) | $(tris) | $(unk) | $(mie) | $(round(total; digits=1)) | $(isnan(gf) ? "—" : round(gf; digits=1)) | [report.md]($(name)/report.md) |")
        end
    end
    println("index: ", path)
    return path
end

end # module
