"""
    MeshMaterialBind.jl

Phase 19.5 — Binding electromagnetic material models to volume mesh elements.

A `BoundMesh` wraps any `AbstractMesh` together with a `Dict{Int, MaterialModel}`
that maps per-element tag IDs (stored in `mesh.tags`) to material models.

Exported types:
  BoundMesh          — mesh + material_map pair

Exported functions:
  bind_materials     — create a BoundMesh from a mesh + tag→model dict
  validate_bindings  — check every unique tag in mesh.tags is covered
  element_material   — retrieve material for a specific element index
"""

export BoundMesh, bind_materials, validate_bindings, element_material

# ─────────────────────────────────────────────────────────────────────────────
# BoundMesh
# ─────────────────────────────────────────────────────────────────────────────

"""
    BoundMesh{MeshT, ModelT}

A volume mesh paired with a material assignment.

# Fields
- `mesh`         : the underlying `AbstractMesh` (e.g. `TetrahedraMesh`)
- `material_map` : `Dict{Int, ModelT}` mapping `mesh.tags[k]` → material

# Construction
Use `bind_materials(mesh, bindings)` rather than constructing directly.
"""
struct BoundMesh{MeshT<:AbstractMesh, ModelT}
    mesh         :: MeshT
    material_map :: Dict{Int,ModelT}
end

# ─────────────────────────────────────────────────────────────────────────────
# bind_materials
# ─────────────────────────────────────────────────────────────────────────────

"""
    bind_materials(mesh, bindings::Dict{Int, <:Any}) → BoundMesh

Associate each tag in `mesh.tags` with a material model from `bindings`.

`bindings` maps integer tag IDs to material objects (any type).  Tags present
in the mesh but absent from `bindings` are silently allowed; use
`validate_bindings` to detect and report such gaps.

# Example
```julia
mesh     = generate_box_tet_mesh(1.0, 1.0, 1.0, 4, 4, 4)
fr4      = Isotropic(4.4 * (1 - 0.02im), 1.0)
bm       = bind_materials(mesh, Dict(0 => fr4))
validate_bindings(bm)  # true
```
"""
function bind_materials(
    mesh::AbstractMesh,
    bindings::Dict{<:Integer, T},
) where {T}
    material_map = Dict{Int, T}(Int(k) => v for (k, v) in bindings)
    return BoundMesh{typeof(mesh), T}(mesh, material_map)
end

# ─────────────────────────────────────────────────────────────────────────────
# validate_bindings
# ─────────────────────────────────────────────────────────────────────────────

"""
    validate_bindings(bm::BoundMesh) → Bool

Return `true` if every unique tag in `bm.mesh.tags` has a corresponding entry
in `bm.material_map`.  Prints a warning listing any unbound tags.

!!! note
    A return value of `false` means at least one element has no material
    assigned, which would cause undefined behaviour in a solver.
"""
function validate_bindings(bm::BoundMesh)
    tags_needed = unique(bm.mesh.tags)
    missing_tags = filter(t -> !haskey(bm.material_map, t), tags_needed)
    if isempty(missing_tags)
        return true
    else
        @warn "validate_bindings: $(length(missing_tags)) unbound tag(s): $missing_tags"
        return false
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# element_material
# ─────────────────────────────────────────────────────────────────────────────

"""
    element_material(bm::BoundMesh, elem_idx::Int)

Return the material assigned to element `elem_idx` (1-based) by looking up
`bm.mesh.tags[elem_idx]` in `bm.material_map`.

Throws `KeyError` if the element's tag has no binding.
"""
function element_material(bm::BoundMesh, elem_idx::Int)
    tag = bm.mesh.tags[elem_idx]
    haskey(bm.material_map, tag) ||
        throw(KeyError("No material bound for tag $tag (element $elem_idx)"))
    return bm.material_map[tag]
end
