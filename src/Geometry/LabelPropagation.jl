"""
    LabelPropagation.jl

Phase 18.4 — Label propagation for B-Rep solids and triangle meshes.

After a boolean operation the resulting BRepSolid needs its `boundary_labels`
and `material_label` metadata re-attached, and surface meshes (TriangleMesh)
need per-triangle label strings derived from those labels.

Exported functions:
  mesh_face_labels   — TriangleMesh + BRepSolid → per-triangle label strings
  label_mesh_tags    — TriangleMesh + BRepSolid → new TriangleMesh with integer IDs
  propagate_labels   — boolean-result BRepSolid ← nearest-centroid label inheritance
"""

using LinearAlgebra
using StaticArrays

export mesh_face_labels, label_mesh_tags, propagate_labels

# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    mesh_face_labels(mesh::TriangleMesh, solid::BRepSolid) → Vector{String}

Return a `Vector{String}` of length `mesh.trinum` where each entry is the
boundary-condition label of the BRepSolid face that triangle `t` belongs to.

The lookup key is `mesh.tags[t]`, which `surface_mesh_gmsh` sets to the
1-based BRepFace index (Phase 18.3).  Triangles on faces that have no entry
in `solid.boundary_labels` receive an empty string `""`.

!!! note
    Requires `length(mesh.tags) == mesh.trinum`; passing a mesh whose `tags`
    field was never populated will throw a `BoundsError`.
"""
function mesh_face_labels(
    mesh::TriangleMesh,
    solid::BRepSolid{<:AbstractFloat},
)
    labels = Vector{String}(undef, mesh.trinum)
    @inbounds for t in 1:mesh.trinum
        labels[t] = get(solid.boundary_labels, mesh.tags[t], "")
    end
    return labels
end

"""
    label_mesh_tags(mesh::TriangleMesh, solid::BRepSolid) → TriangleMesh

Return a new `TriangleMesh` whose `tags` field holds integer label IDs instead
of raw face indices.  IDs are assigned in ascending lexicographic order of the
unique label strings (1-based); triangles on unlabeled faces get tag `0`.

This is the integer-valued counterpart of `mesh_face_labels` and is suitable
for use as a boundary-condition index in solvers.
"""
function label_mesh_tags(
    mesh::TriangleMesh,
    solid::BRepSolid{<:AbstractFloat},
)
    labels   = mesh_face_labels(mesh, solid)
    uniq     = sort(unique(filter(!isempty, labels)))
    label_id = Dict{String,Int}(l => i for (i, l) in enumerate(uniq))
    new_tags = [get(label_id, l, 0) for l in labels]
    return TriangleMesh(mesh.trinum, mesh.node, mesh.triangles, new_tags)
end

"""
    propagate_labels(result::BRepSolid, sources) → BRepSolid

Reconstruct `boundary_labels` for `result` (typically output of
`intersect_solids`) by geometric proximity: for each face of `result`, find the
face in any source solid (that has a non-empty label) whose centroid is closest
and inherit its label.

Returns a new `BRepSolid` identical to `result` but with an updated
`boundary_labels` dict.  Faces with no labelled match in any source receive no
entry (equivalent to empty label).

# Arguments
- `result`  : output of a boolean operation
- `sources` : `AbstractVector` of the operand `BRepSolid` objects
"""
function propagate_labels(
    result::BRepSolid{FT},
    sources::AbstractVector{<:BRepSolid{<:AbstractFloat}},
) where {FT<:AbstractFloat}
    new_labels = Dict{Int,String}()

    # Pre-convert each source's vertices to FT once (avoid O(F_result) re-alloc).
    src_verts_list = [_convert_vertices(FT, src.vertices) for src in sources]

    for (fi, rface) in enumerate(result.faces)
        rc         = _face_centroid(result.vertices, rface)
        best_dist  = Inf
        best_label = ""

        for (k, src) in enumerate(sources)
            src_verts = src_verts_list[k]
            for (si, sface) in enumerate(src.faces)
                label = get(src.boundary_labels, si, "")
                isempty(label) && continue
                sc   = _face_centroid(src_verts, sface)
                dist = norm(rc - sc)
                if dist < best_dist
                    best_dist  = dist
                    best_label = label
                end
            end
        end

        isempty(best_label) || (new_labels[fi] = best_label)
    end

    return BRepSolid{FT}(
        result.vertices,
        result.edges,
        result.faces,
        result.material_label,
        new_labels,
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

# Compute the centroid of a BRepFace given the solid's vertex array.
function _face_centroid(
    vertices::AbstractVector{SVector{3,FT}},
    face::BRepFace,
) where {FT}
    n = length(face.vertex_indices)
    n == 0 && return zero(SVector{3,FT})
    return sum(vertices[i] for i in face.vertex_indices) / FT(n)
end

# Convert vertex array to a different float type (for mixed-type sources).
# Same-type specialisation: return verts unchanged (zero-copy).
function _convert_vertices(
    ::Type{FT},
    verts::AbstractVector{SVector{3,FT}},
) where {FT}
    return verts
end
# General fallback: allocate a converted copy.
function _convert_vertices(
    ::Type{FT},
    verts::AbstractVector,
) where {FT}
    return [SVector{3,FT}(v) for v in verts]
end
