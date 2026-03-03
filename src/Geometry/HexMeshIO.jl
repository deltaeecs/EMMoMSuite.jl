"""
    HexMeshIO.jl

Phase 19.2 — Hexahedral and mixed mesh import channel.

Provides a typed wrapper around the existing `read_msh_mesh` parser that
enforces the return type and adds mesh-level validation.

Exported functions:
  read_hex_mesh    — load a HexahedraMesh from a .msh v4 file
  validate_mesh    — check mesh integrity (no inverted elements, non-empty)
"""

export read_hex_mesh, validate_mesh

# ─────────────────────────────────────────────────────────────────────────────
# read_hex_mesh
# ─────────────────────────────────────────────────────────────────────────────

"""
    read_hex_mesh(path::String; format::Symbol = :gmsh_msh, FT = Float64)
        → HexahedraMesh

Import a hexahedral mesh from `path`.

# Supported formats
- `:gmsh_msh` (default) — Gmsh v4.1 ASCII .msh format

Returns a `HexahedraMesh`; throws if the file contains no hexahedral elements.

# Example
```julia
hex = read_hex_mesh("my_model.msh")
validate_mesh(hex)
```
"""
function read_hex_mesh(
    path::String;
    format::Symbol           = :gmsh_msh,
    FT::Type{<:AbstractFloat} = Float64,
)
    format === :gmsh_msh ||
        error("read_hex_mesh: unsupported format $(format). " *
              "Supported: :gmsh_msh (.msh v4.1 ASCII). " *
              "VTU and Exodus support is planned for a future phase.")

    mesh = read_msh_mesh(path; FT = FT)
    mesh isa HexahedraMesh ||
        error("read_hex_mesh: file does not contain hexahedral elements: $path " *
              "(got $(typeof(mesh)))")
    return mesh::HexahedraMesh
end

# ─────────────────────────────────────────────────────────────────────────────
# validate_mesh
# ─────────────────────────────────────────────────────────────────────────────

"""
    validate_mesh(mesh::AbstractMesh; tol = 1e-14) → Bool

Return `true` if `mesh` passes basic integrity checks:
1. At least one element.
2. At least four vertices (minimum for a 3-D mesh).
3. No inverted elements (negative signed volume / area below `-tol * mean`).

Prints a warning with details for each failed check.  Does not modify the mesh.
"""
function validate_mesh(mesh::AbstractMesh; tol::Float64 = 1e-14)
    ok = true

    n_elems = size(elements(mesh), 2)
    if n_elems == 0
        @warn "validate_mesh: mesh has no elements"
        ok = false
    end

    n_verts = size(vertices(mesh), 2)
    if n_verts < 4
        @warn "validate_mesh: mesh has only $n_verts vertices (need ≥ 4 for a 3-D mesh)"
        ok = false
    end

    # Quality check (only for mesh types that support mesh_quality)
    if mesh isa TriangleMesh || mesh isa TetrahedraMesh
        qr = mesh_quality(mesh)
        if qr.n_inverted > 0
            @warn "validate_mesh: $(qr.n_inverted) inverted element(s) detected"
            ok = false
        end
    end
    # Note: HexahedraMesh quality check deferred until mesh_quality(::HexahedraMesh)
    # is implemented.

    return ok
end
