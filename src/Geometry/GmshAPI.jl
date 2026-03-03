"""
    GmshAPI.jl

Phase 16.7 — Gmsh API integration for programmatic mesh generation.

Uses **lazy loading**: Gmsh.jl is loaded only when a GmshAPI function is first
called, keeping the EMSuite load time unaffected when Gmsh is not used.

Requires `Gmsh.jl` to be a project dependency (`Pkg.add("Gmsh")`).

Exported functions:
  generate_gmsh_sphere     — sphere surface → TriangleMesh
  generate_gmsh_box        — box volume → TetrahedraMesh
  generate_gmsh_from_file  — arbitrary .geo file → TriangleMesh or TetrahedraMesh
"""

export generate_gmsh_sphere, generate_gmsh_box, generate_gmsh_from_file

using LinearAlgebra

# ─────────────────────────────────────────────────────────────────────────────
# Lazy Gmsh loader
# ─────────────────────────────────────────────────────────────────────────────

const _GMSH_UUID = Base.UUID("705231aa-382f-11e9-3f0c-b7cb4346fdeb")
const _gmsh_ref  = Ref{Any}(nothing)

"""
    _gmsh() → gmsh

Return the `gmsh` API object (lazy-loaded from Gmsh.jl on first call).
"""
function _gmsh()
    if _gmsh_ref[] === nothing
        try
            mod = Base.require(Base.PkgId(_GMSH_UUID, "Gmsh"))
            _gmsh_ref[] = mod.gmsh
        catch e
            error("""
GmshAPI: Gmsh.jl could not be loaded. Install with: import Pkg; Pkg.add(\"Gmsh\")
Original error: $e""")
        end
    end
    return _gmsh_ref[]
end

"""
    _with_gmsh(f) → result

Call `f(gmsh)` using `Base.invokelatest` so that Gmsh methods (registered in
the latest world age after lazy-loading) are reachable from functions compiled
before Gmsh was loaded.
"""
function _with_gmsh(f::Function)
    _gmsh()  # ensure loaded
    return Base.invokelatest(f, _gmsh_ref[])
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    generate_gmsh_sphere(radius; mesh_size=0.1, FT=Float64) → TriangleMesh

Generate a triangular surface mesh of a sphere with the given radius using
the Gmsh OCC kernel.

# Arguments
- `radius`    : sphere radius
- `mesh_size` : target element size (controls mesh density)
- `FT`        : floating-point type for node coordinates (default `Float64`)
"""
function generate_gmsh_sphere(
    radius::Real;
    mesh_size::Real = 0.1,
    FT::Type{<:AbstractFloat} = Float64,
)
    return _with_gmsh() do gmsh
        gmsh.initialize(["gmsh", "-nopopup"])
        try
            gmsh.option.setNumber("General.Verbosity", 0)
            gmsh.model.add("sphere")
            gmsh.model.occ.addSphere(0.0, 0.0, 0.0, Float64(radius))
            gmsh.model.occ.synchronize()
            _set_mesh_size!(gmsh, Float64(mesh_size))
            gmsh.model.mesh.generate(2)
            return _extract_triangle_mesh(gmsh, FT)
        finally
            gmsh.finalize()
        end
    end
end

"""
    generate_gmsh_box(Lx, Ly, Lz; mesh_size=0.1, FT=Float64) → TetrahedraMesh

Generate a tetrahedral volume mesh of an axis-aligned box
[0,Lx] × [0,Ly] × [0,Lz] using the Gmsh OCC kernel.

# Arguments
- `Lx, Ly, Lz` : box dimensions
- `mesh_size`   : target element size
- `FT`          : floating-point type for node coordinates
"""
function generate_gmsh_box(
    Lx::Real, Ly::Real, Lz::Real;
    mesh_size::Real = 0.1,
    FT::Type{<:AbstractFloat} = Float64,
)
    return _with_gmsh() do gmsh
        gmsh.initialize(["gmsh", "-nopopup"])
        try
            gmsh.option.setNumber("General.Verbosity", 0)
            gmsh.model.add("box")
            gmsh.model.occ.addBox(0.0, 0.0, 0.0, Float64(Lx), Float64(Ly), Float64(Lz))
            gmsh.model.occ.synchronize()
            _set_mesh_size!(gmsh, Float64(mesh_size))
            gmsh.model.mesh.generate(3)
            return _extract_tet_mesh(gmsh, FT)
        finally
            gmsh.finalize()
        end
    end
end

"""
    generate_gmsh_from_file(geo_file; mesh_size=0.1, FT=Float64)
        → TriangleMesh or TetrahedraMesh

Load a Gmsh `.geo` geometry file, mesh it, and return the resulting mesh.
Returns a `TetrahedraMesh` if 3-D tetrahedra are present, otherwise a
`TriangleMesh` of the surface.

# Arguments
- `geo_file`  : path to a `.geo` file (Gmsh geometry script)
- `mesh_size` : global mesh size override (`0` = use file-defined sizes)
- `FT`        : floating-point type for node coordinates

# Errors
Throws `ErrorException` if the file does not exist.
"""
function generate_gmsh_from_file(
    geo_file::AbstractString;
    mesh_size::Real = 0.1,
    FT::Type{<:AbstractFloat} = Float64,
)
    isfile(geo_file) || error("File not found: $geo_file")
    _ms = Float64(mesh_size)
    return _with_gmsh() do gmsh
        gmsh.initialize(["gmsh", "-nopopup"])
        try
            gmsh.option.setNumber("General.Verbosity", 0)
            gmsh.open(geo_file)          # execute the .geo script
            # Synchronize both factories (safe: no-op if factory not active)
            try gmsh.model.occ.synchronize() catch; end
            try gmsh.model.geo.synchronize() catch; end
            _ms > 0 && _set_mesh_size!(gmsh, _ms)
            dim = _highest_dim_with_entities(gmsh)
            gmsh.model.mesh.generate(dim)
            return dim == 3 ? _extract_tet_mesh(gmsh, FT) : _extract_triangle_mesh(gmsh, FT)
        finally
            gmsh.finalize()
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

# Set both max and min mesh size constraints.
# MeshSizeMin = 0.3h prevents Gmsh from over-refining at geometric singularities
# (e.g., sphere poles) while maintaining overall mesh smoothness.
function _set_mesh_size!(gmsh, h::Float64)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h * 0.3)
end

# Return the highest geometric dimension that has model entities.
function _highest_dim_with_entities(gmsh)
    for d in (3, 2, 1)
        !isempty(gmsh.model.getEntities(d)) && return d
    end
    @warn "GmshAPI: model has no entities in dim 1–3; defaulting to dim=2"
    return 2
end

# ─────────────────────────────────────────────────────────────────────────────
# Mesh element extraction
# ─────────────────────────────────────────────────────────────────────────────

# Build (node_matrix, tag→index map) from a flat node_conn vector.
# Queries Gmsh for all node coordinates in one bulk call.
function _extract_nodes_for_elements(gmsh, node_conn::Vector, FT::Type{<:AbstractFloat})
    unique_tags = unique(node_conn)
    tag2idx     = Dict{Int,Int}(t => i for (i, t) in enumerate(unique_tags))
    all_ntags, coord, _ = gmsh.model.mesh.getNodes()
    tagpos      = Dict{Int,Int}(all_ntags[i] => i for i in eachindex(all_ntags))
    Nv          = length(unique_tags)
    node_mat    = Matrix{FT}(undef, 3, Nv)
    for (i, t) in enumerate(unique_tags)
        p = tagpos[t]
        node_mat[1, i] = FT(coord[3*(p-1)+1])
        node_mat[2, i] = FT(coord[3*(p-1)+2])
        node_mat[3, i] = FT(coord[3*(p-1)+3])
    end
    return node_mat, tag2idx
end

# Build a (npe × n) Int32 connectivity matrix from a flat node_conn vector.
function _build_connectivity(node_conn::Vector, tag2idx::Dict{Int,Int}, npe::Int, n::Int)
    conn_mat = Matrix{Int32}(undef, npe, n)
    @inbounds for t in 1:n, k in 1:npe
        conn_mat[k, t] = tag2idx[node_conn[npe*(t-1)+k]]
    end
    return conn_mat
end

"""
    _extract_triangle_mesh(gmsh, FT) → TriangleMesh

Extract all surface triangle elements (Gmsh element type 2) from the current
Gmsh session and return a `TriangleMesh{Int32,FT}`.

Only the nodes referenced by triangle elements are included in the output mesh.
"""
function _extract_triangle_mesh(gmsh, FT::Type{<:AbstractFloat})
    elem_tags, node_conn = gmsh.model.mesh.getElementsByType(2)
    isempty(elem_tags) && error("GmshAPI: no triangle elements found in current mesh")
    ntri     = length(elem_tags)
    node_mat, tag2idx = _extract_nodes_for_elements(gmsh, node_conn, FT)
    conn_mat = _build_connectivity(node_conn, tag2idx, 3, ntri)
    return TriangleMesh(ntri, node_mat, conn_mat, ones(Int, ntri))
end

"""
    _extract_tet_mesh(gmsh, FT) → TetrahedraMesh

Extract all tetrahedral elements (Gmsh element type 4) from the current
Gmsh session and return a `TetrahedraMesh{Int32,FT}`.
"""
function _extract_tet_mesh(gmsh, FT::Type{<:AbstractFloat})
    elem_tags, node_conn = gmsh.model.mesh.getElementsByType(4)
    isempty(elem_tags) &&
        error("GmshAPI: no tetrahedral elements found; verify 3-D mesh was generated")
    ntet     = length(elem_tags)
    node_mat, tag2idx = _extract_nodes_for_elements(gmsh, node_conn, FT)
    conn_mat = _build_connectivity(node_conn, tag2idx, 4, ntet)
    return TetrahedraMesh(ntet, node_mat, conn_mat, ones(Int, ntet))
end
