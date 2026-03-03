"""
    SurfaceMeshing.jl

Phase 18.3 — Surface triangulation of BRepSolid using the Gmsh geo (built-in) kernel.

Exported functions:
  surface_mesh_gmsh   — BRepSolid → TriangleMesh via Gmsh (explicit mesh_size)
  surface_mesh        — BRepSolid → TriangleMesh via Gmsh (size-field parameters)

Both functions tag each output triangle with the 1-based index of the BRepSolid
face it belongs to, enabling Phase 18.4 label propagation.

Requires `Gmsh.jl` to be a project dependency and the lazy loader from GmshAPI.jl.
"""

export surface_mesh_gmsh, surface_mesh

# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    surface_mesh_gmsh(solid::BRepSolid, mesh_size; FT=Float64, algorithm=6)
        → TriangleMesh

Triangulate the closed polyhedral surface of `solid` using the Gmsh mesher.
Each output triangle is tagged with the 1-based index of the BRepSolid face
from which it was generated.

# Arguments
- `solid`      : closed polyhedral solid (BRepSolid)
- `mesh_size`  : target element edge length
- `FT`         : floating-point type for node coordinates
- `algorithm`  : Gmsh 2D mesh algorithm  (5=Delaunay, 6=Frontal-Delaunay)

# Notes
- Uses the Gmsh `geo` (built-in) kernel directly on the polyhedral faces.
- Faces are assumed to be planar and convex, consistent with BRepSolid semantics.
"""
function surface_mesh_gmsh(
    solid::BRepSolid{<:AbstractFloat},
    mesh_size::Real;
    FT::Type{<:AbstractFloat} = Float64,
    algorithm::Int = 6,
)
    _ms     = Float64(mesh_size)
    _ms > 0 || error("SurfaceMeshing: mesh_size must be positive, got $_ms")
    nfaces  = length(solid.faces)
    nfaces == 0 && error("SurfaceMeshing: BRepSolid has no faces")
    isempty(solid.vertices) && error("SurfaceMeshing: BRepSolid has no vertices")

    # Build an edge lookup: canonical (v_lo, v_hi) → index in solid.edges
    edge_map = Dict{Tuple{Int,Int},Int}()
    for (i, e) in enumerate(solid.edges)
        edge_map[e] = i
    end

    return _with_gmsh() do gmsh
        gmsh.initialize(["gmsh", "-nopopup"])
        try
            gmsh.option.setNumber("General.Verbosity", 0)
            gmsh.model.add("brep_surface")

            # ── 1. Points ──────────────────────────────────────────────────
            for (i, v) in enumerate(solid.vertices)
                gmsh.model.geo.addPoint(Float64(v[1]), Float64(v[2]), Float64(v[3]), _ms, i)
            end

            # ── 2. Lines (edges) ────────────────────────────────────────────
            for (i, (v1, v2)) in enumerate(solid.edges)
                gmsh.model.geo.addLine(v1, v2, i)
            end

            # ── 3. Curve loops + plane surfaces (one per BRepFace) ──────────
            for (fi, face) in enumerate(solid.faces)
                vids = face.vertex_indices
                n    = length(vids)
                line_ids = Vector{Int}(undef, n)
                for k in 1:n
                    v1  = vids[k]
                    v2  = vids[mod1(k + 1, n)]
                    key = (min(v1, v2), max(v1, v2))
                    eidx = edge_map[key]
                    # Sign: +eidx if edge is stored v1→v2, else −eidx (reversed)
                    line_ids[k] = solid.edges[eidx][1] == v1 ? eidx : -eidx
                end
                gmsh.model.geo.addCurveLoop(line_ids, fi)
                gmsh.model.geo.addPlaneSurface([fi], fi)
            end

            gmsh.model.geo.synchronize()
            gmsh.option.setNumber("Mesh.Algorithm", algorithm)
            _set_mesh_size!(gmsh, _ms)
            gmsh.model.mesh.generate(2)

            return _extract_tagged_triangle_mesh(gmsh, nfaces, FT)
        finally
            gmsh.finalize()
        end
    end
end

"""
    surface_mesh(solid::BRepSolid; min_size=0.01, max_size=0.1,
                 curvature_factor=0.2, FT=Float64, algorithm=6) → TriangleMesh

Convenience wrapper around `surface_mesh_gmsh`.  `max_size` is passed as the
target element size; `min_size` and `curvature_factor` are accepted for API
compatibility and reserved for Phase 18.4 adaptive size-field support.
"""
function surface_mesh(
    solid::BRepSolid{<:AbstractFloat};
    min_size::Real         = 0.01,
    max_size::Real         = 0.1,
    curvature_factor::Real = 0.2,
    FT::Type{<:AbstractFloat} = Float64,
    algorithm::Int         = 6,
)
    return surface_mesh_gmsh(solid, Float64(max_size); FT=FT, algorithm=algorithm)
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helper
# ─────────────────────────────────────────────────────────────────────────────

"""
    _extract_tagged_triangle_mesh(gmsh, nfaces, FT) → TriangleMesh

Collect all triangular surface elements from the active Gmsh session, one
Gmsh surface entity at a time (entity tag == BRepFace index).
Node coordinates are extracted once via `_extract_nodes_for_elements`.
The `tags` field of the returned `TriangleMesh` contains the 1-based
BRepFace index for every triangle, enabling Phase 18.4 label propagation.
"""
function _extract_tagged_triangle_mesh(gmsh, nfaces::Int, FT::Type{<:AbstractFloat})
    # Accumulate connectivity and face tags over all BRepFace surfaces.
    all_conn  = Int[]    # flat node-tag list (3 entries per triangle)
    elem_face = Int[]    # face index (1-based) for each triangle

    for fi in 1:nfaces
        etags, conn = gmsh.model.mesh.getElementsByType(2, fi)
        nt = length(etags)
        nt == 0 && continue
        append!(all_conn, conn)
        append!(elem_face, fill(fi, nt))
    end

    ntri = length(elem_face)
    ntri == 0 && error("SurfaceMeshing: no triangles generated for BRepSolid")

    node_mat, tag2idx = _extract_nodes_for_elements(gmsh, all_conn, FT)
    conn_mat = _build_connectivity(all_conn, tag2idx, 3, ntri)

    return TriangleMesh(ntri, node_mat, conn_mat, elem_face)
end
