"""
    TetMeshing.jl

Phase 19.1 — Unstructured tetrahedral volume mesh generation from BRepSolid
using the Gmsh built-in (geo) kernel.

Exported functions:
  tet_mesh_gmsh   — BRepSolid → TetrahedraMesh via Gmsh (explicit mesh_size)
  tet_mesh        — BRepSolid → TetrahedraMesh (convenience keyword API)

Each output tetrahedron inherits the Gmsh entity tag of the volume it belongs
to (stored in `mesh.tags`).

Requires `Gmsh.jl` (via the lazy loader `_with_gmsh` from GmshAPI.jl).
"""

export tet_mesh_gmsh, tet_mesh

# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    tet_mesh_gmsh(solid::BRepSolid, mesh_size; FT, optimize, algorithm3d)
        → TetrahedraMesh

Generate an unstructured tetrahedral volume mesh for `solid` using the Gmsh
mesher.  The polyhedral faces of `solid` define the domain boundary.

# Arguments
- `solid`       : closed polyhedral solid (BRepSolid)
- `mesh_size`   : target element edge length (metres)
- `FT`          : floating-point type for node coordinates (default `Float64`)
- `optimize`    : run `"Netgen"` mesh optimiser after generation (default `true`)
- `algorithm3d` : Gmsh 3-D algorithm id  (1=Delaunay, 4=Frontal, 10=HXT)

# Notes
- The surface is meshed first (2-D), then extended to volume (3-D).
- Tetrahedra with `n_inverted > 0` indicate a geometry problem in `solid`.
"""
function tet_mesh_gmsh(
    solid::BRepSolid{<:AbstractFloat},
    mesh_size::Real;
    FT::Type{<:AbstractFloat} = Float64,
    optimize::Bool             = true,
    algorithm3d::Int           = 1,
)
    _ms = Float64(mesh_size)
    _ms > 0 || error("TetMeshing: mesh_size must be positive, got $_ms")
    nfaces = length(solid.faces)
    nfaces == 0 && error("TetMeshing: BRepSolid has no faces")
    isempty(solid.vertices) && error("TetMeshing: BRepSolid has no vertices")

    # Canonical edge lookup (same as SurfaceMeshing)
    edge_map = Dict{Tuple{Int,Int},Int}()
    for (i, e) in enumerate(solid.edges)
        edge_map[e] = i
    end

    return _with_gmsh() do gmsh
        gmsh.initialize(["gmsh", "-nopopup"])
        try
            gmsh.option.setNumber("General.Verbosity", 0)
            gmsh.model.add("brep_volume")

            # ── 1. Points ──────────────────────────────────────────────────
            for (i, v) in enumerate(solid.vertices)
                gmsh.model.geo.addPoint(
                    Float64(v[1]), Float64(v[2]), Float64(v[3]), _ms, i)
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
                    line_ids[k] = solid.edges[eidx][1] == v1 ? eidx : -eidx
                end
                gmsh.model.geo.addCurveLoop(line_ids, fi)
                gmsh.model.geo.addPlaneSurface([fi], fi)
            end

            # ── 4. Surface loop → Volume ─────────────────────────────────────
            # All BRepFace surfaces form the outer boundary; surface loop tag=1
            surface_ids = collect(1:nfaces)
            gmsh.model.geo.addSurfaceLoop(surface_ids, 1)
            gmsh.model.geo.addVolume([1], 1)

            gmsh.model.geo.synchronize()

            # ── 5. Meshing options ────────────────────────────────────────────
            gmsh.option.setNumber("Mesh.Algorithm",   6)   # 2-D: Frontal-Delaunay
            gmsh.option.setNumber("Mesh.Algorithm3D", algorithm3d)
            _set_mesh_size!(gmsh, _ms)

            gmsh.model.mesh.generate(3)
            optimize && gmsh.model.mesh.optimize("Gmsh")   # "Gmsh" = Laplacian smoother, always available

            return _extract_tet_mesh(gmsh, FT)
        finally
            gmsh.finalize()
        end
    end
end

"""
    tet_mesh(solid::BRepSolid; min_size, max_size, optimize, FT) → TetrahedraMesh

Convenience wrapper around `tet_mesh_gmsh` with keyword-only API.
`max_size` is used as the target element edge length.
"""
function tet_mesh(
    solid::BRepSolid{<:AbstractFloat};
    min_size::Real             = 0.01,
    max_size::Real             = 0.1,
    optimize::Bool             = true,
    FT::Type{<:AbstractFloat}  = Float64,
    algorithm3d::Int           = 1,
)
    return tet_mesh_gmsh(solid, Float64(max_size);
                         FT=FT, optimize=optimize, algorithm3d=algorithm3d)
end
