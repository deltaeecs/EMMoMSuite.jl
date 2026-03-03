using StaticArrays, LinearAlgebra

# ─────────────────────────────────────────────────────────────────────────────
# Phase 16.3 — CompositeMesh 与边界提取
# ─────────────────────────────────────────────────────────────────────────────

"""
    CompositeMesh{IT,FT}

Composite mesh carrying both a surface (`TriangleMesh`) and
a volume (`TetrahedraMesh`) mesh simultaneously.
Used for SCFIE workflows where both are required.
"""
struct CompositeMesh{IT,FT} <: AbstractMesh
    surface::TriangleMesh{IT,FT}
    volume::TetrahedraMesh{IT,FT}
end

CoreModule.vertices(m::CompositeMesh)  = m.surface.node
CoreModule.elements(m::CompositeMesh)  = m.surface.triangles
CoreModule.dimension(m::CompositeMesh) = 3

"""
    extract_surface(vol::TetrahedraMesh) → TriangleMesh

Extract the boundary triangular surface of a tetrahedral mesh.
A face is on the boundary if and only if it is shared by exactly one tetrahedron.
Outward normals are enforced: each triangle's normal points away from the
interior vertex of its parent tetrahedron.
"""
function extract_surface(vol::TetrahedraMesh{IT,FT}) where {IT,FT}
    tetras = vol.tetras       # 4 × N
    ntet   = vol.tetnum

    # For each tet, generate 4 oriented faces.
    # face_opp[i] = face opposite vertex i (0-based), with orientation giving outward normal.
    # For tet [a,b,c,d]:
    #   face 0 (opp a): [b,c,d]  — but we store with outward convention below
    # We store: (sorted_key, original_face_3_nodes, opposite_vertex_index)
    # sorted_key  → for deduplication (canonical order, smallest-first)
    # original    → original ordering for normal computation
    # opp_vnode   → global node index of the vertex NOT on this face

    TupleKey = NTuple{3,IT}
    # face_data: sorted_key → Vector{(face_nodes, opp_vertex)}
    face_map = Dict{TupleKey, Vector{Tuple{NTuple{3,IT}, IT}}}()

    # Each tet contributes 4 faces
    # Face indices (1-based) within [v1,v2,v3,v4]:
    #   face 1: (v2,v3,v4) opp v1
    #   face 2: (v1,v4,v3) opp v2  ← note reverse to keep consistent winding
    #   face 3: (v1,v2,v4) opp v3
    #   face 4: (v1,v3,v2) opp v4  ← note reverse
    face_patterns = (
        (2,3,4,1),   # (f1,f2,f3, opp_local)
        (1,4,3,2),
        (1,2,4,3),
        (1,3,2,4),
    )

    for t in 1:ntet
        v = (tetras[1,t], tetras[2,t], tetras[3,t], tetras[4,t])
        for (i1, i2, i3, iopp) in face_patterns
            fn = (v[i1], v[i2], v[i3])
            key = TupleKey(sort([fn[1], fn[2], fn[3]]))
            entry = (fn, v[iopp])
            if haskey(face_map, key)
                push!(face_map[key], entry)
            else
                face_map[key] = [entry]
            end
        end
    end

    # Collect boundary faces (appear exactly once)
    boundary = [(fn, opp) for entries in values(face_map)
                           if length(entries) == 1
                           for (fn, opp) in entries]

    n_surf = length(boundary)
    surf_nodes = vol.node  # share node array
    surf_elems = zeros(IT, 3, n_surf)

    for (k, (fn, opp_v)) in enumerate(boundary)
        a, b, c = fn
        va = @view surf_nodes[:, a]
        vb = @view surf_nodes[:, b]
        vc = @view surf_nodes[:, c]
        vo = @view surf_nodes[:, opp_v]

        # Outward normal check: (b-a)×(c-a) should point away from opp_v
        normal = cross(vb .- va, vc .- va)
        if dot(normal, va .- vo) < 0
            # flip
            surf_elems[1, k] = a
            surf_elems[2, k] = c
            surf_elems[3, k] = b
        else
            surf_elems[1, k] = a
            surf_elems[2, k] = b
            surf_elems[3, k] = c
        end
    end

    return TriangleMesh(n_surf, surf_nodes, surf_elems)
end
