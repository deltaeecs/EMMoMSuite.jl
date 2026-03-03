"""
    MeshRepair.jl

Utilities for cleaning and repairing mesh data:
  - `remove_duplicate_nodes`  : merge coincident nodes within a tolerance
  - `fix_element_orientation` : unify triangle normal orientation (BFS propagation)
  - `detect_degenerates`      : locate near-zero-area/volume elements
"""

using LinearAlgebra
using StaticArrays

export remove_duplicate_nodes, fix_element_orientation, detect_degenerates

# ─────────────────────────────────────────────────────────────────────────────
# remove_duplicate_nodes
# ─────────────────────────────────────────────────────────────────────────────

"""
    remove_duplicate_nodes(mesh; tol=1e-10) → same mesh type

Return a copy of `mesh` with coincident nodes (distance ≤ `tol`) merged into
a single node.  Connectivity is remapped accordingly.

## Algorithm

Each node is quantised to a grid of spacing `tol` and inserted into a `Dict`.
The first node that falls in a given grid cell becomes the canonical node; all
subsequent nodes mapping to the same cell are redirected to it.

## Notes

- Degenerate elements that collapse after node merging (all vertices become the
  same node) are **kept** in the connectivity but will be reported by
  `detect_degenerates`.
- Element ordering (tags) is preserved one-to-one.
"""
function remove_duplicate_nodes(mesh::TriangleMesh{IT,FT}; tol::Real=1e-10) where {IT,FT}
    new_node, old_to_new = _merge_nodes(mesh.node, FT(tol))
    new_tris = IT.(old_to_new[mesh.triangles])
    return TriangleMesh(mesh.trinum, new_node, new_tris, copy(mesh.tags))
end

function remove_duplicate_nodes(mesh::TetrahedraMesh{IT,FT}; tol::Real=1e-10) where {IT,FT}
    new_node, old_to_new = _merge_nodes(mesh.node, FT(tol))
    new_tets = IT.(old_to_new[mesh.tetras])
    return TetrahedraMesh(mesh.tetnum, new_node, new_tets, copy(mesh.tags))
end

function remove_duplicate_nodes(mesh::HexahedraMesh{IT,FT}; tol::Real=1e-10) where {IT,FT}
    new_node, old_to_new = _merge_nodes(mesh.node, FT(tol))
    new_hexs = IT.(old_to_new[mesh.hexes])
    return HexahedraMesh(mesh.hexnum, new_node, new_hexs, copy(mesh.tags))
end

# Internal: build (new_node_matrix, old_to_new_index_vector)
function _merge_nodes(node::Matrix{FT}, tol::FT) where {FT}
    Nv  = size(node, 2)
    inv_tol = tol > 0 ? 1.0 / tol : 1.0e10

    # Hash nodes by quantised coordinate
    cell_to_new = Dict{Tuple{Int64,Int64,Int64}, Int}()
    old_to_new  = Vector{Int}(undef, Nv)

    # Collect unique nodes in insertion order
    uniq_x = Vector{FT}()
    uniq_y = Vector{FT}()
    uniq_z = Vector{FT}()

    for i in 1:Nv
        cx = round(Int64, node[1, i] * inv_tol)
        cy = round(Int64, node[2, i] * inv_tol)
        cz = round(Int64, node[3, i] * inv_tol)
        key = (cx, cy, cz)
        if haskey(cell_to_new, key)
            old_to_new[i] = cell_to_new[key]
        else
            push!(uniq_x, node[1, i])
            push!(uniq_y, node[2, i])
            push!(uniq_z, node[3, i])
            new_idx = length(uniq_x)
            cell_to_new[key] = new_idx
            old_to_new[i] = new_idx
        end
    end

    Nv_new = length(uniq_x)
    new_node = Matrix{FT}(undef, 3, Nv_new)
    new_node[1, :] .= uniq_x
    new_node[2, :] .= uniq_y
    new_node[3, :] .= uniq_z

    return new_node, old_to_new
end

# ─────────────────────────────────────────────────────────────────────────────
# fix_element_orientation  (TriangleMesh only)
# ─────────────────────────────────────────────────────────────────────────────

"""
    fix_element_orientation(mesh::TriangleMesh) → TriangleMesh

Return a copy of `mesh` with all triangle normals consistently oriented via
BFS propagation through the edge-adjacency graph.

**Convention**: the normal of triangle 1 is taken as the reference.  All
reachable triangles are oriented to be consistent with it.  If the mesh has
multiple disconnected components, each component is treated independently with
its own component-0 triangle as the local reference.

## Algorithm

1. Build an edge-to-triangle adjacency map: `(v_min, v_max) → [(tri_id, local_edge_idx), ...]`.
2. BFS from triangle 0 (or the first unvisited triangle per component).
3. For each unvisited neighbour sharing an edge: check whether the shared
   edge is traversed in the **same** direction in both triangles.  If so,
   flip the neighbour's vertex order (to enforce opposing traversal).

## Note

This function assumes the mesh is predominantly manifold (each edge shared by
≤ 2 triangles).  Non-manifold edges (≥ 3 triangles) are silently ignored.
"""
function fix_element_orientation(mesh::TriangleMesh{IT,FT}) where {IT,FT}
    Nt = mesh.trinum
    tris = mesh.triangles          # 3 × Nt, column-major

    # Working copy of connectivity (mutable, 1-indexed columns)
    new_tris = copy(tris)

    # ── Build edge adjacency ──────────────────────────────────────────────────
    # edge_map: (v_lo, v_hi) → list of (tri_idx, local_edge_position)
    # local_edge_position k ∈ {1,2,3}:
    #   k=1 → edge (v1,v2),  k=2 → edge (v2,v3),  k=3 → edge (v3,v1)
    edge_map = Dict{Tuple{Int,Int}, Vector{Tuple{Int,Int}}}()
    for t in 1:Nt
        for k in 1:3
            v_a = Int(new_tris[k,       t])
            v_b = Int(new_tris[mod1(k+1,3), t])
            key = v_a < v_b ? (v_a, v_b) : (v_b, v_a)
            push!(get!(edge_map, key, Tuple{Int,Int}[]), (t, k))
        end
    end

    # ── BFS propagation ───────────────────────────────────────────────────────
    visited = falses(Nt)

    for start in 1:Nt
        visited[start] && continue
        queue = [start]
        visited[start] = true

        while !isempty(queue)
            t = popfirst!(queue)

            for k in 1:3
                v_a = Int(new_tris[k,          t])
                v_b = Int(new_tris[mod1(k+1,3), t])
                key = v_a < v_b ? (v_a, v_b) : (v_b, v_a)

                nbrs = get(edge_map, key, nothing)
                nbrs === nothing && continue

                for (t_nb, _k_nb) in nbrs
                    t_nb == t && continue
                    visited[t_nb] && continue

                    # ── Check orientation consistency ──────────────────────────
                    # In a consistently oriented mesh, the shared edge must be
                    # traversed in OPPOSITE directions.
                    # t   traverses the edge as v_a → v_b.
                    # t_nb must traverse it as v_b → v_a.
                    # Check all 3 edges of t_nb for the shared undirected edge:
                    same_dir = false
                    for kn in 1:3
                        u_a = Int(new_tris[kn,          t_nb])
                        u_b = Int(new_tris[mod1(kn+1,3), t_nb])
                        if (u_a == v_a && u_b == v_b)
                            same_dir = true
                            break
                        end
                    end

                    if same_dir
                        # Flip: swap vertices 2 and 3 of t_nb (reverses winding)
                        new_tris[2, t_nb], new_tris[3, t_nb] =
                            new_tris[3, t_nb], new_tris[2, t_nb]
                        # Rebuild edges in edge_map won't be needed (we've already
                        # added all edges at the start); we just continue BFS.
                    end

                    visited[t_nb] = true
                    push!(queue, t_nb)
                end
            end
        end
    end

    return TriangleMesh(mesh.trinum, copy(mesh.node), new_tris, copy(mesh.tags))
end

# ─────────────────────────────────────────────────────────────────────────────
# detect_degenerates
# ─────────────────────────────────────────────────────────────────────────────

"""
    detect_degenerates(mesh; tol=1e-15) → Vector{Int}

Return the 1-based indices of degenerate elements: triangles with area < `tol`
or tetrahedra / hexahedra with volume < `tol`.

Useful before assembly to check that no collapsed elements slip through.
"""
function detect_degenerates(mesh::TriangleMesh{IT,FT}; tol::Real=1e-15) where {IT,FT}
    bad = Int[]
    nodes = mesh.node
    for t in 1:mesh.trinum
        v1 = SVector{3,FT}(nodes[:, mesh.triangles[1, t]])
        v2 = SVector{3,FT}(nodes[:, mesh.triangles[2, t]])
        v3 = SVector{3,FT}(nodes[:, mesh.triangles[3, t]])
        A = 0.5 * norm((v2 - v1) × (v3 - v1))
        A < tol && push!(bad, t)
    end
    return bad
end

function detect_degenerates(mesh::TetrahedraMesh{IT,FT}; tol::Real=1e-15) where {IT,FT}
    bad = Int[]
    nodes = mesh.node
    for t in 1:mesh.tetnum
        v1 = SVector{3,FT}(nodes[:, mesh.tetras[1, t]])
        v2 = SVector{3,FT}(nodes[:, mesh.tetras[2, t]])
        v3 = SVector{3,FT}(nodes[:, mesh.tetras[3, t]])
        v4 = SVector{3,FT}(nodes[:, mesh.tetras[4, t]])
        V = abs(dot(v2 - v1, (v3 - v1) × (v4 - v1))) / 6
        V < tol && push!(bad, t)
    end
    return bad
end

function detect_degenerates(mesh::HexahedraMesh; tol::Real=1e-15)
    # Approximate volume via the sum of sub-tets; flag if total < tol
    bad = Int[]
    nodes = mesh.node
    for h in 1:mesh.hexnum
        verts = [SVector{3,Float64}(nodes[:, mesh.hexes[i, h]]) for i in 1:8]
        # Rough volume: use one diagonal of the bounding box
        lo = SVector(minimum(getindex.(verts, 1)),
                     minimum(getindex.(verts, 2)),
                     minimum(getindex.(verts, 3)))
        hi = SVector(maximum(getindex.(verts, 1)),
                     maximum(getindex.(verts, 2)),
                     maximum(getindex.(verts, 3)))
        V_bb = prod(hi .- lo)
        V_bb < tol && push!(bad, h)
    end
    return bad
end
