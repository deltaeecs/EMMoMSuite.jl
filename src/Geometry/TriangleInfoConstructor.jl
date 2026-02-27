using StaticArrays
using LinearAlgebra

"""
    TriangleInfo(mesh, tri_idx, [inBfsID])

Construct TriangleInfo from a mesh and triangle index.
"""
function TriangleInfo(mesh::TriangleMesh{IT, FT}, tri_idx::Int, inBfsID::SVector{3, IT} = SVector{3, IT}(0, 0, 0)) where {IT, FT}
    # Vertices indices
    v_indices = mesh.triangles[:, tri_idx]
    
    # Vertices coordinates
    r1 = mesh.node[:, v_indices[1]]
    r2 = mesh.node[:, v_indices[2]]
    r3 = mesh.node[:, v_indices[3]]
    
    vertices = hcat(r1, r2, r3)
    
    # Edges (Aligned with RWG: Edge n opposite to Vertex n)
    # Edge 1: v2 -> v3 (Opposite v1)
    # Edge 2: v3 -> v1 (Opposite v2)
    # Edge 3: v1 -> v2 (Opposite v3)
    e1 = r3 - r2
    e2 = r1 - r3
    e3 = r2 - r1
    
    l1 = norm(e1)
    l2 = norm(e2)
    l3 = norm(e3)
    
    # Area and Normal
    # Use r2-r1 and r3-r1 for area calculation to be safe
    v12 = r2 - r1
    v13 = r3 - r1
    cross_prod = cross(v12, v13)
    area = 0.5 * norm(cross_prod)
    normal = normalize(cross_prod)
    
    center = (r1 + r2 + r3) / 3
    
    # Edge unit vectors
    ev1 = e1 / l1
    ev2 = e2 / l2
    ev3 = e3 / l3
    
    # Edge normals (in plane, outward)
    # n_edge = ev x normal
    en1 = cross(ev1, normal)
    en2 = cross(ev2, normal)
    en3 = cross(ev3, normal)
    
    return TriangleInfo{IT, FT}(
        IT(tri_idx),
        mesh.tags[tri_idx],
        area,
        SVector{3, IT}(v_indices),
        SMatrix{3, 3, FT, 9}(vertices),
        SVector{3, FT}(center),
        SVector{3, FT}(normal),
        SVector{3, FT}(l1, l2, l3),
        SMatrix{3, 3, FT, 9}(hcat(ev1, ev2, ev3)),
        SMatrix{3, 3, FT, 9}(hcat(en1, en2, en3)),
        inBfsID,
        SVector{3, Int}(0, 0, 0)
    )
end
