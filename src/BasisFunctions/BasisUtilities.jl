using ..CoreModule
using ..Geometry
using StaticArrays
using LinearAlgebra

export count_unknowns, get_triangle_info, get_triangles_info, get_tetrahedra_info

"""
    get_triangles_info(mesh, basis)

Construct a vector of TriangleInfo for the entire mesh.
"""
function get_triangles_info(mesh::TriangleMesh{IT, FT}, basis::RWGBasis{IT, FT}) where {IT, FT}
    nt = mesh.trinum
    infos = Vector{TriangleInfo{IT, FT}}(undef, nt)
    Threads.@threads for t in 1:nt
        infos[t] = get_triangle_info(mesh, basis, t)
    end
    return infos
end

"""
    count_unknowns(basis::AbstractBasisFunction)

Return the number of unknowns (degrees of freedom) in the basis.
"""
function count_unknowns(basis::AbstractBasisFunction)
    return num_basis(basis)
end

"""
    get_triangle_info(mesh::TriangleMesh, basis::RWGBasis, t::Int)

Retrieve detailed geometric and topological information for a specific triangle.

Populates the `TriangleInfo` structure, which includes:
- Vertex coordinates and indices.
- Geometric properties (Area, Normal, Center).
- Edge lengths.
- **Basis Function IDs**: The global IDs of the RWG basis functions associated with the three edges of the triangle.

# Arguments
- `mesh`: The triangular mesh.
- `basis`: The RWG basis set.
- `t`: Index of the triangle.

# Returns
- `TriangleInfo`: A fully populated struct ready for integration routines.
"""
function get_triangle_info(mesh::TriangleMesh{IT, FT}, basis::RWGBasis{IT, FT}, t::Int) where {IT, FT}
    nodes = vertices(mesh)
    tris = elements(mesh)
    
    v_ids = tris[:, t]
    v1 = nodes[:, v_ids[1]]
    v2 = nodes[:, v_ids[2]]
    v3 = nodes[:, v_ids[3]]
    
    # Calculate edge vectors and normals (Legacy Convention: EDGEVpINTriVsID - EDGEVmINTriVsID)
    # Edge 1 (opposite v1): v3 - v2  (direction: v2 → v3)
    # Edge 2 (opposite v2): v1 - v3  (direction: v3 → v1)
    # Edge 3 (opposite v3): v2 - v1  (direction: v1 → v2)
    
    e1 = v3 - v2
    e2 = v1 - v3
    e3 = v2 - v1
    
    l1 = norm(e1)
    l2 = norm(e2)
    l3 = norm(e3)
    
    area = 0.5 * norm(cross(v2 - v1, v3 - v1))
    center = (v1 + v2 + v3) / 3
    n = normalize(cross(v2 - v1, v3 - v1))
    
    ev1 = e1 / l1
    ev2 = e2 / l2
    ev3 = e3 / l3
    
    # Edge normals: outward-pointing, in the triangle plane
    # Legacy convention: edgen̂ = cross(edgev̂, facen̂)
    # With correct edgev̂ direction (e.g., v3-v2 for edge 1), cross(edgev̂, facen̂)
    # points OUTWARD (away from the opposite vertex) for CCW-oriented triangles.
    
    en1 = cross(ev1, n)
    en2 = cross(ev2, n)
    en3 = cross(ev3, n)
    
    # Populate basis IDs and signs
    ids = MVector{3, IT}(0, 0, 0)
    signs = MVector{3, Int}(0, 0, 0)
    
    for k in 1:3
        bid = basis.basis_map[k, t]
        if bid != 0
            bf = basis.functions[bid]
            if bf.support[1] == t
                sign = bf.signs[1]
            else
                sign = bf.signs[2]
            end
            ids[k] = bid
            signs[k] = sign
        end
    end
    bfs_id = SVector{3, IT}(ids)
    bfs_sign = SVector{3, Int}(signs)

    edgel = SVector{3, FT}(l1, l2, l3)
    verts = SMatrix{3, 3, FT}(
        v1[1], v1[2], v1[3],
        v2[1], v2[2], v2[3],
        v3[1], v3[2], v3[3]
    )

    return TriangleInfo(
        IT(t),
        mesh.tags[t],
        area,
        SVector{3, IT}(v_ids...),
        verts,
        SVector{3, FT}(center),
        SVector{3, FT}(n),
        edgel,
        SMatrix{3, 3, FT, 9}(hcat(ev1, ev2, ev3)),
        SMatrix{3, 3, FT, 9}(hcat(en1, en2, en3)),
        bfs_id,
        bfs_sign
    )
end

"""
    get_tetrahedra_info(mesh, basis::PWCBasis, permittivities)

Construct a vector of TetrahedraInfo for the entire mesh using PWC basis functions.
PWC uses 3 DOFs per tetrahedron (x, y, z components), stored in inBfsID[1:3].
The 4th entry of inBfsID is set to 0 (unused).
"""
function get_tetrahedra_info(mesh::TetrahedraMesh{IT, FT}, basis::PWCBasis{IT, FT}, permittivities::Vector{ComplexF64}) where {IT, FT}
    ntet = mesh.tetnum
    infos = Vector{TetrahedraInfo{IT, FT, ComplexF64}}(undef, ntet)
    
    Threads.@threads for i in 1:ntet
        # PWC: 3 basis functions (x,y,z) per tet, 4th = 0
        bf_ids = SVector{4, IT}(3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3, 0)
        bf_signs = SVector{4, Int}(1, 1, 1, 0)
        infos[i] = TetrahedraInfo(mesh, i, bf_ids, bf_signs, permittivities[i])
    end
    
    return infos
end

"""
    get_tetrahedra_info(mesh, basis::SWGBasis, permittivities)

Construct a vector of TetrahedraInfo for the entire mesh.
"""
function get_tetrahedra_info(mesh::TetrahedraMesh{IT, FT}, basis::SWGBasis{IT, FT}, permittivities::Vector{ComplexF64}) where {IT, FT}
    ntet = mesh.tetnum
    infos = Vector{TetrahedraInfo{IT, FT, ComplexF64}}(undef, ntet)
    
    # Initialize basis IDs map: tet_idx -> [bf_id_face1, bf_id_face2, ...]
    tet_bfs = [zeros(IT, 4) for _ in 1:ntet]
    tet_signs = [zeros(Int, 4) for _ in 1:ntet]
    
    # Fill basis IDs and signs
    for bf in basis.functions
        # Support 1
        tet1 = bf.support[1]
        if tet1 > 0
            face1 = bf.local_face_idx[1]
            tet_bfs[tet1][face1] = bf.id
            tet_signs[tet1][face1] = bf.signs[1]
        end
        
        # Support 2
        tet2 = bf.support[2]
        if tet2 > 0
            face2 = bf.local_face_idx[2]
            tet_bfs[tet2][face2] = bf.id
            tet_signs[tet2][face2] = bf.signs[2]
        end
    end
    
    Threads.@threads for i in 1:ntet
        infos[i] = TetrahedraInfo(mesh, i, SVector{4, IT}(tet_bfs[i]...), SVector{4, Int}(tet_signs[i]...), permittivities[i])
    end
    
    return infos
end
