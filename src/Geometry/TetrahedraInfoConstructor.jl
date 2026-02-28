module TetrahedraInfoConstructor

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using StaticArrays
using LinearAlgebra

export get_tetrahedra_info

"""
    get_tetrahedra_info(mesh, basis, permittivities)

Construct a vector of TetrahedraInfo for the entire mesh.
"""
function get_tetrahedra_info(
    mesh::TetrahedraMesh{IT,FT},
    basis::SWGBasis{IT,FT},
    permittivities::Vector{ComplexF64},
) where {IT,FT}
    ntet = mesh.tetnum
    infos = Vector{TetrahedraInfo{IT,FT,ComplexF64}}(undef, ntet)

    # Initialize basis IDs map: tet_idx -> [bf_id_face1, bf_id_face2, ...]
    tet_bfs = [zeros(IT, 4) for _ = 1:ntet]
    tet_signs = [zeros(Int, 4) for _ = 1:ntet]

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

    Threads.@threads for i = 1:ntet
        infos[i] = TetrahedraInfo(
            mesh,
            i,
            SVector{4,IT}(tet_bfs[i]...),
            SVector{4,Int}(tet_signs[i]...),
            permittivities[i],
        )
    end

    return infos
end

"""
    TetrahedraInfo(mesh, tet_idx, inBfsID, bfsSign, permittivity)

Construct TetrahedraInfo from a mesh and tetrahedron index.
"""
function TetrahedraInfo(
    mesh::TetrahedraMesh{IT,FT},
    tet_idx::Int,
    inBfsID::SVector{4,IT},
    bfsSign::SVector{4,Int},
    permittivity::ComplexF64,
) where {IT,FT}
    # Vertices indices
    v_indices = mesh.tetras[:, tet_idx]

    # Vertices coordinates
    r1 = mesh.node[:, v_indices[1]]
    r2 = mesh.node[:, v_indices[2]]
    r3 = mesh.node[:, v_indices[3]]
    r4 = mesh.node[:, v_indices[4]]

    vertices = hcat(r1, r2, r3, r4)

    # Volume
    # V = |(r2-r1) . ((r3-r1) x (r4-r1))| / 6
    v12 = r2 - r1
    v13 = r3 - r1
    v14 = r4 - r1
    volume = abs(dot(v12, cross(v13, v14))) / 6.0

    center = (r1 + r2 + r3 + r4) / 4.0

    # Faces
    # Face 1: v2, v3, v4 (Opposite v1)
    # Face 2: v1, v4, v3 (Opposite v2) - Note orientation for outward normal
    # Face 3: v1, v2, v4 (Opposite v3)
    # Face 4: v1, v3, v2 (Opposite v4)

    # Normals must point OUTWARD.
    # Center of face - Center of tet should be in direction of normal? No.
    # (Center of face - Center of tet) . normal > 0

    # Let's compute raw cross products
    # F1: (r3-r2) x (r4-r2)
    n1_raw = cross(r3 - r2, r4 - r2)
    # Check direction: (r1 - r2) . n1_raw should be < 0 (since r1 is inside relative to face 1)
    if dot(r1 - r2, n1_raw) > 0
        n1_raw = -n1_raw
    end

    # F2: (r4-r1) x (r3-r1)
    n2_raw = cross(r4 - r1, r3 - r1)
    if dot(r2 - r1, n2_raw) > 0
        n2_raw = -n2_raw
    end

    # F3: (r2-r1) x (r4-r1)
    n3_raw = cross(r2 - r1, r4 - r1)
    if dot(r3 - r1, n3_raw) > 0
        n3_raw = -n3_raw
    end

    # F4: (r3-r1) x (r2-r1)
    n4_raw = cross(r3 - r1, r2 - r1)
    if dot(r4 - r1, n4_raw) > 0
        n4_raw = -n4_raw
    end

    area1 = 0.5 * norm(n1_raw)
    area2 = 0.5 * norm(n2_raw)
    area3 = 0.5 * norm(n3_raw)
    area4 = 0.5 * norm(n4_raw)

    n1 = normalize(n1_raw)
    n2 = normalize(n2_raw)
    n3 = normalize(n3_raw)
    n4 = normalize(n4_raw)

    facesArea = SVector(area1, area2, area3, area4)
    facesn̂ = hcat(n1, n2, n3, n4)

    # Material properties
    ε = permittivity
    # Contrast κ = (ε - 1) / ε
    κ = (ε - 1.0) / ε

    return TetrahedraInfo{IT,FT,ComplexF64}(
        tet_idx,
        mesh.tags[tet_idx],
        volume,
        vertices,
        center,
        facesArea,
        facesn̂,
        inBfsID,
        bfsSign,
        κ,
        ε,
    )
end

end
