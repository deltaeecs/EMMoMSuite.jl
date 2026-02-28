using StaticArrays

"""
    TriangleMesh{IT, FT} <: AbstractMesh

Triangle mesh data structure.
"""
struct TriangleMesh{IT, FT} <: AbstractMesh
    trinum      ::Int
    node        ::Matrix{FT}      # 3 x Nv
    triangles   ::Matrix{IT}      # 3 x Nt
    tags        ::Vector{Int}     # Nt, Region/Property IDs
end

# Constructor with default tags (0)
function TriangleMesh(trinum::Int, node::Matrix{FT}, triangles::Matrix{IT}) where {IT, FT}
    tags = zeros(Int, trinum)
    return TriangleMesh(trinum, node, triangles, tags)
end

# Implement AbstractMesh interface
CoreModule.vertices(mesh::TriangleMesh) = mesh.node
CoreModule.elements(mesh::TriangleMesh) = mesh.triangles
CoreModule.dimension(mesh::TriangleMesh) = 2

"""
    TetrahedraMesh{IT, FT} <: AbstractMesh

Tetrahedral mesh data structure.
"""
struct TetrahedraMesh{IT, FT} <: AbstractMesh
    tetnum      ::Int
    node        ::Matrix{FT}      # 3 x Nv
    tetras      ::Matrix{IT}      # 4 x Nt
    tags        ::Vector{Int}     # Nt, Region/Property IDs
end

# Constructor with default tags (0)
function TetrahedraMesh(tetnum::Int, node::Matrix{FT}, tetras::Matrix{IT}) where {IT, FT}
    tags = zeros(Int, tetnum)
    return TetrahedraMesh(tetnum, node, tetras, tags)
end

CoreModule.vertices(mesh::TetrahedraMesh) = mesh.node
CoreModule.elements(mesh::TetrahedraMesh) = mesh.tetras
CoreModule.dimension(mesh::TetrahedraMesh) = 3

"""
    HexahedraMesh{IT, FT} <: AbstractMesh

Hexahedral mesh data structure.
"""
struct HexahedraMesh{IT, FT} <: AbstractMesh
    hexnum      ::Int
    node        ::Matrix{FT}      # 3 x Nv
    hexes       ::Matrix{IT}      # 8 x Nt
    tags        ::Vector{Int}     # Nt, Region/Property IDs
end

# Constructor with default tags (0)
function HexahedraMesh(hexnum::Int, node::Matrix{FT}, hexes::Matrix{IT}) where {IT, FT}
    tags = zeros(Int, hexnum)
    return HexahedraMesh(hexnum, node, hexes, tags)
end

CoreModule.vertices(mesh::HexahedraMesh) = mesh.node
CoreModule.elements(mesh::HexahedraMesh) = mesh.hexes
CoreModule.dimension(mesh::HexahedraMesh) = 3

"""
    TriangleInfo{IT<: Integer, FT<:AbstractFloat}

Detailed geometric information for a single triangle element.

Used extensively in integral equation assembly to avoid recomputing geometric properties.

# Fields
- `triID`: Unique identifier of the triangle.
- `tag`: Region or property tag (e.g., for dielectric regions).
- `area`: Area of the triangle (\$A\$).
- `verticesID`: Global indices of the three vertices.
- `vertices`: Coordinates of the three vertices (columns).
- `center`: Centroid of the triangle (\$(\\mathbf{v}_1 + \\mathbf{v}_2 + \\mathbf{v}_3)/3\$).
- `facen̂`: Unit normal vector to the triangle surface (\$\\hat{n}\$).
- `edgel`: Lengths of the three edges.
- `edgev̂`: Unit vectors along the three edges.
- `edgen̂`: Unit vectors normal to the edges, lying in the triangle plane, pointing outwards.
- `inBfsID`: IDs of the basis functions associated with the three edges (if any).
"""
struct TriangleInfo{IT<: Integer, FT<:AbstractFloat}
    triID       ::IT
    tag         ::Int             # Region/Property ID
    area        ::FT
    verticesID  ::SVector{3, IT}
    vertices    ::SMatrix{3, 3, FT, 9}
    center      ::SVector{3, FT}
    facen̂       ::SVector{3, FT}
    edgel       ::SVector{3, FT}
    edgev̂       ::SMatrix{3, 3, FT, 9}
    edgen̂       ::SMatrix{3, 3, FT, 9}
    inBfsID     ::SVector{3, IT}
    bfsSign     ::SVector{3, Int}
end

function TriangleInfo{IT, FT}(triID::IT = zero(IT)) where {IT <: Integer, FT<:AbstractFloat}
    tag         =    0
    area        =    zero(FT)
    verticesID  =    zero(SVector{3, IT})
    vertices    =    zero(SMatrix{3, 3, FT, 9})
    center      =    zero(SVector{3, FT})
    facen̂       =    zero(SVector{3, FT})
    edgel       =    zero(SVector{3, FT})
    edgev̂       =    zero(SMatrix{3, 3, FT, 9})
    edgen̂       =    zero(SMatrix{3, 3, FT, 9})
    inBfsID     =    zero(SVector{3, IT})
    bfsSign     =    zero(SVector{3, Int})
    return TriangleInfo{IT, FT}(triID,  tag,  area,  verticesID, vertices,
                                center, facen̂, edgel, edgev̂, 
                                edgen̂, inBfsID, bfsSign)
end

"""
    TetrahedraInfo{IT<: Integer, FT<:AbstractFloat, CT<:Complex}

Detailed geometric information for a single tetrahedron element.

# Fields
- `tetID`: Unique identifier.
- `tag`: Region tag.
- `volume`: Volume of the tetrahedron.
- `vertices`: Coordinates of the 4 vertices (3x4).
- `center`: Centroid.
- `facesArea`: Areas of the 4 faces.
- `facesn̂`: Unit normal vectors of the 4 faces (pointing outwards).
- `inBfsID`: IDs of the basis functions associated with the 4 faces.
- `κ`: Contrast or other material property.
- `ε`: Permittivity.
"""
struct TetrahedraInfo{IT<: Integer, FT<:AbstractFloat, CT<:Complex}
    tetID       ::IT
    tag         ::Int
    volume      ::FT
    vertices    ::SMatrix{3, 4, FT, 12}
    center      ::SVector{3, FT}
    facesArea   ::SVector{4, FT}
    facesn̂      ::SMatrix{3, 4, FT, 12}
    inBfsID     ::SVector{4, IT}
    bfsSign     ::SVector{4, Int}
    κ           ::CT
    ε           ::CT
end

function TetrahedraInfo{IT, FT, CT}(tetID::IT = zero(IT)) where {IT <: Integer, FT<:AbstractFloat, CT<:Complex}
    return TetrahedraInfo{IT, FT, CT}(
        tetID, 0, zero(FT),
        zero(SMatrix{3, 4, FT, 12}),
        zero(SVector{3, FT}),
        zero(SVector{4, FT}),
        zero(SMatrix{3, 4, FT, 12}),
        zero(SVector{4, IT}),
        zero(SVector{4, Int}),
        zero(CT), zero(CT)
    )
end

"""
    Quads4Hexa{IT, FT}

Quadrilateral face information for a hexahedron.
Stores geometric data needed for surface integral corrections (Fss terms in RBF).

# Fields
- `isbd`: Whether this face is on the boundary (half-basis function).
- `δκ`: Contrast difference across this face.
- `vertices`: Coordinates of the 4 corner points (3×4).
- `edgel`: Lengths of the 4 edges.
- `edgev̂`: Unit vectors along the 4 edges.
- `edgen̂`: Unit outward normal vectors for the 4 edges (in the face plane).
"""
mutable struct Quads4Hexa{FT<:AbstractFloat}
    isbd        ::Bool
    δκ          ::Complex{FT}
    vertices    ::SMatrix{3, 4, FT, 12}
    edgel       ::SVector{4, FT}
    edgev̂       ::SMatrix{3, 4, FT, 12}
    edgen̂       ::SMatrix{3, 4, FT, 12}
end

function Quads4Hexa{FT}() where {FT<:AbstractFloat}
    Quads4Hexa{FT}(
        true, zero(Complex{FT}),
        zero(SMatrix{3, 4, FT, 12}),
        zero(SVector{4, FT}),
        zero(SMatrix{3, 4, FT, 12}),
        zero(SMatrix{3, 4, FT, 12})
    )
end

"""
    area(q::Quads4Hexa)

Compute area of a quadrilateral face (sum of two triangles).
"""
function area(q::Quads4Hexa)
    v1 = q.vertices[:, 1]
    v2 = q.vertices[:, 2]
    v3 = q.vertices[:, 3]
    v4 = q.vertices[:, 4]
    a1 = 0.5 * norm(cross(v2 - v1, v3 - v1))
    a2 = 0.5 * norm(cross(v3 - v1, v4 - v1))
    return a1 + a2
end

"""
    HexahedraInfo{IT, FT, CT}

Detailed geometric information for a single hexahedron element.
Mirrors Legacy `HexahedraInfo` for use in VEFIE/SCFIE assembly.

# Fields
- `hexaID`: Unique identifier.
- `tag`: Region tag.
- `volume`: Volume (computed from 5 sub-tetrahedra).
- `ε`: Relative permittivity (complex).
- `κ`: Contrast \$\\kappa = (\\varepsilon - 1) / \\varepsilon\$.
- `center`: Centroid.
- `verticesID`: Global indices of 8 vertices.
- `vertices`: Coordinates of 8 vertices (3×8).
- `facesn̂`: Unit outward normal vectors for 6 faces (3×6).
- `facesArea`: **Signed** areas of 6 faces (positive = '+' side, negative = '−' side).
- `faces`: Detailed quadrilateral face information (for Fss surface corrections).
- `inBfsID`: IDs of basis functions (6 for RBF, 3 for PWC with padding).
"""
mutable struct HexahedraInfo{IT<:Integer, FT<:AbstractFloat, CT<:Complex}
    hexaID      ::IT
    tag         ::Int
    volume      ::FT
    ε           ::CT
    κ           ::CT
    center      ::SVector{3, FT}
    verticesID  ::SVector{8, IT}
    vertices    ::SMatrix{3, 8, FT, 24}
    facesn̂      ::MMatrix{3, 6, FT, 18}
    facesArea   ::MVector{6, FT}
    faces       ::Vector{Quads4Hexa{FT}}
    inBfsID     ::Vector{IT}
end

function HexahedraInfo{IT, FT, CT}(hexaID::IT = zero(IT)) where {IT<:Integer, FT<:AbstractFloat, CT<:Complex}
    HexahedraInfo{IT, FT, CT}(
        hexaID, 0, zero(FT),
        one(CT), zero(CT),
        zero(SVector{3, FT}),
        zero(SVector{8, IT}),
        zero(SMatrix{3, 8, FT, 24}),
        zero(MMatrix{3, 6, FT, 18}),
        zero(MVector{6, FT}),
        [Quads4Hexa{FT}() for _ in 1:6],
        zeros(IT, 6)
    )
end

"""
    hex_volume(vertices)

Compute volume of a hexahedron by decomposing into 5 tetrahedra.
`vertices` should be indexable with vertices[:, i] giving 3D coordinates.
"""
function hex_volume(v1, v2, v3, v4, v5, v6, v7, v8)
    # Decompose into 5 tetrahedra (same as Legacy)
    vol  = tet_volume(v1, v2, v3, v6)
    vol += tet_volume(v1, v3, v4, v8)
    vol += tet_volume(v1, v5, v6, v8)
    vol += tet_volume(v1, v3, v6, v8)
    vol += tet_volume(v3, v6, v7, v8)
    return vol
end

"""
    tet_volume(v1, v2, v3, v4)

Signed volume of a tetrahedron from 4 vertices.
Note: Returns signed volume (no abs) to match Legacy convention.
"""
function tet_volume(v1, v2, v3, v4)
    e1 = v2 - v1
    e2 = v3 - v1
    e3 = v4 - v1
    return dot(e1, cross(e2, e3)) / 6
end

# Face ordering for RBF: u=1, u=0, v=1, v=0, w=1, w=0
# Each face is defined by 4 vertex indices (1-based into the 8 hex vertices)
const HEXA_FACE_VERTEX_IDS = SMatrix{4, 6, Int}(
    2, 3, 7, 6,   # face 1: u=1
    1, 4, 8, 5,   # face 2: u=0
    4, 3, 7, 8,   # face 3: v=1
    1, 2, 6, 5,   # face 4: v=0
    5, 6, 7, 8,   # face 5: w=1
    1, 2, 3, 4    # face 6: w=0
)

# Opposite face mapping: face i's opposite is HEXA_OPP_FACE[i]
const HEXA_OPP_FACE = SVector{6, Int}(2, 1, 4, 3, 6, 5)

# Opposite face vertex IDs (for getFreeVns)
const HEXA_OPP_FACE_VERTEX_IDS = SMatrix{4, 6, Int}(
    HEXA_FACE_VERTEX_IDS[:, HEXA_OPP_FACE[1]]...,
    HEXA_FACE_VERTEX_IDS[:, HEXA_OPP_FACE[2]]...,
    HEXA_FACE_VERTEX_IDS[:, HEXA_OPP_FACE[3]]...,
    HEXA_FACE_VERTEX_IDS[:, HEXA_OPP_FACE[4]]...,
    HEXA_FACE_VERTEX_IDS[:, HEXA_OPP_FACE[5]]...,
    HEXA_FACE_VERTEX_IDS[:, HEXA_OPP_FACE[6]]...
)

"""
    get_free_vns(hexa::HexahedraInfo, face_idx, quad_gq_coord)

Compute free-end coordinates for RBF on the opposite face.
`quad_gq_coord` is the 4×N quadrature coordinate matrix (shape function values).
Returns a 3×N matrix of free-end coordinates.

The free-end r₀ of RBF on face `face_idx` lies on the **opposite** face.
"""
function get_free_vns(hexa::HexahedraInfo, face_idx::Integer, quad_gq_coord::AbstractMatrix)
    opp_vids = HEXA_OPP_FACE_VERTEX_IDS[:, face_idx]
    # Opposite face vertices: 3×4
    opp_verts = hexa.vertices[:, opp_vids]
    # Interpolate: (3×4) * (4×N) = 3×N
    return opp_verts * quad_gq_coord
end

"""
    gq3d_to_face2d_idx(gq3d_id::NTuple{3,Int}, face_idx::Integer, n1d::Integer)

Map a 3D Gauss quadrature index (i,j,k) to a 2D face quadrature linear index.

For RBF on face `face_idx`, the free-end lies on the opposite face and is parameterized
by 2 of the 3 parametric coordinates (u,v,w):
- face 1,2 (u-faces): use (j, k) → linear = j + (k-1)*n1d
- face 3,4 (v-faces): use (i, k) → linear = i + (k-1)*n1d
- face 5,6 (w-faces): use (i, j) → linear = i + (j-1)*n1d
"""
function gq3d_to_face2d_idx(gq3d_id::NTuple{3,Int}, face_idx::Integer, n1d::Integer)
    if face_idx == 1 || face_idx == 2
        return gq3d_id[2] + (gq3d_id[3] - 1) * n1d
    elseif face_idx == 3 || face_idx == 4
        return gq3d_id[1] + (gq3d_id[3] - 1) * n1d
    else  # face_idx == 5 || face_idx == 6
        return gq3d_id[1] + (gq3d_id[2] - 1) * n1d
    end
end

"""
    construct_gq3d_index_map(n1d::Integer)

Build a vector mapping linear hex GQ index (1:n1d³) to 3D tuple (i,j,k).
Julia column-major: index = i + (j-1)*n1d + (k-1)*n1d²
"""
function construct_gq3d_index_map(n1d::Integer)
    return vec([(i, j, k) for i in 1:n1d, j in 1:n1d, k in 1:n1d])
end

"""
    set_delta_kappa!(hexas_info::Vector{<:HexahedraInfo})

Set δκ for each face of each hexahedron.
δκ accumulates: +κ for '+' face (facesArea > 0), −κ for '−' face (facesArea < 0).
Shared internal face with uniform medium: δκ = κ − κ = 0.
Boundary face: δκ = ±κ (only one side contributes).
"""
function set_delta_kappa!(hexas_info::Vector{<:HexahedraInfo})
    # Reset
    for hexa in hexas_info
        for f in 1:6
            hexa.faces[f] = Quads4Hexa(
                hexa.faces[f].isbd, zero(hexa.faces[f].δκ),
                hexa.faces[f].vertices, hexa.faces[f].edgel,
                hexa.faces[f].edgev̂, hexa.faces[f].edgen̂
            )
        end
    end
    # Accumulate
    for hexa in hexas_info
        κ = hexa.κ
        for f in 1:6
            temp = hexa.facesArea[f] > 0 ? κ : -κ
            old_face = hexa.faces[f]
            hexa.faces[f] = Quads4Hexa(
                old_face.isbd, old_face.δκ + temp,
                old_face.vertices, old_face.edgel,
                old_face.edgev̂, old_face.edgen̂
            )
        end
    end
end

