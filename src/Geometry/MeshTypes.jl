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
