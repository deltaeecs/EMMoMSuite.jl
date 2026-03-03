"""
    MeshTransforms.jl

Affine transformations (translate, scale, rotate, general affine) and mesh
merging operations for all AbstractMesh subtypes.

Exported:
  translate_mesh, scale_mesh, rotate_mesh, transform_mesh, merge_meshes
"""

# ──────────────────────────────────────────────────────────────────────────────
# Internal: reconstruct a mesh of the same type with a new node matrix
# ──────────────────────────────────────────────────────────────────────────────

function _rebuild(m::TriangleMesh{IT,FT}, new_node::Matrix) where {IT,FT}
    TriangleMesh(m.trinum, convert(Matrix{eltype(new_node)}, new_node),
                 copy(m.triangles), copy(m.tags))
end

function _rebuild(m::TetrahedraMesh{IT,FT}, new_node::Matrix) where {IT,FT}
    TetrahedraMesh(m.tetnum, convert(Matrix{eltype(new_node)}, new_node),
                   copy(m.tetras), copy(m.tags))
end

function _rebuild(m::HexahedraMesh{IT,FT}, new_node::Matrix) where {IT,FT}
    HexahedraMesh(m.hexnum, convert(Matrix{eltype(new_node)}, new_node),
                  copy(m.hexes), copy(m.tags))
end

# ──────────────────────────────────────────────────────────────────────────────
# translate_mesh(mesh, displacement)
# ──────────────────────────────────────────────────────────────────────────────

"""
    translate_mesh(mesh, displacement)

Return a copy of `mesh` with all nodes shifted by `displacement` (length-3 vector).
Connectivity and tags are preserved.
"""
function translate_mesh(m::AbstractMesh, d)
    d3 = reshape(Float64[d[1], d[2], d[3]], 3, 1)
    return _rebuild(m, m.node .+ d3)
end

# ──────────────────────────────────────────────────────────────────────────────
# scale_mesh(mesh, factor)  /  scale_mesh(mesh, sx, sy, sz)
# ──────────────────────────────────────────────────────────────────────────────

"""
    scale_mesh(mesh, factor::Real)

Isotropic scaling: multiply all node coordinates by `factor`.
"""
function scale_mesh(m::AbstractMesh, factor::Real)
    return _rebuild(m, Float64(factor) .* m.node)
end

"""
    scale_mesh(mesh, sx, sy, sz)

Anisotropic scaling: multiply x/y/z coordinates by `sx`, `sy`, `sz` respectively.
"""
function scale_mesh(m::AbstractMesh, sx::Real, sy::Real, sz::Real)
    new_node = similar(m.node, Float64)
    @. new_node[1,:] = sx * m.node[1,:]
    @. new_node[2,:] = sy * m.node[2,:]
    @. new_node[3,:] = sz * m.node[3,:]
    return _rebuild(m, new_node)
end

# ──────────────────────────────────────────────────────────────────────────────
# rotate_mesh(mesh, axis, angle)
# ──────────────────────────────────────────────────────────────────────────────

"""
    _rodrigues(axis, angle) → Matrix{Float64}

Rodrigues rotation formula: 3×3 rotation matrix for rotation of `angle` radians
about the unit vector `axis`.

``R = \\cos(θ) I + (1-\\cos(θ))\\, \\hat{n}\\hat{n}^\\top + \\sin(θ)[\\hat{n}]_×``
"""
function _rodrigues(axis, angle::Real)
    nrm = sqrt(sum(x^2 for x in axis))
    nx, ny, nz = axis[1]/nrm, axis[2]/nrm, axis[3]/nrm
    c = cos(angle)
    s = sin(angle)
    t = 1.0 - c
    return Float64[
        t*nx*nx + c     t*nx*ny - s*nz  t*nx*nz + s*ny ;
        t*nx*ny + s*nz  t*ny*ny + c     t*ny*nz - s*nx ;
        t*nx*nz - s*ny  t*ny*nz + s*nx  t*nz*nz + c
    ]
end

"""
    rotate_mesh(mesh, axis, angle::Real)

Rotate `mesh` by `angle` radians about `axis` (length-3 vector) using the
Rodrigues rotation formula.  The rotation is applied about the origin.
"""
function rotate_mesh(m::AbstractMesh, axis, angle::Real)
    R = _rodrigues(axis, angle)
    return _rebuild(m, R * m.node)
end

# ──────────────────────────────────────────────────────────────────────────────
# transform_mesh(mesh, R, t)  — general affine
# ──────────────────────────────────────────────────────────────────────────────

"""
    transform_mesh(mesh, R, t)

Apply the affine transformation ``x \\mapsto R x + t`` to every node.
`R` must be a 3×3 matrix and `t` a length-3 vector.
"""
function transform_mesh(m::AbstractMesh, R, t)
    t3 = reshape(Float64[t[1], t[2], t[3]], 3, 1)
    return _rebuild(m, R * m.node .+ t3)
end

# ──────────────────────────────────────────────────────────────────────────────
# merge_meshes(meshes::Vector{<:TriangleMesh})   → TriangleMesh
# merge_meshes(meshes::Vector{<:TetrahedraMesh}) → TetrahedraMesh
# merge_meshes(meshes::Vector{<:HexahedraMesh})  → HexahedraMesh
# ──────────────────────────────────────────────────────────────────────────────

"""
    merge_meshes(meshes)

Concatenate a collection of meshes of the **same concrete type** into a single
mesh.  Node arrays are concatenated horizontally; element connectivity indices
are offset so each mesh's elements reference the correct block in the merged
node array.  Tags are also concatenated.
"""
function merge_meshes(meshes::AbstractVector{<:TriangleMesh{IT,FT}}) where {IT,FT}
    total_tri  = sum(m.trinum for m in meshes)
    total_node = sum(size(m.node, 2) for m in meshes)

    out_node  = Matrix{FT}(undef, 3, total_node)
    out_tri   = Matrix{IT}(undef, 3, total_tri)
    out_tags  = Vector{Int}(undef, total_tri)

    noff = 0   # node offset
    eoff = 0   # element offset
    for m in meshes
        nv = size(m.node, 2)
        out_node[:, noff+1:noff+nv] .= m.node
        @. out_tri[:, eoff+1:eoff+m.trinum] = m.triangles + IT(noff)
        out_tags[eoff+1:eoff+m.trinum] .= m.tags
        noff += nv
        eoff += m.trinum
    end

    return TriangleMesh(total_tri, out_node, out_tri, out_tags)
end

function merge_meshes(meshes::AbstractVector{<:TetrahedraMesh{IT,FT}}) where {IT,FT}
    total_tet  = sum(m.tetnum for m in meshes)
    total_node = sum(size(m.node, 2) for m in meshes)

    out_node  = Matrix{FT}(undef, 3, total_node)
    out_tet   = Matrix{IT}(undef, 4, total_tet)
    out_tags  = Vector{Int}(undef, total_tet)

    noff = 0
    eoff = 0
    for m in meshes
        nv = size(m.node, 2)
        out_node[:, noff+1:noff+nv] .= m.node
        @. out_tet[:, eoff+1:eoff+m.tetnum] = m.tetras + IT(noff)
        out_tags[eoff+1:eoff+m.tetnum] .= m.tags
        noff += nv
        eoff += m.tetnum
    end

    return TetrahedraMesh(total_tet, out_node, out_tet, out_tags)
end

function merge_meshes(meshes::AbstractVector{<:HexahedraMesh{IT,FT}}) where {IT,FT}
    total_hex  = sum(m.hexnum for m in meshes)
    total_node = sum(size(m.node, 2) for m in meshes)

    out_node  = Matrix{FT}(undef, 3, total_node)
    out_hex   = Matrix{IT}(undef, 8, total_hex)
    out_tags  = Vector{Int}(undef, total_hex)

    noff = 0
    eoff = 0
    for m in meshes
        nv = size(m.node, 2)
        out_node[:, noff+1:noff+nv] .= m.node
        @. out_hex[:, eoff+1:eoff+m.hexnum] = m.hexes + IT(noff)
        out_tags[eoff+1:eoff+m.hexnum] .= m.tags
        noff += nv
        eoff += m.hexnum
    end

    return HexahedraMesh(total_hex, out_node, out_hex, out_tags)
end
