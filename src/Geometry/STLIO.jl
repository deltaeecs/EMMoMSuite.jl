"""
    STLIO.jl

STL (Stereolithography) format I/O for TriangleMesh.

Exported:
  read_stl_mesh   — read ASCII or binary STL → TriangleMesh (shared-node topology)
  write_stl_mesh  — write TriangleMesh → ASCII STL
"""

using LinearAlgebra

export read_stl_mesh, write_stl_mesh

# ─────────────────────────────────────────────────────────────────────────────
# read_stl_mesh
# ─────────────────────────────────────────────────────────────────────────────

"""
    read_stl_mesh(pathname::String; FT=Float64, tol=1e-10) → TriangleMesh

Read a binary or ASCII STL file and return a `TriangleMesh` with shared-node
topology.

## Details

STL stores **unshared** vertices (each triangle owns 3 independent vertex
records).  After reading the raw triangles, `remove_duplicate_nodes` is called
with the given `tol` to merge coincident nodes and reconstruct the adjacency
graph.  The merged mesh may have fewer nodes than `3 × Ntri`.
"""
function read_stl_mesh(pathname::String; FT::Type{<:AbstractFloat}=Float64, tol::Real=1e-10)
    isfile(pathname) || error("File not found: $pathname")

    # Peek at the first 80 bytes to decide ASCII vs binary
    raw = read(pathname)            # read entire file as bytes
    is_ascii = _is_ascii_stl(raw)

    if is_ascii
        tris_raw = _read_ascii_stl(String(raw), FT)
    else
        tris_raw = _read_binary_stl(raw, FT)
    end

    Ntri = length(tris_raw)
    Ntri == 0 && return TriangleMesh(0, zeros(FT, 3, 0), zeros(Int, 3, 0), Int[])

    # Build unshared-node mesh (node i of tri t → column 3*(t-1)+i)
    node_raw = Matrix{FT}(undef, 3, 3 * Ntri)
    tris     = Matrix{Int}(undef, 3, Ntri)

    for (ti, (v1, v2, v3)) in enumerate(tris_raw)
        base = 3 * (ti - 1)
        node_raw[:, base + 1] .= v1
        node_raw[:, base + 2] .= v2
        node_raw[:, base + 3] .= v3
        tris[:, ti] = [base + 1, base + 2, base + 3]
    end

    mesh_raw = TriangleMesh(Ntri, node_raw, tris, ones(Int, Ntri))

    # Merge coincident nodes
    merged = remove_duplicate_nodes(mesh_raw; tol=FT(tol))
    return merged
end

# ─── ASCII STL parser ─────────────────────────────────────────────────────────

function _is_ascii_stl(raw::Vector{UInt8})
    # Binary STL is ≥ 84 bytes, first 80 bytes are a free-form header.
    # ASCII STL starts with "solid".
    length(raw) < 5 && return false
    header_bytes = String(raw[1:min(80, length(raw))])
    # If the file is a valid binary STL, we can check the size:
    #   84 + 50 * n_tris == length(raw)
    if length(raw) >= 84
        n = reinterpret(UInt32, raw[81:84])[1]
        expected_size = 84 + 50 * Int(n)
        expected_size == length(raw) && return false   # binary signature match
    end
    return startswith(lstrip(header_bytes), "solid")
end

function _read_ascii_stl(content::String, ::Type{FT}) where {FT}
    tris = NTuple{3, Vector{FT}}[]          # each element is (v1, v2, v3)

    v_buf = Vector{Vector{FT}}()

    for line in eachline(IOBuffer(content))
        s = strip(line)
        if startswith(s, "vertex")
            parts = split(s)
            length(parts) >= 4 || continue
            x = parse(FT, parts[2])
            y = parse(FT, parts[3])
            z = parse(FT, parts[4])
            push!(v_buf, [x, y, z])
        elseif startswith(s, "endfacet")
            if length(v_buf) >= 3
                push!(tris, (v_buf[1], v_buf[2], v_buf[3]))
            end
            empty!(v_buf)
        end
    end

    return tris
end

# ─── Binary STL parser ────────────────────────────────────────────────────────

function _read_binary_stl(raw::Vector{UInt8}, ::Type{FT}) where {FT}
    length(raw) < 84 && error("File too short to be a valid binary STL")

    n_tris = Int(reinterpret(UInt32, raw[81:84])[1])
    length(raw) < 84 + 50 * n_tris &&
        error("Binary STL: declared $(n_tris) triangles but file is too short")

    tris = Vector{Tuple{Vector{FT}, Vector{FT}, Vector{FT}}}(undef, n_tris)
    ptr = 85      # 1-based; first triangle starts at byte 85

    for ti in 1:n_tris
        # Skip normal (3 × Float32 = 12 bytes)
        ptr += 12
        # v1, v2, v3 (each 3 × Float32 = 12 bytes); @view avoids per-vertex allocation
        v1 = FT.(reinterpret(Float32, @view raw[ptr:ptr+11]));  ptr += 12
        v2 = FT.(reinterpret(Float32, @view raw[ptr:ptr+11]));  ptr += 12
        v3 = FT.(reinterpret(Float32, @view raw[ptr:ptr+11]));  ptr += 12
        # Attribute byte count (2 bytes, usually 0)
        ptr += 2
        tris[ti] = (v1, v2, v3)
    end

    return tris
end

# ─────────────────────────────────────────────────────────────────────────────
# write_stl_mesh
# ─────────────────────────────────────────────────────────────────────────────

"""
    write_stl_mesh(pathname::String, mesh::TriangleMesh; name="EMSuite")

Write `mesh` to an ASCII STL file.

The facet normal for each triangle is computed from the vertex ordering (`v2-v1`
× `v3-v1`, normalised).  If the triangle is degenerate (zero-area), a zero
normal is written.
"""
function write_stl_mesh(pathname::String, mesh::TriangleMesh; name::String="EMSuite")
    open(pathname, "w") do f
        println(f, "solid $name")

        for t in 1:mesh.trinum
            v1 = mesh.node[:, mesh.triangles[1, t]]
            v2 = mesh.node[:, mesh.triangles[2, t]]
            v3 = mesh.node[:, mesh.triangles[3, t]]

            n_vec = cross(v2 - v1, v3 - v1)
            len   = norm(n_vec)
            n_hat = len > 0 ? n_vec ./ len : zeros(3)

            println(f, "  facet normal $(n_hat[1]) $(n_hat[2]) $(n_hat[3])")
            println(f, "    outer loop")
            println(f, "      vertex $(v1[1]) $(v1[2]) $(v1[3])")
            println(f, "      vertex $(v2[1]) $(v2[2]) $(v2[3])")
            println(f, "      vertex $(v3[1]) $(v3[2]) $(v3[3])")
            println(f, "    endloop")
            println(f, "  endfacet")
        end

        println(f, "endsolid $name")
    end
end
