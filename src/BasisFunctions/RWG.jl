using StaticArrays
using LinearAlgebra

"""
    RWG{IT, FT}

RWG (Rao-Wilton-Glisson) basis function data structure.

Represents a divergence-conforming basis function defined on a pair of triangles sharing a common edge.

# Mathematical Definition

The RWG basis function \$\\mathbf{f}_n(\\mathbf{r})\$ associated with the \$n\$-th edge is defined as:

```math
\\mathbf{f}_n(\\mathbf{r}) = \\begin{cases}
\\frac{l_n}{2A_n^+} \\boldsymbol{\\rho}_n^+ & \\mathbf{r} \\in T_n^+ \\\\
\\frac{l_n}{2A_n^-} \\boldsymbol{\\rho}_n^- & \\mathbf{r} \\in T_n^- \\\\
0 & \\text{otherwise}
\\end{cases}
```

where:
- \$l_n\$ is the length of the common edge.
- \$A_n^\\pm\$ is the area of the triangle \$T_n^\\pm\$.
- \$\\boldsymbol{\\rho}_n^+ = \\mathbf{r} - \\mathbf{v}_n^+\$, where \$\\mathbf{v}_n^+\$ is the free vertex of \$T_n^+\$.
- \$\\boldsymbol{\\rho}_n^- = \\mathbf{v}_n^- - \\mathbf{r}\$, where \$\\mathbf{v}_n^-\$ is the free vertex of \$T_n^-\$.

The surface divergence is constant in each triangle:

```math
\\nabla_s \\cdot \\mathbf{f}_n(\\mathbf{r}) = \\begin{cases}
\\frac{l_n}{A_n^+} & \\mathbf{r} \\in T_n^+ \\\\
-\\frac{l_n}{A_n^-} & \\mathbf{r} \\in T_n^- \\\\
0 & \\text{otherwise}
\\end{cases}
```

# Fields
- `id`: Unique identifier for the basis function.
- `is_boundary`: True if the edge is on the boundary (only one support triangle).
- `edge_length`: Length of the common edge (\$l_n\$).
- `support`: Indices of the two support triangles (or one if boundary).
- `local_edge_idx`: Local index (1-3) of the edge in each support triangle.
- `signs`: Orientation sign relative to the local edge definition (+1 or -1).
- `center`: Midpoint of the common edge.
"""
struct RWG{IT,FT}
    id::IT
    is_boundary::Bool
    edge_length::FT

    # Support triangles (indices in mesh)
    # For boundary edges, the second index is 0 or same as first
    support::SVector{2,IT}

    # Local edge index in the support triangle (1, 2, or 3)
    local_edge_idx::SVector{2,IT}

    # Sign relative to the edge direction in the triangle
    # +1 if edge direction matches basis direction, -1 otherwise
    # Typically RWG is defined as flowing from T+ to T-
    # We store signs to handle orientation
    signs::SVector{2,Int}

    center::SVector{3,FT}
end

"""
    RWGBasis{IT, FT} <: AbstractBasisFunction

Collection of RWG basis functions defined on a triangular mesh.

This structure manages the mapping between mesh edges and basis functions.
It handles the connectivity logic to identify common edges and assign orientations.

# Boundary Edge Design

Unlike `SWGBasis` (which **includes** boundary faces for VEFIE/SCFIE flux continuity),
`RWGBasis` **excludes** boundary edges. This is intentional: for Surface Integral 
Equations on Perfect Electric Conductor (PEC) surfaces, current cannot flow out of 
the closed surface, so boundary edges represent unphysical discontinuities.

Only internal edges (shared between two triangles) are assigned RWG basis functions.
Boundary edges are marked in `basis_map` with ID = 0.

# Fields
- `mesh`: The underlying triangular mesh.
- `functions`: Vector of `RWG` basis function objects (excludes boundary edges).
- `basis_map`: A `3 \\times N_t` matrix mapping (local edge, triangle) pairs to basis function IDs.
  - `basis_map[k, t]` gives the ID of the basis function associated with the \$k\$-th edge of triangle \$t\$.
  - If 0, no basis function is assigned (boundary edge excluded for PEC EFIE).
"""
struct RWGBasis{IT,FT} <: AbstractBasisFunction
    mesh::TriangleMesh{IT,FT}
    functions::Vector{RWG{IT,FT}}
    # Map from (local_edge, triangle) to Basis ID
    # 3 x Nt matrix. 0 indicates no basis (or boundary if not handled)
    basis_map::Matrix{IT}
end

CoreModule.num_basis(basis::RWGBasis) = length(basis.functions)

function CoreModule.support(basis::RWGBasis, i::Int)
    return basis.functions[i].support
end

function CoreModule.evaluate(basis::RWGBasis, i::Int, r::AbstractVector)
    error(
        "RWGBasis.evaluate() is not implemented. " *
        "Use evaluate_rwg(basis.functions[i], supp_idx, r, vertices, area) for direct evaluation, " *
        "or access precomputed quadrature points via TriangleInfoConstructor.",
    )
end

"""
    RWGBasis(mesh::TriangleMesh)

Construct a set of RWG basis functions from a triangular mesh.

# Algorithm
1.  **Edge Extraction**: Iterates through all triangles to identify all edges.
2.  **Edge Matching**: Sorts edges by vertex indices to find pairs of triangles sharing a common edge.
3.  **Basis Creation**:
    - For each internal edge (shared by two triangles), creates an `RWG` basis function.
    - Assigns the "plus" and "minus" triangles based on the global vertex ordering or explicit convention.
    - Computes edge length and orientation signs.
4.  **Boundary Handling**: Boundary edges (shared by only one triangle) are typically ignored for EFIE/MFIE unless specific boundary conditions are applied (e.g., port excitations).

# Arguments
- `mesh`: A `TriangleMesh` object defining the geometry.

# Returns
- An `RWGBasis` object containing the generated basis functions.
"""
function RWGBasis(mesh::TriangleMesh{IT,FT}) where {IT,FT}
    # 1. Extract all edges
    # Each triangle has 3 edges.
    # Edge k connects vertices (k+1)%3 and (k+2)%3 ? 
    # Let's define local edges:
    # Edge 1: v2 -> v3
    # Edge 2: v3 -> v1
    # Edge 3: v1 -> v2

    nt = num_elements(mesh)
    tris = elements(mesh)
    nodes = vertices(mesh)

    # Store edge info: (min_v, max_v, tri_idx, local_edge_idx)
    # We use min/max to identify the edge regardless of direction
    EdgeInfo = Tuple{IT,IT,IT,IT}
    all_edges = Vector{EdgeInfo}(undef, nt * 3)

    idx = 1
    # Match Legacy Order: [All Edge 1s, All Edge 2s, All Edge 3s]
    # This is crucial for stable sort to produce the same T+/T- assignment as legacy code.

    # Edge 1: v2-v3
    for t = 1:nt
        v1, v2, v3 = tris[:, t]
        all_edges[idx] = (min(v2, v3), max(v2, v3), t, 1)
        idx += 1
    end

    # Edge 2: v3-v1
    for t = 1:nt
        v1, v2, v3 = tris[:, t]
        all_edges[idx] = (min(v3, v1), max(v3, v1), t, 2)
        idx += 1
    end

    # Edge 3: v1-v2
    for t = 1:nt
        v1, v2, v3 = tris[:, t]
        all_edges[idx] = (min(v1, v2), max(v1, v2), t, 3)
        idx += 1
    end

    # 2. Sort edges to find pairs
    sort!(all_edges, by = x -> (x[1], x[2]), alg = Base.Sort.MergeSort)

    # 3. Create RWG functions
    functions = Vector{RWG{IT,FT}}()
    basis_map = zeros(IT, 3, nt)

    i = 1
    while i <= length(all_edges)
        e1 = all_edges[i]

        # Check if next edge is the same (internal edge)
        if i < length(all_edges) && all_edges[i+1][1] == e1[1] && all_edges[i+1][2] == e1[2]
            e2 = all_edges[i+1]

            # Internal edge
            # Construct RWG
            # We need to calculate edge length and center
            v_start = nodes[:, e1[1]]
            v_end = nodes[:, e1[2]]
            len = norm(v_start - v_end)
            center = (v_start + v_end) / 2

            # Determine signs/orientation
            # This requires careful definition. 
            # Standard RWG: f(r) = L/2A * rho
            # Current flows from T+ to T- across the edge.
            # Let's assign T+ = e1.tri, T- = e2.tri

            rwg_id = IT(length(functions) + 1)
            rwg = RWG(
                rwg_id,
                false,
                len,
                SVector(e1[3], e2[3]),
                SVector(e1[4], e2[4]),
                SVector(1, -1), # Placeholder signs
                SVector{3,FT}(center),
            )
            push!(functions, rwg)

            # Update map
            basis_map[e1[4], e1[3]] = rwg_id
            basis_map[e2[4], e2[3]] = rwg_id

            i += 2
        else
            # Boundary edge
            # Only one triangle
            # Construct boundary RWG (half basis) if needed, or skip
            # Usually we skip boundary edges for PEC EFIE
            # But let's store it marked as boundary

            v_start = nodes[:, e1[1]]
            v_end = nodes[:, e1[2]]
            len = norm(v_start - v_end)
            center = (v_start + v_end) / 2

            rwg_id = IT(length(functions) + 1)
            rwg = RWG(
                rwg_id,
                true,
                len,
                SVector(e1[3], e1[3]), # Second tri is same (or 0)
                SVector(e1[4], e1[4]),
                SVector(1, 0),
                SVector{3,FT}(center),
            )
            # Boundary edge: skip for PEC EFIE.
            # Unlike SWGBasis (which includes boundary faces for VEFIE flux continuity),
            # RWGBasis excludes boundary edges because there is no current flowing out
            # of a closed PEC surface along a boundary edge.
            # push!(functions, rwg) -- boundary RWG not used; basis_map stays 0

            i += 1
        end
    end

    return RWGBasis(mesh, functions, basis_map)
end

"""
    RWGBasis(comp::CompositeMesh)

Construct an `RWGBasis` from a `CompositeMesh` by delegating to the embedded
surface (`comp.surface`).  This allows SCFIE workflows to pass a single
`CompositeMesh` (which carries both surface and volume information) without
manually extracting the surface first.
"""
RWGBasis(comp::CompositeMesh) = RWGBasis(comp.surface)
