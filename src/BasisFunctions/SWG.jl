using StaticArrays
using LinearAlgebra

"""
    SWG{IT, FT}

SWG (Schaubert-Wilton-Glisson) basis function data structure.

Represents a divergence-conforming basis function defined on a pair of tetrahedra sharing a common face.
Commonly used in Volume Integral Equation (VIE) methods.

# Mathematical Definition

The SWG basis function \$\\mathbf{f}_n(\\mathbf{r})\$ associated with the \$n\$-th face is defined as:

```math
\\mathbf{f}_n(\\mathbf{r}) = \\begin{cases}
\\frac{A_n}{3V_n^+} \\boldsymbol{\\rho}_n^+ & \\mathbf{r} \\in T_n^+ \\\\
\\frac{A_n}{3V_n^-} \\boldsymbol{\\rho}_n^- & \\mathbf{r} \\in T_n^- \\\\
0 & \\text{otherwise}
\\end{cases}
```

where:
- \$A_n\$ is the area of the common face.
- \$V_n^\\pm\$ is the volume of the tetrahedron \$T_n^\\pm\$.
- \$\\boldsymbol{\\rho}_n^+ = \\mathbf{r} - \\mathbf{v}_n^+\$, where \$\\mathbf{v}_n^+\$ is the free vertex of \$T_n^+\$.
- \$\\boldsymbol{\\rho}_n^- = \\mathbf{v}_n^- - \\mathbf{r}\$, where \$\\mathbf{v}_n^-\$ is the free vertex of \$T_n^-\$.

The volume divergence is constant in each tetrahedron:

```math
\\nabla \\cdot \\mathbf{f}_n(\\mathbf{r}) = \\begin{cases}
\\frac{A_n}{V_n^+} & \\mathbf{r} \\in T_n^+ \\\\
-\\frac{A_n}{V_n^-} & \\mathbf{r} \\in T_n^- \\\\
0 & \\text{otherwise}
\\end{cases}
```

# Fields
- `id`: Unique identifier for the basis function.
- `is_boundary`: True if the face is on the boundary (only one support tetrahedron).
- `area`: Area of the common face (\$A_n\$).
- `support`: Indices of the two support tetrahedra (or one if boundary).
- `local_face_idx`: Local index (1-4) of the face in each support tetrahedron.
- `signs`: Orientation sign relative to the local face definition.
- `center`: Centroid of the common face.
"""
struct SWG{IT,FT}
    id::IT
    is_boundary::Bool
    area::FT

    # Support tetrahedra (indices in mesh)
    support::SVector{2,IT}

    # Local face index in the support tetrahedron (1 to 4)
    local_face_idx::SVector{2,IT}

    # Sign relative to the face normal
    signs::SVector{2,Int}

    center::SVector{3,FT}
end

"""
    SWGBasis{IT, FT} <: AbstractBasisFunction

Collection of SWG basis functions defined on a tetrahedral mesh.

This structure manages the mapping between mesh faces and basis functions.

# Fields
- `mesh`: The underlying tetrahedral mesh.
- `functions`: Vector of `SWG` basis function objects.
"""
struct SWGBasis{IT,FT} <: AbstractBasisFunction
    mesh::TetrahedraMesh{IT,FT}
    functions::Vector{SWG{IT,FT}}
end

CoreModule.num_basis(basis::SWGBasis) = length(basis.functions)

function CoreModule.support(basis::SWGBasis, i::Int)
    return basis.functions[i].support
end

function CoreModule.evaluate(basis::SWGBasis, i::Int, r::AbstractVector)
    # TODO: Implement SWG evaluation
    return SVector(0.0, 0.0, 0.0)
end

"""
    SWGBasis(mesh::TetrahedraMesh)

Construct a set of SWG basis functions from a tetrahedral mesh.

# Algorithm
1.  **Face Extraction**: Iterates through all tetrahedra to identify all faces.
2.  **Face Matching**: Sorts faces by vertex indices to find pairs of tetrahedra sharing a common face.
3.  **Basis Creation**:
    - For each internal face (shared by two tetrahedra), creates an `SWG` basis function.
    - Assigns the "plus" and "minus" tetrahedra.
    - Computes face area and orientation signs.
4.  **Boundary Handling**: Boundary faces are typically ignored unless specific boundary conditions are applied.

# Arguments
- `mesh`: A `TetrahedraMesh` object defining the geometry.

# Returns
- An `SWGBasis` object containing the generated basis functions.
"""
function SWGBasis(mesh::TetrahedraMesh{IT,FT}) where {IT,FT}
    nt = num_elements(mesh)
    tets = elements(mesh)
    nodes = vertices(mesh)

    # Store face info: (v1, v2, v3, tet_idx, local_face_idx)
    # v1 < v2 < v3 for unique identification
    FaceInfo = Tuple{IT,IT,IT,IT,IT}
    all_faces = Vector{FaceInfo}(undef, nt * 4)

    idx = 1
    for t = 1:nt
        v1, v2, v3, v4 = tets[:, t]

        # Face 1: v2, v3, v4
        f1 = sort(SVector(v2, v3, v4))
        all_faces[idx] = (f1[1], f1[2], f1[3], t, 1)
        idx += 1

        # Face 2: v1, v4, v3
        f2 = sort(SVector(v1, v4, v3))
        all_faces[idx] = (f2[1], f2[2], f2[3], t, 2)
        idx += 1

        # Face 3: v1, v2, v4
        f3 = sort(SVector(v1, v2, v4))
        all_faces[idx] = (f3[1], f3[2], f3[3], t, 3)
        idx += 1

        # Face 4: v1, v3, v2
        f4 = sort(SVector(v1, v3, v2))
        all_faces[idx] = (f4[1], f4[2], f4[3], t, 4)
        idx += 1
    end

    # Sort faces to find pairs
    sort!(all_faces, by = x -> (x[1], x[2], x[3]))

    functions = Vector{SWG{IT,FT}}()

    i = 1
    while i <= length(all_faces)
        f1 = all_faces[i]

        # Check if next face is the same (internal face)
        if i < length(all_faces) &&
           all_faces[i+1][1] == f1[1] &&
           all_faces[i+1][2] == f1[2] &&
           all_faces[i+1][3] == f1[3]

            f2 = all_faces[i+1]

            # Internal face
            # Calculate area and center
            v_a = nodes[:, f1[1]]
            v_b = nodes[:, f1[2]]
            v_c = nodes[:, f1[3]]

            # Area of triangle
            ab = v_b - v_a
            ac = v_c - v_a
            area = 0.5 * norm(cross(ab, ac))
            center = (v_a + v_b + v_c) / 3

            swg = SWG(
                IT(length(functions) + 1),
                false,
                area,
                SVector(f1[4], f2[4]),
                SVector(f1[5], f2[5]),
                SVector(1, -1), # Placeholder signs
                SVector{3,FT}(center),
            )
            push!(functions, swg)

            i += 2
        else
            # Boundary face
            v_a = nodes[:, f1[1]]
            v_b = nodes[:, f1[2]]
            v_c = nodes[:, f1[3]]

            ab = v_b - v_a
            ac = v_c - v_a
            area = 0.5 * norm(cross(ab, ac))
            center = (v_a + v_b + v_c) / 3

            swg = SWG(
                IT(length(functions) + 1),
                true,
                area,
                SVector(f1[4], 0),
                SVector(f1[5], 0),
                SVector(1, 0),
                SVector{3,FT}(center),
            )
            push!(functions, swg) # Uncomment to include boundary faces

            i += 1
        end
    end

    return SWGBasis(mesh, functions)
end

"""
    evaluate_swg(bf::SWG, i_supp::Int, r::AbstractVector, verts::AbstractMatrix, vol::Real)

Evaluate the SWG basis function and its divergence at point r.
"""
function evaluate_swg(bf::SWG, i_supp::Int, r::AbstractVector, verts::AbstractMatrix, vol::Real)
    local_face = bf.local_face_idx[i_supp]
    v_free = verts[:, local_face]

    sign = bf.signs[i_supp]

    if sign > 0
        rho = r - v_free
    else
        rho = v_free - r
    end

    area = bf.area

    f = (area / (3 * vol)) * rho
    div_f = (area / vol) * sign

    return f, div_f
end
