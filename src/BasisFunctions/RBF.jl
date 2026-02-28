using StaticArrays
using LinearAlgebra

"""
    RBF{IT, FT}

Rooftop Basis Function (RBF).
Defined on faces shared by hexahedra.
"""
struct RBF{IT, FT}
    id          ::IT
    is_boundary ::Bool
    area        ::FT
    
    # Support hexahedra (indices in mesh)
    support     ::SVector{2, IT} 
    
    # Local face index in the support hexahedron (1 to 6)
    local_face_idx ::SVector{2, IT}
    
    # Sign relative to the face normal
    signs       ::SVector{2, Int}
    
    center      ::SVector{3, FT}
end

"""
    RBFBasis{IT, FT} <: AbstractBasisFunction

Collection of RBF basis functions on a hexahedral mesh.
"""
struct RBFBasis{IT, FT} <: AbstractBasisFunction
    mesh::HexahedraMesh{IT, FT}
    functions::Vector{RBF{IT, FT}}
end

CoreModule.num_basis(basis::RBFBasis) = length(basis.functions)

function CoreModule.support(basis::RBFBasis, i::Int)
    return basis.functions[i].support
end

function CoreModule.evaluate(basis::RBFBasis, i::Int, r::AbstractVector)
    # RBF vector ρ = r - r₀
    # For correct evaluation in assembly, use get_free_vns and gq3d_to_face2d_idx
    # This generic evaluate is a fallback using face centroid as free-end approximation
    bf = basis.functions[i]
    # Use the center of the basis function (face center) as an approximation
    # For accurate computation in assembly, the free-end r₀ should be interpolated
    # on the opposite face using parametric coordinates
    return SVector{3, Float64}(r[1] - bf.center[1], r[2] - bf.center[2], r[3] - bf.center[3])
end

"""
    RBFBasis(mesh::HexahedraMesh)

Construct RBF basis functions from a hexahedral mesh.
"""
function RBFBasis(mesh::HexahedraMesh{IT, FT}) where {IT, FT}
    nt = num_elements(mesh)
    hexes = elements(mesh)
    nodes = vertices(mesh)
    
    # Store face info: (v1, v2, v3, v4, hex_idx, local_face_idx)
    # Sorted vertices for unique identification
    FaceInfo = Tuple{IT, IT, IT, IT, IT, IT}
    all_faces = Vector{FaceInfo}(undef, nt * 6)
    
    # Face definitions (indices into 1-8 vertices of hex)
    # 1. (2,3,7,6)
    # 2. (1,4,8,5)
    # 3. (4,3,7,8)
    # 4. (1,2,6,5)
    # 5. (5,6,7,8)
    # 6. (1,2,3,4)
    face_indices = [
        (2,3,7,6), (1,4,8,5), (4,3,7,8),
        (1,2,6,5), (5,6,7,8), (1,2,3,4)
    ]
    
    idx = 1
    for t in 1:nt
        h = hexes[:, t]
        
        for (f_idx, indices) in enumerate(face_indices)
            v1 = h[indices[1]]
            v2 = h[indices[2]]
            v3 = h[indices[3]]
            v4 = h[indices[4]]
            
            # Sort vertices to identify face uniquely
            sorted_v = sort(SVector(v1, v2, v3, v4))
            
            all_faces[idx] = (sorted_v[1], sorted_v[2], sorted_v[3], sorted_v[4], t, f_idx)
            idx += 1
        end
    end
    
    # Sort faces to find pairs
    sort!(all_faces, by = x -> (x[1], x[2], x[3], x[4]))
    
    functions = Vector{RBF{IT, FT}}()
    
    i = 1
    while i <= length(all_faces)
        f1 = all_faces[i]
        
        # Check if next face is the same (internal face)
        if i < length(all_faces) && 
           all_faces[i+1][1] == f1[1] && 
           all_faces[i+1][2] == f1[2] && 
           all_faces[i+1][3] == f1[3] &&
           all_faces[i+1][4] == f1[4]
            
            f2 = all_faces[i+1]
            
            # Internal face
            # Calculate area and center
            # Assuming planar quad for area approximation (split into 2 triangles)
            # Or just use centroid
            v_a = nodes[:, f1[1]]
            v_b = nodes[:, f1[2]]
            v_c = nodes[:, f1[3]]
            v_d = nodes[:, f1[4]]
            
            center = (v_a + v_b + v_c + v_d) / 4
            
            # Area: sum of two triangles (a,b,c) and (a,c,d) ? 
            # Depends on ordering in the quad. 
            # Since we sorted vertices, the order is lost for geometry.
            # But we can use the original vertices from the mesh if we kept them.
            # For now, let's approximate or retrieve original vertices.
            # To retrieve original, we need to look up the element and face index.
            # But we have that in f1[5] (hex_idx) and f1[6] (local_face_idx).
            
            # Let's get original vertices for area calculation
            h_orig = hexes[:, f1[5]]
            indices_orig = face_indices[f1[6]]
            v1_o = nodes[:, h_orig[indices_orig[1]]]
            v2_o = nodes[:, h_orig[indices_orig[2]]]
            v3_o = nodes[:, h_orig[indices_orig[3]]]
            v4_o = nodes[:, h_orig[indices_orig[4]]]
            
            # Area of quad v1-v2-v3-v4
            # Split into v1-v2-v3 and v1-v3-v4
            area1 = 0.5 * norm(cross(v2_o - v1_o, v3_o - v1_o))
            area2 = 0.5 * norm(cross(v3_o - v1_o, v4_o - v1_o))
            area = area1 + area2
            
            rbf = RBF(
                IT(length(functions) + 1),
                false,
                area,
                SVector(f1[5], f2[5]),
                SVector(f1[6], f2[6]),
                SVector(1, -1), # Placeholder signs
                SVector{3, FT}(center)
            )
            push!(functions, rbf)
            
            i += 2
        else
            # Boundary face
            h_orig = hexes[:, f1[5]]
            indices_orig = face_indices[f1[6]]
            v1_o = nodes[:, h_orig[indices_orig[1]]]
            v2_o = nodes[:, h_orig[indices_orig[2]]]
            v3_o = nodes[:, h_orig[indices_orig[3]]]
            v4_o = nodes[:, h_orig[indices_orig[4]]]
            
            area1 = 0.5 * norm(cross(v2_o - v1_o, v3_o - v1_o))
            area2 = 0.5 * norm(cross(v3_o - v1_o, v4_o - v1_o))
            area = area1 + area2
            
            center = (v1_o + v2_o + v3_o + v4_o) / 4
            
            rbf = RBF(
                IT(length(functions) + 1),
                true,
                area,
                SVector(f1[5], f1[5]),
                SVector(f1[6], f1[6]),
                SVector(1, 0),
                SVector{3, FT}(center)
            )
            push!(functions, rbf)
            
            i += 1
        end
    end
    
    return RBFBasis(mesh, functions)
end
