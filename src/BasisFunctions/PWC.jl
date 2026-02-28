using StaticArrays
using LinearAlgebra

"""
    PWC{IT, FT}

Piecewise Constant (PWC) basis function group.

Each tetrahedron contributes 3 scalar unknowns (x, y, z components of the 
piecewise-constant electric current/flux). The `inBfsID` field stores the 
3 global basis function IDs corresponding to the x, y, z components.

# Legacy Parity
Matches `MoM_Basics` PWC definition where `nPWC = 3 * num_tetrahedra` and
`inBfsID = [3(t-1)+1, 3(t-1)+2, 3(t-1)+3]` for tetrahedron `t`.
"""
struct PWC{IT, FT}
    id          ::IT            # Tetrahedron index (group ID)
    support     ::IT            # Support tetrahedron (same as id for PWC)
    center      ::SVector{3, FT}
    volume      ::FT
    inBfsID     ::SVector{3, IT}  # 3 global basis function IDs (x,y,z components)
end

"""
    PWCBasis{IT, FT} <: AbstractBasisFunction

Collection of PWC basis functions on a tetrahedral mesh.

Each tetrahedron has 3 DOFs (x, y, z components), so the total number
of basis functions is `3 * length(functions)`.
"""
struct PWCBasis{IT, FT} <: AbstractBasisFunction
    mesh::TetrahedraMesh{IT, FT}
    functions::Vector{PWC{IT, FT}}
end

"""Total number of unknowns: 3 per tetrahedron."""
CoreModule.num_basis(basis::PWCBasis) = 3 * length(basis.functions)

function CoreModule.support(basis::PWCBasis, i::Int)
    # Map global basis ID to tetrahedron index
    tet_idx = div(i - 1, 3) + 1
    return basis.functions[tet_idx].support
end

"""
    evaluate(basis::PWCBasis, i, r)

Evaluate the i-th PWC basis function at point r.
PWC functions are constant unit vectors (x̂, ŷ, ẑ) over the support tetrahedron.
- i mod 3 == 1 → x̂
- i mod 3 == 2 → ŷ  
- i mod 3 == 0 → ẑ
"""
function CoreModule.evaluate(basis::PWCBasis, i::Int, r::AbstractVector)
    comp = mod1(i, 3)  # 1=x, 2=y, 3=z
    if comp == 1
        return SVector(1.0, 0.0, 0.0)
    elseif comp == 2
        return SVector(0.0, 1.0, 0.0)
    else
        return SVector(0.0, 0.0, 1.0)
    end
end

"""
    PWCBasis(mesh::TetrahedraMesh)

Construct PWC basis functions from a tetrahedral mesh.
Three basis functions per tetrahedron (x, y, z components).
Total unknowns: 3 * num_tetrahedra.
"""
function PWCBasis(mesh::TetrahedraMesh{IT, FT}) where {IT, FT}
    nt = num_elements(mesh)
    tets = elements(mesh)
    nodes = vertices(mesh)
    
    functions = Vector{PWC{IT, FT}}(undef, nt)
    
    for t in 1:nt
        v1, v2, v3, v4 = tets[:, t]
        
        v_a = nodes[:, v1]
        v_b = nodes[:, v2]
        v_c = nodes[:, v3]
        v_d = nodes[:, v4]
        
        # Center
        center = (v_a + v_b + v_c + v_d) / 4
        
        # Volume = |det(v_b-v_a, v_c-v_a, v_d-v_a)| / 6
        vol = abs(det(hcat(v_b-v_a, v_c-v_a, v_d-v_a))) / 6
        
        # 3 global basis function IDs for this tetrahedron
        bf_ids = SVector{3, IT}(3*(t-1)+1, 3*(t-1)+2, 3*(t-1)+3)
        
        functions[t] = PWC(
            IT(t),
            IT(t),
            SVector{3, FT}(center),
            vol,
            bf_ids
        )
    end
    
    return PWCBasis(mesh, functions)
end
