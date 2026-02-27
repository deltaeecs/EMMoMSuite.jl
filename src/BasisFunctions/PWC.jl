using StaticArrays
using LinearAlgebra

"""
    PWC{IT, FT}

Piecewise Constant (PWC) basis function.
Defined on a single tetrahedron.
"""
struct PWC{IT, FT}
    id          ::IT
    
    # Support tetrahedron (index in mesh)
    support     ::IT
    
    center      ::SVector{3, FT}
    volume      ::FT
end

"""
    PWCBasis{IT, FT} <: AbstractBasisFunction

Collection of PWC basis functions on a tetrahedral mesh.
"""
struct PWCBasis{IT, FT} <: AbstractBasisFunction
    mesh::TetrahedraMesh{IT, FT}
    functions::Vector{PWC{IT, FT}}
end

CoreModule.num_basis(basis::PWCBasis) = length(basis.functions)

function CoreModule.support(basis::PWCBasis, i::Int)
    return basis.functions[i].support
end

function CoreModule.evaluate(basis::PWCBasis, i::Int, r::AbstractVector)
    # TODO: Implement PWC evaluation
    return SVector(0.0, 0.0, 0.0)
end

"""
    PWCBasis(mesh::TetrahedraMesh)

Construct PWC basis functions from a tetrahedral mesh.
One basis function per tetrahedron.
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
        
        functions[t] = PWC(
            IT(t),
            IT(t),
            SVector{3, FT}(center),
            vol
        )
    end
    
    return PWCBasis(mesh, functions)
end
