module RadiationIntegral

using StaticArrays
using LinearAlgebra
using ...Geometry
using ...Utilities.Parameters
using ...BasisFunctions
using ...CoreModule: num_elements, vertices, elements

export r̂θϕInfo, raditionalIntegralNθϕCal, radiation_integral_rwg, radiation_integral_swg, radiation_integral_pwc, radiation_integral_pwc_hex, radiation_integral_rbf, ∠Info

struct ∠Info{FT<:Real}
    val::FT
    sin::FT
    cos::FT
end
∠Info(val::FT) where FT = ∠Info(val, sin(val), cos(val))

struct r̂θϕInfo{FT<:Real}
    r̂::SVector{3, FT}
    θhat::SVector{3, FT}
    ϕhat::SVector{3, FT}
    θϕ::SVector{2, FT}
end

function r̂θϕInfo(θ::∠Info{FT}, ϕ::∠Info{FT}) where FT
    st, ct = θ.sin, θ.cos
    sp, cp = ϕ.sin, ϕ.cos
    r̂ = SVector(st*cp, st*sp, ct)
    θhat = SVector(ct*cp, ct*sp, -st)
    ϕhat = SVector(-sp, cp, zero(FT))
    θϕ = SVector(θ.val, ϕ.val)
    return r̂θϕInfo(r̂, θhat, ϕhat, θϕ)
end


"""
    raditionalIntegralNθϕCal(r̂θϕ, trianglesInfo, Jtris)

Calculate the radiation integral N(θ, ϕ) using averaged currents (Legacy/Fast).
"""
function raditionalIntegralNθϕCal(r̂θϕ::r̂θϕInfo{FT}, trianglesInfo::AbstractVector{TriangleInfo{IT, FT}}, 
                                  Jtris::AbstractMatrix{CT}) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    
    Nxyz = zero(MVector{3, CT})
    Nθϕ = zero(MVector{2, CT})
    
    r̂ = r̂θϕ.r̂
    θhat = r̂θϕ.θhat
    ϕhat = r̂θϕ.ϕhat
    
    k0 = get_k0()
    jk0 = im * k0
    
    # Get Quadrature points
    points, weights = Geometry.gaussQuadratureTri(3, FT)
    n_points = length(weights)
    
    for ti in 1:length(trianglesInfo)
        tri = trianglesInfo[ti]
        Jtri = @view Jtris[:, ti]
        
        JSexp = zero(MVector{3, CT})
        
        for gi in 1:n_points
            # Calculate point position
            u = points[1, gi]
            v = points[2, gi]
            w = points[3, gi]
            rgi = u * tri.vertices[:, 1] + v * tri.vertices[:, 2] + w * tri.vertices[:, 3]
            
            phase = exp(jk0 * dot(r̂, rgi))
            JSexp .+= Jtri .* (phase * weights[gi])
        end
        
        JSexp .*= tri.area
        
        Nxyz .+= JSexp
    end
    
    # Debug
    # println("Nxyz: ", Nxyz)
    
    Nθϕ[1] = dot(θhat, Nxyz)
    Nθϕ[2] = dot(ϕhat, Nxyz)
    
    return Nθϕ
end

"""
    radiation_integral_rwg(r̂θϕ, basis, ICoeff)

Calculate the radiation integral N(θ, ϕ) using exact RWG basis functions.
"""
function radiation_integral_rwg(r̂θϕ::r̂θϕInfo{FT}, basis::RWGBasis{IT, FT}, ICoeff::Vector{CT}) where {IT, FT, CT}
    Nxyz = zero(MVector{3, CT})
    Nθϕ = zero(MVector{2, CT})
    
    r̂ = r̂θϕ.r̂
    θhat = r̂θϕ.θhat
    ϕhat = r̂θϕ.ϕhat
    
    k0 = get_k0()
    jk0 = im * k0
    
    # Quadrature
    points, weights = Geometry.gaussQuadratureTri(3, FT)
    n_points = length(weights)
    
    mesh = basis.mesh
    verts = vertices(mesh)
    elems = elements(mesh)
    ntri = num_elements(mesh)
    
    for t in 1:ntri
        # Get triangle vertices
        v_indices = elems[:, t]
        r1 = verts[:, v_indices[1]]
        r2 = verts[:, v_indices[2]]
        r3 = verts[:, v_indices[3]]
        
        # Area
        cross_prod = cross(r2 - r1, r3 - r1)
        area = 0.5 * norm(cross_prod)
        
        JSexp = zero(MVector{3, CT})
        
        for gi in 1:n_points
            u = points[1, gi]
            v = points[2, gi]
            w = points[3, gi]
            rgi = u * r1 + v * r2 + w * r3
            
            # Compute J at rgi
            J_at_r = zero(MVector{3, CT})
            
            for k in 1:3
                bf_id = basis.basis_map[k, t]
                if bf_id == 0 continue end
                
                bf = basis.functions[bf_id]
                
                # Sign
                if bf.support[1] == t
                    sign_val = 1.0
                else
                    sign_val = -1.0
                end
                
                v_op = verts[:, v_indices[k]]
                rho = rgi - v_op
                
                val = ICoeff[bf_id] * (bf.edge_length / (2 * area)) * sign_val * rho
                J_at_r .+= val
            end
            
            phase = exp(jk0 * dot(r̂, rgi))
            JSexp .+= J_at_r .* (phase * weights[gi])
        end
        
        JSexp .*= area
        Nxyz .+= JSexp
    end
    
    Nθϕ[1] = dot(θhat, Nxyz)
    Nθϕ[2] = dot(ϕhat, Nxyz)
    
    return Nθϕ
end

"""
    radiation_integral_swg(r̂θϕ, basis, ICoeff, permittivities)

Calculate the radiation integral N(θ, ϕ) using SWG basis functions for VEFIE.
"""
function radiation_integral_swg(r̂θϕ::r̂θϕInfo{FT}, basis::SWGBasis{IT, FT}, ICoeff::Vector{CT}, permittivities::Vector{CT}) where {IT, FT, CT}
    Nxyz = zero(MVector{3, CT})
    Nθϕ = zero(MVector{2, CT})
    
    r̂ = r̂θϕ.r̂
    θhat = r̂θϕ.θhat
    ϕhat = r̂θϕ.ϕhat
    
    k0 = get_k0()
    jk0 = im * k0
    omega = k0 * 299792458.0 # Approx c0
    
    # Quadrature for Tetrahedron (4 points or 5 points)
    # Use 4 points for efficiency
    points, weights = Geometry.gaussQuadratureTet(4, FT)
    n_points = length(weights)
    
    mesh = basis.mesh
    nodes = mesh.node
    tetras = mesh.tetras
    ntetra = num_elements(mesh)
    
    # Precompute kappa for each tetrahedron
    # kappa = (eps_r - 1) / eps_r
    # J_eq = j * omega * kappa * D
    # So we integrate: j * omega * kappa * D * exp(...)
    
    # Iterate over tetrahedra
    for t in 1:ntetra
        eps_r = permittivities[t]
        kappa = (eps_r - 1.0) / eps_r
        # Legacy Parity: ICoeff are already scaled by j*omega
        # J = kappa * I_coeff * f_n
        factor = kappa
        
        # Get tetra vertices
        v_indices = tetras[:, t]
        r1 = nodes[:, v_indices[1]]
        r2 = nodes[:, v_indices[2]]
        r3 = nodes[:, v_indices[3]]
        r4 = nodes[:, v_indices[4]]
        
        # Volume
        mat = hcat(r2-r1, r3-r1, r4-r1)
        vol = abs(det(mat)) / 6.0
        
        # Accumulate integral for this tetrahedron
        integral_val = zero(MVector{3, CT})
        
        for gi in 1:n_points
            u = points[1, gi]
            v = points[2, gi]
            w = points[3, gi]
            x = points[4, gi] # 1 - u - v - w
            
            rgi = u * r1 + v * r2 + w * r3 + x * r4
            
            # Compute D at rgi
            D_at_r = zero(MVector{3, CT})
            
            # Iterate over 4 faces of the tetrahedron
            # Each face corresponds to a local basis function
            for local_face in 1:4
                # Map local face to global basis function ID
                # We need to know which basis function is defined on this face of this tetrahedron
                # SWGBasis stores this in basis.basis_map? No, SWG is face-based.
                # basis.tet_faces[t, local_face] gives the global face index.
                # Then we check if that face has an assigned basis function.
                
                # Wait, SWGBasis structure:
                # functions::Vector{SWG}
                # Each SWG has plus_tet, minus_tet, common_face.
                
                # We need a map from (tet, local_face) -> basis_id
                # Or iterate over all basis functions and add their contribution.
                # Iterating over basis functions is better if N is large, but here we iterate over elements (integration domain).
                # So we need (tet) -> list of basis functions.
                
                # Let's assume we have a map or we can deduce it.
                # SWGBasis usually has `tet_edges` or similar.
                # Actually, for SWG, the unknowns are on faces.
                # Each tetrahedron has 4 faces.
                # Let's find the global face index.
                # We need to look at `SWGBasis` implementation.
                nothing
            end
        end
    end
    
    # Iterate over Basis Functions (Support)
    for n in 1:length(basis.functions)
        bf = basis.functions[n]
        coeff = ICoeff[n]
        
        # Plus Tetrahedron
        t_plus = bf.support[1]
        if t_plus > 0
            eps_r = permittivities[t_plus]
            kappa = (eps_r - 1.0) / eps_r
            # Legacy Parity: ICoeff are already scaled by j*omega
            factor = kappa * coeff
            
            
            # Vertices of Plus Tet
            v_indices = tetras[:, t_plus]
            r1 = nodes[:, v_indices[1]]
            r2 = nodes[:, v_indices[2]]
            r3 = nodes[:, v_indices[3]]
            r4 = nodes[:, v_indices[4]]
            
            # Volume
            mat = hcat(r2-r1, r3-r1, r4-r1)
            vol = abs(det(mat)) / 6.0
            
            # Free vertex of Plus Tet
            # Assuming local_face_idx points to the vertex opposite to the face
            local_face = bf.local_face_idx[1]
            v_free_idx = v_indices[local_face]
            v_free = nodes[:, v_free_idx]
            
            # Constant part of basis function: A / (3 * V)
            const_bf = bf.area / (3.0 * vol)
            
            # Integration
            integral_val = zero(MVector{3, CT})
            
            for gi in 1:n_points
                u = points[1, gi]
                v = points[2, gi]
                w = points[3, gi]
                x = points[4, gi] # 1 - u - v - w
                
                rgi = u * r1 + v * r2 + w * r3 + x * r4
                
                # f_n = const_bf * (r - v_free)
                # Note: SWG definition for Plus tet is (r - v_free)
                f_val = const_bf * (rgi - v_free)
                
                phase = exp(jk0 * dot(r̂, rgi))
                
                integral_val .+= f_val .* (phase * weights[gi])
            end
            
            # Multiply by volume (Jacobian) and factor
            integral_val .*= vol
            Nxyz .+= factor .* integral_val
        end
        
        # Minus Tetrahedron
        if length(bf.support) > 1
            t_minus = bf.support[2]
            if t_minus > 0
                eps_r = permittivities[t_minus]
                kappa = (eps_r - 1.0) / eps_r
                # Legacy Parity: ICoeff are already scaled by j*omega
                factor = kappa * coeff
                
                # Vertices
                v_indices = tetras[:, t_minus]
                r1 = nodes[:, v_indices[1]]
                r2 = nodes[:, v_indices[2]]
                r3 = nodes[:, v_indices[3]]
                r4 = nodes[:, v_indices[4]]
                
                # Volume
                mat = hcat(r2-r1, r3-r1, r4-r1)
                vol = abs(det(mat)) / 6.0
                
                # Free vertex
                local_face = bf.local_face_idx[2]
                v_free_idx = v_indices[local_face]
                v_free = nodes[:, v_free_idx]
                
                # Constant part: A / (3 * V)
                const_bf = bf.area / (3.0 * vol)
                
                # Integration
                integral_val = zero(MVector{3, CT})
                
                for gi in 1:n_points
                    u = points[1, gi]
                    v = points[2, gi]
                    w = points[3, gi]
                    x = points[4, gi]
                    
                    rgi = u * r1 + v * r2 + w * r3 + x * r4
                    
                    # f_n = const_bf * (v_free - r)
                    f_val = const_bf * (v_free - rgi)
                    
                    phase = exp(jk0 * dot(r̂, rgi))
                    
                    integral_val .+= f_val .* (phase * weights[gi])
                end
                
                integral_val .*= vol
                Nxyz .+= factor .* integral_val
            end
        end
    end
    
    Nθϕ[1] = dot(θhat, Nxyz)
    Nθϕ[2] = dot(ϕhat, Nxyz)
    
    return Nθϕ
end

"""
    radiation_integral_pwc(r̂θϕ, basis, ICoeff, permittivities)

Calculate the radiation integral N(θ, ϕ) using PWC basis functions for VEFIE.

PWC stores 3 DOFs per tetrahedron (x, y, z components). The equivalent current
density is J_eq = jωκD, where D is the unknown. For PWC, D is constant over
each tetrahedron with components [I[3(t-1)+1], I[3(t-1)+2], I[3(t-1)+3]].

N(θ, ϕ) = Σ_t V_t * κ_t * Σ_gq J_t * exp(jk r̂·r_gq) * w_gq

# Legacy Parity
Matches `MoM_Kernels` `raditionalIntegralNθϕCal` for tetrahedra with PWC.
"""
function radiation_integral_pwc(r̂θϕ::r̂θϕInfo{FT}, basis::PWCBasis{IT, FT}, ICoeff::Vector{CT}, permittivities::Vector{CT}) where {IT, FT, CT}
    Nxyz = zero(MVector{3, CT})
    Nθϕ = zero(MVector{2, CT})
    
    r̂ = r̂θϕ.r̂
    θhat = r̂θϕ.θhat
    ϕhat = r̂θϕ.ϕhat
    
    k0 = get_k0()
    jk0 = im * k0
    
    # Quadrature for Tetrahedron
    points, weights = Geometry.gaussQuadratureTet(4, FT)
    n_points = length(weights)
    
    mesh = basis.mesh
    nodes = mesh.node
    tetras_elem = mesh.tetras
    ntet = length(basis.functions)
    
    for t in 1:ntet
        pwc = basis.functions[t]
        vol = pwc.volume
        
        # Material kappa
        eps_r = permittivities[t]
        kappa = (eps_r - 1.0) / eps_r
        
        # Current coefficients for this tetrahedron (3 components)
        Jt = SVector{3, CT}(ICoeff[pwc.inBfsID[1]], ICoeff[pwc.inBfsID[2]], ICoeff[pwc.inBfsID[3]])
        
        # Integrate J * exp(jk r̂·r) over tetrahedron
        v_indices = tetras_elem[:, t]
        r1 = nodes[:, v_indices[1]]
        r2 = nodes[:, v_indices[2]]
        r3 = nodes[:, v_indices[3]]
        r4 = nodes[:, v_indices[4]]
        
        Jtexp = zero(MVector{3, CT})
        
        for gi in 1:n_points
            u = points[1, gi]
            v = points[2, gi]
            w = points[3, gi]
            x = points[4, gi]
            
            rgi = u * r1 + v * r2 + w * r3 + x * r4
            
            phase = exp(jk0 * dot(r̂, rgi))
            Jtexp .+= Jt .* (phase * weights[gi])
        end
        
        Jtexp .*= vol * kappa
        Nxyz .+= Jtexp
    end
    
    Nθϕ[1] = dot(θhat, Nxyz)
    Nθϕ[2] = dot(ϕhat, Nxyz)
    
    return Nθϕ
end

# ============================================================================
# Radiation Integral for PWCHexBasis (hexahedra)
# ============================================================================

"""
    radiation_integral_pwc_hex(r̂θϕ, basis, ICoeff, permittivities)

Calculate the radiation integral N(θ, ϕ) using PWC hexahedra basis functions.

PWC stores 3 DOFs per hexahedron (x, y, z). The equivalent current is:
    J_eq = κ × [I[3(h-1)+1], I[3(h-1)+2], I[3(h-1)+3]]

N(θ, ϕ) = Σ_h V_h × κ_h × Σ_gq J_h × exp(jk r̂·r_gq) × w_gq

# Legacy Parity
Matches `MoM_Kernels` `raditionalIntegralNθϕCal` for hexahedra with PWC.
"""
function radiation_integral_pwc_hex(r̂θϕ::r̂θϕInfo{FT}, basis::PWCHexBasis{IT, FT}, ICoeff::Vector{CT}, permittivities::Vector{CT}) where {IT, FT, CT}
    Nxyz = zero(MVector{3, CT})
    Nθϕ = zero(MVector{2, CT})
    
    r̂ = r̂θϕ.r̂
    θhat = r̂θϕ.θhat
    ϕhat = r̂θϕ.ϕhat
    
    k0 = get_k0()
    jk0 = im * k0
    
    # Hex GQ
    gq_hex = GaussQuadratureInfo(:Hexahedron, 8, FT)
    Nq = length(gq_hex.weight)
    
    mesh = basis.mesh
    nodes = mesh.node
    hexas_elem = mesh.hexas  # 8×nhex
    nhex = length(basis.functions)
    
    for h in 1:nhex
        pwc = basis.functions[h]
        vol = pwc.volume
        
        # Material kappa
        eps_r = permittivities[h]
        kappa = (eps_r - 1.0) / eps_r
        
        # Current coefficients (3 components)
        Jt = SVector{3, CT}(ICoeff[pwc.inBfsID[1]], ICoeff[pwc.inBfsID[2]], ICoeff[pwc.inBfsID[3]])
        
        # Get hex vertices as 3×8
        v_idx = hexas_elem[:, h]
        hex_verts = hcat([nodes[:, v_idx[i]] for i in 1:8]...)
        
        Jtexp = zero(MVector{3, CT})
        
        for gi in 1:Nq
            rgi = hex_verts * gq_hex.coordinate[:, gi]
            phase = exp(jk0 * dot(r̂, rgi))
            Jtexp .+= Jt .* (phase * gq_hex.weight[gi])
        end
        
        Jtexp .*= vol * kappa
        Nxyz .+= Jtexp
    end
    
    Nθϕ[1] = dot(θhat, Nxyz)
    Nθϕ[2] = dot(ϕhat, Nxyz)
    
    return Nθϕ
end

# ============================================================================
# Radiation Integral for RBFBasis (hexahedra)
# ============================================================================

"""
    radiation_integral_rbf(r̂θϕ, basis, ICoeff, permittivities)

Calculate the radiation integral N(θ, ϕ) using RBF basis functions on hexahedra.

RBF has 6 DOFs per hexahedron (one per face). The equivalent current is:
    J(r) = κ × Σ_face I[n] × (A_face / V) × ρ(r)
where ρ(r) = r - r_free(r) with sliding free-end on opposite face.

N(θ, ϕ) = Σ_hex V_h × κ_h × Σ_gq J(r_gq) × exp(jk r̂·r_gq) × w_gq

# Legacy Parity
Matches `MoM_Kernels` `geoElectricJCal` + `raditionalIntegralNθϕCal` for
hexahedra with RBF (`LinearBasisFunction`).
"""
function radiation_integral_rbf(r̂θϕ::r̂θϕInfo{FT}, basis::RBFBasis{IT, FT}, ICoeff::Vector{CT}, permittivities::Vector{CT}) where {IT, FT, CT}
    Nxyz = zero(MVector{3, CT})
    Nθϕ = zero(MVector{2, CT})
    
    r̂ = r̂θϕ.r̂
    θhat = r̂θϕ.θhat
    ϕhat = r̂θϕ.ϕhat
    
    k0 = get_k0()
    jk0 = im * k0
    
    # Hex GQ
    gq_hex = GaussQuadratureInfo(:Hexahedron, 8, FT)
    Nq = length(gq_hex.weight)
    
    mesh = basis.mesh
    nhex = size(mesh.hexas, 2)
    
    # Get hexahedra info
    dummy_perm = permittivities
    hexas = get_hexahedra_info(mesh, basis, dummy_perm)
    
    for (jh, hex) in enumerate(hexas)
        vol = hex.volume
        kappa = hex.κ
        
        # Precompute hex GQ points
        rq_hex = hex.vertices * gq_hex.coordinate
        
        # Precompute free-end points for all 6 faces
        freeVns_all = [get_free_vns(hex, fi, gq_hex.coordinate) for fi in 1:6]
        
        Jtexp = zero(MVector{3, CT})
        
        for gi in 1:Nq
            rgi = @view rq_hex[:, gi]
            
            # Compute J at this GQ point
            J_at_r = zero(MVector{3, CT})
            for mi in 1:6
                n = hex.inBfsID[mi]
                arean = hex.facesArea[mi]
                freeV = @view freeVns_all[mi][:, gi]
                
                # ρ = r - r_free
                ρ_x = rgi[1] - freeV[1]
                ρ_y = rgi[2] - freeV[2]
                ρ_z = rgi[3] - freeV[3]
                
                # J contribution: I[n] × A/V × ρ (Legacy uses just I[n] × A × ρ, then divides by V later)
                coeff = ICoeff[n] * arean
                J_at_r[1] += coeff * ρ_x
                J_at_r[2] += coeff * ρ_y
                J_at_r[3] += coeff * ρ_z
            end
            
            phase = exp(jk0 * dot(r̂, rgi))
            Jtexp .+= J_at_r .* (phase * gq_hex.weight[gi])
        end
        
        # Legacy: Jhexa .*= κ / V, then integral = V × Σ w × Jhexa
        # Combined: integral = κ × Σ w × J_at_r (where J_at_r already has A factors)
        Jtexp .*= kappa
        Nxyz .+= Jtexp
    end
    
    Nθϕ[1] = dot(θhat, Nxyz)
    Nθϕ[2] = dot(ϕhat, Nxyz)
    
    return Nθϕ
end

end

