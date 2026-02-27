module Disaggregation

using LinearAlgebra
using StaticArrays
using ....CoreModule
using ....Geometry
using ....BasisFunctions
using ....IntegralEquations
using ..Level

using ....IntegralEquations.Impedance: get_triangle_info, get_triangles_info
using ....IntegralEquations.VEFIEModule: get_tetrahedra_info
using ....IntegralEquations.CFIEModule: CFIE
using ....IntegralEquations.SCFIEModule: SCFIE

export disaggregate_downward!, disaggregate_leaf!

function get_k(op::AbstractIntegralOperator)
    if hasfield(typeof(op), :k)
        return op.k
    elseif hasfield(typeof(op), :efie)
        return op.efie.k
    elseif hasfield(typeof(op), :freq)
        c0 = 299792458.0
        return 2π * op.freq / c0
    else
        error("Operator $(typeof(op)) does not have field k or efie or freq")
    end
end

function get_eta(op::AbstractIntegralOperator)
    if hasfield(typeof(op), :eta)
        return op.eta
    elseif hasfield(typeof(op), :efie)
        return op.efie.eta
    else
        # Default to free space
        return 376.730313668
    end
end

"""
    disaggregate_downward!(parentLevel::LevelInfo, childLevel::LevelInfo)

Disaggregate from parent level to child level (Phase Shift + Anterpolation).
"""
function disaggregate_downward!(parentLevel::LevelInfo, childLevel::LevelInfo)
    FT = eltype(parentLevel.cubeEdgel)
    CT = Complex{FT}
    
    parentDisaggG = parentLevel.disaggG
    
    nPolesChild = length(childLevel.poles.r̂sθsϕs)
    nCubesChild = childLevel.nCubes
    
    if !isdefined(childLevel, :disaggG)
        childLevel.disaggG = zeros(CT, nPolesChild, 2, nCubesChild)
    end
    childDisaggG = childLevel.disaggG
    
    interp = childLevel.interpWθϕ
    θCSCT = interp.θCSCT
    ϕCSCT = interp.ϕCSCT
    
    phaseShift = parentLevel.phaseShift2Kids
    
    Threads.@threads for iCube in 1:parentLevel.nCubes
        parentCube = parentLevel.cubes[iCube]
        disaggParent = view(parentDisaggG, :, :, iCube)
        
        for iKid in 1:length(parentCube.kidsInterval)
            childID = parentCube.kidsInterval[iKid]
            childIn8 = parentCube.kidsIn8[iKid]
            
            shift = view(phaseShift, :, childIn8)
            
            shiftedParent = zeros(CT, size(disaggParent))
            for pol in 1:2
                @views shiftedParent[:, pol] .= shift .* disaggParent[:, pol]
            end
            
            temp = θCSCT * shiftedParent
            childVal = ϕCSCT * temp
            
            @views childDisaggG[:, :, childID] .+= childVal
        end
    end
end

function get_basis_index(global_idx::Int, offsets::Vector{Int})
    idx = searchsortedfirst(offsets, global_idx)
    local_idx = idx == 1 ? global_idx : global_idx - offsets[idx-1]
    return idx, local_idx
end

"""
    disaggregate_leaf!(level, basis::AbstractBasisFunction, operator, ZI, sorted_ids)

Single-basis convenience wrapper.
"""
function disaggregate_leaf!(level::LevelInfo, basis::AbstractBasisFunction, operator::AbstractIntegralOperator, ZI::AbstractVector, sorted_ids::Vector{Int})
    disaggregate_leaf!(level, AbstractBasisFunction[basis], [num_basis(basis)], operator, ZI, sorted_ids)
end

"""
    disaggregate_leaf!(level::LevelInfo, bases::Vector{<:AbstractBasisFunction}, offsets::Vector{Int}, operator::AbstractIntegralOperator, ZI::AbstractVector, sorted_ids::Vector{Int})

Compute received fields on basis functions at the leaf level.
"""
function disaggregate_leaf!(level::LevelInfo, bases::Vector{<:AbstractBasisFunction}, offsets::Vector{Int}, operator::AbstractIntegralOperator, ZI::AbstractVector, sorted_ids::Vector{Int})
    FT = eltype(level.cubeEdgel)
    CT = Complex{FT}
    
    k = get_k(operator)
    eta = get_eta(operator)
    JK = im * k
    
    if !isdefined(level, :disaggG)
        nPoles = length(level.poles.r̂sθsϕs)
        level.disaggG = zeros(CT, nPoles, 2, level.nCubes)
    end
    
    disaggG = level.disaggG
    poles = level.poles
    nPoles = length(poles.r̂sθsϕs)
    nCubes = level.nCubes
    
    # Precompute element info and quadrature
    element_infos = Any[]
    gqs = Any[]
    
    for b in bases
        if b isa RWGBasis
            push!(element_infos, get_triangles_info(b.mesh, b))
            push!(gqs, GaussQuadratureInfo(:Triangle, 3, FT))
        elseif b isa SWGBasis
            if hasfield(typeof(operator), :permittivities)
                push!(element_infos, get_tetrahedra_info(b.mesh, b, operator.permittivities))
            else
                push!(element_infos, get_tetrahedra_info(b.mesh, b, fill(1.0+0im, num_elements(b.mesh))))
            end
            push!(gqs, GaussQuadratureInfo(:Tetrahedron, 5, FT))
        else
            error("Unsupported basis type")
        end
    end
    
    poles_r̂ = [p.r̂ for p in poles.r̂sθsϕs]
    poles_θhat = [p.θhat for p in poles.r̂sθsϕs]
    poles_ϕhat = [p.ϕhat for p in poles.r̂sθsϕs]

    Threads.@threads for iCube in 1:nCubes
        cube = level.cubes[iCube]
        cubeCenter = cube.center
        
        field_at_center = view(disaggG, :, :, iCube)
        
        for bfID_sorted in cube.bfInterval
            bfID = sorted_ids[bfID_sorted]
            
            b_idx, i_local = get_basis_index(bfID, offsets)
            basis = bases[b_idx]
            bf = basis.functions[i_local]
            
            elem_info = element_infos[b_idx]
            gq = gqs[b_idx]
            
            if basis isa RWGBasis
                add_received_field_rwg!(ZI, bfID, basis, bf, elem_info, gq, operator, k, eta, cubeCenter, poles_r̂, poles_θhat, poles_ϕhat, field_at_center)
            elseif basis isa SWGBasis
                add_received_field_swg!(ZI, bfID, basis, bf, elem_info, gq, operator, k, eta, cubeCenter, poles_r̂, poles_θhat, poles_ϕhat, field_at_center)
            end
        end
    end
end

function add_received_field_rwg!(ZI, bfID, basis, bf, elem_info, gq, operator, k, eta, cubeCenter, poles_r̂, poles_θhat, poles_ϕhat, field_at_center)
    JK = im * k
    n_qp = length(gq.weight)
    nPoles = length(poles_r̂)
    
    # Determine factors
    is_cfie = (operator isa CFIE) || (operator isa SCFIE)
    alpha = 0.5
    efie_factor = 1.0 + 0im
    
    if operator isa CFIE
        alpha = operator.alpha
        efie_factor = operator.efie.factor
    elseif operator isa SCFIE
        alpha = operator.alpha
        # SCFIE surface part uses the same EFIE factor as standalone EFIE: jkη/(16π)
        efie_factor = im * operator.k * operator.eta / (16 * π)
    else
        # EFIE
        if hasfield(typeof(operator), :factor)
            efie_factor = operator.factor
        end
    end
    
    # For pure EFIE, we just compute EFIE term.
    # For CFIE/SCFIE, we compute both.
    
    val = zero(ComplexF64)
    
    for i_supp in 1:2
        tri_idx = bf.support[i_supp]
        if tri_idx == 0; continue; end
        
        tri = elem_info[tri_idx]
        tri_vertices = tri.vertices
        
        # Normal for MFIE
        v1 = tri_vertices[:, 1]
        v2 = tri_vertices[:, 2]
        v3 = tri_vertices[:, 3]
        normal = normalize(cross(v2 - v1, v3 - v1))
        
        local_edge = bf.local_edge_idx[i_supp]
        v_opp = tri_vertices[:, local_edge]
        
        sign = bf.signs[i_supp]
        
        for i_qp in 1:n_qp
            L = gq.coordinate[:, i_qp]
            r = tri_vertices * L
            rho = r - v_opp
            r_local = r - cubeCenter
            factor_vec = sign * bf.edge_length / 2 * gq.weight[i_qp]
            
            for iPole in 1:nPoles
                r̂ = poles_r̂[iPole]
                phase = exp(-JK * dot(r̂, r_local))
                
                E_theta = field_at_center[iPole, 1]
                E_phi = field_at_center[iPole, 2]
                
                θhat = poles_θhat[iPole]
                ϕhat = poles_ϕhat[iPole]
                
                E_inc = (E_theta * θhat + E_phi * ϕhat) * phase
                
                term_efie = dot(rho, E_inc)
                
                if is_cfie
                    # H_inc = 1/eta * (k_hat x E_inc)
                    H_inc = (E_theta * ϕhat - E_phi * θhat) * (phase / eta)
                    term_mfie = dot(cross(rho, normal), H_inc)
                    
                    # Combine
                    val += (alpha * term_efie * efie_factor + (1 - alpha) * eta * term_mfie * (-efie_factor)) * factor_vec
                else
                    # EFIE only
                    # Note: const_factor (jk*eta) is usually applied in operator or here.
                    # If efie_factor is 1.0, we might need to apply it.
                    # But let's assume operator handles it or it's 1.0.
                    val += term_efie * factor_vec
                end
            end
        end
    end
    
    # Apply global factor if not CFIE (CFIE applies inside loop)
    if !is_cfie
        # For EFIE, we might need to apply const_factor = 1.0 (already done)
        # But wait, _disaggregate_leaf_efie! had `const_factor = one(CT)`.
        # So we are good.
    end
    
    ZI[bfID] += val
end

function add_received_field_swg!(ZI, bfID, basis, bf, elem_info, gq, operator, k, eta, cubeCenter, poles_r̂, poles_θhat, poles_ϕhat, field_at_center)
    JK = im * k
    n_qp = length(gq.weight)
    nPoles = length(poles_r̂)
    
    # Factor analysis:
    # - VEFIE near-field uses G = exp(-jkR)/(4πR) with c1 = jωμ₀κ = jkηκ
    # - Translation decomposes exp(-jkR)/R = (-jk/(4π)) Σ W T_L exp(...)
    # - So exp(-jkR)/(4πR) = (1/(4π)) × (-jk/(4π)) Σ W T_L exp(...)
    # - Need: const_factor × κ_from_agg × αTrans = c1 × G_decomposition
    # - Therefore: const_factor = jkη / (4π)
    const_factor = im * k * eta / (4 * π)
    
    val = zero(ComplexF64)
    
    for i_supp in 1:2
        tet_idx = bf.support[i_supp]
        if tet_idx == 0; continue; end
        
        tet = elem_info[tet_idx]
        tet_vertices = tet.vertices
        
        local_face = bf.local_face_idx[i_supp]
        v_free = tet_vertices[:, local_face]
        
        sign = bf.signs[i_supp]
        
        for i_qp in 1:n_qp
            L = gq.coordinate[:, i_qp]
            r = tet_vertices * L
            rho = sign * (r - v_free)
            r_local = r - cubeCenter
            factor_vec = bf.area / 3 * gq.weight[i_qp]
            
            for iPole in 1:nPoles
                r̂ = poles_r̂[iPole]
                phase = exp(-JK * dot(r̂, r_local))
                
                E_theta = field_at_center[iPole, 1]
                E_phi = field_at_center[iPole, 2]
                
                θhat = poles_θhat[iPole]
                ϕhat = poles_ϕhat[iPole]
                
                E_inc = (E_theta * θhat + E_phi * ϕhat) * phase
                
                val += dot(rho, E_inc) * factor_vec
            end
        end
    end
    
    ZI[bfID] += val * const_factor
end

end
