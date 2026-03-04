module Excitation

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..EFIEModule
using ..MFIEModule
using ..CFIEModule
using ..VEFIEModule
using ..PMCHWModule
using LinearAlgebra
using StaticArrays

using ..CoreModule: AbstractSource, PlaneWave, DeltaGapSource

import ..CoreModule: excitation_vector

export excitation_vector

# Default to EFIE if no operator provided (Backward compatibility)
excitation_vector(source::AbstractSource, basis::RWGBasis) =
    excitation_vector(EFIE(0.0), source, basis)

function excitation_vector(op::EFIE, source::PlaneWave, basis::RWGBasis{IT,FT}) where {IT,FT}
    N = num_basis(basis)
    V = zeros(Complex{FT}, N)

    # Use 3-point quadrature for triangles
    quad = GaussQuadratureInfo(:Triangle, 3, FT)
    num_q = length(quad.weight)

    mesh = basis.mesh
    verts = vertices(mesh)
    elems = elements(mesh)

    for n = 1:N
        bf = basis.functions[n]
        val = zero(Complex{FT})

        # Loop over support triangles
        for k = 1:2
            tri_idx = bf.support[k]
            if tri_idx == 0
                continue
            end

            sign = bf.signs[k]
            local_edge = bf.local_edge_idx[k]

            # Get triangle vertices
            v_indices = elems[:, tri_idx]
            r1 = verts[:, v_indices[1]]
            r2 = verts[:, v_indices[2]]
            r3 = verts[:, v_indices[3]]

            # Triangle area
            cross_prod = cross(r2 - r1, r3 - r1)
            area = 0.5 * norm(cross_prod)

            # Opposite vertex
            v_op = verts[:, v_indices[local_edge]]

            # Integration
            for q = 1:num_q
                # Barycentric coordinates
                u = quad.coordinate[1, q]
                v = quad.coordinate[2, q]
                w = quad.coordinate[3, q]

                # Real coordinate
                r = u * r1 + v * r2 + w * r3

                # Basis function value
                # f(r) = sign * (l / 2A) * (r - v_op)
                rho = r - v_op
                f_val = sign * (bf.edge_length / (2 * area)) * rho

                # Incident field
                E_inc = incident_field(source, r)

                # Dot product
                # Note: dot(a, b) conjugates a. We want E_inc . f_val.
                # Since f_val is real, dot(f_val, E_inc) is correct.
                integrand = dot(f_val, E_inc)

                # Accumulate (Area * weight * integrand)
                val += area * quad.weight[q] * integrand
            end
        end
        V[n] = val
    end
    return V
end

function excitation_vector(op::MFIE, source::PlaneWave, basis::RWGBasis{IT,FT}) where {IT,FT}
    N = num_basis(basis)
    V = zeros(Complex{FT}, N)

    quad = GaussQuadratureInfo(:Triangle, 3, FT)
    num_q = length(quad.weight)

    mesh = basis.mesh
    verts = vertices(mesh)
    elems = elements(mesh)

    # H_inc = (k x E_inc) / eta
    # k_hat = (sin t cos p, sin t sin p, cos t)
    st, ct = sincos(source.theta)
    sp, cp = sincos(source.phi)
    k_hat = SVector{3,FT}(st * cp, st * sp, ct)
    eta = op.eta

    for n = 1:N
        bf = basis.functions[n]
        val = zero(Complex{FT})

        for k = 1:2
            tri_idx = bf.support[k]
            if tri_idx == 0
                continue
            end

            sign = bf.signs[k]
            local_edge = bf.local_edge_idx[k]

            v_indices = elems[:, tri_idx]
            r1 = verts[:, v_indices[1]]
            r2 = verts[:, v_indices[2]]
            r3 = verts[:, v_indices[3]]

            cross_prod = cross(r2 - r1, r3 - r1)
            area = 0.5 * norm(cross_prod)
            normal = normalize(cross_prod)

            v_op = verts[:, v_indices[local_edge]]

            for q = 1:num_q
                u = quad.coordinate[1, q]
                v = quad.coordinate[2, q]
                w = quad.coordinate[3, q]
                r = u * r1 + v * r2 + w * r3

                rho = r - v_op
                f_val = sign * (bf.edge_length / (2 * area)) * rho

                E_inc = incident_field(source, r)
                H_inc = cross(k_hat, E_inc) / eta

                # n x H_inc
                nxH = cross(normal, H_inc)

                integrand = dot(f_val, nxH)
                val += area * quad.weight[q] * integrand
            end
        end
        V[n] = val * eta
    end
    return V
end

function excitation_vector(op::CFIE, source::PlaneWave, basis::RWGBasis)
    V_efie = excitation_vector(op.efie, source, basis)
    V_mfie = excitation_vector(op.mfie, source, basis)

    # V_cfie = alpha * V_efie + (1-alpha) * V_mfie
    # Note: V_mfie is already scaled by eta in excitation_vector(MFIE)
    return op.alpha * V_efie + (1 - op.alpha) * V_mfie
end

"""
    excitation_vector(op::PMCHW, source::PlaneWave, basis::RWGBasis) → Vector{Complex}

计算 PMCHWT 激励向量 V（长度 2N）。

```
V = [V_E]   V_E[m] = ∫ f_m(r) · E^inc(r) dS   （电场 RHS, EFIE 公式）
    [V_H]   V_H[m] = ∫ f_m(r) · H^inc(r) dS   （磁场 RHS, 直接点积）
```

注意：V_H 必须直接积分 f · H^inc，不能使用 MFIE 激励接口。
MFIE 激励计算 η₀ × ∫ f · (n̂ × H^inc) dS，与 PMCHW 的 ∫ f · H^inc dS 不同。
"""
function excitation_vector(op::PMCHW{FT,CT}, source::PlaneWave, basis::RWGBasis{IT,FT}) where {IT,FT,CT}
    N = num_basis(basis)
    V = zeros(CT, 2N)

    # V_E：EFIE 激励 ∫ f_m · E^inc dS
    efie_dummy = EFIE(op.freq)
    V[1:N] .= excitation_vector(efie_dummy, source, basis)

    # V_H：PMCHW 正确公式 ∫ f_m · H^inc dS
    # H^inc = (k̂ × E^inc) / η₀（平面波磁场，直接点积，无 n̂ 叉积，无额外 η₀ 因子）
    # 不能用 MFIE 激励（MFIE 计算 η₀ × ∫ f · (n̂ × H^inc) dS，物理含义不同）
    V[N+1:2N] .= _pmchw_excitation_H(source, basis, op.eta0)

    return V
end

"""
    _pmchw_excitation_H(source, basis, eta0) → Vector{Complex}

计算 PMCHW 磁场激励向量：

    V_H[m] = ∫ f_m(r) · H^inc(r) dS

其中 H^inc(r) = (k̂ × E^inc(r)) / η₀。

与 MFIE 激励的区别：
- MFIE：η₀ × ∫ f · (n̂ × H^inc) dS（含面法向叉积与 η₀ 因子，用于 CFIE）
- PMCHW V_H：∫ f · H^inc dS（直接点积，无 n̂ 叉积，无额外因子）
"""
function _pmchw_excitation_H(source::PlaneWave, basis::RWGBasis{IT,FT}, eta0::FT) where {IT,FT}
    N = num_basis(basis)
    CT = Complex{FT}
    V  = zeros(CT, N)

    quad   = GaussQuadratureInfo(:Triangle, 3, FT)
    num_q  = length(quad.weight)

    mesh  = basis.mesh
    verts = vertices(mesh)
    elems = elements(mesh)

    # 波传播方向单位矢量 k̂
    st, ct_val = sincos(source.theta)
    sp, cp     = sincos(source.phi)
    k_hat = SVector{3,FT}(st * cp, st * sp, ct_val)

    for n = 1:N
        bf  = basis.functions[n]
        val = zero(CT)

        for s = 1:2
            tri_idx = bf.support[s]
            if tri_idx == 0
                continue
            end

            sign       = bf.signs[s]
            local_edge = bf.local_edge_idx[s]

            v_indices = elems[:, tri_idx]
            r1 = verts[:, v_indices[1]]
            r2 = verts[:, v_indices[2]]
            r3 = verts[:, v_indices[3]]

            cross_prod = cross(r2 - r1, r3 - r1)
            area       = FT(0.5) * norm(cross_prod)

            v_op = verts[:, v_indices[local_edge]]

            for q = 1:num_q
                u = quad.coordinate[1, q]
                v = quad.coordinate[2, q]
                w = quad.coordinate[3, q]
                r = u * r1 + v * r2 + w * r3

                rho   = r - v_op
                f_val = sign * (bf.edge_length / (2 * area)) * rho

                # H^inc = (k̂ × E^inc) / η₀
                E_inc = incident_field(source, r)
                H_inc = cross(k_hat, E_inc) / eta0

                # 直接点积 f · H^inc（f_val 为实数矢量）
                integrand = dot(f_val, H_inc)
                val += area * quad.weight[q] * integrand
            end
        end
        V[n] = val
    end
    return V
end

# Delta Gap (Same for all operators usually, as it forces V)
excitation_vector(op::AbstractIntegralOperator, source::DeltaGapSource, basis::RWGBasis) =
    excitation_vector(source, basis)

function excitation_vector(source::DeltaGapSource, basis::RWGBasis{IT,FT}) where {IT,FT}
    N = num_basis(basis)
    V = zeros(Complex{FT}, N)

    # Delta gap source is applied directly to the edge
    # V_n = V_source * L_n

    for idx in source.edge_indices
        if 1 <= idx <= N
            # Access edge length from the basis function struct
            edge_len = basis.functions[idx].edge_length
            V[idx] = source.voltage * edge_len
        else
            @warn "Delta gap source edge index $idx out of bounds (1:$N)"
        end
    end

    return V
end

function excitation_vector(
    op::VEFIE,
    source::PlaneWave,
    basis::SWGBasis{IT,FT},
    permittivities::Vector{ComplexF64},
) where {IT,FT}
    N = num_basis(basis)
    V = zeros(Complex{FT}, N)

    quad = op.gq_info
    num_q = length(quad.weight)

    mesh = basis.mesh
    verts = vertices(mesh)
    elems = elements(mesh)

    for n = 1:N
        bf = basis.functions[n]
        val = zero(Complex{FT})

        for k = 1:2
            tet_idx = bf.support[k]
            if tet_idx == 0
                continue
            end

            sign = bf.signs[k]
            local_face = bf.local_face_idx[k]

            v_indices = elems[:, tet_idx]
            v1 = verts[:, v_indices[1]]
            v2 = verts[:, v_indices[2]]
            v3 = verts[:, v_indices[3]]
            v4 = verts[:, v_indices[4]]

            # Volume
            vol = abs(dot(v2 - v1, cross(v3 - v1, v4 - v1))) / 6.0

            # Opposite vertex (assuming Face i is opposite to Vertex i)
            v_op = verts[:, v_indices[local_face]]

            # Factor A / 3V
            # bf.area is A.
            factor = bf.area / (3.0 * vol)

            for q = 1:num_q
                u = quad.coordinate[1, q]
                v = quad.coordinate[2, q]
                w = quad.coordinate[3, q]
                t = quad.coordinate[4, q]

                r = u * v1 + v * v2 + w * v3 + t * v4

                # f(r) = sign * factor * (r - v_op)
                f_val = sign * factor * (r - v_op)

                E_inc = incident_field(source, r)

                integrand = dot(f_val, E_inc)

                val += vol * quad.weight[q] * integrand
            end
        end
        V[n] = val
    end
    return V
end

function excitation_vector(source::PlaneWave, basis::SWGBasis{IT,FT}) where {IT,FT}
    N = num_basis(basis)
    V = zeros(Complex{FT}, N)

    # Use 4-point quadrature for tetrahedra
    quad = GaussQuadratureInfo(:Tetrahedron, 4, FT)
    num_q = length(quad.weight)

    mesh = basis.mesh
    verts = vertices(mesh)
    elems = elements(mesh)

    for n = 1:N
        bf = basis.functions[n]
        val = zero(Complex{FT})

        for k = 1:2
            tet_idx = bf.support[k]
            if tet_idx == 0
                continue
            end

            sign = bf.signs[k]
            local_face = bf.local_face_idx[k]

            v_indices = elems[:, tet_idx]
            v1 = verts[:, v_indices[1]]
            v2 = verts[:, v_indices[2]]
            v3 = verts[:, v_indices[3]]
            v4 = verts[:, v_indices[4]]

            # Volume
            vol = abs(dot(v2 - v1, cross(v3 - v1, v4 - v1))) / 6.0

            # Opposite vertex (assuming Face i is opposite to Vertex i)
            v_op = verts[:, v_indices[local_face]]

            # Factor A / 3V
            # bf.area is A.
            factor = bf.area / (3.0 * vol)

            for q = 1:num_q
                u = quad.coordinate[1, q]
                v = quad.coordinate[2, q]
                w = quad.coordinate[3, q]
                t = quad.coordinate[4, q]

                r = u * v1 + v * v2 + w * v3 + t * v4

                # f(r) = sign * factor * (r - v_op)
                f_val = sign * factor * (r - v_op)

                E_inc = incident_field(source, r)

                integrand = dot(f_val, E_inc)

                val += vol * quad.weight[q] * integrand
            end
        end
        V[n] = val
    end
    return V
end

function excitation_vector(source::AbstractSource, surface_basis::RWGBasis, volume_basis::SWGBasis)
    V_surf = excitation_vector(source, surface_basis)
    V_vol = excitation_vector(source, volume_basis)
    return [V_surf; V_vol]
end

"""
    excitation_vector(op::VEFIE, source::PlaneWave, basis::PWCBasis)

Compute the VEFIE excitation vector for PWC basis functions.

For each tetrahedron t, the 3 excitation components (x, y, z) are:
    V[3(t-1)+i] = V_t * Σ_gq w_gq * E_inc(r_gq)[i]

# Legacy Parity
Matches `MoM_Kernels` `excitationVectorEFIE` for `ConstBasisFunction` (tetrahedra).
"""
function excitation_vector(op::VEFIE, source::PlaneWave, basis::PWCBasis{IT,FT}) where {IT,FT}
    CT = Complex{FT}
    N = num_basis(basis)
    V = zeros(CT, N)

    quad = op.gq_info
    num_q = length(quad.weight)

    mesh = basis.mesh
    verts = vertices(mesh)
    elems = elements(mesh)
    ntet = length(basis.functions)

    for t = 1:ntet
        pwc = basis.functions[t]

        # Get tetra vertices
        v_indices = elems[:, t]
        v1 = verts[:, v_indices[1]]
        v2 = verts[:, v_indices[2]]
        v3 = verts[:, v_indices[3]]
        v4 = verts[:, v_indices[4]]

        vol = pwc.volume

        # Integrate E_inc over tetrahedron
        E_sum = zero(MVector{3,CT})
        for q = 1:num_q
            u = quad.coordinate[1, q]
            v = quad.coordinate[2, q]
            w = quad.coordinate[3, q]
            s = quad.coordinate[4, q]

            r = u * v1 + v * v2 + w * v3 + s * v4

            E_inc = incident_field(source, r)
            E_sum .+= E_inc .* quad.weight[q]
        end
        E_sum .*= vol

        # Scatter to global indices
        for i = 1:3
            V[pwc.inBfsID[i]] = E_sum[i]
        end
    end

    return V
end

"""
    excitation_vector(source::PlaneWave, basis::PWCBasis)

Compute the VEFIE excitation vector for PWC basis functions (without VEFIE operator).
Uses 4-point GQ by default.
"""
function excitation_vector(source::PlaneWave, basis::PWCBasis{IT,FT}) where {IT,FT}
    CT = Complex{FT}
    N = num_basis(basis)
    V = zeros(CT, N)

    quad = GaussQuadratureInfo(:Tetrahedron, 4, FT)
    num_q = length(quad.weight)

    mesh = basis.mesh
    verts = vertices(mesh)
    elems = elements(mesh)
    ntet = length(basis.functions)

    for t = 1:ntet
        pwc = basis.functions[t]

        v_indices = elems[:, t]
        v1 = verts[:, v_indices[1]]
        v2 = verts[:, v_indices[2]]
        v3 = verts[:, v_indices[3]]
        v4 = verts[:, v_indices[4]]

        vol = pwc.volume

        E_sum = zero(MVector{3,CT})
        for q = 1:num_q
            u = quad.coordinate[1, q]
            v = quad.coordinate[2, q]
            w = quad.coordinate[3, q]
            s = quad.coordinate[4, q]

            r = u * v1 + v * v2 + w * v3 + s * v4
            E_inc = incident_field(source, r)
            E_sum .+= E_inc .* quad.weight[q]
        end
        E_sum .*= vol

        for i = 1:3
            V[pwc.inBfsID[i]] = E_sum[i]
        end
    end

    return V
end

"""
    excitation_vector(source, surface_basis::RWGBasis, volume_basis::PWCBasis)

Compute combined excitation vector for SCFIE (RWG + PWC).
"""
function excitation_vector(source::AbstractSource, surface_basis::RWGBasis, volume_basis::PWCBasis)
    V_surf = excitation_vector(source, surface_basis)
    V_vol = excitation_vector(source, volume_basis)
    return [V_surf; V_vol]
end

# ============================================================================
# Hex Volume Excitation: PWCHexBasis
# ============================================================================

"""
    excitation_vector(op::VEFIE, source::PlaneWave, basis::PWCHexBasis)

Compute the VEFIE excitation vector for PWC hexahedra basis functions.

For each hexahedron h, the 3 excitation components (x, y, z) are:
    V[3(h-1)+i] = V_h × Σ_gq w_gq × E_inc(r_gq)[i]

# Legacy Parity
Matches `MoM_Kernels` `excitationVectorEFIE` for `ConstBasisFunction` (hexahedra).
"""
function excitation_vector(source::PlaneWave, basis::PWCHexBasis{IT,FT}) where {IT,FT}
    CT = Complex{FT}
    N = num_basis(basis)
    V = zeros(CT, N)

    # Use hex GQ
    quad = GaussQuadratureInfo(:Hexahedron, 8, FT)
    num_q = length(quad.weight)

    mesh = basis.mesh
    verts = vertices(mesh)
    elems = mesh.hexes  # 8×nhex
    nhex = length(basis.functions)

    for h = 1:nhex
        pwc = basis.functions[h]
        vol = pwc.volume

        # Get hex vertices as 3×8 matrix
        v_idx = elems[:, h]
        hex_verts = hcat([verts[:, v_idx[i]] for i = 1:8]...)

        # Integrate E_inc over hexahedron
        E_sum = zero(MVector{3,CT})
        for q = 1:num_q
            # Hex shape function: r = Σ_i N_i(u,v,w) × v_i 
            # quad.coordinate is 8×Nq (8-node shape function values)
            r = hex_verts * quad.coordinate[:, q]
            E_inc = incident_field(source, r)
            E_sum .+= E_inc .* quad.weight[q]
        end
        E_sum .*= vol

        # Scatter to global indices
        for i = 1:3
            V[pwc.inBfsID[i]] = E_sum[i]
        end
    end

    return V
end

"""
    excitation_vector(source, surface_basis::RWGBasis, volume_basis::PWCHexBasis)

Compute combined excitation vector for SCFIE (RWG + PWCHex).
"""
function excitation_vector(
    source::AbstractSource,
    surface_basis::RWGBasis,
    volume_basis::PWCHexBasis,
)
    V_surf = excitation_vector(source, surface_basis)
    V_vol = excitation_vector(source, volume_basis)
    return [V_surf; V_vol]
end

# ============================================================================
# Hex Volume Excitation: RBFBasis
# ============================================================================

"""
    excitation_vector(op::VEFIE, source::PlaneWave, basis::RBFBasis)

Compute the VEFIE excitation vector for RBF (rooftop) basis functions.

For each face m on hexahedron h:
    V_m = A_m × Σ_gq w_gq × ρ_m(r_gq) · E_inc(r_gq)

where ρ_m(r) = r - r_free(gq) and r_free is the corresponding point on the
opposite face (varies with quadrature point, unlike SWG's fixed opposite vertex).

# Legacy Parity
Matches `MoM_Kernels` `excitationVectorEFIE` for `LinearBasisFunction` (hexahedra).
"""
function excitation_vector(source::PlaneWave, basis::RBFBasis{IT,FT}) where {IT,FT}
    CT = Complex{FT}
    N = num_basis(basis)
    V = zeros(CT, N)

    # Use hex GQ for volume integration
    quad = GaussQuadratureInfo(:Hexahedron, 8, FT)
    quad_face = GaussQuadratureInfo(:Quadrangle, 4, FT)
    num_q = length(quad.weight)

    # Build 3D -> 2D GQ index map for free-end lookup (n1d=2 for 8-pt hex)
    gq3d_map = construct_gq3d_index_map(2)

    mesh = basis.mesh
    elems = mesh.hexes  # 8×nhex
    nhex = size(elems, 2)

    # Get hexahedra info (needed for face areas and free-end computation)
    # We create dummy permittivities for geometry-only purposes
    dummy_perm = fill(Complex{FT}(1.0), nhex)
    hexas = get_hexahedra_info(mesh, basis, dummy_perm)

    for (jh, hex) in enumerate(hexas)
        # Precompute hex GQ points: 3×Nq
        rq_hex = hex.vertices * quad.coordinate

        # For each RBF face
        for mi = 1:6
            n = hex.inBfsID[mi]

            arean = hex.facesArea[mi]
            # Free-end coords on opposite face (3×Nq_quad, using face GQ)
            freeVns = get_free_vns(hex, mi, quad_face.coordinate)

            val = zero(CT)
            @inbounds for q = 1:num_q
                rq = @view rq_hex[:, q]
                # Map 3D hex GQ index to 2D face index for free-end lookup
                id_face = gq3d_to_face2d_idx(gq3d_map[q], mi, 2)
                freeV = @view freeVns[:, id_face]

                # ρ_m = r - r_free
                ρ_x = rq[1] - freeV[1]
                ρ_y = rq[2] - freeV[2]
                ρ_z = rq[3] - freeV[3]

                E_inc = incident_field(source, rq)

                # ρ · E_inc
                val += (ρ_x * E_inc[1] + ρ_y * E_inc[2] + ρ_z * E_inc[3]) * quad.weight[q]
            end
            val *= arean

            V[n] += val
        end
    end

    return V
end

"""
    excitation_vector(source, surface_basis::RWGBasis, volume_basis::RBFBasis)

Compute combined excitation vector for SCFIE (RWG + RBF).
"""
function excitation_vector(source::AbstractSource, surface_basis::RWGBasis, volume_basis::RBFBasis)
    V_surf = excitation_vector(source, surface_basis)
    V_vol = excitation_vector(source, volume_basis)
    return [V_surf; V_vol]
end

end
