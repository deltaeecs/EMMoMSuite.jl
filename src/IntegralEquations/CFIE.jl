module CFIEModule

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..EFIEModule
using ..MFIEModule
using ..Impedance
using ..Kernels
import ..EFIEModule: calc_self_interaction!, calc_near_interaction!, is_adjacent
import ..MFIEModule: calc_self_term!, calc_k_term_fast!
using StaticArrays
using LinearAlgebra
using Base.Threads

import ..CoreModule: assemble_impedance_matrix

export CFIE, assemble_impedance_matrix

"""
    CFIE{FT, CT} <: AbstractIntegralOperator

Combined Field Integral Equation (CFIE) operator.

A linear combination of EFIE and MFIE used to eliminate internal resonance problems that plague EFIE and MFIE for closed structures at specific frequencies.

# Mathematical Formulation

The CFIE is defined as:
```math
\\text{CFIE} = \\alpha \\cdot \\text{EFIE} + (1-\\alpha) \\eta \\cdot \\text{MFIE}
```

The impedance matrix is constructed as:
```math
\\mathbf{Z}_{CFIE} = \\alpha \\mathbf{Z}_{EFIE} + (1-\\alpha) \\eta \\mathbf{Z}_{MFIE}
```

where:
- \$\\alpha\$ is the weighting factor (typically 0.5).
- \$\\eta\$ is the intrinsic impedance of the medium (used to balance units).

# Fields
- `freq`: Operating frequency.
- `alpha`: Weighting factor \$\\alpha\$ (0 to 1).
- `efie`: Underlying `EFIE` operator.
- `mfie`: Underlying `MFIE` operator.
"""
struct CFIE{FT<:AbstractFloat, CT<:Complex} <: AbstractIntegralOperator
    freq::FT
    alpha::FT
    efie::EFIE{FT, CT}
    mfie::MFIE{FT, CT}
end

function CFIE(freq::FT, alpha::FT = 0.5) where {FT}
    efie = EFIE(freq)
    mfie = MFIE(freq)
    return CFIE{FT, Complex{FT}}(freq, alpha, efie, mfie)
end

"""
    assemble_impedance_matrix(cfie::CFIE, basis::RWGBasis)

Assemble the impedance matrix Z for the CFIE using merged single-pass assembly.

## Optimisation: Merged EFIE+MFIE Single-Pass Assembly

Instead of assembling Z_efie and Z_mfie separately (2 full traversals + matrix add),
this performs a single symmetric traversal over triangle pairs, computing both EFIE and MFIE
contributions simultaneously with shared precomputation:

- **Far pairs (>99%)**: Shared Green's function evaluation, EFIE symmetric + MFIE both directions
- **Near pairs**: EFIE 7pt GQ + MFIE 4pt GQ (separate but same traversal)
- **Self pairs**: EFIE analytical + MFIE mass matrix

Saves ~30-40% vs separate assembly by eliminating:
1. Duplicate Green's function evaluations (~2× saved)
2. Duplicate TriangleInfo precomputation
3. N×N matrix addition and extra allocation
"""
function assemble_impedance_matrix(cfie::CFIE{FT, CT}, basis::RWGBasis{IT, FT}) where {IT, FT, CT}
    N = num_basis(basis)
    Z = zeros(CT, N, N)
    
    mesh = basis.mesh
    nt = num_elements(mesh)
    
    efie = cfie.efie
    mfie = cfie.mfie
    alpha = FT(cfie.alpha)
    mfie_factor = FT(1.0) - alpha
    
    # Precompute TriangleInfo (shared between EFIE and MFIE — done once instead of twice)
    tris_info = Vector{TriangleInfo{IT, FT}}(undef, nt)
    Threads.@threads for t in 1:nt
        tris_info[t] = get_triangle_info(mesh, basis, t)
    end
    
    # Precompute 4-point quad points (shared between EFIE far-field and MFIE)
    gq_far = efie.gq_far
    N_far = length(gq_far.weight)
    quad_points = Vector{SVector{N_far, SVector{3, FT}}}(undef, nt)
    Threads.@threads for t in 1:nt
        v_indices = mesh.triangles[:, t]
        v1 = SVector{3, FT}(mesh.node[:, v_indices[1]])
        v2 = SVector{3, FT}(mesh.node[:, v_indices[2]])
        v3 = SVector{3, FT}(mesh.node[:, v_indices[3]])
        quad_points[t] = SVector{N_far, SVector{3, FT}}(
            v1 * gq_far.coordinate[1, i] + v2 * gq_far.coordinate[2, i] + v3 * gq_far.coordinate[3, i]
            for i in 1:N_far
        )
    end
    
    n_threads = Threads.nthreads()
    row_locks = [SpinLock() for _ in 1:N]
    
    # Precomputed constants
    C4divk2 = efie.C4divk2
    k = efie.k
    efie_factor_val = efie.factor  # jkη/16π
    JK_0 = im * mfie.k
    eta_div_16pi = mfie.eta / (16 * FT(π))
    
    Threads.@threads for tid in 1:n_threads
        # Thread-local buffers
        Z_efie_local = zeros(CT, 3, 3)
        Z_mfie_fwd = zeros(CT, 3, 3)
        Z_mfie_bwd = zeros(CT, 3, 3)
        Z_combined_fwd = zeros(CT, 3, 3)
        Z_combined_bwd = zeros(CT, 3, 3)
        
        for t_test in tid:n_threads:nt
            tri_test = tris_info[t_test]
            
            for t_source in t_test:nt  # symmetric traversal
                tri_source = tris_info[t_source]
                
                fill!(Z_efie_local, zero(CT))
                fill!(Z_mfie_fwd, zero(CT))
                fill!(Z_mfie_bwd, zero(CT))
                
                if t_test == t_source
                    # === Self-pair ===
                    # EFIE: analytical singular integration
                    calc_self_interaction!(Z_efie_local, efie, tri_test)
                    Z_efie_local .*= efie_factor_val
                    
                    # MFIE: mass matrix (symmetric → fwd only)
                    calc_self_term!(Z_mfie_fwd, mfie, tri_test)
                    
                    # Combined forward (self-pair, no backward needed)
                    @inbounds for n in 1:3, m in 1:3
                        Z_combined_fwd[m, n] = alpha * Z_efie_local[m, n] + mfie_factor * Z_mfie_fwd[m, n]
                    end
                    
                    # Write self-pair to Z
                    @inbounds for i in 1:3
                        row_idx = tri_test.inBfsID[i]
                        row_idx == 0 && continue
                        sign_test = tri_test.bfsSign[i]
                        lock(row_locks[row_idx])
                        for j in 1:3
                            col_idx = tri_test.inBfsID[j]
                            if col_idx != 0
                                Z[row_idx, col_idx] += Z_combined_fwd[i, j] * sign_test * tri_test.bfsSign[j]
                            end
                        end
                        unlock(row_locks[row_idx])
                    end
                    
                elseif is_adjacent(tri_test, tri_source)
                    # === Near-pair (adjacent triangles) ===
                    # EFIE: 7-point singular integration (asymmetric — enforce ordering)
                    if tri_test.triID < tri_source.triID
                        calc_near_interaction!(Z_efie_local, efie, tri_test, tri_source)
                    else
                        Z_temp = zeros(CT, 3, 3)
                        calc_near_interaction!(Z_temp, efie, tri_source, tri_test)
                        for m in 1:3, n in 1:3
                            Z_efie_local[m, n] += Z_temp[n, m]
                        end
                    end
                    Z_efie_local .*= efie_factor_val
                    
                    # MFIE: K-term (4-point, both directions)
                    r_test_pts = quad_points[tri_test.triID]
                    r_src_pts = quad_points[tri_source.triID]
                    calc_k_term_fast!(Z_mfie_fwd, mfie, tri_test, tri_source, r_test_pts, r_src_pts)
                    calc_k_term_fast!(Z_mfie_bwd, mfie, tri_source, tri_test, r_src_pts, r_test_pts)
                    
                    # Combined forward: Z[test_bf, source_bf]
                    @inbounds for n in 1:3, m in 1:3
                        Z_combined_fwd[m, n] = alpha * Z_efie_local[m, n] + mfie_factor * Z_mfie_fwd[m, n]
                    end
                    
                    # Combined backward: Z[source_bf, test_bf]
                    # EFIE symmetric: Z_efie_backward[m,n] = Z_efie_local[n,m]
                    @inbounds for n in 1:3, m in 1:3
                        Z_combined_bwd[m, n] = alpha * Z_efie_local[n, m] + mfie_factor * Z_mfie_bwd[m, n]
                    end
                    
                    # Write forward
                    _write_to_Z!(Z, row_locks, Z_combined_fwd, tri_test, tri_source)
                    # Write backward
                    _write_to_Z!(Z, row_locks, Z_combined_bwd, tri_source, tri_test)
                    
                else
                    # === Far-pair (merged kernel with shared Green's function) ===
                    _cfie_far_merged!(Z_efie_local, Z_mfie_fwd, Z_mfie_bwd,
                                      efie, mfie, tri_test, tri_source,
                                      quad_points, C4divk2, k, efie_factor_val,
                                      JK_0, eta_div_16pi)
                    
                    # Combined forward
                    @inbounds for n in 1:3, m in 1:3
                        Z_combined_fwd[m, n] = alpha * Z_efie_local[m, n] + mfie_factor * Z_mfie_fwd[m, n]
                    end
                    
                    # Combined backward (EFIE symmetric → transpose)
                    @inbounds for n in 1:3, m in 1:3
                        Z_combined_bwd[m, n] = alpha * Z_efie_local[n, m] + mfie_factor * Z_mfie_bwd[m, n]
                    end
                    
                    # Write forward and backward
                    _write_to_Z!(Z, row_locks, Z_combined_fwd, tri_test, tri_source)
                    _write_to_Z!(Z, row_locks, Z_combined_bwd, tri_source, tri_test)
                end
            end
        end
    end
    
    return Z
end

"""Write 3×3 Z_local contributions to global Z matrix using row locks."""
@inline function _write_to_Z!(Z, row_locks, Z_local, tri_row, tri_col)
    @inbounds for i in 1:3
        row_idx = tri_row.inBfsID[i]
        row_idx == 0 && continue
        sign_row = tri_row.bfsSign[i]
        lock(row_locks[row_idx])
        for j in 1:3
            col_idx = tri_col.inBfsID[j]
            if col_idx != 0
                Z[row_idx, col_idx] += Z_local[i, j] * sign_row * tri_col.bfsSign[j]
            end
        end
        unlock(row_locks[row_idx])
    end
end

"""
Merged EFIE+MFIE far-field kernel with shared Green's function evaluation.

For each (i,j) quadrature point pair, computes G and ∇G once, then accumulates:
- EFIE: ρ·ρ' and scalar divergence terms (symmetric)
- MFIE forward: K-operator with n̂_test
- MFIE backward: K-operator with n̂_source (swapped roles)
"""
function _cfie_far_merged!(Z_efie::Matrix{CT}, Z_mfie_fwd::Matrix{CT}, Z_mfie_bwd::Matrix{CT},
                            efie, mfie,
                            tri_test::TriangleInfo{IT, FT}, tri_source::TriangleInfo{IT, FT},
                            quad_points, C4divk2, k, efie_factor_val,
                            JK_0, eta_div_16pi) where {IT, FT, CT}
    gq = efie.gq_far
    w = gq.weight
    n_pts = length(w)
    
    r_test = quad_points[tri_test.triID]
    r_src = quad_points[tri_source.triID]
    
    # EFIE vertex precomputation
    v_src_1 = tri_source.vertices[:, 1]
    v_src_2 = tri_source.vertices[:, 2]
    v_src_3 = tri_source.vertices[:, 3]
    v_test_1 = tri_test.vertices[:, 1]
    v_test_2 = tri_test.vertices[:, 2]
    v_test_3 = tri_test.vertices[:, 3]
    
    # MFIE precomputation
    nt_hat = tri_test.facen̂    # test normal
    ns_hat = tri_source.facen̂  # source normal (for backward MFIE)
    
    v_test_sv = SVector{3, SVector{3, FT}}(SVector{3,FT}(v_test_1), SVector{3,FT}(v_test_2), SVector{3,FT}(v_test_3))
    v_src_sv = SVector{3, SVector{3, FT}}(SVector{3,FT}(v_src_1), SVector{3,FT}(v_src_2), SVector{3,FT}(v_src_3))
    
    @inbounds for j in 1:n_pts
        rgj = r_src[j]
        wj = w[j]
        
        # Precompute source rho vectors (shared)
        rho_src = (rgj - v_src_sv[1], rgj - v_src_sv[2], rgj - v_src_sv[3])
        
        # MFIE precomputation for source: n̂_test · ρ_src
        nt_dot_rho_src = (dot(nt_hat, rho_src[1]), dot(nt_hat, rho_src[2]), dot(nt_hat, rho_src[3]))
        
        for i in 1:n_pts
            rgi = r_test[i]
            wi = w[i]
            
            # === Shared: Green's function and gradient ===
            rvec = rgi - rgj
            diff1 = rvec[1]; diff2 = rvec[2]; diff3 = rvec[3]
            R = sqrt(diff1^2 + diff2^2 + diff3^2)
            R < FT(1e-12) && continue
            divr = one(FT) / R
            
            @fastmath ejkR = exp(-JK_0 * R)
            G_over_R = ejkR * divr         # exp(-jkR)/R
            
            # === EFIE accumulation ===
            efie_val_common = G_over_R * (wi * wj * tri_test.area * tri_source.area)
            
            rho_test = (rgi - v_test_sv[1], rgi - v_test_sv[2], rgi - v_test_sv[3])
            
            for n in 1:3
                rho_n = rho_src[n]
                for m in 1:3
                    Z_efie[m, n] += (dot(rho_test[m], rho_n) - C4divk2) * efie_val_common
                end
            end
            
            # === MFIE accumulation (forward: test→source) ===
            gw_ij = G_over_R * divr * wi * wj    # exp(-jkR)/R² * wi * wj
            temp = (JK_0 + divr) * gw_ij          # (jk + 1/R) * exp(-jkR)/R² * wi * wj
            
            nt_dot_rvec = dot(nt_hat, rvec)
            
            for n in 1:3
                for m in 1:3
                    term1 = dot(rho_test[m], rvec) * nt_dot_rho_src[n]
                    term2 = nt_dot_rvec * dot(rho_test[m], rho_src[n])
                    Z_mfie_fwd[m, n] += (term1 - term2) * temp
                end
            end
            
            # === MFIE accumulation (backward: source→test, swapped roles) ===
            # Now source triangle is "test", test triangle is "source"
            # rvec_bwd = rgj - rgi = -rvec
            # n̂ for backward = n̂_source
            ns_dot_rvec_neg = dot(ns_hat, -rvec)  # n̂_src · (-rvec)
            
            # Backward rho_test = ρ of source BFs at r_src (rgj)
            # Backward rho_src = ρ of test BFs at r_test (rgi)
            ns_dot_rho_test_bwd = (dot(ns_hat, rho_test[1]), dot(ns_hat, rho_test[2]), dot(ns_hat, rho_test[3]))
            
            for n in 1:3  # n indexes test BFs (backward "source")
                for m in 1:3  # m indexes source BFs (backward "test")
                    # Backward: ρ_m_bwd = rho_src[m], rvec_bwd = -rvec, n̂_bwd = ns_hat
                    # ρ_n_bwd = rho_test[n] (test BFs evaluated at rgi)
                    term1_bwd = dot(rho_src[m], -rvec) * ns_dot_rho_test_bwd[n]
                    term2_bwd = ns_dot_rvec_neg * dot(rho_src[m], rho_test[n])
                    Z_mfie_bwd[m, n] += (term1_bwd - term2_bwd) * temp
                end
            end
        end
    end
    
    # Apply EFIE scaling: edge lengths, area normalization, factor
    inv_areas_efie = one(FT) / (tri_test.area * tri_source.area)
    @inbounds for n in 1:3
        ln = tri_source.edgel[n]
        for m in 1:3
            lm = tri_test.edgel[m]
            Z_efie[m, n] *= lm * ln * inv_areas_efie * efie_factor_val
        end
    end
    
    # Apply MFIE forward scaling: edge lengths and constant
    @inbounds for n in 1:3
        ln = tri_source.edgel[n]
        for m in 1:3
            lm = tri_test.edgel[m]
            Z_mfie_fwd[m, n] *= lm * ln * eta_div_16pi
        end
    end
    
    # Apply MFIE backward scaling: edge lengths (source/test swapped) and constant
    @inbounds for n in 1:3
        ln = tri_test.edgel[n]  # backward "source" = original test
        for m in 1:3
            lm = tri_source.edgel[m]  # backward "test" = original source
            Z_mfie_bwd[m, n] *= lm * ln * eta_div_16pi
        end
    end
    
    return nothing
end

end
