module SCFIEModule

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..Kernels
using ..CFIEModule
using ..VEFIEModule
using StaticArrays
using LinearAlgebra
using SparseArrays
using Base.Threads

import ..CoreModule: assemble_impedance_matrix

export SCFIE, assemble_impedance_matrix

"""
    SCFIE{FT, CT, N_GQ_S, N_GQ_V} <: AbstractIntegralOperator

Surface-Volume Combined Field Integral Equation (SCFIE) operator.
Couples Surface EFIE/MFIE with Volume EFIE.

# Fields
- `freq`: Operating frequency.
- `k`: Wavenumber (background).
- `eta`: Intrinsic impedance (background).
- `alpha`: Coupling parameter (0.0 = EFIE, 1.0 = MFIE).
- `gq_surf`: Gauss quadrature info for triangle.
- `gq_vol`: Gauss quadrature info for tetrahedron.
- `permittivities`: Vector of permittivities for volume elements.
"""
struct SCFIE{FT<:AbstractFloat, CT<:Complex, N_GQ_S, N_GQ_V} <: AbstractIntegralOperator
    freq::FT
    k::FT
    eta::FT
    alpha::FT # For CFIE part on surface
    gq_surf::GaussQuadratureInfoStruct{FT, N_GQ_S, 3}
    gq_vol::GaussQuadratureInfoStruct{FT, N_GQ_V, 4}
    permittivities::Vector{CT}
end

function SCFIE(freq::FT, permittivities::Vector{Complex{FT}}; alpha::FT=0.5) where {FT}
    c0 = 299792458.0
    mu0 = 4π * 1e-7
    eps0 = 1.0 / (c0^2 * mu0)
    k = 2π * freq / c0
    eta = sqrt(mu0 / eps0)
    
    # Quadrature rules
    gq_surf = GaussQuadratureInfo(:Triangle, 7, FT) # 7-point for surface
    gq_vol = GaussQuadratureInfo(:Tetrahedron, 5, FT) # 5-point for volume
    
    return SCFIE{FT, Complex{FT}, 7, 5}(freq, k, eta, alpha, gq_surf, gq_vol, permittivities)
end

"""
    assemble_impedance_matrix(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::SWGBasis)

Assemble the full coupled impedance matrix for SCFIE.
Structure:
[ Z_SS  Z_SV ]
[ Z_VS  Z_VV ]
"""
function assemble_impedance_matrix(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::SWGBasis)
    FT = eltype(scfie.freq)
    CT = Complex{FT}
    
    n_surf = num_basis(surf_basis)
    n_vol = num_basis(vol_basis)
    n_total = n_surf + n_vol
    
    Z = zeros(CT, n_total, n_total)
    
    # 1. Surface-Surface Block (Z_SS)
    # We can reuse EFIE/MFIE/CFIE logic here, but we need to be careful about the implementation.
    # For now, let's implement a dedicated loop or reuse existing if possible.
    # Ideally, we should call `assemble_impedance_matrix(cfie, surf_basis)`
    # But we need to construct a CFIE object.
    
    # Construct temporary CFIE object
    # Note: CFIE in EMSuite might not be fully ready or we want to be explicit.
    # Let's assume we can use the existing CFIE implementation.
    # But wait, CFIE.jl exists.
    
    # 2. Volume-Volume Block (Z_VV)
    # Reuse VEFIE logic.
    
    # 3. Surface-Volume (Z_SV) and Volume-Surface (Z_VS) Blocks
    # These are the new parts.
    
    # Let's implement the full assembly loop to handle everything consistently.
    # Or better, assemble blocks and concatenate.
    
    # Block 1: Z_SS
    # We need a CFIE object.
    cfie = CFIE(scfie.freq, scfie.alpha)
    Z_SS = assemble_impedance_matrix(cfie, surf_basis)
    Z[1:n_surf, 1:n_surf] = Z_SS
    
    # Block 4: Z_VV
    vefie = VEFIE(scfie.freq, scfie.permittivities)
    Z_VV = assemble_impedance_matrix(vefie, vol_basis)
    Z[n_surf+1:end, n_surf+1:end] = Z_VV
    
    # Block 2 & 3: Coupling
    assemble_coupling_blocks!(Z, scfie, surf_basis, vol_basis)
    
    return Z
end

function assemble_coupling_blocks!(Z::Matrix{CT}, scfie::SCFIE, surf_basis::RWGBasis, vol_basis::SWGBasis) where {CT}
    # Implement Z_SV and Z_VS
    # Loop over surface triangles and volume tetrahedra
    
    # Precompute geometry
    tris = get_triangles_info(surf_basis.mesh, surf_basis)
    tetras = get_tetrahedra_info(vol_basis.mesh, vol_basis, scfie.permittivities)
    
    ntri = length(tris)
    ntet = length(tetras)
    n_surf = num_basis(surf_basis)
    n_total = size(Z, 1)
    
    # Locks for thread safety (one per row)
    row_locks = [SpinLock() for _ in 1:n_total]
    
    println("SCFIE Coupling Assembly: $ntri triangles x $ntet tetrahedra.")
    
    Threads.@threads for it in 1:ntri
        tri = tris[it]
        
        for js in 1:ntet
            tet = tetras[js]
            
            # Compute interaction
            # Z_sv_elem: Surface Test (tri), Volume Source (tet)
            # Z_vs_elem: Volume Test (tet), Surface Source (tri)
            Z_sv_elem, Z_vs_elem = scfie_coupling_interaction(scfie, tri, tet)
            
            # Fill Z_SV (Top Right)
            # Rows: Surface (m), Cols: Volume (n)
            for i in 1:3
                m = tri.inBfsID[i]
                if m == 0; continue; end
                
                for j in 1:4
                    n = tet.inBfsID[j]
                    if n == 0; continue; end
                    
                    # Z[m, n_surf + n] += Z_sv_elem[i, j]
                    lock(row_locks[m])
                    try
                        Z[m, n_surf + n] += Z_sv_elem[i, j]
                    finally
                        unlock(row_locks[m])
                    end
                    
                    # Fill Z_VS (Bottom Left)
                    # Rows: Volume (n), Cols: Surface (m)
                    # Z[n_surf + n, m] += Z_vs_elem[j, i]
                    row_idx = n_surf + n
                    lock(row_locks[row_idx])
                    try
                        Z[row_idx, m] += Z_vs_elem[j, i]
                    finally
                        unlock(row_locks[row_idx])
                    end
                end
            end
        end
    end
    println("SCFIE Coupling Assembly Completed.")
end

function scfie_coupling_interaction(scfie::SCFIE{FT, CT, N_GQ_S, N_GQ_V}, tri::TriangleInfo, tet::TetrahedraInfo) where {FT, CT, N_GQ_S, N_GQ_V}
    # FT = eltype(scfie.freq)
    # CT = Complex{FT}
    
    Z_sv = @MMatrix zeros(CT, 3, 4) # Surface Test, Volume Source
    Z_vs = @MMatrix zeros(CT, 4, 3) # Volume Test, Surface Source
    
    k = scfie.k
    omega = 2π * scfie.freq
    mu0 = 4π * 1e-7
    eps0 = 8.854187817e-12
    
    # Material properties of volume source
    κ_vol = tet.κ
    
    # Quadrature
    gq_s = scfie.gq_surf
    gq_v = scfie.gq_vol
    
    Nq_s = length(gq_s.weight)
    Nq_v = length(gq_v.weight)
    
    r_q_s = tri.vertices * gq_s.coordinate
    r_q_v = tet.vertices * gq_v.coordinate
    
    # Constants
    # Both Z_SV and Z_VS use the same EFIE L-operator:
    #   L(f_test, f_src) = jωμ₀ <f_test, f_src>_G + 1/(jωε₀) <∇·f_test, ∇'·f_src>_G
    # This comes from testing E^scat = -(jωA + ∇Φ) with integration by parts:
    #   <f, -E^scat> = jωμ₀ <f, f'>_G + 1/(jωε₀) <∇·f, ∇'·f'>_G
    # Note: 1/(jωε₀) = -j/(ωε₀) = -jωμ₀/k² (negative imaginary)
    #
    # Z_SV = κ_vol × L(f_surf, f_vol)    [Surface Test, Volume Source]
    # Z_VS = L(f_vol, f_surf)             [Volume Test, Surface Source]
    #
    # Legacy uses: jωμ₀/(4π) × (ρ·ρ/6 - 1/k²) × exp(-jkR)/R
    # which is equivalent to: jωμ₀ × f·f × G₄π + 1/(jωε₀) × div·div × G₄π
    
    # Z_SV coefficients (with κ)
    c1_sv = im * omega * mu0 * κ_vol
    c2_sv = 1.0 / (im * omega * eps0) * κ_vol
    
    # Z_VS coefficients (without κ)
    c1_vs = im * omega * mu0
    c2_vs = 1.0 / (im * omega * eps0)
    
    vol_factor = tri.area * tet.volume
    
    # Double loop over quadrature points
    for j in 1:Nq_v # Volume (Tet)
        w_j = gq_v.weight[j]
        r_j = r_q_v[:, j]
        
        for i in 1:Nq_s # Surface (Tri)
            w_i = gq_s.weight[i]
            r_i = r_q_s[:, i]
            
            # Green's function
            R_vec = r_i - r_j
            R = norm(R_vec)
            
            if R < 1e-10
                G = zero(CT)
            else
                G = exp(-im * k * R) / (4π * R)
            end
            
            factor = w_i * w_j * vol_factor * G
            
            # Accumulate
            for m in 1:3 # Surface (RWG)
                # Compute f_m on the fly
                div_val_s = (tri.edgel[m] / tri.area) * tri.bfsSign[m]
                v_free_s = tri.vertices[:, m]
                const_val_s = (tri.edgel[m] / (2 * tri.area)) * tri.bfsSign[m]
                f_m = const_val_s * (r_i - v_free_s)
                d_m = div_val_s
                
                for n in 1:4 # Volume (SWG)
                    # Compute f_n on the fly
                    div_val_v = (tet.facesArea[n] / tet.volume) * tet.bfsSign[n]
                    v_free_v = tet.vertices[:, n]
                    const_val_v = (tet.facesArea[n] / (3 * tet.volume)) * tet.bfsSign[n]
                    f_n = const_val_v * (r_j - v_free_v)
                    d_n = div_val_v
                    
                    # Z_SV (Test=m, Source=n): κ × L(f_surf, f_vol)
                    term1 = c1_sv * dot(f_m, f_n)
                    term2 = c2_sv * d_m * d_n
                    Z_sv[m, n] += (term1 + term2) * factor
                    
                    # Z_VS (Test=n, Source=m): L(f_vol, f_surf)
                    term1_vs = c1_vs * dot(f_n, f_m)
                    term2_vs = c2_vs * d_n * d_m
                    Z_vs[n, m] += (term1_vs + term2_vs) * factor
                end
            end
        end
    end
    
    return SMatrix(Z_sv), SMatrix(Z_vs)
end
    
end
