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

export SCFIE, assemble_impedance_matrix, assemble_fss_boundary_correction_sparse,
       scfie_sv_only_interaction

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
struct SCFIE{FT<:AbstractFloat,CT<:Complex,N_GQ_S,N_GQ_V} <: AbstractIntegralOperator
    freq::FT
    k::FT
    eta::FT
    alpha::FT # For CFIE part on surface
    gq_surf::GaussQuadratureInfoStruct{FT,N_GQ_S,3}
    gq_vol::GaussQuadratureInfoStruct{FT,N_GQ_V,4}
    permittivities::Vector{CT}
end

function SCFIE(freq::FT, permittivities::Vector{Complex{FT}}; alpha::FT = 0.5) where {FT}
    c0 = 299792458.0
    mu0 = 4π * 1e-7
    eps0 = 1.0 / (c0^2 * mu0)
    k = 2π * freq / c0
    eta = sqrt(mu0 / eps0)

    # Quadrature rules
    gq_surf = GaussQuadratureInfo(:Triangle, 7, FT) # 7-point for surface
    gq_vol = GaussQuadratureInfo(:Tetrahedron, 5, FT) # 5-point for volume

    return SCFIE{FT,Complex{FT},7,5}(freq, k, eta, alpha, gq_surf, gq_vol, permittivities)
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

function assemble_coupling_blocks!(
    Z::Matrix{CT},
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::SWGBasis,
) where {CT}
    # Precompute geometry
    tris = get_triangles_info(surf_basis.mesh, surf_basis)
    tetras = get_tetrahedra_info(vol_basis.mesh, vol_basis, scfie.permittivities)

    ntri = length(tris)
    ntet = length(tetras)
    n_surf = num_basis(surf_basis)
    n_total = size(Z, 1)
    n_threads = Threads.nthreads()

    # Row SpinLocks — plain lock/unlock (no try/finally overhead for SpinLock)
    row_locks = [SpinLock() for _ = 1:n_total]

    # Cyclic scheduling over test triangles (same pattern as EFIE/VEFIE).
    # Each thread owns surface rows for its assigned triangles → minimal Z_SV contention.
    # Z_VS uses the identity Z_vs[n,m] = Z_sv[m,n] / κ_vol (proven from c1_vs = c1_sv/κ).
    Threads.@threads for tid = 1:n_threads
        for it = tid:n_threads:ntri
            tri = tris[it]
            for js = 1:ntet
                tet = tetras[js]

                # Compute Z_sv only (3×4); Z_vs derived via κ reciprocity
                Z_sv_local = scfie_sv_only_interaction(scfie, tri, tet)
                κ_inv = iszero(tet.κ) ? zero(CT) : CT(1.0) / tet.κ

                # ── Z_SV (surface rows) ─────────────────────────────────────────
                @inbounds for i = 1:3
                    m = tri.inBfsID[i]
                    m == 0 && continue
                    lock(row_locks[m])
                    @inbounds for j = 1:4
                        n = tet.inBfsID[j]
                        n == 0 && continue
                        Z[m, n_surf+n] += Z_sv_local[i, j]
                    end
                    unlock(row_locks[m])
                end

                # ── Z_VS (volume rows) via reciprocity: Z_vs[n,m] = Z_sv[m,n]/κ ─
                @inbounds for j = 1:4
                    n = tet.inBfsID[j]
                    n == 0 && continue
                    row_vs = n_surf + n
                    lock(row_locks[row_vs])
                    @inbounds for i = 1:3
                        m = tri.inBfsID[i]
                        m == 0 && continue
                        Z[row_vs, m] += Z_sv_local[i, j] * κ_inv
                    end
                    unlock(row_locks[row_vs])
                end
            end
        end
    end

    # ── Fss boundary correction ─────────────────────────────────────────────────
    assemble_fss_boundary_correction!(Z, scfie, surf_basis, vol_basis, tris)
end

"""
    assemble_fss_boundary_correction!(Z, scfie, surf_basis, vol_basis, tris)

Add boundary surface integral corrections (Fss terms) for half-SWG basis functions.

For SWG basis functions on the volume mesh boundary, the scalar potential
requires an additional surface integral on the boundary face. This corrects
for the surface divergence contribution:

    ΔZ_VS[n,m] = jωμ₀/(4πk²) × l_m × |A_n| × ∫∫ G(r,r') dS_tri dS_face
    ΔZ_SV[m,n] = κ × ΔZ_VS[n,m]

where G(r,r') = exp(-jkR)/R (without 1/4π factor).
"""
function assemble_fss_boundary_correction!(
    Z::Matrix{CT},
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::SWGBasis,
    tris,
) where {CT}

    k = scfie.k
    omega = 2π * scfie.freq
    mu0 = 4π * 1e-7

    # Coefficient: jωμ₀/(4π) × 1/k²
    coeff = im * omega * mu0 / (4π * k^2)

    n_surf = num_basis(surf_basis)
    N_total = size(Z, 1)

    # Triangle quadrature for the face-to-face integral
    FT = eltype(scfie.freq)
    gq = scfie.gq_surf  # 7-point triangle quadrature
    Nq = length(gq.weight)

    # Volume mesh data
    vol_mesh = vol_basis.mesh
    vol_nodes = vol_mesh.node
    vol_elems = vol_mesh.tetras

    # Collect boundary SWG indices for parallel iteration
    boundary_indices = Int[]
    for n = 1:num_basis(vol_basis)
        if vol_basis.functions[n].is_boundary
            push!(boundary_indices, n)
        end
    end
    n_boundary = length(boundary_indices)

    if n_boundary == 0
        println("Fss Boundary Correction: 0 boundary SWG functions (skipped).")
        return
    end

    # Per-row locks for thread-safe Z updates
    row_locks = [Base.Threads.SpinLock() for _ = 1:N_total]

    # Parallel over boundary SWG functions
    Threads.@threads for idx = 1:n_boundary
        n = boundary_indices[idx]
        bf = vol_basis.functions[n]

        # Get boundary face vertices
        tet_idx = bf.support[1]
        local_face = bf.local_face_idx[1]

        v_indices = vol_elems[:, tet_idx]

        # Face i is opposite to vertex i — get the other 3 vertices
        face_v = Vector{SVector{3,FT}}(undef, 3)
        fi = 1
        for kv = 1:4
            if kv != local_face
                face_v[fi] = SVector{3,FT}(vol_nodes[:, v_indices[kv]])
                fi += 1
            end
        end

        # Material contrast
        perm_idx = tet_idx
        eps_r = scfie.permittivities[perm_idx]
        κ_tet = (eps_r - 1.0) / eps_r
        δκ = κ_tet

        # Face area
        abs_arean = bf.area

        # Precompute quadrature points on the boundary face
        r_face = Matrix{FT}(undef, 3, Nq)
        for q = 1:Nq
            u = gq.coordinate[1, q]
            v = gq.coordinate[2, q]
            w = gq.coordinate[3, q]
            r_face[:, q] = u * face_v[1] + v * face_v[2] + w * face_v[3]
        end

        # Iterate over test triangles
        for it = 1:length(tris)
            tri = tris[it]

            # Quadrature points on test triangle 
            r_tri = tri.vertices * gq.coordinate

            # Compute Fss = ∫∫ G(r_tri, r_face) dS_tri dS_face
            Fss = zero(CT)
            for gi = 1:Nq
                @views rgi = r_tri[:, gi]
                for gj = 1:Nq
                    @views rgj = r_face[:, gj]
                    R_vec = rgi - rgj
                    R = norm(R_vec)
                    if R < 1e-10
                        continue
                    end
                    @fastmath G = exp(-im * k * R) / R
                    Fss += G * gq.weight[gi] * gq.weight[gj]
                end
            end

            # For each surface basis function on this triangle
            for mi = 1:3
                m = tri.inBfsID[mi]
                if m == 0
                    continue
                end

                lm = tri.edgel[mi]
                temp = coeff * lm * abs_arean * Fss

                # Z_VS correction: Z[n_surf + n, m] += temp (thread-safe)
                vol_row = n_surf + n
                lock(row_locks[vol_row])
                Z[vol_row, m] += temp
                unlock(row_locks[vol_row])

                # Z_SV correction: Z[m, n_surf + n] += δκ × temp (thread-safe)
                lock(row_locks[m])
                Z[m, vol_row] += δκ * temp
                unlock(row_locks[m])
            end
        end
    end

    println("Fss Boundary Correction: $n_boundary boundary SWG functions processed (parallel).")
end

"""
    scfie_sv_only_interaction(scfie, tri, tet)

Compute only the Surface-Test / Volume-Source block Z_sv (3×4).

The Volume-Test / Surface-Source block Z_vs is derived from Z_sv via the
identity Z_vs[n,m] = Z_sv[m,n] / κ_vol, which follows from:
  c1_sv = im*ω*μ₀*κ,   c2_sv = κ/(im*ω*ε₀)
  c1_vs = im*ω*μ₀,      c2_vs = 1/(im*ω*ε₀)
Since the Green's function and quadrature are shared, Z_vs[n,m] = Z_sv[m,n]/κ.
Using this identity avoids computing 12 extra inner products per quadrature pair.
"""
function scfie_sv_only_interaction(
    scfie::SCFIE{FT,CT,N_GQ_S,N_GQ_V},
    tri::TriangleInfo,
    tet::TetrahedraInfo,
) where {FT,CT,N_GQ_S,N_GQ_V}
    Z_sv = @MMatrix zeros(CT, 3, 4)

    k = scfie.k
    omega = 2π * scfie.freq
    mu0 = 4π * 1e-7
    eps0 = 8.854187817e-12

    κ_vol = tet.κ
    c1_sv = im * omega * mu0 * κ_vol
    c2_sv = κ_vol / (im * omega * eps0)
    vol_factor = tri.area * tet.volume

    gq_s = scfie.gq_surf
    gq_v = scfie.gq_vol
    Nq_s = length(gq_s.weight)
    Nq_v = length(gq_v.weight)
    r_q_s = tri.vertices * gq_s.coordinate
    r_q_v = tet.vertices * gq_v.coordinate

    for j = 1:Nq_v
        w_j = gq_v.weight[j]
        r_j = r_q_v[:, j]

        for i = 1:Nq_s
            w_i = gq_s.weight[i]
            r_i = r_q_s[:, i]

            R = norm(r_i - r_j)
            R < 1e-10 && continue
            G = exp(-im * k * R) / (4π * R)
            factor = w_i * w_j * vol_factor * G

            @inbounds for m = 1:3
                c_sm = (tri.edgel[m] / (2 * tri.area)) * tri.bfsSign[m]
                f_m = c_sm * (r_i - tri.vertices[:, m])
                d_m = (tri.edgel[m] / tri.area) * tri.bfsSign[m]
                @inbounds for n = 1:4
                    c_vn = (tet.facesArea[n] / (3 * tet.volume)) * tet.bfsSign[n]
                    f_n = c_vn * (r_j - tet.vertices[:, n])
                    d_n = (tet.facesArea[n] / tet.volume) * tet.bfsSign[n]
                    Z_sv[m, n] += (c1_sv * dot(f_m, f_n) + c2_sv * d_m * d_n) * factor
                end
            end
        end
    end

    return SMatrix(Z_sv)
end

function scfie_coupling_interaction(
    scfie::SCFIE{FT,CT,N_GQ_S,N_GQ_V},
    tri::TriangleInfo,
    tet::TetrahedraInfo,
) where {FT,CT,N_GQ_S,N_GQ_V}
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
    for j = 1:Nq_v # Volume (Tet)
        w_j = gq_v.weight[j]
        r_j = r_q_v[:, j]

        for i = 1:Nq_s # Surface (Tri)
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
            for m = 1:3 # Surface (RWG)
                # Compute f_m on the fly
                div_val_s = (tri.edgel[m] / tri.area) * tri.bfsSign[m]
                v_free_s = tri.vertices[:, m]
                const_val_s = (tri.edgel[m] / (2 * tri.area)) * tri.bfsSign[m]
                f_m = const_val_s * (r_i - v_free_s)
                d_m = div_val_s

                for n = 1:4 # Volume (SWG)
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

"""
    assemble_fss_boundary_correction_sparse(scfie, surf_basis, vol_basis)

Return the Fss boundary correction as a sparse matrix (original indexing).
Used by MLFMAOperator to add the correction to Z_near.
"""
function assemble_fss_boundary_correction_sparse(
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::SWGBasis,
)
    CT = Complex{eltype(scfie.freq)}

    k = scfie.k
    omega = 2π * scfie.freq
    mu0 = 4π * 1e-7
    FT = eltype(scfie.freq)

    coeff = im * omega * mu0 / (4π * k^2)

    n_surf = num_basis(surf_basis)
    n_total = n_surf + num_basis(vol_basis)

    gq = scfie.gq_surf
    Nq = length(gq.weight)

    tris = get_triangles_info(surf_basis.mesh, surf_basis)

    vol_mesh = vol_basis.mesh
    vol_nodes = vol_mesh.node
    vol_elems = vol_mesh.tetras

    Is = Int[]
    Js = Int[]
    Vs = CT[]

    for n = 1:num_basis(vol_basis)
        bf = vol_basis.functions[n]
        if !bf.is_boundary
            continue
        end

        tet_idx = bf.support[1]
        local_face = bf.local_face_idx[1]
        v_indices = vol_elems[:, tet_idx]

        face_v = Vector{SVector{3,FT}}(undef, 3)
        fi = 1
        for kv = 1:4
            if kv != local_face
                face_v[fi] = SVector{3,FT}(vol_nodes[:, v_indices[kv]])
                fi += 1
            end
        end

        eps_r = scfie.permittivities[tet_idx]
        κ_tet = (eps_r - 1.0) / eps_r
        δκ = κ_tet

        abs_arean = bf.area

        r_face = Matrix{FT}(undef, 3, Nq)
        for q = 1:Nq
            u = gq.coordinate[1, q]
            v = gq.coordinate[2, q]
            w = gq.coordinate[3, q]
            r_face[:, q] = u * face_v[1] + v * face_v[2] + w * face_v[3]
        end

        for it = 1:length(tris)
            tri = tris[it]
            r_tri = tri.vertices * gq.coordinate

            Fss = zero(CT)
            for gi = 1:Nq
                @views rgi = r_tri[:, gi]
                for gj = 1:Nq
                    @views rgj = r_face[:, gj]
                    R_vec = rgi - rgj
                    R = norm(R_vec)
                    if R < 1e-10
                        continue
                    end
                    G = exp(-im * k * R) / R
                    Fss += G * gq.weight[gi] * gq.weight[gj]
                end
            end

            for mi = 1:3
                m = tri.inBfsID[mi]
                if m == 0
                    continue
                end

                lm = tri.edgel[mi]
                temp = coeff * lm * abs_arean * Fss

                # Z_VS: row = n_surf + n, col = m
                push!(Is, n_surf + n)
                push!(Js, m)
                push!(Vs, temp)

                # Z_SV: row = m, col = n_surf + n
                push!(Is, m)
                push!(Js, n_surf + n)
                push!(Vs, δκ * temp)
            end
        end
    end

    return sparse(Is, Js, Vs, n_total, n_total)
end

# ============================================================================
# SCFIE for RWG + PWC (Surface-Volume Coupled IE with PWC volume basis)
# ============================================================================

"""
    assemble_impedance_matrix(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::PWCBasis)

Assemble the full coupled impedance matrix for SCFIE using PWC volume basis.

Structure:
    [ Z_SS  Z_SV ]
    [ Z_VS  Z_VV ]

where Z_SS is the CFIE surface block (RWG), Z_VV is the VEFIE volume block (PWC),
and Z_SV/Z_VS are the dyadic coupling blocks.

# Legacy Parity
Matches `MoM_Kernels` `impedancemat4RWGPWC!` coupling with `EFIEOnRWGPWC` kernel.
"""
function assemble_impedance_matrix(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::PWCBasis)
    FT = eltype(scfie.freq)
    CT = Complex{FT}

    n_surf = num_basis(surf_basis)
    n_vol = num_basis(vol_basis)
    n_total = n_surf + n_vol

    Z = zeros(CT, n_total, n_total)

    # Block 1: Z_SS (Surface CFIE)
    cfie = CFIE(scfie.freq, scfie.alpha)
    Z_SS = assemble_impedance_matrix(cfie, surf_basis)
    Z[1:n_surf, 1:n_surf] = Z_SS

    # Block 4: Z_VV (Volume VEFIE with PWC)
    vefie = VEFIE(scfie.freq, scfie.permittivities)
    Z_VV = assemble_impedance_matrix(vefie, vol_basis)
    Z[n_surf+1:end, n_surf+1:end] = Z_VV

    # Block 2 & 3: Coupling (RWG-PWC dyadic coupling)
    assemble_coupling_blocks_pwc!(Z, scfie, surf_basis, vol_basis)

    return Z
end

"""
    assemble_coupling_blocks_pwc!(Z, scfie, surf_basis, vol_basis)

Assemble the surface-volume coupling blocks Z_SV and Z_VS for RWG+PWC.

Uses the dyadic L operator: (k²I + ∇∇) G(R) to couple:
- Surface RWG test functions with volume PWC source unknowns (Z_SV, with κ)
- Volume PWC test functions with surface RWG source current (Z_VS, no κ)

No Fss boundary correction is needed for PWC (PWC has no half-basis functions).
"""
function assemble_coupling_blocks_pwc!(
    Z::Matrix{CT},
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::PWCBasis,
) where {CT}
    # Precompute geometry
    tris = get_triangles_info(surf_basis.mesh, surf_basis)
    tetras = get_tetrahedra_info(vol_basis.mesh, vol_basis, scfie.permittivities)

    ntri = length(tris)
    ntet = length(tetras)
    n_surf = num_basis(surf_basis)
    n_total = size(Z, 1)

    # Constants
    FT = eltype(scfie.freq)
    k = scfie.k
    k² = k^2
    jk = im * k
    omega = 2π * scfie.freq
    mu0 = 4π * 1e-7
    eps0 = 8.854187817e-12
    eta0 = sqrt(mu0 / eps0)
    Jη₀divK = im * eta0 / k  # = j/(ωε₀)
    div4π = 1.0 / (4π)

    # Quadrature
    gq_s = scfie.gq_surf
    gq_v = scfie.gq_vol
    Nq_s = length(gq_s.weight)
    Nq_v = length(gq_v.weight)

    # Thread safety locks
    row_locks = [SpinLock() for _ = 1:n_total]

    println("SCFIE-PWC Coupling Assembly: $ntri triangles x $ntet tetrahedra.")

    Threads.@threads for it = 1:ntri
        tri = tris[it]

        # Precompute triangle quadrature points
        r_q_tri = tri.vertices * gq_s.coordinate

        for js = 1:ntet
            tet = tetras[js]
            κs = tet.κ
            Vs = tet.volume

            # Precompute tetrahedron quadrature points
            r_q_tet = tet.vertices * gq_v.coordinate

            # Compute dyadic L integral for each surface GQ point
            # Z_SV[mi, ni]: surfTest mi, volSource ni
            # Z_VS[ni, mi]: volTest ni, surfSource mi

            for gi = 1:Nq_s
                rgi = @view r_q_tri[:, gi]

                # Compute dyadic Green's function integral over the volume tetrahedron
                # dyadG (3×3) = Σ_gj w_gj * GR * L_dyad
                dyadG = zeros(CT, 3, 3)

                for gj = 1:Nq_v
                    rgj = @view r_q_tet[:, gj]

                    Rx = rgi[1] - rgj[1]
                    Ry = rgi[2] - rgj[2]
                    Rz = rgi[3] - rgj[3]
                    R = sqrt(Rx^2 + Ry^2 + Rz^2)

                    if R < 1e-10
                        continue
                    end

                    divR = 1.0 / R
                    jkplusR1stdivR1st = (jk + divR) * divR

                    R̂x = Rx * divR
                    R̂y = Ry * divR
                    R̂z = Rz * divR

                    GR = exp(-jk * R) * div4π * divR * gq_v.weight[gj]

                    # R̂R̂ dyad
                    RR11 = R̂x * R̂x
                    RR12 = R̂x * R̂y
                    RR13 = R̂x * R̂z
                    RR22 = R̂y * R̂y
                    RR23 = R̂y * R̂z
                    RR33 = R̂z * R̂z

                    # L_dyad diagonal: (1-R̂ᵢR̂ⱼ)*k² - (1-3R̂ᵢR̂ⱼ)*(jk+1/R)/R
                    dyadG[1, 1] += GR * ((1 - RR11) * k² - (1 - 3RR11) * jkplusR1stdivR1st)
                    dyadG[2, 2] += GR * ((1 - RR22) * k² - (1 - 3RR22) * jkplusR1stdivR1st)
                    dyadG[3, 3] += GR * ((1 - RR33) * k² - (1 - 3RR33) * jkplusR1stdivR1st)

                    # L_dyad off-diagonal: -R̂ᵢR̂ⱼ*k² + 3R̂ᵢR̂ⱼ*(jk+1/R)/R
                    od12 = GR * (-RR12 * k² + 3RR12 * jkplusR1stdivR1st)
                    od13 = GR * (-RR13 * k² + 3RR13 * jkplusR1stdivR1st)
                    od23 = GR * (-RR23 * k² + 3RR23 * jkplusR1stdivR1st)

                    dyadG[1, 2] += od12
                    dyadG[2, 1] += od12
                    dyadG[1, 3] += od13
                    dyadG[3, 1] += od13
                    dyadG[2, 3] += od23
                    dyadG[3, 2] += od23
                end

                # Contract with surface RWG basis functions
                for mi = 1:3
                    m = tri.inBfsID[mi]
                    if m == 0
                        continue
                    end

                    lm = tri.edgel[mi]
                    freeVm = tri.vertices[:, mi]
                    ρmi = SVector(rgi[1] - freeVm[1], rgi[2] - freeVm[2], rgi[3] - freeVm[3])

                    temp = gq_s.weight[gi] * lm / 2

                    for ni = 1:3
                        n = tet.inBfsID[ni]

                        # Z_SV: ρ · dyadG[:, ni] — surface test, volume source
                        # In Z_SV, κ is included
                        dot_sv =
                            ρmi[1] * dyadG[1, ni] + ρmi[2] * dyadG[2, ni] + ρmi[3] * dyadG[3, ni]
                        z_sv = temp * dot_sv * Jη₀divK * Vs * κs

                        # Z_VS: ρ · dyadG[ni, :] — volume test, surface source
                        # In Z_VS, no κ
                        dot_vs =
                            ρmi[1] * dyadG[ni, 1] + ρmi[2] * dyadG[ni, 2] + ρmi[3] * dyadG[ni, 3]
                        z_vs = temp * dot_vs * Jη₀divK * Vs

                        # Fill Z_SV (row: m, col: n_surf + n)
                        lock(row_locks[m])
                        try
                            Z[m, n_surf+n] += z_sv
                        finally
                            unlock(row_locks[m])
                        end

                        # Fill Z_VS (row: n_surf + n, col: m)
                        row_idx = n_surf + n
                        lock(row_locks[row_idx])
                        try
                            Z[row_idx, m] += z_vs
                        finally
                            unlock(row_locks[row_idx])
                        end
                    end
                end
            end
        end
    end

    println("SCFIE-PWC Coupling Assembly Completed.")
end

# ============================================================================
# SCFIE for RWG + PWCHexBasis (Surface-Volume with PWC on hexahedra)
# ============================================================================

"""
    assemble_impedance_matrix(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::PWCHexBasis)

Assemble the full coupled impedance matrix for SCFIE using PWC hexahedra volume basis.

Structure:
    [ Z_SS  Z_SV ]
    [ Z_VS  Z_VV ]

# Legacy Parity
Matches `MoM_Kernels` `impedancemat4RWGPWC!` coupling with `EFIEOnRWGPWC` kernel (hexa version).
"""
function assemble_impedance_matrix(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::PWCHexBasis)
    FT = eltype(scfie.freq)
    CT = Complex{FT}

    n_surf = num_basis(surf_basis)
    n_vol = num_basis(vol_basis)
    n_total = n_surf + n_vol

    Z = zeros(CT, n_total, n_total)

    # Block 1: Z_SS (Surface CFIE)
    cfie = CFIE(scfie.freq, scfie.alpha)
    Z_SS = assemble_impedance_matrix(cfie, surf_basis)
    Z[1:n_surf, 1:n_surf] = Z_SS

    # Block 4: Z_VV (Volume VEFIE with PWCHex)
    vefie = VEFIE(scfie.freq, scfie.permittivities)
    Z_VV = assemble_impedance_matrix(vefie, vol_basis)
    Z[n_surf+1:end, n_surf+1:end] = Z_VV

    # Block 2 & 3: Coupling (RWG-PWCHex dyadic coupling)
    assemble_coupling_blocks_pwc_hex!(Z, scfie, surf_basis, vol_basis)

    return Z
end

"""
    assemble_coupling_blocks_pwc_hex!(Z, scfie, surf_basis, vol_basis)

Assemble surface-volume coupling blocks Z_SV and Z_VS for RWG + PWC(Hexa).

Uses the dyadic L operator: (k²I + ∇∇) G(R) to couple surface RWG with volume PWC.
Same formula as RWG+PWC(Tetra) but with hexahedron GQ points.

No Fss boundary correction needed (PWC has no half-basis functions).
"""
function assemble_coupling_blocks_pwc_hex!(
    Z::Matrix{CT},
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::PWCHexBasis,
) where {CT}
    # Precompute geometry
    tris = get_triangles_info(surf_basis.mesh, surf_basis)
    hexas = get_hexahedra_info(vol_basis.mesh, vol_basis, scfie.permittivities)

    ntri = length(tris)
    nhex = length(hexas)
    n_surf = num_basis(surf_basis)
    n_total = size(Z, 1)

    # Constants
    FT = eltype(scfie.freq)
    k = scfie.k
    k² = k^2
    jk = im * k
    eta0 = scfie.eta
    Jη₀divK = im * eta0 / k  # = j/(ωε₀)
    div4π = FT(1.0 / (4π))

    # Hex GQ (created on the fly since SCFIE struct only has tet GQ)
    gq_hex = GaussQuadratureInfo(:Hexahedron, 8, FT)
    Nq_hex = length(gq_hex.weight)

    # Surface triangle GQ
    gq_s = scfie.gq_surf
    Nq_s = length(gq_s.weight)

    # Thread safety locks
    row_locks = [SpinLock() for _ = 1:n_total]

    println("SCFIE-PWCHex Coupling Assembly: $ntri triangles x $nhex hexahedra.")

    Threads.@threads for it = 1:ntri
        tri = tris[it]

        # Precompute triangle quadrature points
        r_q_tri = tri.vertices * gq_s.coordinate

        for js = 1:nhex
            hex = hexas[js]
            κs = hex.κ
            Vs = hex.volume

            # Precompute hexahedron quadrature points
            r_q_hex = hex.vertices * gq_hex.coordinate

            # Compute dyadic L integral for each surface GQ point
            for gi = 1:Nq_s
                rgi = @view r_q_tri[:, gi]

                # Accumulate dyadic Green's function over hex volume
                dyadG = zeros(CT, 3, 3)

                for gj = 1:Nq_hex
                    rgj = @view r_q_hex[:, gj]

                    Rx = rgi[1] - rgj[1]
                    Ry = rgi[2] - rgj[2]
                    Rz = rgi[3] - rgj[3]
                    R = sqrt(Rx^2 + Ry^2 + Rz^2)

                    if R < 1e-10
                        continue
                    end

                    divR = 1.0 / R
                    jkplusR1stdivR1st = (jk + divR) * divR

                    R̂x = Rx * divR
                    R̂y = Ry * divR
                    R̂z = Rz * divR

                    GR = exp(-jk * R) * div4π * divR * gq_hex.weight[gj]

                    # R̂R̂ dyad components
                    RR11 = R̂x * R̂x
                    RR12 = R̂x * R̂y
                    RR13 = R̂x * R̂z
                    RR22 = R̂y * R̂y
                    RR23 = R̂y * R̂z
                    RR33 = R̂z * R̂z

                    # L_dyad: (I-R̂R̂)k² - (I/3-R̂R̂)×3(jk+1/R)/R
                    coeff3 = 3 * jkplusR1stdivR1st
                    dyadG[1, 1] += GR * ((1 - RR11) * k² - (1 / 3 - RR11) * coeff3)
                    dyadG[2, 2] += GR * ((1 - RR22) * k² - (1 / 3 - RR22) * coeff3)
                    dyadG[3, 3] += GR * ((1 - RR33) * k² - (1 / 3 - RR33) * coeff3)

                    od12 = GR * (-RR12 * k² + RR12 * coeff3)
                    od13 = GR * (-RR13 * k² + RR13 * coeff3)
                    od23 = GR * (-RR23 * k² + RR23 * coeff3)

                    dyadG[1, 2] += od12
                    dyadG[2, 1] += od12
                    dyadG[1, 3] += od13
                    dyadG[3, 1] += od13
                    dyadG[2, 3] += od23
                    dyadG[3, 2] += od23
                end

                # Contract with surface RWG basis functions
                for mi = 1:3
                    m = tri.inBfsID[mi]
                    if m == 0
                        continue
                    end

                    lm = tri.edgel[mi]
                    freeVm = tri.vertices[:, mi]
                    ρmi = SVector(rgi[1] - freeVm[1], rgi[2] - freeVm[2], rgi[3] - freeVm[3])

                    temp = gq_s.weight[gi] * lm / 2

                    for ni = 1:3
                        n = hex.inBfsID[ni]

                        # Z_SV: ρ · dyadG[:, ni], with κ
                        dot_sv =
                            ρmi[1] * dyadG[1, ni] + ρmi[2] * dyadG[2, ni] + ρmi[3] * dyadG[3, ni]
                        z_sv = temp * dot_sv * Jη₀divK * Vs * κs

                        # Z_VS: ρ · dyadG[ni, :], no κ
                        dot_vs =
                            ρmi[1] * dyadG[ni, 1] + ρmi[2] * dyadG[ni, 2] + ρmi[3] * dyadG[ni, 3]
                        z_vs = temp * dot_vs * Jη₀divK * Vs

                        lock(row_locks[m])
                        try
                            Z[m, n_surf+n] += z_sv
                        finally
                            unlock(row_locks[m])
                        end

                        row_idx = n_surf + n
                        lock(row_locks[row_idx])
                        try
                            Z[row_idx, m] += z_vs
                        finally
                            unlock(row_locks[row_idx])
                        end
                    end
                end
            end
        end
    end

    println("SCFIE-PWCHex Coupling Assembly Completed.")
end

# ============================================================================
# SCFIE for RWG + RBF (Surface-Volume with RBF hexahedra basis)
# ============================================================================

"""
    assemble_impedance_matrix(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::RBFBasis)

Assemble the full coupled impedance matrix for SCFIE using RBF hexahedra volume basis.

Structure:
    [ Z_SS  Z_SV ]
    [ Z_VS  Z_VV ]

Z_SV and Z_VS use the scalar potential form:
    Fsv = Σ_gi Σ_gj (ρm·ρn/2 - 1/k²) G(rgi, rgj) × w_gi × w_gj
    Z_SV = κs × jkη₀/(4π) × lm × An × Fsv
    Z_VS = jkη₀/(4π) × lm × An × Fsv

Plus Fss boundary correction for half-RBF basis functions on boundary faces.

# Legacy Parity
Matches `MoM_Kernels` `EFIEVSIERWGRBF.jl`.
"""
function assemble_impedance_matrix(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::RBFBasis)
    FT = eltype(scfie.freq)
    CT = Complex{FT}

    n_surf = num_basis(surf_basis)
    n_vol = num_basis(vol_basis)
    n_total = n_surf + n_vol

    Z = zeros(CT, n_total, n_total)

    # Block 1: Z_SS (Surface CFIE)
    cfie = CFIE(scfie.freq, scfie.alpha)
    Z_SS = assemble_impedance_matrix(cfie, surf_basis)
    Z[1:n_surf, 1:n_surf] = Z_SS

    # Block 4: Z_VV (Volume VEFIE with RBF)
    vefie = VEFIE(scfie.freq, scfie.permittivities)
    Z_VV = assemble_impedance_matrix(vefie, vol_basis)
    Z[n_surf+1:end, n_surf+1:end] = Z_VV

    # Block 2 & 3: Coupling (RWG-RBF)
    assemble_coupling_blocks_rbf!(Z, scfie, surf_basis, vol_basis)

    return Z
end

"""
    assemble_coupling_blocks_rbf!(Z, scfie, surf_basis, vol_basis)

Assemble surface-volume coupling blocks Z_SV and Z_VS for RWG + RBF.

Uses scalar potential form (not dyadic L for better numerical stability with RBF):
    Fsv = Σ Σ (ρm·ρn/2 - 1/k²) G w_gi w_gj
    Fsv *= lm × An
    Z_SV[m,n] = κs × jkη₀/(4π) × Fsv
    Z_VS[n,m] = jkη₀/(4π) × Fsv

Plus Fss boundary correction for half-RBF basis functions:
    Fss = ∫∫ G(r_tri, r_face) dS_tri dS_face
    Fss *= lm × |An|
    temp = jkη₀/(4π) × (1/k²) × Fss
    Z_SV[m,n] += δκn × temp
    Z_VS[n,m] += temp (if is_boundary)
"""
function assemble_coupling_blocks_rbf!(
    Z::Matrix{CT},
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::RBFBasis,
) where {CT}
    # Precompute geometry
    tris = get_triangles_info(surf_basis.mesh, surf_basis)
    hexas = get_hexahedra_info(vol_basis.mesh, vol_basis, scfie.permittivities)

    ntri = length(tris)
    nhex = length(hexas)
    n_surf = num_basis(surf_basis)
    n_total = size(Z, 1)

    # Constants
    FT = eltype(scfie.freq)
    k = scfie.k
    k² = k^2
    jk = im * k
    eta0 = scfie.eta
    JKη₀div4π = im * k * eta0 / (4 * FT(π))  # jkη₀/(4π)
    divk² = 1.0 / k²
    div4π = FT(1.0 / (4π))

    # Hex GQ (on the fly)
    gq_hex = GaussQuadratureInfo(:Hexahedron, 8, FT)
    Nq_hex = length(gq_hex.weight)

    # Quad GQ for Fss (on the fly)
    gq_quad = GaussQuadratureInfo(:Quadrangle, 4, FT)
    Nq_quad = length(gq_quad.weight)

    # Surface triangle GQ
    gq_s = scfie.gq_surf
    Nq_s = length(gq_s.weight)

    # Thread safety locks
    row_locks = [SpinLock() for _ = 1:n_total]

    println("SCFIE-RBF Coupling Assembly: $ntri triangles x $nhex hexahedra.")

    Threads.@threads for it = 1:ntri
        tri = tris[it]

        # Precompute triangle quadrature points
        r_q_tri = tri.vertices * gq_s.coordinate

        # Precompute Green's function × weight for reuse across basis functions
        for js = 1:nhex
            hex = hexas[js]
            κs = hex.κ

            # Precompute hex quadrature points
            r_q_hex = hex.vertices * gq_hex.coordinate

            # Precompute gw matrix: G(rgi, rgj) × w_gi × w_gj
            gw = zeros(CT, Nq_s, Nq_hex)
            @inbounds for gj = 1:Nq_hex
                rgj = @view r_q_hex[:, gj]
                for gi = 1:Nq_s
                    rgi = @view r_q_tri[:, gi]
                    Rx = rgi[1] - rgj[1]
                    Ry = rgi[2] - rgj[2]
                    Rz = rgi[3] - rgj[3]
                    R = sqrt(Rx^2 + Ry^2 + Rz^2)
                    if R < 1e-10
                        continue
                    end
                    gw[gi, gj] = exp(-jk * R) * div4π / R * gq_s.weight[gi] * gq_hex.weight[gj]
                end
            end

            # Precompute free-end points for all 6 RBF functions
            freeVns_all = [get_free_vns(hex, fi, gq_hex.coordinate) for fi = 1:6]

            # Loop over RBF source basis functions (6 per hex)
            for ni = 1:6
                arean = hex.facesArea[ni]
                face = hex.faces[ni]
                δκn = face.δκ
                isbdn = face.isbd
                freeVns = freeVns_all[ni]

                # Loop over RWG test basis functions (3 per triangle)
                for mi = 1:3
                    m = tri.inBfsID[mi]
                    m == 0 && continue

                    lm = tri.edgel[mi]
                    freeVm = tri.vertices[:, mi]
                    lman = lm * arean

                    # --- Fsv term ---
                    Fsv = zero(CT)
                    @inbounds for gj = 1:Nq_hex
                        rgj = @view r_q_hex[:, gj]
                        # ρn: source RBF vector at gj
                        # free end for this GQ point (interpolated on opposite face)
                        freeVn = @view freeVns[:, gj]
                        ρnj_x = rgj[1] - freeVn[1]
                        ρnj_y = rgj[2] - freeVn[2]
                        ρnj_z = rgj[3] - freeVn[3]

                        for gi = 1:Nq_s
                            rgi = @view r_q_tri[:, gi]
                            # ρm: test RWG vector at gi
                            ρmi_x = rgi[1] - freeVm[1]
                            ρmi_y = rgi[2] - freeVm[2]
                            ρmi_z = rgi[3] - freeVm[3]

                            dot_val = ρmi_x * ρnj_x + ρmi_y * ρnj_y + ρmi_z * ρnj_z
                            Fsv += (dot_val / 2 - divk²) * gw[gi, gj]
                        end
                    end
                    Fsv *= lman

                    # Compute Z contributions
                    JKη₀div4πFsv = JKη₀div4π * Fsv
                    Zmn = κs * JKη₀div4πFsv  # Z_SV
                    Znm = JKη₀div4πFsv        # Z_VS

                    # --- Fss term (boundary correction) ---
                    if isbdn || !iszero(δκn)
                        # Compute face-to-face integral
                        rq_face = hex.faces[ni].vertices * gq_quad.coordinate

                        Fss = zero(CT)
                        @inbounds for gj = 1:Nq_quad
                            rgj = @view rq_face[:, gj]
                            for gi = 1:Nq_s
                                rgi = @view r_q_tri[:, gi]
                                Rx = rgi[1] - rgj[1]
                                Ry = rgi[2] - rgj[2]
                                Rz = rgi[3] - rgj[3]
                                R = sqrt(Rx^2 + Ry^2 + Rz^2)
                                if R < 1e-10
                                    continue
                                end
                                G = exp(-jk * R) * div4π / R
                                Fss += G * gq_s.weight[gi] * gq_quad.weight[gj]
                            end
                        end
                        Fss *= lm * abs(arean)

                        temp = JKη₀div4π * divk² * Fss
                        if !iszero(δκn)
                            Zmn += δκn * temp
                        end
                        if isbdn
                            Znm += temp
                        end
                    end

                    # Fill Z_SV
                    n = hex.inBfsID[ni]
                    lock(row_locks[m])
                    try
                        Z[m, n_surf+n] += Zmn
                    finally
                        unlock(row_locks[m])
                    end

                    # Fill Z_VS
                    row_idx = n_surf + n
                    lock(row_locks[row_idx])
                    try
                        Z[row_idx, m] += Znm
                    finally
                        unlock(row_locks[row_idx])
                    end
                end
            end
        end
    end

    println("SCFIE-RBF Coupling Assembly Completed.")
end


end
