module VEFIEModule

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..Kernels
using StaticArrays
using LinearAlgebra
using SparseArrays
using Base.Threads

import ..CoreModule: assemble_impedance_matrix

export VEFIE, assemble_impedance_matrix

"""
    VEFIE{FT, CT, N_GQ} <: AbstractIntegralOperator

Volume Electric Field Integral Equation (VEFIE) operator.

Solves for the electric flux density \$\\mathbf{D}\$ in a dielectric volume.

# Fields
- `freq`: Operating frequency.
- `k`: Wavenumber (background).
- `eta`: Intrinsic impedance (background).
- `gq_info`: Gauss quadrature info for tetrahedron.
"""
struct VEFIE{FT<:AbstractFloat, CT<:Complex, N_GQ, N_GQ_FAR} <: AbstractIntegralOperator
    freq::FT
    k::FT
    eta::FT
    gq_info::GaussQuadratureInfoStruct{FT, N_GQ, 4}
    gq_far::GaussQuadratureInfoStruct{FT, N_GQ_FAR, 4}
    permittivities::Vector{CT}
end

struct TetBasisCache{CT, NQ, NQ_FAR}
    r_q    ::SMatrix{3, NQ, Float64}
    f_vals ::SMatrix{4, NQ, SVector{3, CT}}
    div_f  ::SVector{4, CT}
    
    r_q_far    ::SMatrix{3, NQ_FAR, Float64}
    f_vals_far ::SMatrix{4, NQ_FAR, SVector{3, CT}}
end

function VEFIE(freq::FT, permittivities::Vector{Complex{FT}}) where {FT}
    c0 = 299792458.0
    mu0 = 4π * 1e-7
    eps0 = 1.0 / (c0^2 * mu0)
    k = 2π * freq / c0
    eta = sqrt(mu0 / eps0)
    
    # Use 5-point rule for near field
    gq_info = GaussQuadratureInfo(:Tetrahedron, 5, FT)
    # Use 1-point rule for far field
    gq_far = GaussQuadratureInfo(:Tetrahedron, 1, FT)
    
    return VEFIE{FT, Complex{FT}, 5, 1}(freq, k, eta, gq_info, gq_far, permittivities)
end

"""
    assemble_impedance_matrix(vefie::VEFIE, basis::SWGBasis)

Assemble the impedance matrix Z for the VEFIE using stored permittivities.
"""
function assemble_impedance_matrix(vefie::VEFIE, basis::SWGBasis)
    return assemble_impedance_matrix(vefie, basis, vefie.permittivities)
end

"""
    assemble_impedance_matrix(vefie::VEFIE, basis::SWGBasis, permittivities::Vector{ComplexF64})

Assemble the impedance matrix Z for the VEFIE using optimized element-based loop.
"""
function assemble_impedance_matrix(vefie::VEFIE, basis::SWGBasis, permittivities::Vector{ComplexF64})
    FT = eltype(vefie.freq)
    CT = Complex{FT}
    
    nbf = num_basis(basis)
    Z = zeros(CT, nbf, nbf)
    
    # Precompute geometry
    tetras = get_tetrahedra_info(basis.mesh, basis, permittivities)
    ntet = length(tetras)
    
    # Precompute Basis Cache
    println("VEFIE: Precomputing basis functions...")
    basis_cache = precompute_vefie_basis(vefie, tetras)
    
    # Locks for thread safety (one per column to minimize contention)
    col_locks = [SpinLock() for _ in 1:nbf]
    
    # Progress tracking
    progress_counter = Threads.Atomic{Int}(0)
    next_idx = Threads.Atomic{Int}(1) # For dynamic scheduling
    
    println("VEFIE Assembly: $ntet tetrahedra, $nbf unknowns. (Column-Parallel, No Symmetry)")
    
    # Dynamic scheduling loop over Source Tetrahedra (js)
    Threads.@threads for _ in 1:Threads.nthreads()
        # Allocate thread-local buffer for columns
        # We accumulate contributions to the 4 columns of the current source tet
        # Size: nbf x 4
        local_cols = zeros(CT, nbf, 4)
        
        while true
            js = Threads.atomic_add!(next_idx, 1)
            if js > ntet; break; end
            
            tet_s = tetras[js]
            cache_s = basis_cache[js]
            
            # Reset local buffer
            fill!(local_cols, zero(CT))
            
            # Loop over Test Tetrahedra (it) - Full loop (No Symmetry)
            for it in 1:ntet
                tet_t = tetras[it]
                cache_t = basis_cache[it]
                
                # Compute interaction matrix (4x4)
                # We only need Z_ts (Test t, Source s)
                Z_ts = vefie_element_interaction_kernel(vefie, tet_t, tet_s, cache_t, cache_s)
                
                # Accumulate to local_cols
                # local_cols[m, j] += Z_ts[i, j]
                # m is global row index (from tet_t)
                # j is local column index (1..4)
                
                for i in 1:4
                    m = tet_t.inBfsID[i]
                    if m == 0; continue; end
                    
                    for j in 1:4
                        # We don't check n here, we just fill the 4 columns
                        local_cols[m, j] += Z_ts[i, j]
                    end
                end
                
                # Add Mass Matrix if it == js (Self-term)
                if it == js
                    M_t = vefie_mass_matrix_cached(vefie, tet_t, cache_t)
                    for i in 1:4
                        m = tet_t.inBfsID[i]
                        if m == 0; continue; end
                        for j in 1:4
                            local_cols[m, j] += M_t[i, j]
                        end
                    end
                end
            end
            
            # Scatter local_cols to global Z
            # We lock the columns n corresponding to tet_s
            for j in 1:4
                n = tet_s.inBfsID[j]
                if n == 0; continue; end
                
                lock(col_locks[n])
                try
                    # Z[:, n] += local_cols[:, j]
                    # Manual loop for speed? Julia's broadcast is fast.
                    @views Z[:, n] .+= local_cols[:, j]
                finally
                    unlock(col_locks[n])
                end
            end
            
            # Print progress every 10 elements
            c = Threads.atomic_add!(progress_counter, 1)
            if c % 10 == 0
                print("\rVEFIE Assembly: $c / $ntet source elements processed.")
            end
        end
    end
    println("\nVEFIE Assembly Completed.")
    
    return Z
end

function vefie_mass_matrix(vefie::VEFIE, tet::TetrahedraInfo)
    FT = eltype(vefie.freq)
    CT = Complex{FT}
    M = @MMatrix zeros(CT, 4, 4)
    
    eps0 = 8.854187817e-12
    inv_eps = (1.0 - tet.κ) / eps0
    
    gq = vefie.gq_info
    r_q = tet.vertices * gq.coordinate
    
    omega = 2π * vefie.freq
    factor_base = inv_eps / (im * omega) * tet.volume
    
    # Precompute basis values
    f_vals = MVector{4, SVector{3, CT}}(undef)
    
    for k in 1:length(gq.weight)
        w = gq.weight[k]
        r = r_q[:, k]
        factor = w * factor_base
        
        for m in 1:4
            v_free = tet.vertices[:, m]
            rho = r - v_free
            val = (tet.facesArea[m] / (3 * tet.volume)) * rho
            val *= tet.bfsSign[m]
            f_vals[m] = val
        end
        
        for n in 1:4
            f_n = f_vals[n]
            for m in 1:4
                M[m, n] += dot(f_vals[m], f_n) * factor
            end
        end
    end
    
    return SMatrix(M)
end

function vefie_element_interaction(vefie::VEFIE{FT, CT, N_GQ, N_GQ_FAR}, tet_t::TetrahedraInfo, tet_s::TetrahedraInfo) where {FT, CT, N_GQ, N_GQ_FAR}
    # FT = eltype(vefie.freq)
    # CT = Complex{FT}
    
    Z_ts = @MMatrix zeros(CT, 4, 4)
    Z_st = @MMatrix zeros(CT, 4, 4)
    
    # Constants
    k = vefie.k
    omega = 2π * vefie.freq
    mu0 = 4π * 1e-7
    eps0 = 8.854187817e-12
    
    # Material properties
    κs = tet_s.κ
    κt = tet_t.κ
    
    # Quadrature
    gq = vefie.gq_info
    Nq = length(gq.weight)
    
    # Precompute quadrature points
    r_q_t = tet_t.vertices * gq.coordinate
    r_q_s = tet_s.vertices * gq.coordinate
    
    # Constants for terms
    c1_ts = -im * omega * mu0 * κs
    c2_ts = (1/(im * omega * eps0)) * κs
    
    c1_st = -im * omega * mu0 * κt
    c2_st = (1/(im * omega * eps0)) * κt
    
    vol_factor = tet_t.volume * tet_s.volume
    
    # Double loop over quadrature points
    for j in 1:Nq # Source
        w_j = gq.weight[j]
        r_j = r_q_s[:, j]
        
        for i in 1:Nq # Test
            w_i = gq.weight[i]
            r_i = r_q_t[:, i]
            
            # Green's function
            R_vec = r_i - r_j
            R = norm(R_vec)
            
            if R < 1e-10
                G = zero(CT)
            else
                G = exp(-im * k * R) / (4π * R)
            end
            
            factor = w_i * w_j * vol_factor * G
            
            # Vectorized Accumulation
            # Construct F_t (4x3) and D_t (4x1)
            # Rows of F_t are f_m vectors
            
            # Test Basis
            # m=1
            c_t1 = (tet_t.facesArea[1] / (3 * tet_t.volume)) * tet_t.bfsSign[1]
            v_t1 = tet_t.vertices[:, 1]
            f_t1 = c_t1 * (r_i - v_t1)
            d_t1 = (tet_t.facesArea[1] / tet_t.volume) * tet_t.bfsSign[1]
            
            # m=2
            c_t2 = (tet_t.facesArea[2] / (3 * tet_t.volume)) * tet_t.bfsSign[2]
            v_t2 = tet_t.vertices[:, 2]
            f_t2 = c_t2 * (r_i - v_t2)
            d_t2 = (tet_t.facesArea[2] / tet_t.volume) * tet_t.bfsSign[2]
            
            # m=3
            c_t3 = (tet_t.facesArea[3] / (3 * tet_t.volume)) * tet_t.bfsSign[3]
            v_t3 = tet_t.vertices[:, 3]
            f_t3 = c_t3 * (r_i - v_t3)
            d_t3 = (tet_t.facesArea[3] / tet_t.volume) * tet_t.bfsSign[3]
            
            # m=4
            c_t4 = (tet_t.facesArea[4] / (3 * tet_t.volume)) * tet_t.bfsSign[4]
            v_t4 = tet_t.vertices[:, 4]
            f_t4 = c_t4 * (r_i - v_t4)
            d_t4 = (tet_t.facesArea[4] / tet_t.volume) * tet_t.bfsSign[4]
            
            # Source Basis
            # n=1
            c_s1 = (tet_s.facesArea[1] / (3 * tet_s.volume)) * tet_s.bfsSign[1]
            v_s1 = tet_s.vertices[:, 1]
            f_s1 = c_s1 * (r_j - v_s1)
            d_s1 = (tet_s.facesArea[1] / tet_s.volume) * tet_s.bfsSign[1]
            
            # n=2
            c_s2 = (tet_s.facesArea[2] / (3 * tet_s.volume)) * tet_s.bfsSign[2]
            v_s2 = tet_s.vertices[:, 2]
            f_s2 = c_s2 * (r_j - v_s2)
            d_s2 = (tet_s.facesArea[2] / tet_s.volume) * tet_s.bfsSign[2]
            
            # n=3
            c_s3 = (tet_s.facesArea[3] / (3 * tet_s.volume)) * tet_s.bfsSign[3]
            v_s3 = tet_s.vertices[:, 3]
            f_s3 = c_s3 * (r_j - v_s3)
            d_s3 = (tet_s.facesArea[3] / tet_s.volume) * tet_s.bfsSign[3]
            
            # n=4
            c_s4 = (tet_s.facesArea[4] / (3 * tet_s.volume)) * tet_s.bfsSign[4]
            v_s4 = tet_s.vertices[:, 4]
            f_s4 = c_s4 * (r_j - v_s4)
            d_s4 = (tet_s.facesArea[4] / tet_s.volume) * tet_s.bfsSign[4]
            
            # Dot Products Matrix (4x4)
            # Z_dot[m, n] = dot(f_tm, f_sn)
            # Manual unroll or use SMatrix?
            # SMatrix construction is cleaner but maybe verbose.
            # Let's do manual accumulation to avoid intermediate matrix allocation if possible,
            # or use SMatrix if it stays in registers.
            
            # Z_ts += (c1_ts * dot(f_t, f_s) + c2_ts * d_t * d_s) * factor
            
            # Precompute scalar factors
            fac_ts_1 = c1_ts * factor
            fac_ts_2 = c2_ts * factor
            
            fac_st_1 = c1_st * factor
            fac_st_2 = c2_st * factor
            
            # Unrolled Loop (4x4 = 16 ops)
            # m=1
            dot_11 = dot(f_t1, f_s1); div_11 = d_t1 * d_s1
            Z_ts[1, 1] += fac_ts_1 * dot_11 + fac_ts_2 * div_11
            Z_st[1, 1] += fac_st_1 * dot_11 + fac_st_2 * div_11
            
            dot_12 = dot(f_t1, f_s2); div_12 = d_t1 * d_s2
            Z_ts[1, 2] += fac_ts_1 * dot_12 + fac_ts_2 * div_12
            Z_st[2, 1] += fac_st_1 * dot_12 + fac_st_2 * div_12 # Note indices for Z_st
            
            dot_13 = dot(f_t1, f_s3); div_13 = d_t1 * d_s3
            Z_ts[1, 3] += fac_ts_1 * dot_13 + fac_ts_2 * div_13
            Z_st[3, 1] += fac_st_1 * dot_13 + fac_st_2 * div_13
            
            dot_14 = dot(f_t1, f_s4); div_14 = d_t1 * d_s4
            Z_ts[1, 4] += fac_ts_1 * dot_14 + fac_ts_2 * div_14
            Z_st[4, 1] += fac_st_1 * dot_14 + fac_st_2 * div_14
            
            # m=2
            dot_21 = dot(f_t2, f_s1); div_21 = d_t2 * d_s1
            Z_ts[2, 1] += fac_ts_1 * dot_21 + fac_ts_2 * div_21
            Z_st[1, 2] += fac_st_1 * dot_21 + fac_st_2 * div_21
            
            dot_22 = dot(f_t2, f_s2); div_22 = d_t2 * d_s2
            Z_ts[2, 2] += fac_ts_1 * dot_22 + fac_ts_2 * div_22
            Z_st[2, 2] += fac_st_1 * dot_22 + fac_st_2 * div_22
            
            dot_23 = dot(f_t2, f_s3); div_23 = d_t2 * d_s3
            Z_ts[2, 3] += fac_ts_1 * dot_23 + fac_ts_2 * div_23
            Z_st[3, 2] += fac_st_1 * dot_23 + fac_st_2 * div_23
            
            dot_24 = dot(f_t2, f_s4); div_24 = d_t2 * d_s4
            Z_ts[2, 4] += fac_ts_1 * dot_24 + fac_ts_2 * div_24
            Z_st[4, 2] += fac_st_1 * dot_24 + fac_st_2 * div_24
            
            # m=3
            dot_31 = dot(f_t3, f_s1); div_31 = d_t3 * d_s1
            Z_ts[3, 1] += fac_ts_1 * dot_31 + fac_ts_2 * div_31
            Z_st[1, 3] += fac_st_1 * dot_31 + fac_st_2 * div_31
            
            dot_32 = dot(f_t3, f_s2); div_32 = d_t3 * d_s2
            Z_ts[3, 2] += fac_ts_1 * dot_32 + fac_ts_2 * div_32
            Z_st[2, 3] += fac_st_1 * dot_32 + fac_st_2 * div_32
            
            dot_33 = dot(f_t3, f_s3); div_33 = d_t3 * d_s3
            Z_ts[3, 3] += fac_ts_1 * dot_33 + fac_ts_2 * div_33
            Z_st[3, 3] += fac_st_1 * dot_33 + fac_st_2 * div_33
            
            dot_34 = dot(f_t3, f_s4); div_34 = d_t3 * d_s4
            Z_ts[3, 4] += fac_ts_1 * dot_34 + fac_ts_2 * div_34
            Z_st[4, 3] += fac_st_1 * dot_34 + fac_st_2 * div_34
            
            # m=4
            dot_41 = dot(f_t4, f_s1); div_41 = d_t4 * d_s1
            Z_ts[4, 1] += fac_ts_1 * dot_41 + fac_ts_2 * div_41
            Z_st[1, 4] += fac_st_1 * dot_41 + fac_st_2 * div_41
            
            dot_42 = dot(f_t4, f_s2); div_42 = d_t4 * d_s2
            Z_ts[4, 2] += fac_ts_1 * dot_42 + fac_ts_2 * div_42
            Z_st[2, 4] += fac_st_1 * dot_42 + fac_st_2 * div_42
            
            dot_43 = dot(f_t4, f_s3); div_43 = d_t4 * d_s3
            Z_ts[4, 3] += fac_ts_1 * dot_43 + fac_ts_2 * div_43
            Z_st[3, 4] += fac_st_1 * dot_43 + fac_st_2 * div_43
            
            dot_44 = dot(f_t4, f_s4); div_44 = d_t4 * d_s4
            Z_ts[4, 4] += fac_ts_1 * dot_44 + fac_ts_2 * div_44
            Z_st[4, 4] += fac_st_1 * dot_44 + fac_st_2 * div_44

        end
    end
    
    return SMatrix(Z_ts), SMatrix(Z_st)
end

function precompute_vefie_basis(vefie::VEFIE, tetras::Vector{TetrahedraInfo{IT, FT, CT}}) where {IT, FT, CT}
    gq = vefie.gq_info
    gq_far = vefie.gq_far
    Nq = length(gq.weight)
    Nq_far = length(gq_far.weight)
    ntet = length(tetras)
    
    cache = Vector{TetBasisCache{CT, Nq, Nq_far}}(undef, ntet)
    
    Threads.@threads for i in 1:ntet
        tet = tetras[i]
        
        # Near (5-point)
        r_q = tet.vertices * gq.coordinate
        f_vals = MMatrix{4, Nq, SVector{3, CT}}(undef)
        div_f = MVector{4, CT}(undef)
        
        # Far (1-point)
        r_q_far = tet.vertices * gq_far.coordinate
        f_vals_far = MMatrix{4, Nq_far, SVector{3, CT}}(undef)
        
        for m in 1:4
            # Divergence (constant)
            div_val = (tet.facesArea[m] / tet.volume) * tet.bfsSign[m]
            div_f[m] = div_val
            
            # Basis function
            v_free = tet.vertices[:, m]
            const_val = (tet.facesArea[m] / (3 * tet.volume)) * tet.bfsSign[m]
            
            # Near
            for k in 1:Nq
                r = r_q[:, k]
                f_vals[m, k] = const_val * (r - v_free)
            end
            
            # Far
            for k in 1:Nq_far
                r = r_q_far[:, k]
                f_vals_far[m, k] = const_val * (r - v_free)
            end
        end
        
        cache[i] = TetBasisCache{CT, Nq, Nq_far}(SMatrix(r_q), SMatrix(f_vals), SVector(div_f), SMatrix(r_q_far), SMatrix(f_vals_far))
    end
    
    return cache
end

function vefie_mass_matrix_cached(vefie::VEFIE, tet::TetrahedraInfo, cache::TetBasisCache)
    FT = eltype(vefie.freq)
    CT = Complex{FT}
    M = @MMatrix zeros(CT, 4, 4)
    
    eps0 = 8.854187817e-12
    inv_eps = (1.0 - tet.κ) / eps0
    
    gq = vefie.gq_info
    omega = 2π * vefie.freq
    # Corrected: Added division by (im * omega) to match MoM_Kernels scaling
    factor_base = inv_eps / (im * omega) * tet.volume
    
    f_vals = cache.f_vals
    
    for k in 1:length(gq.weight)
        w = gq.weight[k]
        factor = w * factor_base
        
        for n in 1:4
            f_n = f_vals[n, k]
            for m in 1:4
                M[m, n] += dot(f_vals[m, k], f_n) * factor
            end
        end
    end
    
    return SMatrix(M)
end

function vefie_element_interaction_kernel(vefie::VEFIE, tet_t::TetrahedraInfo, tet_s::TetrahedraInfo, cache_t::TetBasisCache, cache_s::TetBasisCache)
    FT = eltype(vefie.freq)
    CT = Complex{FT}
    
    Z_ts = @MMatrix zeros(CT, 4, 4)
    
    # Constants
    k = vefie.k
    omega = 2π * vefie.freq
    mu0 = 4π * 1e-7
    eps0 = 8.854187817e-12
    
    # Material properties
    κs = tet_s.κ
    
    # Constants for terms (Scaled by 1/jw)
    c1_ts = im * omega * mu0 * κs
    c2_ts = 1.0 / (im * omega * eps0) * κs
    
    vol_factor = tet_t.volume * tet_s.volume
    
    # Adaptive Quadrature Check
    # Use centroids (1-point quadrature)
    C_t = cache_t.r_q_far[:, 1]
    C_s = cache_s.r_q_far[:, 1]
    dist = norm(C_t - C_s)
    
    # Threshold: 3.0 * (radius_t + radius_s)
    # radius approx cbrt(vol)
    rad_t = cbrt(tet_t.volume)
    rad_s = cbrt(tet_s.volume)
    threshold = 3.0 * (rad_t + rad_s)
    
    if dist > threshold
        # Use Far (1-point)
        gq = vefie.gq_far
        Nq = length(gq.weight)
        r_q_t = cache_t.r_q_far
        r_q_s = cache_s.r_q_far
        f_vals_t = cache_t.f_vals_far
        f_vals_s = cache_s.f_vals_far
    else
        # Use Near (5-point)
        gq = vefie.gq_info
        Nq = length(gq.weight)
        r_q_t = cache_t.r_q
        r_q_s = cache_s.r_q
        f_vals_t = cache_t.f_vals
        f_vals_s = cache_s.f_vals
    end
    
    div_f_t = cache_t.div_f
    div_f_s = cache_s.div_f
    
    # Double loop over quadrature points
    for j in 1:Nq # Source
        w_j = gq.weight[j]
        r_j = r_q_s[:, j]
        
        for i in 1:Nq # Test
            w_i = gq.weight[i]
            r_i = r_q_t[:, i]
            
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
            for m in 1:4
                f_m = f_vals_t[m, i]
                d_m = div_f_t[m]
                
                for n in 1:4
                    f_n = f_vals_s[n, j]
                    d_n = div_f_s[n]
                    
                    # Z_ts
                    term1 = c1_ts * dot(f_m, f_n)
                    term2 = c2_ts * d_m * d_n
                    Z_ts[m, n] += (term1 + term2) * factor
                end
            end
        end
    end
    
    return Z_ts
end

end
