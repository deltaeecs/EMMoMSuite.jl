module VEFIEModule

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..Kernels
using StaticArrays
using LinearAlgebra
using SparseArrays
using Base.Threads

include("FastExp.jl")
using .FastExpModule

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
struct VEFIE{FT<:AbstractFloat,CT<:Complex,N_GQ,N_GQ_FAR} <: AbstractIntegralOperator
    freq::FT
    k::FT
    eta::FT
    gq_info::GaussQuadratureInfoStruct{FT,N_GQ,4}
    gq_far::GaussQuadratureInfoStruct{FT,N_GQ_FAR,4}
    permittivities::Vector{CT}
    exp_table::FastExpTable{FT}  # Fast exponential lookup table
end

struct TetBasisCache{CT,NQ,NQ_FAR}
    r_q::SMatrix{3,NQ,Float64}
    f_vals::SMatrix{4,NQ,SVector{3,CT}}
    div_f::SVector{4,CT}

    r_q_far::SMatrix{3,NQ_FAR,Float64}
    f_vals_far::SMatrix{4,NQ_FAR,SVector{3,CT}}
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
    
    # Create fast exponential lookup table (20λ range, 10000 entries ≈ 80KB)
    exp_table = FastExpTable(k)

    return VEFIE{FT,Complex{FT},5,1}(freq, k, eta, gq_info, gq_far, permittivities, exp_table)
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
function assemble_impedance_matrix(
    vefie::VEFIE,
    basis::SWGBasis,
    permittivities::Vector{ComplexF64},
)
    FT = eltype(vefie.freq)
    CT = Complex{FT}

    nbf = num_basis(basis)
    Z = zeros(CT, nbf, nbf)

    # Precompute geometry and basis cache
    tetras = get_tetrahedra_info(basis.mesh, basis, permittivities)
    ntet = length(tetras)
    basis_cache = precompute_vefie_basis(vefie, tetras)

    # Symmetry exploitation (Legacy-style):
    #   Outer loop on test tet `it`, inner loop on source `js` from it to ntet.
    #   For js > it: compute Z_ts once, derive Z_st[j,i] = (κ_t/κ_s) * Z_ts[i,j]
    #   (exact for homogeneous; correct scaling for inhomogeneous media).
    #   Write BOTH Z[m,n] and Z[n,m] in one pass → ~2× speedup.
    #
    # Global SpinLock: held for ≤32 scalar writes (~32 ns) with 5 μs compute between
    # lock events → <1% contention probability with 4 threads.
    lockZ = SpinLock()
    next_it = Threads.Atomic{Int}(1)
    n_threads = Threads.nthreads()

    Threads.@threads for _ = 1:n_threads
        while true
            it = Threads.atomic_add!(next_it, 1)
            it > ntet && break

            tet_t = tetras[it]
            cache_t = basis_cache[it]

            # ── Self term (it == it) ────────────────────────────────────────────
            Z_self = vefie_element_interaction_kernel(vefie, tet_t, tet_t, cache_t, cache_t)
            M = vefie_mass_matrix_cached(vefie, tet_t, cache_t)
            lock(lockZ)
            @inbounds for i = 1:4
                m = tet_t.inBfsID[i]
                m == 0 && continue
                @inbounds for j = 1:4
                    n = tet_t.inBfsID[j]
                    n == 0 && continue
                    Z[m, n] += Z_self[i, j] + M[i, j]
                end
            end
            unlock(lockZ)

            # ── Upper tet triangle: js > it ─────────────────────────────────────
            # Z_st[j,i] = (κ_t / κ_s) * Z_ts[i,j]  (from κ-weighted Green's function)
            kappa_t = tet_t.κ
            for js = it+1:ntet
                tet_s = tetras[js]
                cache_s = basis_cache[js]
                Z_ts = vefie_element_interaction_kernel(vefie, tet_t, tet_s, cache_t, cache_s)
                kappa_ratio = kappa_t / tet_s.κ   # scalar (CT), precomputed outside lock

                lock(lockZ)
                @inbounds for i = 1:4
                    m = tet_t.inBfsID[i]
                    m == 0 && continue
                    @inbounds for j = 1:4
                        n = tet_s.inBfsID[j]
                        n == 0 && continue
                        zval = Z_ts[i, j]
                        Z[m, n] += zval                          # Z_ts: test=t, source=s
                        Z[n, m] += kappa_ratio * zval            # Z_st[j,i] by symmetry
                    end
                end
                unlock(lockZ)
            end
        end
    end

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
    f_vals = MVector{4,SVector{3,CT}}(undef)

    for k = 1:length(gq.weight)
        w = gq.weight[k]
        r = r_q[:, k]
        factor = w * factor_base

        for m = 1:4
            v_free = tet.vertices[:, m]
            rho = r - v_free
            val = (tet.facesArea[m] / (3 * tet.volume)) * rho
            val *= tet.bfsSign[m]
            f_vals[m] = val
        end

        for n = 1:4
            f_n = f_vals[n]
            for m = 1:4
                M[m, n] += dot(f_vals[m], f_n) * factor
            end
        end
    end

    return SMatrix(M)
end

function vefie_element_interaction(
    vefie::VEFIE{FT,CT,N_GQ,N_GQ_FAR},
    tet_t::TetrahedraInfo,
    tet_s::TetrahedraInfo,
) where {FT,CT,N_GQ,N_GQ_FAR}
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
    c2_ts = (1 / (im * omega * eps0)) * κs

    c1_st = -im * omega * mu0 * κt
    c2_st = (1 / (im * omega * eps0)) * κt

    vol_factor = tet_t.volume * tet_s.volume

    # Double loop over quadrature points
    for j = 1:Nq # Source
        w_j = gq.weight[j]
        r_j = r_q_s[:, j]

        for i = 1:Nq # Test
            w_i = gq.weight[i]
            r_i = r_q_t[:, i]

            # Green's function (using FastExp lookup table)
            R_vec = r_i - r_j
            R = norm(R_vec)
            
            G = fast_green_func(vefie.exp_table, R)

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
            dot_11 = dot(f_t1, f_s1)
            div_11 = d_t1 * d_s1
            Z_ts[1, 1] += fac_ts_1 * dot_11 + fac_ts_2 * div_11
            Z_st[1, 1] += fac_st_1 * dot_11 + fac_st_2 * div_11

            dot_12 = dot(f_t1, f_s2)
            div_12 = d_t1 * d_s2
            Z_ts[1, 2] += fac_ts_1 * dot_12 + fac_ts_2 * div_12
            Z_st[2, 1] += fac_st_1 * dot_12 + fac_st_2 * div_12 # Note indices for Z_st

            dot_13 = dot(f_t1, f_s3)
            div_13 = d_t1 * d_s3
            Z_ts[1, 3] += fac_ts_1 * dot_13 + fac_ts_2 * div_13
            Z_st[3, 1] += fac_st_1 * dot_13 + fac_st_2 * div_13

            dot_14 = dot(f_t1, f_s4)
            div_14 = d_t1 * d_s4
            Z_ts[1, 4] += fac_ts_1 * dot_14 + fac_ts_2 * div_14
            Z_st[4, 1] += fac_st_1 * dot_14 + fac_st_2 * div_14

            # m=2
            dot_21 = dot(f_t2, f_s1)
            div_21 = d_t2 * d_s1
            Z_ts[2, 1] += fac_ts_1 * dot_21 + fac_ts_2 * div_21
            Z_st[1, 2] += fac_st_1 * dot_21 + fac_st_2 * div_21

            dot_22 = dot(f_t2, f_s2)
            div_22 = d_t2 * d_s2
            Z_ts[2, 2] += fac_ts_1 * dot_22 + fac_ts_2 * div_22
            Z_st[2, 2] += fac_st_1 * dot_22 + fac_st_2 * div_22

            dot_23 = dot(f_t2, f_s3)
            div_23 = d_t2 * d_s3
            Z_ts[2, 3] += fac_ts_1 * dot_23 + fac_ts_2 * div_23
            Z_st[3, 2] += fac_st_1 * dot_23 + fac_st_2 * div_23

            dot_24 = dot(f_t2, f_s4)
            div_24 = d_t2 * d_s4
            Z_ts[2, 4] += fac_ts_1 * dot_24 + fac_ts_2 * div_24
            Z_st[4, 2] += fac_st_1 * dot_24 + fac_st_2 * div_24

            # m=3
            dot_31 = dot(f_t3, f_s1)
            div_31 = d_t3 * d_s1
            Z_ts[3, 1] += fac_ts_1 * dot_31 + fac_ts_2 * div_31
            Z_st[1, 3] += fac_st_1 * dot_31 + fac_st_2 * div_31

            dot_32 = dot(f_t3, f_s2)
            div_32 = d_t3 * d_s2
            Z_ts[3, 2] += fac_ts_1 * dot_32 + fac_ts_2 * div_32
            Z_st[2, 3] += fac_st_1 * dot_32 + fac_st_2 * div_32

            dot_33 = dot(f_t3, f_s3)
            div_33 = d_t3 * d_s3
            Z_ts[3, 3] += fac_ts_1 * dot_33 + fac_ts_2 * div_33
            Z_st[3, 3] += fac_st_1 * dot_33 + fac_st_2 * div_33

            dot_34 = dot(f_t3, f_s4)
            div_34 = d_t3 * d_s4
            Z_ts[3, 4] += fac_ts_1 * dot_34 + fac_ts_2 * div_34
            Z_st[4, 3] += fac_st_1 * dot_34 + fac_st_2 * div_34

            # m=4
            dot_41 = dot(f_t4, f_s1)
            div_41 = d_t4 * d_s1
            Z_ts[4, 1] += fac_ts_1 * dot_41 + fac_ts_2 * div_41
            Z_st[1, 4] += fac_st_1 * dot_41 + fac_st_2 * div_41

            dot_42 = dot(f_t4, f_s2)
            div_42 = d_t4 * d_s2
            Z_ts[4, 2] += fac_ts_1 * dot_42 + fac_ts_2 * div_42
            Z_st[2, 4] += fac_st_1 * dot_42 + fac_st_2 * div_42

            dot_43 = dot(f_t4, f_s3)
            div_43 = d_t4 * d_s3
            Z_ts[4, 3] += fac_ts_1 * dot_43 + fac_ts_2 * div_43
            Z_st[3, 4] += fac_st_1 * dot_43 + fac_st_2 * div_43

            dot_44 = dot(f_t4, f_s4)
            div_44 = d_t4 * d_s4
            Z_ts[4, 4] += fac_ts_1 * dot_44 + fac_ts_2 * div_44
            Z_st[4, 4] += fac_st_1 * dot_44 + fac_st_2 * div_44

        end
    end

    return SMatrix(Z_ts), SMatrix(Z_st)
end

function precompute_vefie_basis(
    vefie::VEFIE,
    tetras::Vector{TetrahedraInfo{IT,FT,CT}},
) where {IT,FT,CT}
    gq = vefie.gq_info
    gq_far = vefie.gq_far
    Nq = length(gq.weight)
    Nq_far = length(gq_far.weight)
    ntet = length(tetras)

    cache = Vector{TetBasisCache{CT,Nq,Nq_far}}(undef, ntet)

    Threads.@threads for i = 1:ntet
        tet = tetras[i]

        # Near (5-point)
        r_q = tet.vertices * gq.coordinate
        f_vals = MMatrix{4,Nq,SVector{3,CT}}(undef)
        div_f = MVector{4,CT}(undef)

        # Far (1-point)
        r_q_far = tet.vertices * gq_far.coordinate
        f_vals_far = MMatrix{4,Nq_far,SVector{3,CT}}(undef)

        for m = 1:4
            # Divergence (constant)
            div_val = (tet.facesArea[m] / tet.volume) * tet.bfsSign[m]
            div_f[m] = div_val

            # Basis function
            v_free = tet.vertices[:, m]
            const_val = (tet.facesArea[m] / (3 * tet.volume)) * tet.bfsSign[m]

            # Near
            for k = 1:Nq
                r = r_q[:, k]
                f_vals[m, k] = const_val * (r - v_free)
            end

            # Far
            for k = 1:Nq_far
                r = r_q_far[:, k]
                f_vals_far[m, k] = const_val * (r - v_free)
            end
        end

        cache[i] = TetBasisCache{CT,Nq,Nq_far}(
            SMatrix(r_q),
            SMatrix(f_vals),
            SVector(div_f),
            SMatrix(r_q_far),
            SMatrix(f_vals_far),
        )
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

    for k = 1:length(gq.weight)
        w = gq.weight[k]
        factor = w * factor_base

        for n = 1:4
            f_n = f_vals[n, k]
            for m = 1:4
                M[m, n] += dot(f_vals[m, k], f_n) * factor
            end
        end
    end

    return SMatrix(M)
end

function vefie_element_interaction_kernel(
    vefie::VEFIE,
    tet_t::TetrahedraInfo,
    tet_s::TetrahedraInfo,
    cache_t::TetBasisCache,
    cache_s::TetBasisCache,
)
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
    for j = 1:Nq # Source
        w_j = gq.weight[j]
        r_j = r_q_s[:, j]

        for i = 1:Nq # Test
            w_i = gq.weight[i]
            r_i = r_q_t[:, i]

            # Green's function (using FastExp lookup table)
            R_vec = r_i - r_j
            R = norm(R_vec)
            
            G = fast_green_func(vefie.exp_table, R)

            factor = w_i * w_j * vol_factor * G

            # Accumulate
            for m = 1:4
                f_m = f_vals_t[m, i]
                d_m = div_f_t[m]

                for n = 1:4
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

# ============================================================================
# PWC (Piecewise Constant) Assembly for VEFIE
# ============================================================================

"""
    assemble_impedance_matrix(vefie::VEFIE, basis::PWCBasis)

Assemble the impedance matrix Z for the VEFIE using PWC basis functions.
Each tetrahedron contributes 3 DOFs (x, y, z components).
The interaction kernel is the dyadic Green's function L operator:
    (k²I + ∇∇) G(R) 
"""
function assemble_impedance_matrix(vefie::VEFIE, basis::PWCBasis)
    return assemble_impedance_matrix(vefie, basis, vefie.permittivities)
end

"""
    assemble_impedance_matrix(vefie::VEFIE, basis::PWCBasis, permittivities)

Assemble the VEFIE impedance matrix for PWC basis functions.

The matrix element between test tet t (component i) and source tet s (component j) is:

    Z[3(t-1)+i, 3(s-1)+j] = (jη₀/k) * κₛ * V_t * V_s * 
        Σ_gi Σ_gj w_i w_j * G₀(R) * L_dyad[i,j]  +  δ_ts * δ_ij * V/(jωε)

where L_dyad is the dyadic L operator:
    L[i,j] = (δ_ij - R̂_i R̂_j) k² - (δ_ij - 3R̂_i R̂_j)(jk + 1/R)/R

# Legacy Parity
Matches `MoM_Kernels` `EFIEPWCTetra.jl` `impedancemat4VIE!` with `discreteVar = "D"`.
"""
function assemble_impedance_matrix(
    vefie::VEFIE,
    basis::PWCBasis,
    permittivities::Vector{ComplexF64},
)
    FT = eltype(vefie.freq)
    CT = Complex{FT}

    nbf = num_basis(basis)
    Z = zeros(CT, nbf, nbf)

    # Precompute geometry
    tetras = get_tetrahedra_info(basis.mesh, basis, permittivities)
    ntet = length(tetras)

    # Constants
    k = vefie.k
    k² = k^2
    jk = im * k
    omega = 2π * vefie.freq
    mu0 = 4π * 1e-7
    eps0 = 8.854187817e-12
    eta0 = sqrt(mu0 / eps0)
    # jη₀/k = j/(ωε₀) 
    Jη₀divK = im * eta0 / k
    div4π = 1.0 / (4π)

    # Progress tracking
    progress_counter = Threads.Atomic{Int}(0)
    next_idx = Threads.Atomic{Int}(1)

    # Locks for thread safety (one per column block)
    col_locks = [SpinLock() for _ = 1:ntet]

    println("VEFIE-PWC Assembly: $ntet tetrahedra, $nbf unknowns. (Symmetric + Threading)")

    # Precompute quadrature points for all tetrahedra
    gq = vefie.gq_info
    gq_far = vefie.gq_far
    Nq = length(gq.weight)
    Nq_far = length(gq_far.weight)

    # Precompute GQ points
    rq_near = Vector{Matrix{FT}}(undef, ntet)
    rq_far_pts = Vector{Matrix{FT}}(undef, ntet)
    Threads.@threads for i = 1:ntet
        rq_near[i] = tetras[i].vertices * gq.coordinate
        rq_far_pts[i] = tetras[i].vertices * gq_far.coordinate
    end

    # Dynamic scheduling: loop over source tets
    Threads.@threads for _ = 1:Threads.nthreads()
        # Thread-local 3×3 buffer for element interaction
        Z_ts_buf = zeros(CT, 3, 3)

        while true
            ti = Threads.atomic_add!(next_idx, 1)
            if ti > ntet
                break
            end

            tet_t = tetras[ti]
            κₜ = tet_t.κ

            for sj = ti:ntet
                tet_s = tetras[sj]
                κₛ = tet_s.κ

                # Compute distance between centers
                dist_ts = norm(tet_t.center - tet_s.center)

                # Adaptive quadrature
                rad_t = cbrt(tet_t.volume)
                rad_s = cbrt(tet_s.volume)
                threshold = 3.0 * (rad_t + rad_s)

                if dist_ts > threshold
                    # Far field: 1-point rule
                    _pwc_dyad_kernel!(
                        Z_ts_buf,
                        tet_t,
                        tet_s,
                        rq_far_pts[ti],
                        rq_far_pts[sj],
                        gq_far,
                        gq_far,
                        Nq_far,
                        Nq_far,
                        k,
                        k²,
                        jk,
                        Jη₀divK,
                        div4π,
                    )
                else
                    # Near field: 5-point rule
                    _pwc_dyad_kernel!(
                        Z_ts_buf,
                        tet_t,
                        tet_s,
                        rq_near[ti],
                        rq_near[sj],
                        gq,
                        gq,
                        Nq,
                        Nq,
                        k,
                        k²,
                        jk,
                        Jη₀divK,
                        div4π,
                    )
                end

                # Fill matrix using symmetry
                if ti == sj
                    # Self-term: Z_ts * κₜ + mass matrix
                    for ni = 1:3
                        n = tet_s.inBfsID[ni]
                        for mi = 1:3
                            m = tet_t.inBfsID[mi]
                            Z[m, n] = Z_ts_buf[mi, ni] * κₜ
                        end
                    end
                    # Add mass matrix diagonal: V/(jωε)
                    selfImp = 1.0 / (im * omega) / tet_t.ε * tet_t.volume
                    for ni = 1:3
                        n = tet_t.inBfsID[ni]
                        Z[n, n] += selfImp
                    end
                else
                    # Off-diagonal: use symmetry
                    for ni = 1:3
                        n = tet_s.inBfsID[ni]
                        for mi = 1:3
                            m = tet_t.inBfsID[mi]
                            Z[m, n] = Z_ts_buf[mi, ni] * κₛ
                            Z[n, m] = Z_ts_buf[mi, ni] * κₜ
                        end
                    end
                end
            end

            # Print progress
            c = Threads.atomic_add!(progress_counter, 1)
            if c % 10 == 0
                print("\rVEFIE-PWC Assembly: $c / $ntet source elements processed.")
            end
        end
    end
    println("\nVEFIE-PWC Assembly Completed.")

    return Z
end

"""
    _pwc_dyad_kernel!(Z_ts, vol_t, vol_s, rq_t, rq_s, gq, Nq, k, k², jk, Jη₀divK, div4π)

Compute the 3×3 dyadic interaction matrix between test volume `vol_t` and source volume `vol_s`
using the L operator: (k²I + ∇∇) G(R).

`vol_t` and `vol_s` can be any struct with a `.volume` field (TetrahedraInfo, HexahedraInfo).
Result stored in Z_ts (3×3 mutable buffer).
Does NOT include κ — caller must multiply by appropriate κ.

# Legacy Parity
Matches `EFIEOnTetrasPWC` / `EFIEOnHexasPWC` in `MoM_Kernels`.
"""
function _pwc_dyad_kernel!(
    Z_ts::Matrix{CT},
    vol_t,
    vol_s,
    rq_t::Matrix{FT},
    rq_s::Matrix{FT},
    gq_t,
    gq_s,
    Nq_t::Int,
    Nq_s::Int,
    k::FT,
    k²::FT,
    jk::CT,
    Jη₀divK::CT,
    div4π::FT,
) where {FT,CT}
    # Reset buffer
    fill!(Z_ts, zero(CT))

    dVtdVs = vol_t.volume * vol_s.volume
    Jη₀divKdVtdVs = Jη₀divK * dVtdVs

    # Double loop over quadrature points
    @inbounds for gj = 1:Nq_s
        rgj = @view rq_s[:, gj]

        for gi = 1:Nq_t
            rgi = @view rq_t[:, gi]

            # Distance vector
            Rx = rgi[1] - rgj[1]
            Ry = rgi[2] - rgj[2]
            Rz = rgi[3] - rgj[3]
            R = sqrt(Rx^2 + Ry^2 + Rz^2)

            if R < 1e-10
                continue  # Skip self-coincident points
            end

            divR = 1.0 / R
            # (jk + 1/R) / R
            jkplusR1stdivR1st = (jk + divR) * divR

            # R̂ components
            R̂x = Rx * divR
            R̂y = Ry * divR
            R̂z = Rz * divR

            # Green's function × quadrature weights
            GR = exp(-jk * R) * div4π * divR * gq_t.weight[gi] * gq_s.weight[gj]

            # Combined constant
            fac = Jη₀divKdVtdVs * GR

            # R̂R̂ dyad components (symmetric)
            RR11 = R̂x * R̂x
            RR12 = R̂x * R̂y
            RR13 = R̂x * R̂z
            RR22 = R̂y * R̂y
            RR23 = R̂y * R̂z
            RR33 = R̂z * R̂z

            # Diagonal terms (i == j): (1 - R̂ᵢR̂ⱼ)*k² - (1 - 3R̂ᵢR̂ⱼ)*(jk+1/R)/R
            Z_ts[1, 1] += fac * ((1 - RR11) * k² - (1 - 3RR11) * jkplusR1stdivR1st)
            Z_ts[2, 2] += fac * ((1 - RR22) * k² - (1 - 3RR22) * jkplusR1stdivR1st)
            Z_ts[3, 3] += fac * ((1 - RR33) * k² - (1 - 3RR33) * jkplusR1stdivR1st)

            # Off-diagonal terms (i ≠ j): -R̂ᵢR̂ⱼ*k² + 3R̂ᵢR̂ⱼ*(jk+1/R)/R
            offdiag12 = fac * (-RR12 * k² + 3RR12 * jkplusR1stdivR1st)
            offdiag13 = fac * (-RR13 * k² + 3RR13 * jkplusR1stdivR1st)
            offdiag23 = fac * (-RR23 * k² + 3RR23 * jkplusR1stdivR1st)

            Z_ts[1, 2] += offdiag12
            Z_ts[2, 1] += offdiag12
            Z_ts[1, 3] += offdiag13
            Z_ts[3, 1] += offdiag13
            Z_ts[2, 3] += offdiag23
            Z_ts[3, 2] += offdiag23
        end
    end

    return nothing
end

# ============================================================================
# PWC Hexahedra Assembly for VEFIE
# ============================================================================

"""
    assemble_impedance_matrix(vefie::VEFIE, basis::PWCHexBasis)

Assemble the impedance matrix Z for the VEFIE using PWC basis functions on hexahedra.
Each hexahedron contributes 3 DOFs (x, y, z components).
The interaction kernel is the same dyadic L operator as PWC on tetra.

# Legacy Parity
Matches `MoM_Kernels` `EFIEPWCHexa.jl` `impedancemat4VIE!` with `discreteVar = "D"`.
"""
function assemble_impedance_matrix(vefie::VEFIE, basis::PWCHexBasis)
    return assemble_impedance_matrix(vefie, basis, vefie.permittivities)
end

function assemble_impedance_matrix(
    vefie::VEFIE,
    basis::PWCHexBasis,
    permittivities::Vector{ComplexF64},
)
    FT = eltype(vefie.freq)
    CT = Complex{FT}

    nbf = num_basis(basis)
    Z = zeros(CT, nbf, nbf)

    # Precompute geometry
    hexas = get_hexahedra_info(basis.mesh, basis, permittivities)
    nhex = length(hexas)

    # Constants
    k = vefie.k
    k² = k^2
    jk = im * k
    omega = 2π * vefie.freq
    mu0 = 4π * 1e-7
    eps0 = 8.854187817e-12
    eta0 = sqrt(mu0 / eps0)
    Jη₀divK = im * eta0 / k
    div4π = 1.0 / (4π)

    # GQ info for hexahedra (8-point near, 1-point far)
    gq_hex = GaussQuadratureInfo(:Hexahedron, 8, FT)
    gq_hex_far = GaussQuadratureInfo(:Hexahedron, 1, FT)
    Nq_hex = length(gq_hex.weight)
    Nq_hex_far = length(gq_hex_far.weight)

    # Progress tracking
    progress_counter = Threads.Atomic{Int}(0)
    next_idx = Threads.Atomic{Int}(1)

    println("VEFIE-PWC-Hexa Assembly: $nhex hexahedra, $nbf unknowns.")

    # Precompute GQ points for all hexahedra (coordinates × shape functions)
    rq_near = Vector{Matrix{FT}}(undef, nhex)
    rq_far_pts = Vector{Matrix{FT}}(undef, nhex)
    Threads.@threads for i = 1:nhex
        rq_near[i] = hexas[i].vertices * gq_hex.coordinate
        rq_far_pts[i] = hexas[i].vertices * gq_hex_far.coordinate
    end

    # Dynamic scheduling
    Threads.@threads for _ = 1:Threads.nthreads()
        Z_ts_buf = zeros(CT, 3, 3)

        while true
            ti = Threads.atomic_add!(next_idx, 1)
            if ti > nhex
                break
            end

            hex_t = hexas[ti]
            κₜ = hex_t.κ

            for sj = ti:nhex
                hex_s = hexas[sj]
                κₛ = hex_s.κ

                # Distance between centers
                dist_ts = norm(hex_t.center - hex_s.center)

                # Adaptive quadrature threshold
                rad_t = cbrt(hex_t.volume)
                rad_s = cbrt(hex_s.volume)
                threshold = 3.0 * (rad_t + rad_s)

                if dist_ts > threshold
                    _pwc_dyad_kernel!(
                        Z_ts_buf,
                        hex_t,
                        hex_s,
                        rq_far_pts[ti],
                        rq_far_pts[sj],
                        gq_hex_far,
                        gq_hex_far,
                        Nq_hex_far,
                        Nq_hex_far,
                        k,
                        k²,
                        jk,
                        Jη₀divK,
                        div4π,
                    )
                else
                    _pwc_dyad_kernel!(
                        Z_ts_buf,
                        hex_t,
                        hex_s,
                        rq_near[ti],
                        rq_near[sj],
                        gq_hex,
                        gq_hex,
                        Nq_hex,
                        Nq_hex,
                        k,
                        k²,
                        jk,
                        Jη₀divK,
                        div4π,
                    )
                end

                if ti == sj
                    # Self-term: Zts * κₜ + mass diagonal
                    for ni = 1:3, mi = 1:3
                        m = hex_t.inBfsID[mi]
                        n = hex_s.inBfsID[ni]
                        Z[m, n] = Z_ts_buf[mi, ni] * κₜ
                    end
                    selfImp = 1.0 / (im * omega) / hex_t.ε * hex_t.volume
                    for ni = 1:3
                        n = hex_t.inBfsID[ni]
                        Z[n, n] += selfImp
                    end
                else
                    for ni = 1:3, mi = 1:3
                        m = hex_t.inBfsID[mi]
                        n = hex_s.inBfsID[ni]
                        Z[m, n] = Z_ts_buf[mi, ni] * κₛ
                        Z[n, m] = Z_ts_buf[mi, ni] * κₜ
                    end
                end
            end

            c = Threads.atomic_add!(progress_counter, 1)
            if c % 10 == 0
                print("\rVEFIE-PWC-Hexa: $c / $nhex processed.")
            end
        end
    end
    println("\nVEFIE-PWC-Hexa Assembly Completed.")

    return Z
end

# ============================================================================
# RBF (Rooftop Basis Function) Hexahedra Assembly for VEFIE
# ============================================================================

"""
    assemble_impedance_matrix(vefie::VEFIE, basis::RBFBasis)

Assemble the impedance matrix Z for the VEFIE using RBF basis functions on hexahedra.
Each hexahedron contributes 6 DOFs (one per face).

The matrix element Z_mn involves 6 integral terms (F₁-F₆):
  Z_mn = C₁F₁ + κₛC₃(-k²F₂ + F₃) - δκₙC₃F₄ - κₛC₃F₅ + δκₙC₃F₆

where:
  C₁ = (1/jω)(1/Vε)AₘAₙ  (mass matrix, self only)
  C₃ = -jη₀/(4πk) AₘAₙ   (Green's function coupling)

# Legacy Parity
Matches `MoM_Kernels` `EFIERBFHexa.jl`.
"""
function assemble_impedance_matrix(vefie::VEFIE, basis::RBFBasis)
    return assemble_impedance_matrix(vefie, basis, vefie.permittivities)
end

function assemble_impedance_matrix(
    vefie::VEFIE,
    basis::RBFBasis,
    permittivities::Vector{ComplexF64},
)
    FT = eltype(vefie.freq)
    CT = Complex{FT}

    nbf = num_basis(basis)
    Z = zeros(CT, nbf, nbf)

    # Precompute geometry
    hexas = get_hexahedra_info(basis.mesh, basis, permittivities)
    nhex = length(hexas)

    # Constants
    k = vefie.k
    k² = k^2
    omega = 2π * vefie.freq
    mu0 = 4π * 1e-7
    eps0 = 8.854187817e-12
    eta0 = sqrt(mu0 / eps0)
    mJη₀div4πK = -im * eta0 / (4π * k)
    divJω = 1.0 / (im * omega)
    div4π = 1.0 / (4π)
    jk = im * k

    # GQ info for hexahedra and quadrangles
    gq_hex = GaussQuadratureInfo(:Hexahedron, 8, FT)
    gq_hex_far = GaussQuadratureInfo(:Hexahedron, 1, FT)
    gq_quad = GaussQuadratureInfo(:Quadrangle, 4, FT)
    Nq_hex = length(gq_hex.weight)
    Nq_hex_far = length(gq_hex_far.weight)
    Nq_quad = length(gq_quad.weight)

    # Precompute GQ points for all hexahedra
    rq_near = Vector{Matrix{FT}}(undef, nhex)
    rq_far = Vector{Matrix{FT}}(undef, nhex)
    Threads.@threads for i = 1:nhex
        rq_near[i] = hexas[i].vertices * gq_hex.coordinate
        rq_far[i] = hexas[i].vertices * gq_hex_far.coordinate
    end

    # Build GQ 3D → 2D index map for free-end computation
    gq3d_map = construct_gq3d_index_map(2)  # n1d=2 for 8-point hex GQ

    # Progress tracking
    progress_counter = Threads.Atomic{Int}(0)
    next_idx = Threads.Atomic{Int}(1)
    lockZ = SpinLock()

    println("VEFIE-RBF Assembly: $nhex hexahedra, $nbf unknowns.")

    # Dynamic scheduling
    Threads.@threads for _ = 1:Threads.nthreads()
        Zts = zeros(CT, 6, 6)
        Zst = zeros(CT, 6, 6)

        while true
            it = Threads.atomic_add!(next_idx, 1)
            if it > nhex
                break
            end

            hex_t = hexas[it]
            κt = hex_t.κ

            for js = it:nhex
                hex_s = hexas[js]
                κs = hex_s.κ

                dist_ts = norm(hex_t.center - hex_s.center)
                rad_t = cbrt(hex_t.volume)
                rad_s = cbrt(hex_s.volume)
                threshold = 3.0 * (rad_t + rad_s)

                if it == js
                    _rbf_self_kernel!(
                        Zts,
                        hex_t,
                        rq_near[it],
                        gq_hex,
                        Nq_hex,
                        gq_quad,
                        Nq_quad,
                        gq3d_map,
                        k²,
                        jk,
                        mJη₀div4πK,
                        divJω,
                        div4π,
                    )
                    lock(lockZ)
                    for ni = 1:6, mi = 1:6
                        m = hex_t.inBfsID[mi]
                        n = hex_s.inBfsID[ni]
                        Z[m, n] += Zts[mi, ni]
                    end
                    unlock(lockZ)
                elseif dist_ts > threshold
                    _rbf_far_kernel!(
                        Zts,
                        Zst,
                        hex_t,
                        hex_s,
                        rq_far[it],
                        rq_far[js],
                        gq_hex_far,
                        Nq_hex_far,
                        gq_quad,
                        Nq_quad,
                        k²,
                        jk,
                        mJη₀div4πK,
                        div4π,
                    )
                    lock(lockZ)
                    for ni = 1:6, mi = 1:6
                        m = hex_t.inBfsID[mi]
                        n = hex_s.inBfsID[ni]
                        Z[m, n] += Zts[mi, ni]
                        Z[n, m] += Zst[ni, mi]
                    end
                    unlock(lockZ)
                else
                    _rbf_near_kernel!(
                        Zts,
                        Zst,
                        hex_t,
                        hex_s,
                        rq_near[it],
                        rq_near[js],
                        gq_hex,
                        Nq_hex,
                        gq_quad,
                        Nq_quad,
                        k²,
                        jk,
                        mJη₀div4πK,
                        div4π,
                    )
                    lock(lockZ)
                    for ni = 1:6, mi = 1:6
                        m = hex_t.inBfsID[mi]
                        n = hex_s.inBfsID[ni]
                        Z[m, n] += Zts[mi, ni]
                        Z[n, m] += Zst[ni, mi]
                    end
                    unlock(lockZ)
                end
            end

            c = Threads.atomic_add!(progress_counter, 1)
            if c % 10 == 0
                print("\rVEFIE-RBF: $c / $nhex processed.")
            end
        end
    end
    println("\nVEFIE-RBF Assembly Completed.")

    return Z
end

"""
Green's function helper: exp(-jkR)/(4πR).
"""
@inline function _greenfunc(rgi, rgj, jk::CT, div4π::FT) where {CT,FT}
    Rx = rgi[1] - rgj[1]
    Ry = rgi[2] - rgj[2]
    Rz = rgi[3] - rgj[3]
    R = sqrt(Rx^2 + Ry^2 + Rz^2)
    if R < 1e-10
        return zero(CT)
    end
    return exp(-jk * R) * div4π / R
end

"""
Get GQ points on a quadrilateral face of a hexahedron.
Returns 3×Nq matrix of physical coordinates.
"""
function _get_quad_gq_points(hex_info::HexahedraInfo, face_idx::Int, gq_quad, Nq_quad)
    face = hex_info.faces[face_idx]
    # face.vertices is 3×4 matrix of face corner coordinates
    # gq_quad.coordinate is 4×Nq_quad (bilinear shape function values)
    return face.vertices * gq_quad.coordinate
end

"""
RBF far-field kernel: compute 6×6 Zts and Zst for two non-overlapping distant hexahedra.
"""
function _rbf_far_kernel!(
    Zts::Matrix{CT},
    Zst::Matrix{CT},
    hex_t::HexahedraInfo,
    hex_s::HexahedraInfo,
    rq_t::Matrix{FT},
    rq_s::Matrix{FT},
    gq_hex,
    Nq_hex::Int,
    gq_quad,
    Nq_quad::Int,
    k²::FT,
    jk::CT,
    mJη₀div4πK::CT,
    div4π::FT,
) where {FT,CT}
    fill!(Zts, zero(CT))
    fill!(Zst, zero(CT))

    κt = hex_t.κ
    κs = hex_s.κ
    facest = hex_t.faces
    facess = hex_s.faces

    # --- Precompute gw matrix (Green's function × weights) ---
    gw = zeros(CT, Nq_hex, Nq_hex)
    @inbounds for gj = 1:Nq_hex
        rgj = @view rq_s[:, gj]
        for gi = 1:Nq_hex
            rgi = @view rq_t[:, gi]
            gw[gi, gj] = _greenfunc(rgi, rgj, jk, div4π) * gq_hex.weight[gi] * gq_hex.weight[gj]
        end
    end

    # F₃ (scalar Green integral, independent of basis functions)
    F₃ = sum(gw)

    # --- F₄ (source face − test volume) ---
    F₄s = zeros(CT, 6)
    @inbounds for ni = 1:6
        arean = hex_s.facesArea[ni]
        δκn = facess[ni].δκ
        isbdn = facess[ni].isbd
        if isbdn || ((δκn != 0) && (arean > 0))
            rq_face_s = _get_quad_gq_points(hex_s, ni, gq_quad, Nq_quad)
            gtemp = zero(CT)
            for gj = 1:Nq_quad
                rgj = @view rq_face_s[:, gj]
                for gi = 1:Nq_hex
                    rgi = @view rq_t[:, gi]
                    gtemp +=
                        _greenfunc(rgi, rgj, jk, div4π) * gq_hex.weight[gi] * gq_quad.weight[gj]
                end
            end
            F₄s[ni] = gtemp
        end
    end

    # --- F₅ (test face − source volume) ---
    F₅t = zeros(CT, 6)
    @inbounds for mi = 1:6
        aream = hex_t.facesArea[mi]
        δκm = facest[mi].δκ
        isbdm = facest[mi].isbd
        if isbdm || ((δκm != 0) && (aream > 0))
            rq_face_t = _get_quad_gq_points(hex_t, mi, gq_quad, Nq_quad)
            gtemp = zero(CT)
            for gj = 1:Nq_hex
                rgj = @view rq_s[:, gj]
                for gi = 1:Nq_quad
                    rgi = @view rq_face_t[:, gi]
                    gtemp +=
                        _greenfunc(rgi, rgj, jk, div4π) * gq_quad.weight[gi] * gq_hex.weight[gj]
                end
            end
            F₅t[mi] = gtemp
        end
    end

    # --- Loop over basis function pairs to compute F₂ and accumulate ---
    @inbounds for ni = 1:6
        arean = hex_s.facesArea[ni]
        δκn = facess[ni].δκ
        isbdn = facess[ni].isbd
        # Precompute free-end coords for source basis ni
        freeVns_ni = get_free_vns(hex_s, ni, gq_hex.coordinate)

        for mi = 1:6
            aream = hex_t.facesArea[mi]
            δκm = facest[mi].δκ
            isbdm = facest[mi].isbd
            aman = aream * arean
            C₃ = mJη₀div4πK * aman
            freeVms_mi = get_free_vns(hex_t, mi, gq_hex.coordinate)

            # F₂ (ρₘ·ρₙ × G)
            F₂ = zero(CT)
            for gj = 1:Nq_hex
                rgj = @view rq_s[:, gj]
                idn = gq3d_to_face2d_idx(gj, ni, 2)
                ρnx = rgj[1] - freeVns_ni[1, idn]
                ρny = rgj[2] - freeVns_ni[2, idn]
                ρnz = rgj[3] - freeVns_ni[3, idn]
                for gi = 1:Nq_hex
                    rgi = @view rq_t[:, gi]
                    idm = gq3d_to_face2d_idx(gi, mi, 2)
                    ρmx = rgi[1] - freeVms_mi[1, idm]
                    ρmy = rgi[2] - freeVms_mi[2, idm]
                    ρmz = rgi[3] - freeVms_mi[3, idm]
                    ρdot = ρmx * ρnx + ρmy * ρny + ρmz * ρnz
                    F₂ += ρdot * gw[gi, gj]
                end
            end

            # CF23 (symmetric part)
            CF23 = C₃ * (-k² * F₂ + F₃)
            Zmn = κs * CF23
            Znm = κt * CF23

            # F₄ term
            (δκn != 0) && (arean > 0) && (Zmn -= δκn * C₃ * F₄s[ni])
            isbdn && (Znm -= κt * C₃ * F₄s[ni])

            # F₅ term
            isbdm && (Zmn -= κs * C₃ * F₅t[mi])
            (δκm != 0) && (aream > 0) && (Znm -= δκm * C₃ * F₅t[mi])

            # F₆ term
            n_global = hex_s.inBfsID[ni]
            m_global = hex_t.inBfsID[mi]
            statem = isbdm && (δκn != 0) && (arean > 0)
            staten = isbdn && (δκm != 0) && (aream > 0)
            if statem || staten
                F₆ = zero(CT)
                rq_face_m = _get_quad_gq_points(hex_t, mi, gq_quad, Nq_quad)
                rq_face_n = _get_quad_gq_points(hex_s, ni, gq_quad, Nq_quad)
                for gj = 1:Nq_quad
                    rgj = @view rq_face_n[:, gj]
                    for gi = 1:Nq_quad
                        rgi = @view rq_face_m[:, gi]
                        F₆ +=
                            _greenfunc(rgi, rgj, jk, div4π) *
                            gq_quad.weight[gi] *
                            gq_quad.weight[gj]
                    end
                end
                C₃F₆ = C₃ * F₆
                statem && (Zmn += δκn * C₃F₆)
                staten && (Znm += δκm * C₃F₆)
            end

            Zts[mi, ni] = Zmn
            Zst[ni, mi] = Znm
        end
    end

    return nothing
end

"""
RBF near-field kernel: same as far but using higher-order GQ.
For simplicity, follows same structure as far kernel.
"""
function _rbf_near_kernel!(
    Zts::Matrix{CT},
    Zst::Matrix{CT},
    hex_t::HexahedraInfo,
    hex_s::HexahedraInfo,
    rq_t::Matrix{FT},
    rq_s::Matrix{FT},
    gq_hex,
    Nq_hex::Int,
    gq_quad,
    Nq_quad::Int,
    k²::FT,
    jk::CT,
    mJη₀div4πK::CT,
    div4π::FT,
) where {FT,CT}
    # Near field uses the same structure as far field but with higher-order GQ
    # (rq_t and rq_s already have the near-field GQ points)
    _rbf_far_kernel!(
        Zts,
        Zst,
        hex_t,
        hex_s,
        rq_t,
        rq_s,
        gq_hex,
        Nq_hex,
        gq_quad,
        Nq_quad,
        k²,
        jk,
        mJη₀div4πK,
        div4π,
    )
end

"""
RBF self kernel: compute 6×6 Ztt for a hexahedron with itself.
Includes F₁ (mass matrix), F₂-F₆ terms, and self-impedance C₁F₁.
"""
function _rbf_self_kernel!(
    Ztt::Matrix{CT},
    hex_t::HexahedraInfo,
    rq_t::Matrix{FT},
    gq_hex,
    Nq_hex::Int,
    gq_quad,
    Nq_quad::Int,
    gq3d_map,
    k²::FT,
    jk::CT,
    mJη₀div4πK::CT,
    divJω::CT,
    div4π::FT,
) where {FT,CT}
    fill!(Ztt, zero(CT))

    κt = hex_t.κ
    faces = hex_t.faces
    divVε = 1.0 / (hex_t.volume * hex_t.ε)

    # gw matrix (self: same hex for test and source)
    gw = zeros(CT, Nq_hex, Nq_hex)
    @inbounds for gj = 1:Nq_hex
        rgj = @view rq_t[:, gj]
        for gi = 1:Nq_hex
            rgi = @view rq_t[:, gi]
            gw[gi, gj] = _greenfunc(rgi, rgj, jk, div4π) * gq_hex.weight[gi] * gq_hex.weight[gj]
        end
    end
    F₃ = sum(gw)

    # F₄/F₅ (same for self-term since hex_t == hex_s)
    F₄s = zeros(CT, 6)
    @inbounds for ni = 1:6
        arean = hex_t.facesArea[ni]
        δκn = faces[ni].δκ
        isbdn = faces[ni].isbd
        if isbdn || ((δκn != 0) && (arean > 0))
            rq_face = _get_quad_gq_points(hex_t, ni, gq_quad, Nq_quad)
            gtemp = zero(CT)
            for gi = 1:Nq_hex
                rgi = @view rq_t[:, gi]
                for gj = 1:Nq_quad
                    rgj = @view rq_face[:, gj]
                    gtemp +=
                        _greenfunc(rgi, rgj, jk, div4π) * gq_hex.weight[gi] * gq_quad.weight[gj]
                end
            end
            F₄s[ni] = gtemp
        end
    end
    F₅t = F₄s  # Same for self-term

    # Precompute free-end coords for all faces
    freeVns = Vector{Matrix{FT}}(undef, 6)
    for fi = 1:6
        freeVns[fi] = get_free_vns(hex_t, fi, gq_hex.coordinate)
    end

    # Loop over basis function pairs (use upper triangle symmetry mi >= ni)
    @inbounds for ni = 1:6
        arean = hex_t.facesArea[ni]
        δκn = faces[ni].δκ
        isbdn = faces[ni].isbd

        for mi = ni:6
            aream = hex_t.facesArea[mi]
            δκm = faces[mi].δκ
            isbdm = faces[mi].isbd
            aman = aream * arean
            C₁ = divJω * divVε * aman
            C₃ = mJη₀div4πK * aman

            # F₁ (mass matrix) and F₂ (ρ·ρ' × G)
            F₁ = zero(FT)
            F₂ = zero(CT)
            for gi = 1:Nq_hex
                rgi = @view rq_t[:, gi]
                idm = gq3d_to_face2d_idx(gi, mi, 2)
                idn_self = gq3d_to_face2d_idx(gi, ni, 2)
                ρmx = rgi[1] - freeVns[mi][1, idm]
                ρmy = rgi[2] - freeVns[mi][2, idm]
                ρmz = rgi[3] - freeVns[mi][3, idm]
                ρnx_self = rgi[1] - freeVns[ni][1, idn_self]
                ρny_self = rgi[2] - freeVns[ni][2, idn_self]
                ρnz_self = rgi[3] - freeVns[ni][3, idn_self]
                ρmρn_self = ρmx * ρnx_self + ρmy * ρny_self + ρmz * ρnz_self
                F₁ += gq_hex.weight[gi] * ρmρn_self

                for gj = 1:Nq_hex
                    rgj = @view rq_t[:, gj]
                    idn = gq3d_to_face2d_idx(gj, ni, 2)
                    ρnx = rgj[1] - freeVns[ni][1, idn]
                    ρny = rgj[2] - freeVns[ni][2, idn]
                    ρnz = rgj[3] - freeVns[ni][3, idn]
                    ρdot = ρmx * ρnx + ρmy * ρny + ρmz * ρnz
                    F₂ += ρdot * gw[gi, gj]
                end
            end

            C₁F₁ = C₁ * F₁
            CF23 = C₃ * (-k² * F₂ + F₃)

            Zmn = C₁F₁ + κt * CF23
            Znm = C₁F₁ + κt * CF23

            # F₄/F₅ terms
            (δκn != 0) && (arean > 0) && (Zmn -= δκn * C₃ * F₄s[ni])
            isbdn && (Znm -= κt * C₃ * F₄s[ni])
            isbdm && (Zmn -= κt * C₃ * F₅t[mi])
            (δκm != 0) && (aream > 0) && (Znm -= δκm * C₃ * F₅t[mi])

            # F₆ terms
            m_global = hex_t.inBfsID[mi]
            n_global = hex_t.inBfsID[ni]
            statem = isbdm && (δκn != 0) && (arean > 0)
            staten = isbdn && (δκm != 0) && (aream > 0)
            if statem || staten
                F₆ = zero(CT)
                if m_global != n_global
                    rq_face_m = _get_quad_gq_points(hex_t, mi, gq_quad, Nq_quad)
                    rq_face_n = _get_quad_gq_points(hex_t, ni, gq_quad, Nq_quad)
                    for gj = 1:Nq_quad
                        rgj = @view rq_face_n[:, gj]
                        for gi = 1:Nq_quad
                            rgi = @view rq_face_m[:, gi]
                            F₆ +=
                                _greenfunc(rgi, rgj, jk, div4π) *
                                gq_quad.weight[gi] *
                                gq_quad.weight[gj]
                        end
                    end
                else
                    # Same face: use self-GQ (no singularity extraction, just skip R=0)
                    rq_face = _get_quad_gq_points(hex_t, mi, gq_quad, Nq_quad)
                    for gj = 1:Nq_quad
                        rgj = @view rq_face[:, gj]
                        for gi = 1:Nq_quad
                            rgi = @view rq_face[:, gi]
                            F₆ +=
                                _greenfunc(rgi, rgj, jk, div4π) *
                                gq_quad.weight[gi] *
                                gq_quad.weight[gj]
                        end
                    end
                end
                C₃F₆ = C₃ * F₆
                statem && (Zmn += δκn * C₃F₆)
                staten && (Znm += δκm * C₃F₆)
            end

            Ztt[mi, ni] = Zmn
            Ztt[ni, mi] = Znm
        end
    end

    return nothing
end

# ============================================================================
# Mixed Tetra-Hexa PWC Assembly for VEFIE
# ============================================================================

"""
    assemble_impedance_matrix(vefie::VEFIE, basis_tet::PWCBasis, basis_hex::PWCHexBasis)

Assemble the coupling impedance matrix between PWC on tetrahedra and PWC on hexahedra.
Returns the full matrix including both self-blocks and cross-blocks.

# Legacy Parity
Matches `MoM_Kernels` `EFIEPWCTetraHexa.jl`.
"""
function assemble_impedance_matrix(vefie::VEFIE, basis_tet::PWCBasis, basis_hex::PWCHexBasis)
    return assemble_impedance_matrix(vefie, basis_tet, basis_hex, vefie.permittivities)
end

function assemble_impedance_matrix(
    vefie::VEFIE,
    basis_tet::PWCBasis,
    basis_hex::PWCHexBasis,
    permittivities::Vector{ComplexF64},
)
    FT = eltype(vefie.freq)
    CT = Complex{FT}

    nbf_tet = num_basis(basis_tet)
    nbf_hex = num_basis(basis_hex)
    nbf = nbf_tet + nbf_hex
    Z = zeros(CT, nbf, nbf)

    # Constants
    k = vefie.k
    k² = k^2
    jk = im * k
    omega = 2π * vefie.freq
    mu0 = 4π * 1e-7
    eps0 = 8.854187817e-12
    eta0 = sqrt(mu0 / eps0)
    Jη₀divK = im * eta0 / k
    div4π = 1.0 / (4π)

    # Precompute geometry for tetrahedra
    tetras = get_tetrahedra_info(basis_tet.mesh, basis_tet, permittivities)
    ntet = length(tetras)

    # Precompute geometry for hexahedra
    hexas = get_hexahedra_info(basis_hex.mesh, basis_hex, permittivities)
    nhex = length(hexas)

    # GQ info
    gq_tet = vefie.gq_info
    gq_tet_far = vefie.gq_far
    gq_hex = GaussQuadratureInfo(:Hexahedron, 8, FT)
    gq_hex_far = GaussQuadratureInfo(:Hexahedron, 1, FT)
    Nq_tet = length(gq_tet.weight)
    Nq_tet_far = length(gq_tet_far.weight)
    Nq_hex = length(gq_hex.weight)
    Nq_hex_far = length(gq_hex_far.weight)

    # Precompute GQ points
    rq_tet_near = Vector{Matrix{FT}}(undef, ntet)
    rq_tet_far_pts = Vector{Matrix{FT}}(undef, ntet)
    Threads.@threads for i = 1:ntet
        rq_tet_near[i] = tetras[i].vertices * gq_tet.coordinate
        rq_tet_far_pts[i] = tetras[i].vertices * gq_tet_far.coordinate
    end

    rq_hex_near = Vector{Matrix{FT}}(undef, nhex)
    rq_hex_far_pts = Vector{Matrix{FT}}(undef, nhex)
    Threads.@threads for i = 1:nhex
        rq_hex_near[i] = hexas[i].vertices * gq_hex.coordinate
        rq_hex_far_pts[i] = hexas[i].vertices * gq_hex_far.coordinate
    end

    println("VEFIE-PWC Mixed Assembly: $ntet tetra + $nhex hexa, $nbf unknowns.")

    # 1. Assemble tetra-tetra block
    Z_tt = assemble_impedance_matrix(vefie, basis_tet, permittivities)
    Z[1:nbf_tet, 1:nbf_tet] .= Z_tt

    # 2. Assemble hexa-hexa block
    Z_hh = assemble_impedance_matrix(vefie, basis_hex, permittivities)
    Z[nbf_tet+1:nbf, nbf_tet+1:nbf] .= Z_hh

    # 3. Assemble cross-coupling blocks (hexa-tetra and tetra-hexa)
    Z_ts_buf = zeros(CT, 3, 3)

    for ti = 1:nhex
        hex_t = hexas[ti]
        κₜ = hex_t.κ

        for sj = 1:ntet
            tet_s = tetras[sj]
            κₛ = tet_s.κ

            dist_ts = norm(hex_t.center - tet_s.center)
            rad_t = cbrt(hex_t.volume)
            rad_s = cbrt(tet_s.volume)
            threshold = 3.0 * (rad_t + rad_s)

            if dist_ts > threshold
                _pwc_dyad_kernel!(
                    Z_ts_buf,
                    hex_t,
                    tet_s,
                    rq_hex_far_pts[ti],
                    rq_tet_far_pts[sj],
                    gq_hex_far,
                    gq_tet_far,
                    Nq_hex_far,
                    Nq_tet_far,
                    k,
                    k²,
                    jk,
                    Jη₀divK,
                    div4π,
                )
            else
                _pwc_dyad_kernel!(
                    Z_ts_buf,
                    hex_t,
                    tet_s,
                    rq_hex_near[ti],
                    rq_tet_near[sj],
                    gq_hex,
                    gq_tet,
                    Nq_hex,
                    Nq_tet,
                    k,
                    k²,
                    jk,
                    Jη₀divK,
                    div4π,
                )
            end

            for ni = 1:3, mi = 1:3
                m = hex_t.inBfsID[mi]
                n = tet_s.inBfsID[ni]
                Z[m, n] = Z_ts_buf[mi, ni] * κₛ
                Z[n, m] = Z_ts_buf[mi, ni] * κₜ
            end
        end
    end

    println("VEFIE-PWC Mixed Assembly Completed.")

    return Z
end

end
