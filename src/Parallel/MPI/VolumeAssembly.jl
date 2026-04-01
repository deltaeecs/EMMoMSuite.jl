# VolumeAssembly.jl
# MPI 骞惰浣撶Н鍒嗘柟绋嬬粍瑁咃細鍦?Parallel 鍛藉悕绌洪棿鎵╁睍
# `Assembly.assemble_impedance_matrix_parallel`銆?
#
# 褰撳墠瀹炵幇閲囩敤 MPI 鍒楀垎鍖猴紙姣忎釜 rank 鍙寔鏈夊叏灞€鐭╅樀鐨勮繛缁垪鑼冨洿锛夛紝
# 骞跺湪姣忎釜 rank 鍐呬娇鐢ㄥ绾跨▼骞惰閬嶅巻婧?娴嬭瘯鍑犱綍鍏冪礌銆?

import .Assembly: assemble_impedance_matrix_parallel

using ..IntegralEquations.VEFIEModule:
    VEFIE,
    get_tetrahedra_info,
    precompute_vefie_basis,
    vefie_element_interaction_kernel,
    vefie_mass_matrix_cached,
    _pwc_dyad_kernel!

using ..IntegralEquations.SCFIEModule:
    SCFIE,
    scfie_sv_only_interaction,
    assemble_fss_boundary_correction_sparse

using ..BasisFunctions: get_triangles_info, SWGBasis, RWGBasis,
    PWCBasis, PWCHexBasis, get_hexahedra_info

using ..Geometry: GaussQuadratureInfo

using ..IntegralEquations.EFIEModule: efie_interaction!
using ..IntegralEquations.MFIEModule: mfie_interaction!
using ..IntegralEquations.CFIEModule: CFIE

using SparseArrays
using LinearAlgebra
using StaticArrays: SVector
using Logging
using ...CoreModule: Constants

# =============================================================================
# Common MPI column-partition setup helpers
# =============================================================================

"""
    _mpi_surf_vol_setup(CT, n_surf, n_vol; comm, label) 鈫?NamedTuple

MPI 鍒楀垎鍖哄垵濮嬪寲杈呭姪锛堣〃闈?浣撶Н娣峰悎鐭╅樀锛夈€?

杩斿洖锛?
- `comm, rank, n_procs` 鈥?MPI 鐜
- `Z`                   鈥?鍒嗗竷寮忕煩闃碉紝澶у皬 (n_surf+n_vol) 脳 (n_surf+n_vol)
- `local_surf_cols`     鈥?鏈?rank 璐熻矗鐨勫垪涓睘浜庤〃闈?DOF锛?:n_surf锛夐儴鍒?
- `local_vol_cols`      鈥?鏈?rank 璐熻矗鐨勫垪涓睘浜庝綋绉?DOF锛坣_surf+1:n_total锛夐儴鍒?
- `row_locks`           鈥?姣忚涓€鎶?SpinLock锛岄槻姝㈣绾у埆鍐?鍐欑珵浜?
- `n_threads`           鈥?Threads.nthreads()

# 浣跨敤鑰呰礋璐ｏ細鎵撳嵃鐘舵€併€佸厓绱犲姞杞姐€佸叿浣撶Н鍒嗗惊鐜€?
"""
function _mpi_surf_vol_setup(CT::Type, n_surf::Int, n_vol::Int; comm = MPI.COMM_WORLD, label::String = "")
    rank    = MPI.Comm_rank(comm)
    n_procs = MPI.Comm_size(comm)
    n_total = n_surf + n_vol

    Z = mpiarray(CT, n_total, n_total; comm = comm)
    local_cols = Z.indices[2]

    surf_lo = max(first(local_cols), 1)
    surf_hi = min(last(local_cols), n_surf)
    vol_lo  = max(first(local_cols), n_surf + 1)
    vol_hi  = min(last(local_cols), n_total)

    local_surf_cols = surf_lo <= surf_hi ? (surf_lo:surf_hi) : (1:0)
    local_vol_cols  = vol_lo  <= vol_hi  ? (vol_lo:vol_hi)   : (1:0)

    if !isempty(label)
        rank == 0 && @info(
            "$label (n_surf=$n_surf, n_vol=$n_vol, n_total=$n_total, procs=$n_procs, threads=$(Threads.nthreads()))",
        )
        MPI.Barrier(comm)
        @info "  Rank $rank: surf_cols=$(length(local_surf_cols)), vol_cols=$(length(local_vol_cols))"
        MPI.Barrier(comm)
        flush(stdout)
    end

    row_locks = [SpinLock() for _ = 1:n_total]
    n_threads = Threads.nthreads()

    return (; comm, rank, n_procs, Z, local_cols = Z.indices[2], local_surf_cols, local_vol_cols, row_locks, n_threads)
end

# =============================================================================
# VEFIE + SWGBasis (MPI)
# =============================================================================
"""
    assemble_impedance_matrix_parallel(vefie, basis::SWGBasis, permittivities)

MPI 骞惰缁勮 VEFIE 闃绘姉鐭╅樀锛堝垪鍒嗗尯锛夈€?

姣忎釜 rank 浠呰礋璐ｆ湰鍦板垪 `local_cols`锛屽苟鍦ㄦ湰鍦扮嚎绋嬩腑瀹屾垚鐩稿叧婧愬洓闈綋璐＄尞鍐欏叆銆?
"""
function assemble_impedance_matrix_parallel(
    vefie::VEFIE,
    basis::SWGBasis,
    permittivities::AbstractVector,
)
    comm    = MPI.COMM_WORLD
    rank    = MPI.Comm_rank(comm)
    n_procs = MPI.Comm_size(comm)

    FT = typeof(real(vefie.freq))
    CT = Complex{FT}
    N  = num_basis(basis)

    Z = mpiarray(CT, N, N; comm = comm)
    local_cols = Z.indices[2]

    rank == 0 && @info(
        "VEFIE-SWG MPI Assembly (column partition): N=$N, procs=$n_procs, threads=$(Threads.nthreads())",
    )
    MPI.Barrier(comm)
    @info "  Rank $rank owns columns $(first(local_cols)):$(last(local_cols))"
    MPI.Barrier(comm)
    flush(stdout)

    tetras      = get_tetrahedra_info(basis.mesh, basis, permittivities)
    ntet        = length(tetras)
    basis_cache = precompute_vefie_basis(vefie, tetras)

    src_tets = Int[]
    for js = 1:ntet
        tet = tetras[js]
        for j = 1:4
            n = tet.inBfsID[j]
            if n != 0 && n in local_cols
                push!(src_tets, js)
                break
            end
        end
    end

    # Lockfree threading: partition local_cols into T exclusive sub-ranges.
    # Each thread owns col_lo..col_hi and writes Z[*, col] only for those cols.
    # Since col ranges are disjoint, no (row, col) pair is written concurrently.
    n_threads    = Threads.nthreads()
    n_local_cols = length(local_cols)
    col_chunk    = max(1, cld(n_local_cols, n_threads))

    Threads.@threads for tid = 1:n_threads
        col_lo = first(local_cols) + (tid - 1) * col_chunk
        col_hi = min(first(local_cols) + tid * col_chunk - 1, last(local_cols))

        for js in src_tets
            tet_s   = tetras[js]
            cache_s = basis_cache[js]

            # Quick skip if tet_s has no DOF in this thread's col range
            has_col = false
            @inbounds for j = 1:4
                n = tet_s.inBfsID[j]
                if n != 0 && col_lo <= n <= col_hi
                    has_col = true; break
                end
            end
            has_col || continue

            # Self-interaction (test == source)
            Z_self = vefie_element_interaction_kernel(vefie, tet_s, tet_s, cache_s, cache_s)
            M      = vefie_mass_matrix_cached(vefie, tet_s, cache_s)
            @inbounds for i = 1:4
                m = tet_s.inBfsID[i]
                m == 0 && continue
                @inbounds for j = 1:4
                    n = tet_s.inBfsID[j]
                    (n == 0 || !(col_lo <= n <= col_hi)) && continue
                    Z[m, n] += Z_self[i, j] + M[i, j]   # no lock: col exclusive
                end
            end

            # Off-diagonal interactions
            for it = 1:ntet
                it == js && continue
                tet_t   = tetras[it]
                cache_t = basis_cache[it]
                Z_ts    = vefie_element_interaction_kernel(vefie, tet_t, tet_s, cache_t, cache_s)
                @inbounds for i = 1:4
                    m = tet_t.inBfsID[i]
                    m == 0 && continue
                    @inbounds for j = 1:4
                        n = tet_s.inBfsID[j]
                        (n == 0 || !(col_lo <= n <= col_hi)) && continue
                        Z[m, n] += Z_ts[i, j]   # no lock: col exclusive
                    end
                end
            end
        end
    end

    sync!(Z)

    rank == 0 && println("VEFIE-SWG MPI Assembly Completed.")
    rank == 0 && flush(stdout)

    return Z
end


# =============================================================================
# SCFIE + RWGBasis + SWGBasis (MPI)
# =============================================================================
"""
    assemble_impedance_matrix_parallel(scfie, surf_basis::RWGBasis, vol_basis::SWGBasis)

MPI 骞惰缁勮 SCFIE 鑰﹀悎鐭╅樀锛堝垪鍒嗗尯锛夈€?

鍏ㄥ眬绱㈠紩甯冨眬:
- `1:n_surf`锛氳〃闈?RWG 鑷敱搴?
- `n_surf+1:n_total`锛氫綋 SWG 鑷敱搴?

鍚勫瓙鍧楁寜鈥滄簮鍒楁槸鍚﹀綊鏈?rank鈥濆啓鍏ユ湰鍦板垪鍧楋紝鏈€缁?`sync!(Z)` 鍚屾銆?
"""
function assemble_impedance_matrix_parallel(
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::SWGBasis,
)
    FT = typeof(real(scfie.freq))
    CT = Complex{FT}

    n_surf  = num_basis(surf_basis)
    n_vol   = num_basis(vol_basis)
    n_total = n_surf + n_vol

    (; comm, rank, n_procs, Z, local_cols, local_surf_cols, local_vol_cols, row_locks, n_threads) =
        _mpi_surf_vol_setup(CT, n_surf, n_vol; label = "SCFIE-SWG MPI Assembly (column partition)")

    tris        = get_triangles_info(surf_basis.mesh, surf_basis)
    ntri        = length(tris)
    tetras      = get_tetrahedra_info(vol_basis.mesh, vol_basis, scfie.permittivities)
    ntet        = length(tetras)
    vefie_inner = VEFIE(scfie.freq, scfie.permittivities)
    basis_cache = precompute_vefie_basis(vefie_inner, tetras)
    cfie        = CFIE(scfie.freq, scfie.alpha)

    src_tris_surf = Int[]
    for it = 1:ntri
        tri = tris[it]
        for j = 1:3
            m = tri.inBfsID[j]
            if m != 0 && m in local_surf_cols
                push!(src_tris_surf, it)
                break
            end
        end
    end

    src_tets_vol = Int[]
    for js = 1:ntet
        tet = tetras[js]
        for j = 1:4
            n = tet.inBfsID[j]
            if n != 0 && (n_surf + n) in local_vol_cols
                push!(src_tets_vol, js)
                break
            end
        end
    end

    # row_locks 鐢?_mpi_surf_vol_setup 鎻愪緵锛堝凡涓嶅啀鐢ㄤ簬鐑惊鐜紝淇濈暀渚?Fss 淇锛?

    rank == 0 && println("  [SCFIE 1+3b] Z_SS + Z_VS (surf cols)...")
    rank == 0 && flush(stdout)

    # Phase 1 lockfree: partition local_surf_cols into T exclusive sub-ranges.
    # Each thread writes only to Z[*, col] for col in its own surf col range.
    # Different threads never write the same (row,col) 鈫?no race, no lock.
    if !isempty(local_surf_cols)
        n_surf_local  = length(local_surf_cols)
        surf_col_chunk = max(1, cld(n_surf_local, n_threads))

        Threads.@threads for tid = 1:n_threads
            col_lo = first(local_surf_cols) + (tid - 1) * surf_col_chunk
            col_hi = min(first(local_surf_cols) + tid * surf_col_chunk - 1, last(local_surf_cols))
            l_efie = zeros(CT, 3, 3)
            l_mfie = zeros(CT, 3, 3)
            l_loc  = zeros(CT, 3, 3)

            for t_src in src_tris_surf
                tri_s = tris[t_src]

                # Quick skip if tri_s has no surf DOF in this thread's col range
                has_col = false
                @inbounds for j = 1:3
                    c = tri_s.inBfsID[j]
                    if c != 0 && col_lo <= c <= col_hi
                        has_col = true; break
                    end
                end
                has_col || continue

                # Z_SS: CFIE surf-surf
                for t_tst = 1:ntri
                    tri_t = tris[t_tst]
                    fill!(l_efie, zero(CT)); fill!(l_mfie, zero(CT)); fill!(l_loc, zero(CT))
                    efie_interaction!(l_efie, cfie.efie, tri_t, tri_s)
                    mfie_interaction!(l_mfie, cfie.mfie, tri_t, tri_s)
                    alpha_val = scfie.alpha
                    @. l_loc = alpha_val * l_efie + (1 - alpha_val) * l_mfie
                    @inbounds for j = 1:3
                        col = tri_s.inBfsID[j]
                        (col == 0 || !(col_lo <= col <= col_hi)) && continue
                        bf_s   = surf_basis.functions[col]
                        sign_s = (bf_s.support[1] == t_src) ? bf_s.signs[1] : bf_s.signs[2]
                        @inbounds for i = 1:3
                            row = tri_t.inBfsID[i]
                            row == 0 && continue
                            bf_t   = surf_basis.functions[row]
                            sign_t = (bf_t.support[1] == t_tst) ? bf_t.signs[1] : bf_t.signs[2]
                            Z[row, col] += l_loc[i, j] * sign_t * sign_s   # no lock: col exclusive
                        end
                    end
                end

                # Z_VS: vol rows, surf cols
                for js = 1:ntet
                    tet       = tetras[js]
                    Z_sv      = scfie_sv_only_interaction(scfie, tri_s, tet)
                    kappa_inv = iszero(tet.κ) ? zero(CT) : CT(1) / tet.κ
                    @inbounds for i = 1:3
                        col = tri_s.inBfsID[i]
                        (col == 0 || !(col_lo <= col <= col_hi)) && continue
                        @inbounds for j = 1:4
                            n = tet.inBfsID[j]
                            n == 0 && continue
                            row = n_surf + n
                            Z[row, col] += Z_sv[i, j] * kappa_inv   # no lock: col exclusive
                        end
                    end
                end
            end
        end
    end

    rank == 0 && println("  [SCFIE 2+3a] Z_VV + Z_SV (vol cols)...")
    rank == 0 && flush(stdout)

    # Phase 2 lockfree: partition local_vol_cols into T exclusive sub-ranges.
    if !isempty(local_vol_cols)
        n_vol_local  = length(local_vol_cols)
        vol_col_chunk = max(1, cld(n_vol_local, n_threads))

        Threads.@threads for tid = 1:n_threads
            col_lo = first(local_vol_cols) + (tid - 1) * vol_col_chunk
            col_hi = min(first(local_vol_cols) + tid * vol_col_chunk - 1, last(local_vol_cols))

            for js in src_tets_vol
                tet_s   = tetras[js]
                cache_s = basis_cache[js]

                # Quick skip if tet_s has no vol col in this thread's range
                has_col = false
                @inbounds for j = 1:4
                    n = tet_s.inBfsID[j]
                    c = n_surf + n
                    if n != 0 && col_lo <= c <= col_hi
                        has_col = true; break
                    end
                end
                has_col || continue

                # Z_VV: vol-vol (SWG)
                for it = 1:ntet
                    tet_t   = tetras[it]
                    cache_t = basis_cache[it]
                    if it == js
                        Z_self = vefie_element_interaction_kernel(
                            vefie_inner, tet_s, tet_s, cache_s, cache_s,
                        )
                        M_mat = vefie_mass_matrix_cached(vefie_inner, tet_s, cache_s)
                        @inbounds for j = 1:4
                            n = tet_s.inBfsID[j]
                            n == 0 && continue
                            col = n_surf + n
                            col_lo <= col <= col_hi || continue
                            @inbounds for i = 1:4
                                m = tet_s.inBfsID[i]
                                m == 0 && continue
                                row = n_surf + m
                                Z[row, col] += Z_self[i, j] + M_mat[i, j]   # no lock
                            end
                        end
                    else
                        Z_ts = vefie_element_interaction_kernel(
                            vefie_inner, tet_t, tet_s, cache_t, cache_s,
                        )
                        @inbounds for j = 1:4
                            n = tet_s.inBfsID[j]
                            n == 0 && continue
                            col = n_surf + n
                            col_lo <= col <= col_hi || continue
                            @inbounds for i = 1:4
                                m = tet_t.inBfsID[i]
                                m == 0 && continue
                                row = n_surf + m
                                Z[row, n_surf + n] += Z_ts[i, j]   # no lock
                            end
                        end
                    end
                end

                # Z_SV: surf rows, vol cols
                for it_tri = 1:ntri
                    tri  = tris[it_tri]
                    Z_sv = scfie_sv_only_interaction(scfie, tri, tet_s)
                    @inbounds for i = 1:3
                        m = tri.inBfsID[i]
                        m == 0 && continue
                        @inbounds for j = 1:4
                            n = tet_s.inBfsID[j]
                            n == 0 && continue
                            col = n_surf + n
                            col_lo <= col <= col_hi || continue
                            Z[m, col] += Z_sv[i, j]   # no lock: col exclusive
                        end
                    end
                end
            end
        end
    end

    rank == 0 && println("  [SCFIE Fss] boundary correction...")
    rank == 0 && flush(stdout)

    Z_fss = assemble_fss_boundary_correction_sparse(scfie, surf_basis, vol_basis)
    rows_f, cols_f, vals_f = findnz(Z_fss)
    @inbounds for k = 1:length(vals_f)
        c = cols_f[k]
        c in local_cols || continue
        Z[rows_f[k], c] += vals_f[k]
    end

    sync!(Z)

    rank == 0 && println("SCFIE MPI Assembly Completed.")
    rank == 0 && flush(stdout)

    return Z
end

# =============================================================================
# VEFIE + PWCBasis (MPI)
# =============================================================================
"""
    assemble_impedance_matrix_parallel(vefie::VEFIE, basis::PWCBasis, permittivities)

MPI 骞惰缁勮 VEFIE 闃绘姉鐭╅樀锛圥WC 鍥涢潰浣擄紝鍒楀垎鍖猴級銆?

姣忎釜鍥涢潰浣撹础鐚?3 涓?DOF锛坸/y/z 鍒嗛噺锛夈€?
鏃犲绉版€у埄鐢細婧愬垪杩囨护鍐冲畾姣忎釜 rank 澶勭悊鐨勬簮鍥涢潰浣擄紝
娴嬭瘯鍥涢潰浣撻亶鍘嗗叏閮紙鍚嚜椤瑰拰浜掗」锛夈€?

# Legacy Parity
涓庝覆琛?`assemble_impedance_matrix(vefie, basis::PWCBasis)` 缁撴灉涓€鑷淬€?
"""
function assemble_impedance_matrix_parallel(
    vefie::VEFIE,
    basis::PWCBasis,
    permittivities::AbstractVector,
)
    comm    = MPI.COMM_WORLD
    rank    = MPI.Comm_rank(comm)
    n_procs = MPI.Comm_size(comm)

    FT = typeof(real(vefie.freq))
    CT = Complex{FT}
    N  = num_basis(basis)

    Z = mpiarray(CT, N, N; comm = comm)
    local_cols = Z.indices[2]

    rank == 0 && @info(
        "VEFIE-PWC MPI Assembly (column partition): N=$N, procs=$n_procs, threads=$(Threads.nthreads())",
    )
    MPI.Barrier(comm)
    @info "  Rank $rank owns columns $(first(local_cols)):$(last(local_cols))"
    MPI.Barrier(comm)
    flush(stdout)

    tetras = get_tetrahedra_info(basis.mesh, basis, permittivities)
    ntet   = length(tetras)

    # Constants
    k          = vefie.k
    k虏         = k^2
    jk         = im * k
    omega      = 2π * vefie.freq
    mu0        = Constants.mu0
    eps0       = Constants.eps0
    eta0       = sqrt(mu0 / eps0)
    Jη₀divK   = im * eta0 / k
    div4π      = FT(1) / (4 * FT(π))

    # Quadrature from VEFIE struct
    gq     = vefie.gq_info
    gq_far = vefie.gq_far
    Nq     = length(gq.weight)
    Nq_far = length(gq_far.weight)

    # Precompute GQ points for all tetrahedra
    rq_near    = Vector{Matrix{FT}}(undef, ntet)
    rq_far_pts = Vector{Matrix{FT}}(undef, ntet)
    Threads.@threads for i = 1:ntet
        rq_near[i]    = tetras[i].vertices * gq.coordinate
        rq_far_pts[i] = tetras[i].vertices * gq_far.coordinate
    end

    # Source tets whose DOF columns fall in local_cols
    src_tets = Int[]
    for js = 1:ntet
        tet = tetras[js]
        for ni = 1:3
            n = tet.inBfsID[ni]
            if n != 0 && n in local_cols
                push!(src_tets, js)
                break
            end
        end
    end

    n_threads = Threads.nthreads()

    # Lockfree: each thread owns an exclusive sub-range of local_cols.
    # col (source DOF) is exclusive per thread 鈫?no (row,col) write conflict.
    n_local_cols = length(local_cols)
    col_chunk    = max(1, cld(n_local_cols, n_threads))

    Threads.@threads for tid = 1:n_threads
        col_lo   = first(local_cols) + (tid - 1) * col_chunk
        col_hi   = min(first(local_cols) + tid * col_chunk - 1, last(local_cols))
        Z_ts_buf = zeros(CT, 3, 3)

        for js in src_tets
            tet_s = tetras[js]
            κₛ   = tet_s.κ

            # Quick skip if tet_s has no source DOF in this thread's col range
            has_col = false
            @inbounds for ni = 1:3
                n = tet_s.inBfsID[ni]
                if n != 0 && col_lo <= n <= col_hi
                    has_col = true; break
                end
            end
            has_col || continue

            rad_s = cbrt(tet_s.volume)

            # Self-term (ti == js)
            _pwc_dyad_kernel!(
                Z_ts_buf,
                tet_s, tet_s,
                rq_near[js], rq_near[js],
                gq, gq,
                Nq, Nq,
                k, k虏, jk, Jη₀divK, div4π,
            )
            selfImp = CT(1) / (im * omega) / tet_s.ε * tet_s.volume
            @inbounds for ni = 1:3
                n = tet_s.inBfsID[ni]
                (n == 0 || !(col_lo <= n <= col_hi)) && continue
                for mi = 1:3
                    m = tet_s.inBfsID[mi]
                    m == 0 && continue
                    val = Z_ts_buf[mi, ni] * κₛ
                    mi == ni && (val += selfImp)
                    Z[m, n] += val   # no lock: n exclusive to tid
                end
            end

            # Off-diagonal terms (ti != js)
            for ti = 1:ntet
                ti == js && continue
                tet_t = tetras[ti]

                dist_ts   = norm(tet_t.center - tet_s.center)
                rad_t     = cbrt(tet_t.volume)
                threshold = FT(3) * (rad_t + rad_s)

                if dist_ts > threshold
                    _pwc_dyad_kernel!(
                        Z_ts_buf,
                        tet_t, tet_s,
                        rq_far_pts[ti], rq_far_pts[js],
                        gq_far, gq_far,
                        Nq_far, Nq_far,
                        k, k虏, jk, Jη₀divK, div4π,
                    )
                else
                    _pwc_dyad_kernel!(
                        Z_ts_buf,
                        tet_t, tet_s,
                        rq_near[ti], rq_near[js],
                        gq, gq,
                        Nq, Nq,
                        k, k虏, jk, Jη₀divK, div4π,
                    )
                end

                @inbounds for ni = 1:3
                    n = tet_s.inBfsID[ni]
                    (n == 0 || !(col_lo <= n <= col_hi)) && continue
                    for mi = 1:3
                        m = tet_t.inBfsID[mi]
                        m == 0 && continue
                        Z[m, n] += Z_ts_buf[mi, ni] * κₛ  # no lock: n exclusive to tid
                    end
                end
            end
        end
    end

    sync!(Z)

    rank == 0 && println("VEFIE-PWC MPI Assembly Completed.")
    rank == 0 && flush(stdout)

    return Z
end

"""
    assemble_impedance_matrix_parallel(vefie::VEFIE, basis::PWCBasis)

浣跨敤 `vefie.permittivities` 鐨勪究鎹烽噸杞姐€?
"""
function assemble_impedance_matrix_parallel(vefie::VEFIE, basis::PWCBasis)
    return assemble_impedance_matrix_parallel(vefie, basis, vefie.permittivities)
end

# =============================================================================
# VEFIE + PWCHexBasis (MPI)
# =============================================================================
"""
    assemble_impedance_matrix_parallel(vefie::VEFIE, basis::PWCHexBasis)

MPI 骞惰缁勮 VEFIE 闃绘姉鐭╅樀锛圥WC 鍏潰浣擄紝鍒楀垎鍖猴級銆?

姣忎釜鍏潰浣撹础鐚?3 涓?DOF锛坸/y/z锛夈€備娇鐢?Hexahedron 8 鐐癸紙杩戝満锛?1 鐐癸紙杩滃満锛夐珮鏂Н鍒嗐€?

# Legacy Parity
涓庝覆琛?`assemble_impedance_matrix(vefie, basis::PWCHexBasis)` 缁撴灉涓€鑷淬€?
"""
function assemble_impedance_matrix_parallel(vefie::VEFIE, basis::PWCHexBasis)
    comm    = MPI.COMM_WORLD
    rank    = MPI.Comm_rank(comm)
    n_procs = MPI.Comm_size(comm)

    FT = typeof(real(vefie.freq))
    CT = Complex{FT}
    N  = num_basis(basis)

    Z = mpiarray(CT, N, N; comm = comm)
    local_cols = Z.indices[2]

    rank == 0 && @info(
        "VEFIE-PWC-Hex MPI Assembly (column partition): N=$N, procs=$n_procs, threads=$(Threads.nthreads())",
    )
    MPI.Barrier(comm)
    @info "  Rank $rank owns columns $(first(local_cols)):$(last(local_cols))"
    MPI.Barrier(comm)
    flush(stdout)

    hexas = get_hexahedra_info(basis.mesh, basis, vefie.permittivities)
    nhex  = length(hexas)

    # Constants
    k         = vefie.k
    k虏        = k^2
    jk        = im * k
    omega     = 2π * vefie.freq
    mu0       = Constants.mu0
    eps0      = Constants.eps0
    eta0      = sqrt(mu0 / eps0)
    Jη₀divK  = im * eta0 / k
    div4π     = FT(1) / (4 * FT(π))

    # Local GQ for hexahedra
    gq_hex     = GaussQuadratureInfo(:Hexahedron, 8, FT)
    gq_hex_far = GaussQuadratureInfo(:Hexahedron, 1, FT)
    Nq_hex     = length(gq_hex.weight)
    Nq_hex_far = length(gq_hex_far.weight)

    # Precompute GQ points
    rq_near    = Vector{Matrix{FT}}(undef, nhex)
    rq_far_pts = Vector{Matrix{FT}}(undef, nhex)
    Threads.@threads for i = 1:nhex
        rq_near[i]    = hexas[i].vertices * gq_hex.coordinate
        rq_far_pts[i] = hexas[i].vertices * gq_hex_far.coordinate
    end

    # Source hexahedra whose DOF columns fall in local_cols
    src_hexas = Int[]
    for js = 1:nhex
        hex = hexas[js]
        for ni = 1:3
            n = hex.inBfsID[ni]
            if n != 0 && n in local_cols
                push!(src_hexas, js)
                break
            end
        end
    end

    n_threads = Threads.nthreads()

    # Lockfree: each thread owns an exclusive sub-range of local_cols.
    # col (source DOF) is exclusive per thread 鈫?no (row,col) write conflict.
    n_local_cols = length(local_cols)
    col_chunk    = max(1, cld(n_local_cols, n_threads))

    Threads.@threads for tid = 1:n_threads
        col_lo   = first(local_cols) + (tid - 1) * col_chunk
        col_hi   = min(first(local_cols) + tid * col_chunk - 1, last(local_cols))
        Z_ts_buf = zeros(CT, 3, 3)

        for js in src_hexas
            hex_s = hexas[js]
            κₛ   = hex_s.κ

            # Quick skip if hex_s has no source DOF in this thread's col range
            has_col = false
            @inbounds for ni = 1:3
                n = hex_s.inBfsID[ni]
                if n != 0 && col_lo <= n <= col_hi
                    has_col = true; break
                end
            end
            has_col || continue

            rad_s = cbrt(hex_s.volume)

            # Self-term
            _pwc_dyad_kernel!(
                Z_ts_buf,
                hex_s, hex_s,
                rq_near[js], rq_near[js],
                gq_hex, gq_hex,
                Nq_hex, Nq_hex,
                k, k虏, jk, Jη₀divK, div4π,
            )
            selfImp = CT(1) / (im * omega) / hex_s.ε * hex_s.volume
            @inbounds for ni = 1:3
                n = hex_s.inBfsID[ni]
                (n == 0 || !(col_lo <= n <= col_hi)) && continue
                for mi = 1:3
                    m = hex_s.inBfsID[mi]
                    m == 0 && continue
                    val = Z_ts_buf[mi, ni] * κₛ
                    mi == ni && (val += selfImp)
                    Z[m, n] += val   # no lock: n exclusive to tid
                end
            end

            # Off-diagonal terms
            for ti = 1:nhex
                ti == js && continue
                hex_t = hexas[ti]

                dist_ts   = norm(hex_t.center - hex_s.center)
                rad_t     = cbrt(hex_t.volume)
                threshold = FT(3) * (rad_t + rad_s)

                if dist_ts > threshold
                    _pwc_dyad_kernel!(
                        Z_ts_buf,
                        hex_t, hex_s,
                        rq_far_pts[ti], rq_far_pts[js],
                        gq_hex_far, gq_hex_far,
                        Nq_hex_far, Nq_hex_far,
                        k, k虏, jk, Jη₀divK, div4π,
                    )
                else
                    _pwc_dyad_kernel!(
                        Z_ts_buf,
                        hex_t, hex_s,
                        rq_near[ti], rq_near[js],
                        gq_hex, gq_hex,
                        Nq_hex, Nq_hex,
                        k, k虏, jk, Jη₀divK, div4π,
                    )
                end

                @inbounds for ni = 1:3
                    n = hex_s.inBfsID[ni]
                    (n == 0 || !(col_lo <= n <= col_hi)) && continue
                    for mi = 1:3
                        m = hex_t.inBfsID[mi]
                        m == 0 && continue
                        Z[m, n] += Z_ts_buf[mi, ni] * κₛ  # no lock: n exclusive to tid
                    end
                end
            end
        end
    end

    sync!(Z)

    rank == 0 && println("VEFIE-PWC-Hex MPI Assembly Completed.")
    rank == 0 && flush(stdout)

    return Z
end

# =============================================================================
# SCFIE + RWGBasis + PWCBasis (MPI)
# =============================================================================
"""
    assemble_impedance_matrix_parallel(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::PWCBasis)

MPI 骞惰缁勮 SCFIE 鑰﹀悎鐭╅樀锛堝垪鍒嗗尯锛夛紝浣撶Н鍒嗘柟绋嬪熀鍑芥暟涓?PWCBasis锛堝洓闈綋锛夈€?

鍏ㄥ眬绱㈠紩甯冨眬:
- `1:n_surf`  : 琛ㄩ潰 RWG 鑷敱搴?
- `n_surf+1:n_total` : 浣?PWC 鑷敱搴︼紙姣忓洓闈綋 3 涓?x/y/z 鍒嗛噺锛?

缁勮鍒嗕袱闃舵锛?
1. **surf cols**锛堟簮鍦ㄨ〃闈級: 濉厖 Z_SS + Z_VS
2. **vol cols** 锛堟簮鍦ㄤ綋绉級: 濉厖 Z_VV + Z_SV

鏃?Fss 杈圭晫淇锛圥WC 鏃犲崐鍩哄嚱鏁帮級銆?

# Legacy Parity
涓庝覆琛?`assemble_impedance_matrix(scfie, surf_basis::RWGBasis, vol_basis::PWCBasis)` 涓€鑷淬€?
"""
function assemble_impedance_matrix_parallel(
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::PWCBasis,
)
    FT = typeof(real(scfie.freq))
    CT = Complex{FT}

    n_surf  = num_basis(surf_basis)
    n_vol   = num_basis(vol_basis)
    n_total = n_surf + n_vol

    (; comm, rank, n_procs, Z, local_surf_cols, local_vol_cols, row_locks, n_threads) =
        _mpi_surf_vol_setup(CT, n_surf, n_vol; label = "SCFIE-PWC MPI Assembly (column partition)")

    tris   = get_triangles_info(surf_basis.mesh, surf_basis)
    ntri   = length(tris)
    tetras = get_tetrahedra_info(vol_basis.mesh, vol_basis, scfie.permittivities)
    ntet   = length(tetras)

    cfie = CFIE(scfie.freq, scfie.alpha)

    # Physical constants (for coupling dyad and Z_VV kernel)
    k         = scfie.k
    k虏        = k^2
    jk        = im * k
    omega     = FT(2π) * scfie.freq
    eta0      = scfie.eta
    Jη₀divK  = im * eta0 / k
    div4π     = FT(1) / (4 * FT(π))

    # Quadrature
    gq_s       = scfie.gq_surf                            # Triangle 7-pt
    gq_v       = scfie.gq_vol                             # Tet 5-pt (near)
    gq_far_tet = GaussQuadratureInfo(:Tetrahedron, 1, FT) # Tet 1-pt (far Z_VV)
    Nq_s       = length(gq_s.weight)
    Nq_v       = length(gq_v.weight)
    Nq_far     = length(gq_far_tet.weight)

    # Precompute GQ points for all elements
    rq_tri     = [tris[i].vertices * gq_s.coordinate   for i = 1:ntri]
    rq_tet_near = [tetras[i].vertices * gq_v.coordinate     for i = 1:ntet]
    rq_tet_far  = [tetras[i].vertices * gq_far_tet.coordinate for i = 1:ntet]

    # Source element sets
    src_tris_surf = Int[]
    for it = 1:ntri
        tri = tris[it]
        for j = 1:3
            m = tri.inBfsID[j]
            if m != 0 && m in local_surf_cols
                push!(src_tris_surf, it)
                break
            end
        end
    end

    src_tets_vol = Int[]
    for js = 1:ntet
        tet = tetras[js]
        for j = 1:3  # PWC: 3 DOFs per tet
            n = tet.inBfsID[j]
            if n != 0 && (n_surf + n) in local_vol_cols
                push!(src_tets_vol, js)
                break
            end
        end
    end


    # 鈹€鈹€鈹€ Phase 1: Z_SS + Z_VS (surf cols) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    rank == 0 && println("  [SCFIE-PWC 1] Z_SS + Z_VS (surf cols)...")
    rank == 0 && flush(stdout)

    # Lockfree: each thread owns exclusive surf_col sub-range 鈫?no (row,col) race.
    if !isempty(local_surf_cols)
        n_surf_local   = length(local_surf_cols)
        surf_col_chunk = max(1, cld(n_surf_local, n_threads))

    Threads.@threads for tid = 1:n_threads
        col_lo = first(local_surf_cols) + (tid - 1) * surf_col_chunk
        col_hi = min(first(local_surf_cols) + tid * surf_col_chunk - 1, last(local_surf_cols))
        l_efie = zeros(CT, 3, 3)
        l_mfie = zeros(CT, 3, 3)
        l_loc  = zeros(CT, 3, 3)
        dyadG  = zeros(CT, 3, 3)

        for t_src in src_tris_surf
            tri_s = tris[t_src]

            # Quick skip if tri_s has no surf DOF in this thread's col range
            has_col = false
            @inbounds for j = 1:3
                c = tri_s.inBfsID[j]
                if c != 0 && col_lo <= c <= col_hi; has_col = true; break; end
            end
            has_col || continue

            # Z_SS: CFIE surface-surface interactions
            for t_tst = 1:ntri
                tri_t = tris[t_tst]
                fill!(l_efie, zero(CT)); fill!(l_mfie, zero(CT)); fill!(l_loc, zero(CT))
                efie_interaction!(l_efie, cfie.efie, tri_t, tri_s)
                mfie_interaction!(l_mfie, cfie.mfie, tri_t, tri_s)
                alpha_val = scfie.alpha
                @. l_loc = alpha_val * l_efie + (1 - alpha_val) * l_mfie
                @inbounds for j = 1:3
                    col = tri_s.inBfsID[j]
                    (col == 0 || !(col_lo <= col <= col_hi)) && continue
                    bf_s   = surf_basis.functions[col]
                    sign_s = (bf_s.support[1] == t_src) ? bf_s.signs[1] : bf_s.signs[2]
                    @inbounds for i = 1:3
                        row = tri_t.inBfsID[i]
                        row == 0 && continue
                        bf_t   = surf_basis.functions[row]
                        sign_t = (bf_t.support[1] == t_tst) ? bf_t.signs[1] : bf_t.signs[2]
                        Z[row, col] += l_loc[i, j] * sign_t * sign_s   # no lock: col exclusive
                    end
                end
            end

            # Z_VS: source = tri_s (col = m 鈭?local_surf_cols), test = all tets
            r_q_tri = rq_tri[t_src]
            for js = 1:ntet
                tet   = tetras[js]
                Vs    = tet.volume
                r_q_v = rq_tet_near[js]

                for gi = 1:Nq_s
                    rgi = @view r_q_tri[:, gi]
                    fill!(dyadG, zero(CT))

                    for gj = 1:Nq_v
                        rgj = @view r_q_v[:, gj]
                        Rx = rgi[1] - rgj[1]
                        Ry = rgi[2] - rgj[2]
                        Rz = rgi[3] - rgj[3]
                        R  = sqrt(Rx^2 + Ry^2 + Rz^2)
                        R < FT(1e-10) && continue

                        divR  = FT(1) / R
                        jkpR  = (jk + divR) * divR
                        R虃x = Rx * divR; R虃y = Ry * divR; R虃z = Rz * divR
                        GR   = exp(-jk * R) * div4π * divR * gq_v.weight[gj]

                        RR11 = R虃x*R虃x; RR12 = R虃x*R虃y; RR13 = R虃x*R虃z
                        RR22 = R虃y*R虃y; RR23 = R虃y*R虃z; RR33 = R虃z*R虃z

                        dyadG[1,1] += GR * ((1 - RR11) * k虏 - (1 - 3RR11) * jkpR)
                        dyadG[2,2] += GR * ((1 - RR22) * k虏 - (1 - 3RR22) * jkpR)
                        dyadG[3,3] += GR * ((1 - RR33) * k虏 - (1 - 3RR33) * jkpR)
                        od12 = GR * (-RR12 * k虏 + 3RR12 * jkpR)
                        od13 = GR * (-RR13 * k虏 + 3RR13 * jkpR)
                        od23 = GR * (-RR23 * k虏 + 3RR23 * jkpR)
                        dyadG[1,2] += od12; dyadG[2,1] += od12
                        dyadG[1,3] += od13; dyadG[3,1] += od13
                        dyadG[2,3] += od23; dyadG[3,2] += od23
                    end

                    for mi = 1:3
                        m = tri_s.inBfsID[mi]
                        m == 0 && continue
                        col_lo <= m <= col_hi || continue   # this thread's surf cols only
                        lm    = tri_s.edgel[mi]
                        freeV = tri_s.vertices[:, mi]
                        蟻mi   = SVector(rgi[1] - freeV[1], rgi[2] - freeV[2], rgi[3] - freeV[3])
                        temp  = gq_s.weight[gi] * lm / 2

                        for ni = 1:3
                            n     = tet.inBfsID[ni]
                            dot_vs = 蟻mi[1]*dyadG[ni,1] + 蟻mi[2]*dyadG[ni,2] + 蟻mi[3]*dyadG[ni,3]
                            z_vs   = temp * dot_vs * Jη₀divK * Vs
                            row_vs = n_surf + n
                            Z[row_vs, m] += z_vs   # no lock: col=m exclusive to tid
                        end
                    end
                end
            end
        end
    end   # threads
    end   # if !isempty(local_surf_cols)

    # 鈹€鈹€鈹€ Phase 2: Z_VV + Z_SV (vol cols) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    rank == 0 && println("  [SCFIE-PWC 2] Z_VV + Z_SV (vol cols)...")
    rank == 0 && flush(stdout)

    # Lockfree: each thread owns exclusive vol_col sub-range.
    if !isempty(local_vol_cols)
        n_vol_local   = length(local_vol_cols)
        vol_col_chunk = max(1, cld(n_vol_local, n_threads))

    Threads.@threads for tid = 1:n_threads
        col_lo   = first(local_vol_cols) + (tid - 1) * vol_col_chunk
        col_hi   = min(first(local_vol_cols) + tid * vol_col_chunk - 1, last(local_vol_cols))
        Z_ts_buf = zeros(CT, 3, 3)
        dyadG    = zeros(CT, 3, 3)

        for js in src_tets_vol
            tet_s = tetras[js]
            κₛ   = tet_s.κ
            rad_s = cbrt(tet_s.volume)

            # Quick skip if tet_s has no vol col in this thread's range
            has_col = false
            @inbounds for ni = 1:3
                n = tet_s.inBfsID[ni]
                c = n_surf + n
                if n != 0 && col_lo <= c <= col_hi; has_col = true; break; end
            end
            has_col || continue

            # Z_VV: _pwc_dyad_kernel! for (tet_s, all tets)
            for ti = 1:ntet
                tet_t = tetras[ti]

                dist_ts   = norm(tet_t.center - tet_s.center)
                rad_t     = cbrt(tet_t.volume)
                threshold = FT(3) * (rad_t + rad_s)
                is_far    = dist_ts > threshold

                _pwc_dyad_kernel!(
                    Z_ts_buf,
                    tet_t, tet_s,
                    is_far ? rq_tet_far[ti] : rq_tet_near[ti],
                    is_far ? rq_tet_far[js] : rq_tet_near[js],
                    is_far ? gq_far_tet : gq_v,
                    is_far ? gq_far_tet : gq_v,
                    is_far ? Nq_far : Nq_v,
                    is_far ? Nq_far : Nq_v,
                    k, k虏, jk, Jη₀divK, div4π,
                )

                @inbounds for ni = 1:3
                    n = tet_s.inBfsID[ni]
                    col = n_surf + n
                    col_lo <= col <= col_hi || continue
                    for mi = 1:3
                        m   = tet_t.inBfsID[mi]
                        row = n_surf + m
                        if ti == js
                            val = Z_ts_buf[mi, ni] * κₛ
                            if mi == ni
                                selfImp = CT(1) / (im * omega) / tet_s.ε * tet_s.volume
                                val += selfImp
                            end
                            Z[row, col] += val   # no lock: col exclusive
                        else
                            Z[row, col] += Z_ts_buf[mi, ni] * κₛ  # no lock
                        end
                    end
                end
            end

            # Z_SV: source = tet_s (col 鈭?thread's vol range), test = all tris
            r_q_v = rq_tet_near[js]
            for it = 1:ntri
                tri   = tris[it]
                r_q_t = rq_tri[it]

                for gi = 1:Nq_s
                    rgi = @view r_q_t[:, gi]
                    fill!(dyadG, zero(CT))

                    for gj = 1:Nq_v
                        rgj = @view r_q_v[:, gj]
                        Rx = rgi[1] - rgj[1]
                        Ry = rgi[2] - rgj[2]
                        Rz = rgi[3] - rgj[3]
                        R  = sqrt(Rx^2 + Ry^2 + Rz^2)
                        R < FT(1e-10) && continue

                        divR  = FT(1) / R
                        jkpR  = (jk + divR) * divR
                        R虃x = Rx * divR; R虃y = Ry * divR; R虃z = Rz * divR
                        GR   = exp(-jk * R) * div4π * divR * gq_v.weight[gj]

                        RR11 = R虃x*R虃x; RR12 = R虃x*R虃y; RR13 = R虃x*R虃z
                        RR22 = R虃y*R虃y; RR23 = R虃y*R虃z; RR33 = R虃z*R虃z

                        dyadG[1,1] += GR * ((1 - RR11) * k虏 - (1 - 3RR11) * jkpR)
                        dyadG[2,2] += GR * ((1 - RR22) * k虏 - (1 - 3RR22) * jkpR)
                        dyadG[3,3] += GR * ((1 - RR33) * k虏 - (1 - 3RR33) * jkpR)
                        od12 = GR * (-RR12 * k虏 + 3RR12 * jkpR)
                        od13 = GR * (-RR13 * k虏 + 3RR13 * jkpR)
                        od23 = GR * (-RR23 * k虏 + 3RR23 * jkpR)
                        dyadG[1,2] += od12; dyadG[2,1] += od12
                        dyadG[1,3] += od13; dyadG[3,1] += od13
                        dyadG[2,3] += od23; dyadG[3,2] += od23
                    end

                    for mi = 1:3
                        m = tri.inBfsID[mi]
                        m == 0 && continue
                        lm    = tri.edgel[mi]
                        freeV = tri.vertices[:, mi]
                        蟻mi   = SVector(rgi[1] - freeV[1], rgi[2] - freeV[2], rgi[3] - freeV[3])
                        temp  = gq_s.weight[gi] * lm / 2

                        for ni = 1:3
                            n   = tet_s.inBfsID[ni]
                            col = n_surf + n
                            col_lo <= col <= col_hi || continue
                            dot_sv = 蟻mi[1]*dyadG[1,ni] + 蟻mi[2]*dyadG[2,ni] + 蟻mi[3]*dyadG[3,ni]
                            z_sv   = temp * dot_sv * Jη₀divK * tet_s.volume * κₛ
                            Z[m, col] += z_sv   # no lock: col exclusive
                        end
                    end
                end
            end
        end
    end   # threads
    end   # if !isempty(local_vol_cols)

    sync!(Z)

    rank == 0 && println("SCFIE-PWC MPI Assembly Completed.")
    rank == 0 && flush(stdout)

    return Z
end

# =============================================================================
# SCFIE + RWGBasis + PWCHexBasis (MPI)
# =============================================================================
"""
    assemble_impedance_matrix_parallel(scfie::SCFIE, surf_basis::RWGBasis, vol_basis::PWCHexBasis)

MPI 骞惰缁勮 SCFIE 鑰﹀悎鐭╅樀锛堝垪鍒嗗尯锛夛紝浣撶Н鍒嗘柟绋嬪熀鍑芥暟涓?PWCHexBasis锛堝叚闈綋锛夈€?

涓?PWCBasis 鐗堟湰缁撴瀯鐩稿悓锛屽尯鍒湪浜庯細鍏潰浣?GQ锛?鐐硅繎鍦?1鐐硅繙鍦猴級鍦ㄦ湰鍦板垱寤恒€?

# Legacy Parity
涓庝覆琛?`assemble_impedance_matrix(scfie, surf_basis::RWGBasis, vol_basis::PWCHexBasis)` 涓€鑷淬€?
"""
function assemble_impedance_matrix_parallel(
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::PWCHexBasis,
)
    FT = typeof(real(scfie.freq))
    CT = Complex{FT}

    n_surf  = num_basis(surf_basis)
    n_vol   = num_basis(vol_basis)
    n_total = n_surf + n_vol

    (; comm, rank, n_procs, Z, local_surf_cols, local_vol_cols, row_locks, n_threads) =
        _mpi_surf_vol_setup(CT, n_surf, n_vol; label = "SCFIE-PWCHex MPI Assembly (column partition)")

    tris  = get_triangles_info(surf_basis.mesh, surf_basis)
    ntri  = length(tris)
    hexas = get_hexahedra_info(vol_basis.mesh, vol_basis, scfie.permittivities)
    nhex  = length(hexas)

    cfie = CFIE(scfie.freq, scfie.alpha)

    # Physical constants
    k         = scfie.k
    k虏        = k^2
    jk        = im * k
    omega     = FT(2π) * scfie.freq
    eta0      = scfie.eta
    Jη₀divK  = im * eta0 / k
    div4π     = FT(1) / (4 * FT(π))

    # Quadrature
    gq_s        = scfie.gq_surf                              # Triangle 7-pt
    gq_hex_near = GaussQuadratureInfo(:Hexahedron, 8, FT)    # Hex 8-pt (near)
    gq_hex_far  = GaussQuadratureInfo(:Hexahedron, 1, FT)    # Hex 1-pt (far)
    Nq_s        = length(gq_s.weight)
    Nq_hex      = length(gq_hex_near.weight)
    Nq_hex_far  = length(gq_hex_far.weight)

    # Precompute GQ points
    rq_tri      = [tris[i].vertices * gq_s.coordinate          for i = 1:ntri]
    rq_hex_near = [hexas[i].vertices * gq_hex_near.coordinate   for i = 1:nhex]
    rq_hex_far  = [hexas[i].vertices * gq_hex_far.coordinate    for i = 1:nhex]

    # Source element sets
    src_tris_surf = Int[]
    for it = 1:ntri
        tri = tris[it]
        for j = 1:3
            m = tri.inBfsID[j]
            if m != 0 && m in local_surf_cols
                push!(src_tris_surf, it)
                break
            end
        end
    end

    src_hexas_vol = Int[]
    for js = 1:nhex
        hex = hexas[js]
        for j = 1:3
            n = hex.inBfsID[j]
            if n != 0 && (n_surf + n) in local_vol_cols
                push!(src_hexas_vol, js)
                break
            end
        end
    end


    # 鈹€鈹€鈹€ Phase 1: Z_SS + Z_VS (surf cols) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    rank == 0 && println("  [SCFIE-PWCHex 1] Z_SS + Z_VS (surf cols)...")
    rank == 0 && flush(stdout)

    # Lockfree: each thread owns exclusive surf_col sub-range.
    if !isempty(local_surf_cols)
        n_surf_local   = length(local_surf_cols)
        surf_col_chunk = max(1, cld(n_surf_local, n_threads))

    Threads.@threads for tid = 1:n_threads
        col_lo = first(local_surf_cols) + (tid - 1) * surf_col_chunk
        col_hi = min(first(local_surf_cols) + tid * surf_col_chunk - 1, last(local_surf_cols))
        l_efie = zeros(CT, 3, 3)
        l_mfie = zeros(CT, 3, 3)
        l_loc  = zeros(CT, 3, 3)
        dyadG  = zeros(CT, 3, 3)

        for t_src in src_tris_surf
            tri_s = tris[t_src]

            # Quick skip: any surf col in this thread's range?
            has_col = false
            @inbounds for j = 1:3
                c = tri_s.inBfsID[j]
                if c != 0 && col_lo <= c <= col_hi; has_col = true; break; end
            end
            has_col || continue

            # Z_SS
            for t_tst = 1:ntri
                tri_t = tris[t_tst]
                fill!(l_efie, zero(CT)); fill!(l_mfie, zero(CT)); fill!(l_loc, zero(CT))
                efie_interaction!(l_efie, cfie.efie, tri_t, tri_s)
                mfie_interaction!(l_mfie, cfie.mfie, tri_t, tri_s)
                alpha_val = scfie.alpha
                @. l_loc = alpha_val * l_efie + (1 - alpha_val) * l_mfie
                @inbounds for j = 1:3
                    col = tri_s.inBfsID[j]
                    (col == 0 || !(col_lo <= col <= col_hi)) && continue
                    bf_s   = surf_basis.functions[col]
                    sign_s = (bf_s.support[1] == t_src) ? bf_s.signs[1] : bf_s.signs[2]
                    @inbounds for i = 1:3
                        row = tri_t.inBfsID[i]
                        row == 0 && continue
                        bf_t   = surf_basis.functions[row]
                        sign_t = (bf_t.support[1] == t_tst) ? bf_t.signs[1] : bf_t.signs[2]
                        Z[row, col] += l_loc[i, j] * sign_t * sign_s   # no lock: col exclusive
                    end
                end
            end

            # Z_VS: source = tri_s, test = all hexas
            r_q_tri = rq_tri[t_src]
            for js = 1:nhex
                hex   = hexas[js]
                Vs    = hex.volume
                r_q_v = rq_hex_near[js]

                for gi = 1:Nq_s
                    rgi = @view r_q_tri[:, gi]
                    fill!(dyadG, zero(CT))

                    for gj = 1:Nq_hex
                        rgj = @view r_q_v[:, gj]
                        Rx = rgi[1] - rgj[1]
                        Ry = rgi[2] - rgj[2]
                        Rz = rgi[3] - rgj[3]
                        R  = sqrt(Rx^2 + Ry^2 + Rz^2)
                        R < FT(1e-10) && continue

                        divR  = FT(1) / R
                        jkpR  = (jk + divR) * divR
                        R虃x = Rx * divR; R虃y = Ry * divR; R虃z = Rz * divR
                        GR   = exp(-jk * R) * div4π * divR * gq_hex_near.weight[gj]

                        RR11 = R虃x*R虃x; RR12 = R虃x*R虃y; RR13 = R虃x*R虃z
                        RR22 = R虃y*R虃y; RR23 = R虃y*R虃z; RR33 = R虃z*R虃z

                        dyadG[1,1] += GR * ((1 - RR11) * k虏 - (1 - 3RR11) * jkpR)
                        dyadG[2,2] += GR * ((1 - RR22) * k虏 - (1 - 3RR22) * jkpR)
                        dyadG[3,3] += GR * ((1 - RR33) * k虏 - (1 - 3RR33) * jkpR)
                        od12 = GR * (-RR12 * k虏 + 3RR12 * jkpR)
                        od13 = GR * (-RR13 * k虏 + 3RR13 * jkpR)
                        od23 = GR * (-RR23 * k虏 + 3RR23 * jkpR)
                        dyadG[1,2] += od12; dyadG[2,1] += od12
                        dyadG[1,3] += od13; dyadG[3,1] += od13
                        dyadG[2,3] += od23; dyadG[3,2] += od23
                    end

                    for mi = 1:3
                        m = tri_s.inBfsID[mi]
                        m == 0 && continue
                        col_lo <= m <= col_hi || continue   # this thread's surf cols only
                        lm    = tri_s.edgel[mi]
                        freeV = tri_s.vertices[:, mi]
                        蟻mi   = SVector(rgi[1]-freeV[1], rgi[2]-freeV[2], rgi[3]-freeV[3])
                        temp  = gq_s.weight[gi] * lm / 2

                        for ni = 1:3
                            n     = hex.inBfsID[ni]
                            dot_vs = 蟻mi[1]*dyadG[ni,1] + 蟻mi[2]*dyadG[ni,2] + 蟻mi[3]*dyadG[ni,3]
                            z_vs   = temp * dot_vs * Jη₀divK * Vs
                            row_vs = n_surf + n
                            Z[row_vs, m] += z_vs   # no lock: col=m exclusive to tid
                        end
                    end
                end
            end
        end
    end   # threads
    end   # if !isempty(local_surf_cols)

    # 鈹€鈹€鈹€ Phase 2: Z_VV + Z_SV (vol cols) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    rank == 0 && println("  [SCFIE-PWCHex 2] Z_VV + Z_SV (vol cols)...")
    rank == 0 && flush(stdout)

    # Lockfree: each thread owns exclusive vol_col sub-range.
    if !isempty(local_vol_cols)
        n_vol_local   = length(local_vol_cols)
        vol_col_chunk = max(1, cld(n_vol_local, n_threads))

    Threads.@threads for tid = 1:n_threads
        col_lo   = first(local_vol_cols) + (tid - 1) * vol_col_chunk
        col_hi   = min(first(local_vol_cols) + tid * vol_col_chunk - 1, last(local_vol_cols))
        Z_ts_buf = zeros(CT, 3, 3)
        dyadG    = zeros(CT, 3, 3)

        for js in src_hexas_vol
            hex_s = hexas[js]
            κₛ   = hex_s.κ
            rad_s = cbrt(hex_s.volume)

            # Quick skip if hex_s has no vol col in this thread's range
            has_col = false
            @inbounds for ni = 1:3
                n = hex_s.inBfsID[ni]
                c = n_surf + n
                if n != 0 && col_lo <= c <= col_hi; has_col = true; break; end
            end
            has_col || continue

            # Z_VV: _pwc_dyad_kernel! for (hex_s, all hexas)
            for ti = 1:nhex
                hex_t = hexas[ti]

                dist_ts   = norm(hex_t.center - hex_s.center)
                rad_t     = cbrt(hex_t.volume)
                threshold = FT(3) * (rad_t + rad_s)
                is_far    = dist_ts > threshold

                _pwc_dyad_kernel!(
                    Z_ts_buf,
                    hex_t, hex_s,
                    is_far ? rq_hex_far[ti] : rq_hex_near[ti],
                    is_far ? rq_hex_far[js] : rq_hex_near[js],
                    is_far ? gq_hex_far : gq_hex_near,
                    is_far ? gq_hex_far : gq_hex_near,
                    is_far ? Nq_hex_far : Nq_hex,
                    is_far ? Nq_hex_far : Nq_hex,
                    k, k虏, jk, Jη₀divK, div4π,
                )

                @inbounds for ni = 1:3
                    n   = hex_s.inBfsID[ni]
                    col = n_surf + n
                    col_lo <= col <= col_hi || continue
                    for mi = 1:3
                        m   = hex_t.inBfsID[mi]
                        row = n_surf + m
                        if ti == js
                            val = Z_ts_buf[mi, ni] * κₛ
                            if mi == ni
                                selfImp = CT(1) / (im * omega) / hex_s.ε * hex_s.volume
                                val += selfImp
                            end
                            Z[row, col] += val   # no lock: col exclusive
                        else
                            Z[row, col] += Z_ts_buf[mi, ni] * κₛ  # no lock
                        end
                    end
                end
            end

            # Z_SV: source = hex_s, test = all tris
            r_q_v = rq_hex_near[js]
            for it = 1:ntri
                tri   = tris[it]
                r_q_t = rq_tri[it]

                for gi = 1:Nq_s
                    rgi = @view r_q_t[:, gi]
                    fill!(dyadG, zero(CT))

                    for gj = 1:Nq_hex
                        rgj = @view r_q_v[:, gj]
                        Rx = rgi[1] - rgj[1]
                        Ry = rgi[2] - rgj[2]
                        Rz = rgi[3] - rgj[3]
                        R  = sqrt(Rx^2 + Ry^2 + Rz^2)
                        R < FT(1e-10) && continue

                        divR  = FT(1) / R
                        jkpR  = (jk + divR) * divR
                        R虃x = Rx * divR; R虃y = Ry * divR; R虃z = Rz * divR
                        GR   = exp(-jk * R) * div4π * divR * gq_hex_near.weight[gj]

                        RR11 = R虃x*R虃x; RR12 = R虃x*R虃y; RR13 = R虃x*R虃z
                        RR22 = R虃y*R虃y; RR23 = R虃y*R虃z; RR33 = R虃z*R虃z

                        dyadG[1,1] += GR * ((1 - RR11) * k虏 - (1 - 3RR11) * jkpR)
                        dyadG[2,2] += GR * ((1 - RR22) * k虏 - (1 - 3RR22) * jkpR)
                        dyadG[3,3] += GR * ((1 - RR33) * k虏 - (1 - 3RR33) * jkpR)
                        od12 = GR * (-RR12 * k虏 + 3RR12 * jkpR)
                        od13 = GR * (-RR13 * k虏 + 3RR13 * jkpR)
                        od23 = GR * (-RR23 * k虏 + 3RR23 * jkpR)
                        dyadG[1,2] += od12; dyadG[2,1] += od12
                        dyadG[1,3] += od13; dyadG[3,1] += od13
                        dyadG[2,3] += od23; dyadG[3,2] += od23
                    end

                    for mi = 1:3
                        m = tri.inBfsID[mi]
                        m == 0 && continue
                        lm    = tri.edgel[mi]
                        freeV = tri.vertices[:, mi]
                        蟻mi   = SVector(rgi[1]-freeV[1], rgi[2]-freeV[2], rgi[3]-freeV[3])
                        temp  = gq_s.weight[gi] * lm / 2

                        for ni = 1:3
                            n   = hex_s.inBfsID[ni]
                            col = n_surf + n
                            col_lo <= col <= col_hi || continue
                            dot_sv = 蟻mi[1]*dyadG[1,ni] + 蟻mi[2]*dyadG[2,ni] + 蟻mi[3]*dyadG[3,ni]
                            z_sv   = temp * dot_sv * Jη₀divK * hex_s.volume * κₛ
                            Z[m, col] += z_sv   # no lock: col exclusive
                        end
                    end
                end
            end
        end
    end   # threads
    end   # if !isempty(local_vol_cols)

    sync!(Z)

    rank == 0 && println("SCFIE-PWCHex MPI Assembly Completed.")
    rank == 0 && flush(stdout)

    return Z
end
