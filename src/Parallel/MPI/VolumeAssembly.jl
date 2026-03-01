# VolumeAssembly.jl
# MPI 骞惰浣撶Н绉垎鏂圭▼缁勮 鈥?鐩存帴鎵╁睍 Assembly.assemble_impedance_matrix_parallel
#
# 绛栫暐: 寰幆鍒嗛厤澶栧眰娴嬭瘯鍥涢潰浣?it, 姣忎釜 MPI 杩涚▼澶勭悊 (it-1) % n_procs == rank 鐨?it.
# 鍒╃敤 Z 瀵圭О鎬?(js > it): 姣忔璁＄畻 Z_ts, 鍚屾椂濉厖 Z[m,n] 鍜?Z[n,m].
# 鏈€鍚?MPI.Allreduce! 鍦ㄦ墍鏈夎繘绋嬮棿瀵瑰眬閮ㄧ煩闃垫眰鍜?
#
# File included by Parallel.jl (no module wrapper: methods live in Parallel namespace)

import .Assembly: assemble_impedance_matrix_parallel

using ..IntegralEquations.VEFIEModule:
    VEFIE,
    get_tetrahedra_info,
    precompute_vefie_basis,
    vefie_element_interaction_kernel,
    vefie_mass_matrix_cached

using ..IntegralEquations.SCFIEModule:
    SCFIE,
    scfie_sv_only_interaction,
    assemble_fss_boundary_correction_sparse

using ..BasisFunctions: get_triangles_info, SWGBasis, RWGBasis

using ..IntegralEquations.EFIEModule: efie_interaction!
using ..IntegralEquations.MFIEModule: mfie_interaction!
using ..IntegralEquations.CFIEModule: CFIE

using SparseArrays

# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
#  VEFIE + SWGBasis  (MPI 骞惰)
# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
"""
    assemble_impedance_matrix_parallel(vefie, basis::SWGBasis, permittivities)

MPI 骞惰缁勮 VEFIE 闃绘姉鐭╅樀.

**鍒嗛厤绛栫暐**: 娴嬭瘯鍥涢潰浣?`it` 鎸夊惊鐜垎閰嶇粰鍚勮繘绋?
瀵圭О鍒╃敤: 姣忎釜 (it, js>it) 瀵圭敱 `it` 鐨勬墍鏈夎€呯嫭绔嬭绠?
鏈€鍚?`MPI.Allreduce!` 璺ㄨ繘绋嬫眰鍜屽緱鍒板畬鏁?Z.
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

    # Phase 14.1: 列分区策略 — 每进程仅分配 N × (N/P) 列，消除 Allreduce 全复制
    Z = mpiarray(CT, N, N; comm = comm)
    local_cols    = Z.indices[2]          # 本进程负责的列范围 (UnitRange)
    local_col_set = Set{Int}(local_cols)  # O(1) 成员查询

    rank == 0 && println("VEFIE-SWG MPI Assembly (column partition): N=$N, procs=$n_procs, threads=$(Threads.nthreads())")
    MPI.Barrier(comm)
    println("  Rank $rank owns columns $(first(local_cols)):$(last(local_cols))")
    MPI.Barrier(comm)
    flush(stdout)

    tetras      = get_tetrahedra_info(basis.mesh, basis, permittivities)
    ntet        = length(tetras)
    basis_cache = precompute_vefie_basis(vefie, tetras)

    # 找出源四面体：其 SWG 基函数 ID 至少有一个落在 local_cols 内
    src_tets = Int[]
    for js = 1:ntet
        tet = tetras[js]
        for j = 1:4
            n = tet.inBfsID[j]
            if n != 0 && n in local_col_set
                push!(src_tets, js)
                break
            end
        end
    end

    lockZ     = SpinLock()
    next_idx  = Threads.Atomic{Int}(1)
    n_threads = Threads.nthreads()

    Threads.@threads for _ = 1:n_threads
        while true
            tidx = Threads.atomic_add!(next_idx, 1)
            tidx > length(src_tets) && break

            js      = src_tets[tidx]
            tet_s   = tetras[js]
            cache_s = basis_cache[js]

            # 自项 (test == source == tet_s)
            Z_self = vefie_element_interaction_kernel(vefie, tet_s, tet_s, cache_s, cache_s)
            M      = vefie_mass_matrix_cached(vefie, tet_s, cache_s)
            lock(lockZ)
            @inbounds for j = 1:4
                n = tet_s.inBfsID[j]; n == 0 && continue
                n in local_col_set || continue
                @inbounds for i = 1:4
                    m = tet_s.inBfsID[i]; m == 0 && continue
                    Z[m, n] += Z_self[i, j] + M[i, j]
                end
            end
            unlock(lockZ)

            # 交叉项：所有 it != js 测试四面体
            for it = 1:ntet
                it == js && continue
                tet_t   = tetras[it]
                cache_t = basis_cache[it]
                # Z_ts[i,j] = 测试 tet_t 行 i、源 tet_s 列 j 的贡献
                Z_ts = vefie_element_interaction_kernel(vefie, tet_t, tet_s, cache_t, cache_s)
                lock(lockZ)
                @inbounds for j = 1:4
                    n = tet_s.inBfsID[j]; n == 0 && continue
                    n in local_col_set || continue
                    @inbounds for i = 1:4
                        m = tet_t.inBfsID[i]; m == 0 && continue
                        Z[m, n] += Z_ts[i, j]
                    end
                end
                unlock(lockZ)
            end
        end
    end

    sync!(Z)  # 同步 ghost 数据（列分区下通常为空操作）

    rank == 0 && println("VEFIE-SWG MPI Assembly Completed.")
    rank == 0 && flush(stdout)

    return Z  # 返回 MPIMatrix（列分布式）
end


# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
#  SCFIE + RWGBasis + SWGBasis  (MPI 骞惰)
# 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
"""
    assemble_impedance_matrix_parallel(scfie, surf_basis::RWGBasis, vol_basis::SWGBasis)

MPI 骞惰缁勮 SCFIE 鑰﹀悎闃绘姉鐭╅樀.

鐭╅樀绱㈠紩甯冨眬:
  琛?鍒?1..n_surf         鈫?闈?DOFs (RWG)
  琛?鍒?n_surf+1..n_total 鈫?浣?DOFs (SWG)

鍒嗛厤绛栫暐:
  Z_SS: 娴嬭瘯闈笁瑙?it 鎸夊惊鐜垎閰?
  Z_VV: 娴嬭瘯鍥涢潰浣?it 鎸夊惊鐜垎閰?(瀵圭О)
  Z_SV/Z_VS: 闈?it 鎸夊惊鐜垎閰? 鍚屾椂鐢ㄤ簰鏄撴€у～ Z_VS
  Fss: 姣忚繘绋嬬嫭绔嬭绠? 闄や互 n_procs 鍚?Allreduce 鑷劧姹囧悎
"""
function assemble_impedance_matrix_parallel(
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::SWGBasis,
)
    comm    = MPI.COMM_WORLD
    rank    = MPI.Comm_rank(comm)
    n_procs = MPI.Comm_size(comm)

    FT = typeof(real(scfie.freq))
    CT = Complex{FT}

    n_surf  = num_basis(surf_basis)
    n_vol   = num_basis(vol_basis)
    n_total = n_surf + n_vol

    # Phase 14.1: 列分区策略 — 混合空间全局矩阵 n_total × n_total，列分布
    Z = mpiarray(CT, n_total, n_total; comm = comm)
    local_cols        = Z.indices[2]
    local_col_set     = Set{Int}(local_cols)
    # 拆分为面列 (1..n_surf) 和体列 (n_surf+1..n_total)
    local_surf_cols   = Set{Int}(c for c in local_cols if c <= n_surf)
    local_vol_cols    = Set{Int}(c for c in local_cols if c > n_surf)

    rank == 0 && println("SCFIE MPI Assembly (column partition): n_surf=$n_surf, n_vol=$n_vol, n_total=$n_total, procs=$n_procs, threads=$(Threads.nthreads())")
    MPI.Barrier(comm)
    println("  Rank $rank: surf_cols=$(length(local_surf_cols)), vol_cols=$(length(local_vol_cols))")
    MPI.Barrier(comm)
    flush(stdout)

    tris        = get_triangles_info(surf_basis.mesh, surf_basis)
    ntri        = length(tris)
    tetras      = get_tetrahedra_info(vol_basis.mesh, vol_basis, scfie.permittivities)
    ntet        = length(tetras)
    vefie_inner = VEFIE(scfie.freq, scfie.permittivities)
    basis_cache = precompute_vefie_basis(vefie_inner, tetras)
    cfie        = CFIE(scfie.freq, scfie.alpha)  # 用于 Z_SS 块的 EFIE+MFIE 交互

    # 找本进程负责的源三角形（surf 列）和源四面体（vol 列）
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

    lockZ     = SpinLock()
    n_threads = Threads.nthreads()

    # ── [1] Z_SS (surf行 surf列) + [3b] Z_VS (vol行 surf列) ───────────────────
    # 迭代源三角形 → 填写 col ∈ local_surf_cols 的全部行
    rank == 0 && println("  [SCFIE 1+3b] Z_SS + Z_VS (surf cols)..."); rank == 0 && flush(stdout)

    next_ss = Threads.Atomic{Int}(1)
    Threads.@threads for _ = 1:n_threads
        l_efie = zeros(CT, 3, 3)
        l_mfie = zeros(CT, 3, 3)
        l_loc  = zeros(CT, 3, 3)
        while true
            tidx = Threads.atomic_add!(next_ss, 1)
            tidx > length(src_tris_surf) && break
            t_src = src_tris_surf[tidx]
            tri_s = tris[t_src]

            # Z_SS: 对所有测试三角形
            for t_tst = 1:ntri
                tri_t = tris[t_tst]
                fill!(l_efie, zero(CT)); fill!(l_mfie, zero(CT)); fill!(l_loc, zero(CT))
                efie_interaction!(l_efie, cfie.efie, tri_t, tri_s)
                mfie_interaction!(l_mfie, cfie.mfie, tri_t, tri_s)
                alpha_val = scfie.alpha
                @. l_loc = alpha_val * l_efie + (1 - alpha_val) * l_mfie
                lock(lockZ)
                @inbounds for j = 1:3
                    col = tri_s.inBfsID[j]; col == 0 && continue
                    col in local_surf_cols || continue
                    bf_s   = surf_basis.functions[col]
                    sign_s = (bf_s.support[1] == t_src) ? bf_s.signs[1] : bf_s.signs[2]
                    @inbounds for i = 1:3
                        row = tri_t.inBfsID[i]; row == 0 && continue
                        bf_t   = surf_basis.functions[row]
                        sign_t = (bf_t.support[1] == t_tst) ? bf_t.signs[1] : bf_t.signs[2]
                        Z[row, col] += l_loc[i, j] * sign_t * sign_s
                    end
                end
                unlock(lockZ)
            end

            # Z_VS: row = n_surf+n (vol), col = m (surf)
            # scfie_sv_only_interaction(scfie, tri_test, tet_source) → Z_SV_elem[tri_i, tet_j]
            # Z_VS[n_surf+n, m] = Z_SV[m, n_surf+n] * κ_inv
            # 这里 tri_s 扮演 test（行），col=tri_s.inBfsID[i]=m 扮演 surf col
            for js = 1:ntet
                tet   = tetras[js]
                Z_sv  = scfie_sv_only_interaction(scfie, tri_s, tet)
                kappa_inv = iszero(tet.κ) ? zero(CT) : CT(1) / tet.κ
                lock(lockZ)
                @inbounds for i = 1:3
                    col = tri_s.inBfsID[i]; col == 0 && continue
                    col in local_surf_cols || continue
                    @inbounds for j = 1:4
                        n = tet.inBfsID[j]; n == 0 && continue
                        Z[n_surf + n, col] += Z_sv[i, j] * kappa_inv
                    end
                end
                unlock(lockZ)
            end
        end
    end

    # ── [2] Z_VV (vol行 vol列) + [3a] Z_SV (surf行 vol列) ────────────────────
    # 迭代源四面体 → 填写 col ∈ local_vol_cols 的全部行
    rank == 0 && println("  [SCFIE 2+3a] Z_VV + Z_SV (vol cols)..."); rank == 0 && flush(stdout)

    next_vv = Threads.Atomic{Int}(1)
    Threads.@threads for _ = 1:n_threads
        while true
            tidx = Threads.atomic_add!(next_vv, 1)
            tidx > length(src_tets_vol) && break
            js      = src_tets_vol[tidx]
            tet_s   = tetras[js]
            cache_s = basis_cache[js]

            # Z_VV: 对所有测试四面体（列分区下无对称利用）
            for it = 1:ntet
                tet_t   = tetras[it]
                cache_t = basis_cache[it]
                if it == js
                    Z_self = vefie_element_interaction_kernel(vefie_inner, tet_s, tet_s, cache_s, cache_s)
                    M_mat  = vefie_mass_matrix_cached(vefie_inner, tet_s, cache_s)
                    lock(lockZ)
                    @inbounds for j = 1:4
                        n = tet_s.inBfsID[j]; n == 0 && continue
                        (n_surf + n) in local_vol_cols || continue
                        @inbounds for i = 1:4
                            m = tet_s.inBfsID[i]; m == 0 && continue
                            Z[n_surf + m, n_surf + n] += Z_self[i, j] + M_mat[i, j]
                        end
                    end
                    unlock(lockZ)
                else
                    Z_ts = vefie_element_interaction_kernel(vefie_inner, tet_t, tet_s, cache_t, cache_s)
                    lock(lockZ)
                    @inbounds for j = 1:4
                        n = tet_s.inBfsID[j]; n == 0 && continue
                        (n_surf + n) in local_vol_cols || continue
                        @inbounds for i = 1:4
                            m = tet_t.inBfsID[i]; m == 0 && continue
                            Z[n_surf + m, n_surf + n] += Z_ts[i, j]
                        end
                    end
                    unlock(lockZ)
                end
            end

            # Z_SV: row = m (surf), col = n_surf+n (vol)
            for it_tri = 1:ntri
                tri   = tris[it_tri]
                Z_sv  = scfie_sv_only_interaction(scfie, tri, tet_s)
                lock(lockZ)
                @inbounds for i = 1:3
                    m = tri.inBfsID[i]; m == 0 && continue
                    @inbounds for j = 1:4
                        n = tet_s.inBfsID[j]; n == 0 && continue
                        (n_surf + n) in local_vol_cols || continue
                        Z[m, n_surf + n] += Z_sv[i, j]
                    end
                end
                unlock(lockZ)
            end
        end
    end

    # ── [4] Fss 边界修正（稀疏，仅写本地列）──────────────────────────────────
    rank == 0 && println("  [SCFIE Fss] boundary correction..."); rank == 0 && flush(stdout)

    Z_fss = assemble_fss_boundary_correction_sparse(scfie, surf_basis, vol_basis)
    rows_f, cols_f, vals_f = findnz(Z_fss)
    @inbounds for k = 1:length(vals_f)
        c = cols_f[k]
        c in local_col_set || continue
        Z[rows_f[k], c] += vals_f[k]  # 不再除以 n_procs（列分区下每列仅一个进程负责）
    end

    sync!(Z)

    rank == 0 && println("SCFIE MPI Assembly Completed.")
    rank == 0 && flush(stdout)

    return Z  # MPIMatrix（列分布式）
end

