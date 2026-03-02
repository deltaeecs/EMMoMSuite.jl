# VolumeAssembly.jl
# MPI 并行体积分方程组装：在 Parallel 命名空间扩展
# `Assembly.assemble_impedance_matrix_parallel`。
#
# 当前实现采用 MPI 列分区（每个 rank 只持有全局矩阵的连续列范围），
# 并在每个 rank 内使用多线程并行遍历源/测试几何元素。

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

# =============================================================================
# VEFIE + SWGBasis (MPI)
# =============================================================================
"""
    assemble_impedance_matrix_parallel(vefie, basis::SWGBasis, permittivities)

MPI 并行组装 VEFIE 阻抗矩阵（列分区）。

每个 rank 仅负责本地列 `local_cols`，并在本地线程中完成相关源四面体贡献写入。
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

    rank == 0 && println(
        "VEFIE-SWG MPI Assembly (column partition): N=$N, procs=$n_procs, threads=$(Threads.nthreads())",
    )
    MPI.Barrier(comm)
    println("  Rank $rank owns columns $(first(local_cols)):$(last(local_cols))")
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

    row_locks = [SpinLock() for _ = 1:N]
    next_idx  = Threads.Atomic{Int}(1)
    n_threads = Threads.nthreads()

    Threads.@threads for _ = 1:n_threads
        while true
            tidx = Threads.atomic_add!(next_idx, 1)
            tidx > length(src_tets) && break

            js      = src_tets[tidx]
            tet_s   = tetras[js]
            cache_s = basis_cache[js]

            Z_self = vefie_element_interaction_kernel(vefie, tet_s, tet_s, cache_s, cache_s)
            M      = vefie_mass_matrix_cached(vefie, tet_s, cache_s)
            @inbounds for i = 1:4
                m = tet_s.inBfsID[i]
                m == 0 && continue
                lock(row_locks[m])
                @inbounds for j = 1:4
                    n = tet_s.inBfsID[j]
                    n == 0 && continue
                    n in local_cols || continue
                    Z[m, n] += Z_self[i, j] + M[i, j]
                end
                unlock(row_locks[m])
            end

            for it = 1:ntet
                it == js && continue
                tet_t   = tetras[it]
                cache_t = basis_cache[it]
                Z_ts    = vefie_element_interaction_kernel(vefie, tet_t, tet_s, cache_t, cache_s)
                @inbounds for i = 1:4
                    m = tet_t.inBfsID[i]
                    m == 0 && continue
                    lock(row_locks[m])
                    @inbounds for j = 1:4
                        n = tet_s.inBfsID[j]
                        n == 0 && continue
                        n in local_cols || continue
                        Z[m, n] += Z_ts[i, j]
                    end
                    unlock(row_locks[m])
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

MPI 并行组装 SCFIE 耦合矩阵（列分区）。

全局索引布局:
- `1:n_surf`：表面 RWG 自由度
- `n_surf+1:n_total`：体 SWG 自由度

各子块按“源列是否归本 rank”写入本地列块，最终 `sync!(Z)` 同步。
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

    Z = mpiarray(CT, n_total, n_total; comm = comm)
    local_cols = Z.indices[2]

    surf_lo = max(first(local_cols), 1)
    surf_hi = min(last(local_cols), n_surf)
    vol_lo  = max(first(local_cols), n_surf + 1)
    vol_hi  = min(last(local_cols), n_total)

    local_surf_cols = surf_lo <= surf_hi ? (surf_lo:surf_hi) : (1:0)
    local_vol_cols  = vol_lo <= vol_hi ? (vol_lo:vol_hi) : (1:0)

    rank == 0 && println(
        "SCFIE MPI Assembly (column partition): n_surf=$n_surf, n_vol=$n_vol, n_total=$n_total, procs=$n_procs, threads=$(Threads.nthreads())",
    )
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

    # 行级锁数组：n_total 把锁，每输出行一把，竞争范围从「全局/子块」降至「单行」
    row_locks = [SpinLock() for _ = 1:n_total]
    n_threads = Threads.nthreads()

    rank == 0 && println("  [SCFIE 1+3b] Z_SS + Z_VS (surf cols)...")
    rank == 0 && flush(stdout)

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

            for t_tst = 1:ntri
                tri_t = tris[t_tst]
                fill!(l_efie, zero(CT))
                fill!(l_mfie, zero(CT))
                fill!(l_loc, zero(CT))
                efie_interaction!(l_efie, cfie.efie, tri_t, tri_s)
                mfie_interaction!(l_mfie, cfie.mfie, tri_t, tri_s)
                alpha_val = scfie.alpha
                @. l_loc = alpha_val * l_efie + (1 - alpha_val) * l_mfie
                @inbounds for j = 1:3
                    col = tri_s.inBfsID[j]
                    col == 0 && continue
                    col in local_surf_cols || continue
                    bf_s      = surf_basis.functions[col]
                    sign_s    = (bf_s.support[1] == t_src) ? bf_s.signs[1] : bf_s.signs[2]
                    @inbounds for i = 1:3
                        row = tri_t.inBfsID[i]
                        row == 0 && continue
                        bf_t      = surf_basis.functions[row]
                        sign_t    = (bf_t.support[1] == t_tst) ? bf_t.signs[1] : bf_t.signs[2]
                        lock(row_locks[row])
                        Z[row, col] += l_loc[i, j] * sign_t * sign_s
                        unlock(row_locks[row])
                    end
                end
            end

            for js = 1:ntet
                tet       = tetras[js]
                Z_sv      = scfie_sv_only_interaction(scfie, tri_s, tet)
                kappa_inv = iszero(tet.κ) ? zero(CT) : CT(1) / tet.κ
                @inbounds for i = 1:3
                    col = tri_s.inBfsID[i]
                    col == 0 && continue
                    col in local_surf_cols || continue
                    @inbounds for j = 1:4
                        n = tet.inBfsID[j]
                        n == 0 && continue
                        row = n_surf + n
                        lock(row_locks[row])
                        Z[row, col] += Z_sv[i, j] * kappa_inv
                        unlock(row_locks[row])
                    end
                end
            end
        end
    end

    rank == 0 && println("  [SCFIE 2+3a] Z_VV + Z_SV (vol cols)...")
    rank == 0 && flush(stdout)

    next_vv = Threads.Atomic{Int}(1)
    Threads.@threads for _ = 1:n_threads
        while true
            tidx = Threads.atomic_add!(next_vv, 1)
            tidx > length(src_tets_vol) && break
            js      = src_tets_vol[tidx]
            tet_s   = tetras[js]
            cache_s = basis_cache[js]

            for it = 1:ntet
                tet_t   = tetras[it]
                cache_t = basis_cache[it]
                if it == js
                    Z_self = vefie_element_interaction_kernel(
                        vefie_inner,
                        tet_s,
                        tet_s,
                        cache_s,
                        cache_s,
                    )
                    M_mat = vefie_mass_matrix_cached(vefie_inner, tet_s, cache_s)
                    @inbounds for j = 1:4
                        n = tet_s.inBfsID[j]
                        n == 0 && continue
                        (n_surf + n) in local_vol_cols || continue
                        @inbounds for i = 1:4
                            m = tet_s.inBfsID[i]
                            m == 0 && continue
                            row = n_surf + m
                            lock(row_locks[row])
                            Z[row, n_surf + n] += Z_self[i, j] + M_mat[i, j]
                            unlock(row_locks[row])
                        end
                    end
                else
                    Z_ts = vefie_element_interaction_kernel(
                        vefie_inner,
                        tet_t,
                        tet_s,
                        cache_t,
                        cache_s,
                    )
                    @inbounds for j = 1:4
                        n = tet_s.inBfsID[j]
                        n == 0 && continue
                        (n_surf + n) in local_vol_cols || continue
                        @inbounds for i = 1:4
                            m = tet_t.inBfsID[i]
                            m == 0 && continue
                            row = n_surf + m
                            lock(row_locks[row])
                            Z[row, n_surf + n] += Z_ts[i, j]
                            unlock(row_locks[row])
                        end
                    end
                end
            end

            for it_tri = 1:ntri
                tri  = tris[it_tri]
                Z_sv = scfie_sv_only_interaction(scfie, tri, tet_s)
                @inbounds for i = 1:3
                    m = tri.inBfsID[i]
                    m == 0 && continue
                    @inbounds for j = 1:4
                        n = tet_s.inBfsID[j]
                        n == 0 && continue
                        (n_surf + n) in local_vol_cols || continue
                        lock(row_locks[m])
                        Z[m, n_surf + n] += Z_sv[i, j]
                        unlock(row_locks[m])
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
