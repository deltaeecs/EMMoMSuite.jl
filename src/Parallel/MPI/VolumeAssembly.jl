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

# =============================================================================
# VEFIE + PWCBasis (MPI)
# =============================================================================
"""
    assemble_impedance_matrix_parallel(vefie::VEFIE, basis::PWCBasis, permittivities)

MPI 并行组装 VEFIE 阻抗矩阵（PWC 四面体，列分区）。

每个四面体贡献 3 个 DOF（x/y/z 分量）。
无对称性利用：源列过滤决定每个 rank 处理的源四面体，
测试四面体遍历全部（含自项和互项）。

# Legacy Parity
与串行 `assemble_impedance_matrix(vefie, basis::PWCBasis)` 结果一致。
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

    rank == 0 && println(
        "VEFIE-PWC MPI Assembly (column partition): N=$N, procs=$n_procs, threads=$(Threads.nthreads())",
    )
    MPI.Barrier(comm)
    println("  Rank $rank owns columns $(first(local_cols)):$(last(local_cols))")
    MPI.Barrier(comm)
    flush(stdout)

    tetras = get_tetrahedra_info(basis.mesh, basis, permittivities)
    ntet   = length(tetras)

    # Constants
    k          = vefie.k
    k²         = k^2
    jk         = im * k
    omega      = 2π * vefie.freq
    mu0        = FT(4π * 1e-7)
    eps0       = FT(8.854187817e-12)
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

    row_locks = [SpinLock() for _ = 1:N]
    next_idx  = Threads.Atomic{Int}(1)
    n_threads = Threads.nthreads()

    Threads.@threads for _ = 1:n_threads
        Z_ts_buf = zeros(CT, 3, 3)

        while true
            tidx = Threads.atomic_add!(next_idx, 1)
            tidx > length(src_tets) && break

            js    = src_tets[tidx]
            tet_s = tetras[js]
            κₛ    = tet_s.κ

            # --- Self-term (ti == js) ---
            rad_s = cbrt(tet_s.volume)
            # Self: always near-field
            _pwc_dyad_kernel!(
                Z_ts_buf,
                tet_s, tet_s,
                rq_near[js], rq_near[js],
                gq, gq,
                Nq, Nq,
                k, k², jk, Jη₀divK, div4π,
            )
            selfImp = CT(1) / (im * omega) / tet_s.ε * tet_s.volume
            @inbounds for ni = 1:3
                n = tet_s.inBfsID[ni]
                (n == 0 || !(n in local_cols)) && continue
                lock(row_locks[n])
                for mi = 1:3
                    m = tet_s.inBfsID[mi]
                    m == 0 && continue
                    val = Z_ts_buf[mi, ni] * κₛ
                    if mi == ni
                        val += selfImp
                    end
                    Z[m, n] += val
                end
                unlock(row_locks[n])
            end

            # --- Off-diagonal terms (ti != js) ---
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
                        k, k², jk, Jη₀divK, div4π,
                    )
                else
                    _pwc_dyad_kernel!(
                        Z_ts_buf,
                        tet_t, tet_s,
                        rq_near[ti], rq_near[js],
                        gq, gq,
                        Nq, Nq,
                        k, k², jk, Jη₀divK, div4π,
                    )
                end

                # Write: test row m, source col n ∈ local_cols
                @inbounds for ni = 1:3
                    n = tet_s.inBfsID[ni]
                    (n == 0 || !(n in local_cols)) && continue
                    lock(row_locks[n])
                    for mi = 1:3
                        m = tet_t.inBfsID[mi]
                        m == 0 && continue
                        Z[m, n] += Z_ts_buf[mi, ni] * κₛ
                    end
                    unlock(row_locks[n])
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

使用 `vefie.permittivities` 的便捷重载。
"""
function assemble_impedance_matrix_parallel(vefie::VEFIE, basis::PWCBasis)
    return assemble_impedance_matrix_parallel(vefie, basis, vefie.permittivities)
end

# =============================================================================
# VEFIE + PWCHexBasis (MPI)
# =============================================================================
"""
    assemble_impedance_matrix_parallel(vefie::VEFIE, basis::PWCHexBasis)

MPI 并行组装 VEFIE 阻抗矩阵（PWC 六面体，列分区）。

每个六面体贡献 3 个 DOF（x/y/z）。使用 Hexahedron 8 点（近场）/1 点（远场）高斯积分。

# Legacy Parity
与串行 `assemble_impedance_matrix(vefie, basis::PWCHexBasis)` 结果一致。
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

    rank == 0 && println(
        "VEFIE-PWC-Hex MPI Assembly (column partition): N=$N, procs=$n_procs, threads=$(Threads.nthreads())",
    )
    MPI.Barrier(comm)
    println("  Rank $rank owns columns $(first(local_cols)):$(last(local_cols))")
    MPI.Barrier(comm)
    flush(stdout)

    hexas = get_hexahedra_info(basis.mesh, basis, vefie.permittivities)
    nhex  = length(hexas)

    # Constants
    k         = vefie.k
    k²        = k^2
    jk        = im * k
    omega     = 2π * vefie.freq
    mu0       = FT(4π * 1e-7)
    eps0      = FT(8.854187817e-12)
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

    row_locks = [SpinLock() for _ = 1:N]
    next_idx  = Threads.Atomic{Int}(1)
    n_threads = Threads.nthreads()

    Threads.@threads for _ = 1:n_threads
        Z_ts_buf = zeros(CT, 3, 3)

        while true
            tidx = Threads.atomic_add!(next_idx, 1)
            tidx > length(src_hexas) && break

            js    = src_hexas[tidx]
            hex_s = hexas[js]
            κₛ    = hex_s.κ

            # --- Self-term ---
            rad_s = cbrt(hex_s.volume)
            _pwc_dyad_kernel!(
                Z_ts_buf,
                hex_s, hex_s,
                rq_near[js], rq_near[js],
                gq_hex, gq_hex,
                Nq_hex, Nq_hex,
                k, k², jk, Jη₀divK, div4π,
            )
            selfImp = CT(1) / (im * omega) / hex_s.ε * hex_s.volume
            @inbounds for ni = 1:3
                n = hex_s.inBfsID[ni]
                (n == 0 || !(n in local_cols)) && continue
                lock(row_locks[n])
                for mi = 1:3
                    m = hex_s.inBfsID[mi]
                    m == 0 && continue
                    val = Z_ts_buf[mi, ni] * κₛ
                    if mi == ni
                        val += selfImp
                    end
                    Z[m, n] += val
                end
                unlock(row_locks[n])
            end

            # --- Off-diagonal terms ---
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
                        k, k², jk, Jη₀divK, div4π,
                    )
                else
                    _pwc_dyad_kernel!(
                        Z_ts_buf,
                        hex_t, hex_s,
                        rq_near[ti], rq_near[js],
                        gq_hex, gq_hex,
                        Nq_hex, Nq_hex,
                        k, k², jk, Jη₀divK, div4π,
                    )
                end

                @inbounds for ni = 1:3
                    n = hex_s.inBfsID[ni]
                    (n == 0 || !(n in local_cols)) && continue
                    lock(row_locks[n])
                    for mi = 1:3
                        m = hex_t.inBfsID[mi]
                        m == 0 && continue
                        Z[m, n] += Z_ts_buf[mi, ni] * κₛ
                    end
                    unlock(row_locks[n])
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

MPI 并行组装 SCFIE 耦合矩阵（列分区），体积分方程基函数为 PWCBasis（四面体）。

全局索引布局:
- `1:n_surf`  : 表面 RWG 自由度
- `n_surf+1:n_total` : 体 PWC 自由度（每四面体 3 个 x/y/z 分量）

组装分两阶段：
1. **surf cols**（源在表面）: 填充 Z_SS + Z_VS
2. **vol cols** （源在体积）: 填充 Z_VV + Z_SV

无 Fss 边界修正（PWC 无半基函数）。

# Legacy Parity
与串行 `assemble_impedance_matrix(scfie, surf_basis::RWGBasis, vol_basis::PWCBasis)` 一致。
"""
function assemble_impedance_matrix_parallel(
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::PWCBasis,
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
    local_vol_cols  = vol_lo  <= vol_hi  ? (vol_lo:vol_hi)   : (1:0)

    rank == 0 && println(
        "SCFIE-PWC MPI Assembly (column partition): n_surf=$n_surf, n_vol=$n_vol, n_total=$n_total, procs=$n_procs, threads=$(Threads.nthreads())",
    )
    MPI.Barrier(comm)
    println("  Rank $rank: surf_cols=$(length(local_surf_cols)), vol_cols=$(length(local_vol_cols))")
    MPI.Barrier(comm)
    flush(stdout)

    tris   = get_triangles_info(surf_basis.mesh, surf_basis)
    ntri   = length(tris)
    tetras = get_tetrahedra_info(vol_basis.mesh, vol_basis, scfie.permittivities)
    ntet   = length(tetras)

    cfie = CFIE(scfie.freq, scfie.alpha)

    # Physical constants (for coupling dyad and Z_VV kernel)
    k         = scfie.k
    k²        = k^2
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

    row_locks = [SpinLock() for _ = 1:n_total]
    n_threads = Threads.nthreads()

    # ─── Phase 1: Z_SS + Z_VS (surf cols) ────────────────────────────────────
    rank == 0 && println("  [SCFIE-PWC 1] Z_SS + Z_VS (surf cols)...")
    rank == 0 && flush(stdout)

    next_ss = Threads.Atomic{Int}(1)
    Threads.@threads for _ = 1:n_threads
        l_efie = zeros(CT, 3, 3)
        l_mfie = zeros(CT, 3, 3)
        l_loc  = zeros(CT, 3, 3)
        dyadG  = zeros(CT, 3, 3)

        while true
            tidx = Threads.atomic_add!(next_ss, 1)
            tidx > length(src_tris_surf) && break

            t_src = src_tris_surf[tidx]
            tri_s = tris[t_src]

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
                    col == 0 && continue
                    col in local_surf_cols || continue
                    bf_s   = surf_basis.functions[col]
                    sign_s = (bf_s.support[1] == t_src) ? bf_s.signs[1] : bf_s.signs[2]
                    @inbounds for i = 1:3
                        row = tri_t.inBfsID[i]
                        row == 0 && continue
                        bf_t   = surf_basis.functions[row]
                        sign_t = (bf_t.support[1] == t_tst) ? bf_t.signs[1] : bf_t.signs[2]
                        lock(row_locks[row])
                        Z[row, col] += l_loc[i, j] * sign_t * sign_s
                        unlock(row_locks[row])
                    end
                end
            end

            # Z_VS: source = tri_s (col = m ∈ local_surf_cols), test = all tets
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
                        R̂x = Rx * divR; R̂y = Ry * divR; R̂z = Rz * divR
                        GR   = exp(-jk * R) * div4π * divR * gq_v.weight[gj]

                        RR11 = R̂x*R̂x; RR12 = R̂x*R̂y; RR13 = R̂x*R̂z
                        RR22 = R̂y*R̂y; RR23 = R̂y*R̂z; RR33 = R̂z*R̂z

                        dyadG[1,1] += GR * ((1 - RR11) * k² - (1 - 3RR11) * jkpR)
                        dyadG[2,2] += GR * ((1 - RR22) * k² - (1 - 3RR22) * jkpR)
                        dyadG[3,3] += GR * ((1 - RR33) * k² - (1 - 3RR33) * jkpR)
                        od12 = GR * (-RR12 * k² + 3RR12 * jkpR)
                        od13 = GR * (-RR13 * k² + 3RR13 * jkpR)
                        od23 = GR * (-RR23 * k² + 3RR23 * jkpR)
                        dyadG[1,2] += od12; dyadG[2,1] += od12
                        dyadG[1,3] += od13; dyadG[3,1] += od13
                        dyadG[2,3] += od23; dyadG[3,2] += od23
                    end

                    for mi = 1:3
                        m = tri_s.inBfsID[mi]
                        m == 0 && continue
                        m in local_surf_cols || continue
                        lm    = tri_s.edgel[mi]
                        freeV = tri_s.vertices[:, mi]
                        ρmi   = SVector(rgi[1] - freeV[1], rgi[2] - freeV[2], rgi[3] - freeV[3])
                        temp  = gq_s.weight[gi] * lm / 2

                        for ni = 1:3
                            n     = tet.inBfsID[ni]
                            # Z_VS: test=vol, source=surf → row = n_surf+n, col = m
                            dot_vs = ρmi[1]*dyadG[ni,1] + ρmi[2]*dyadG[ni,2] + ρmi[3]*dyadG[ni,3]
                            z_vs   = temp * dot_vs * Jη₀divK * Vs
                            row_vs = n_surf + n
                            lock(row_locks[row_vs])
                            Z[row_vs, m] += z_vs
                            unlock(row_locks[row_vs])
                        end
                    end
                end
            end
        end
    end

    # ─── Phase 2: Z_VV + Z_SV (vol cols) ─────────────────────────────────────
    rank == 0 && println("  [SCFIE-PWC 2] Z_VV + Z_SV (vol cols)...")
    rank == 0 && flush(stdout)

    next_vv = Threads.Atomic{Int}(1)
    Threads.@threads for _ = 1:n_threads
        Z_ts_buf = zeros(CT, 3, 3)
        dyadG    = zeros(CT, 3, 3)

        while true
            tidx = Threads.atomic_add!(next_vv, 1)
            tidx > length(src_tets_vol) && break

            js    = src_tets_vol[tidx]
            tet_s = tetras[js]
            κₛ    = tet_s.κ
            rad_s = cbrt(tet_s.volume)

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
                    k, k², jk, Jη₀divK, div4π,
                )

                @inbounds for ni = 1:3
                    n = tet_s.inBfsID[ni]
                    col = n_surf + n
                    col in local_vol_cols || continue
                    lock(row_locks[col])
                    for mi = 1:3
                        m   = tet_t.inBfsID[mi]
                        row = n_surf + m
                        if ti == js
                            val = Z_ts_buf[mi, ni] * κₛ
                            if mi == ni
                                selfImp = CT(1) / (im * omega) / tet_s.ε * tet_s.volume
                                val += selfImp
                            end
                            Z[row, col] += val
                        else
                            Z[row, col] += Z_ts_buf[mi, ni] * κₛ
                        end
                    end
                    unlock(row_locks[col])
                end
            end

            # Z_SV: source = tet_s (col = n_surf+n ∈ local_vol_cols), test = all tris
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
                        R̂x = Rx * divR; R̂y = Ry * divR; R̂z = Rz * divR
                        GR   = exp(-jk * R) * div4π * divR * gq_v.weight[gj]

                        RR11 = R̂x*R̂x; RR12 = R̂x*R̂y; RR13 = R̂x*R̂z
                        RR22 = R̂y*R̂y; RR23 = R̂y*R̂z; RR33 = R̂z*R̂z

                        dyadG[1,1] += GR * ((1 - RR11) * k² - (1 - 3RR11) * jkpR)
                        dyadG[2,2] += GR * ((1 - RR22) * k² - (1 - 3RR22) * jkpR)
                        dyadG[3,3] += GR * ((1 - RR33) * k² - (1 - 3RR33) * jkpR)
                        od12 = GR * (-RR12 * k² + 3RR12 * jkpR)
                        od13 = GR * (-RR13 * k² + 3RR13 * jkpR)
                        od23 = GR * (-RR23 * k² + 3RR23 * jkpR)
                        dyadG[1,2] += od12; dyadG[2,1] += od12
                        dyadG[1,3] += od13; dyadG[3,1] += od13
                        dyadG[2,3] += od23; dyadG[3,2] += od23
                    end

                    for mi = 1:3
                        m = tri.inBfsID[mi]
                        m == 0 && continue
                        lm    = tri.edgel[mi]
                        freeV = tri.vertices[:, mi]
                        ρmi   = SVector(rgi[1] - freeV[1], rgi[2] - freeV[2], rgi[3] - freeV[3])
                        temp  = gq_s.weight[gi] * lm / 2

                        for ni = 1:3
                            n   = tet_s.inBfsID[ni]
                            col = n_surf + n
                            col in local_vol_cols || continue
                            # Z_SV: test=surf, source=vol → row = m, col = n_surf+n
                            dot_sv = ρmi[1]*dyadG[1,ni] + ρmi[2]*dyadG[2,ni] + ρmi[3]*dyadG[3,ni]
                            z_sv   = temp * dot_sv * Jη₀divK * tet_s.volume * κₛ
                            lock(row_locks[m])
                            Z[m, col] += z_sv
                            unlock(row_locks[m])
                        end
                    end
                end
            end
        end
    end

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

MPI 并行组装 SCFIE 耦合矩阵（列分区），体积分方程基函数为 PWCHexBasis（六面体）。

与 PWCBasis 版本结构相同，区别在于：六面体 GQ（8点近场/1点远场）在本地创建。

# Legacy Parity
与串行 `assemble_impedance_matrix(scfie, surf_basis::RWGBasis, vol_basis::PWCHexBasis)` 一致。
"""
function assemble_impedance_matrix_parallel(
    scfie::SCFIE,
    surf_basis::RWGBasis,
    vol_basis::PWCHexBasis,
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
    local_vol_cols  = vol_lo  <= vol_hi  ? (vol_lo:vol_hi)   : (1:0)

    rank == 0 && println(
        "SCFIE-PWCHex MPI Assembly (column partition): n_surf=$n_surf, n_vol=$n_vol, n_total=$n_total, procs=$n_procs, threads=$(Threads.nthreads())",
    )
    MPI.Barrier(comm)
    println("  Rank $rank: surf_cols=$(length(local_surf_cols)), vol_cols=$(length(local_vol_cols))")
    MPI.Barrier(comm)
    flush(stdout)

    tris  = get_triangles_info(surf_basis.mesh, surf_basis)
    ntri  = length(tris)
    hexas = get_hexahedra_info(vol_basis.mesh, vol_basis, scfie.permittivities)
    nhex  = length(hexas)

    cfie = CFIE(scfie.freq, scfie.alpha)

    # Physical constants
    k         = scfie.k
    k²        = k^2
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

    row_locks = [SpinLock() for _ = 1:n_total]
    n_threads = Threads.nthreads()

    # ─── Phase 1: Z_SS + Z_VS (surf cols) ────────────────────────────────────
    rank == 0 && println("  [SCFIE-PWCHex 1] Z_SS + Z_VS (surf cols)...")
    rank == 0 && flush(stdout)

    next_ss = Threads.Atomic{Int}(1)
    Threads.@threads for _ = 1:n_threads
        l_efie = zeros(CT, 3, 3)
        l_mfie = zeros(CT, 3, 3)
        l_loc  = zeros(CT, 3, 3)
        dyadG  = zeros(CT, 3, 3)

        while true
            tidx = Threads.atomic_add!(next_ss, 1)
            tidx > length(src_tris_surf) && break

            t_src = src_tris_surf[tidx]
            tri_s = tris[t_src]

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
                    col == 0 && continue
                    col in local_surf_cols || continue
                    bf_s   = surf_basis.functions[col]
                    sign_s = (bf_s.support[1] == t_src) ? bf_s.signs[1] : bf_s.signs[2]
                    @inbounds for i = 1:3
                        row = tri_t.inBfsID[i]
                        row == 0 && continue
                        bf_t   = surf_basis.functions[row]
                        sign_t = (bf_t.support[1] == t_tst) ? bf_t.signs[1] : bf_t.signs[2]
                        lock(row_locks[row])
                        Z[row, col] += l_loc[i, j] * sign_t * sign_s
                        unlock(row_locks[row])
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
                        R̂x = Rx * divR; R̂y = Ry * divR; R̂z = Rz * divR
                        GR   = exp(-jk * R) * div4π * divR * gq_hex_near.weight[gj]

                        RR11 = R̂x*R̂x; RR12 = R̂x*R̂y; RR13 = R̂x*R̂z
                        RR22 = R̂y*R̂y; RR23 = R̂y*R̂z; RR33 = R̂z*R̂z

                        dyadG[1,1] += GR * ((1 - RR11) * k² - (1 - 3RR11) * jkpR)
                        dyadG[2,2] += GR * ((1 - RR22) * k² - (1 - 3RR22) * jkpR)
                        dyadG[3,3] += GR * ((1 - RR33) * k² - (1 - 3RR33) * jkpR)
                        od12 = GR * (-RR12 * k² + 3RR12 * jkpR)
                        od13 = GR * (-RR13 * k² + 3RR13 * jkpR)
                        od23 = GR * (-RR23 * k² + 3RR23 * jkpR)
                        dyadG[1,2] += od12; dyadG[2,1] += od12
                        dyadG[1,3] += od13; dyadG[3,1] += od13
                        dyadG[2,3] += od23; dyadG[3,2] += od23
                    end

                    for mi = 1:3
                        m = tri_s.inBfsID[mi]
                        m == 0 && continue
                        m in local_surf_cols || continue
                        lm    = tri_s.edgel[mi]
                        freeV = tri_s.vertices[:, mi]
                        ρmi   = SVector(rgi[1]-freeV[1], rgi[2]-freeV[2], rgi[3]-freeV[3])
                        temp  = gq_s.weight[gi] * lm / 2

                        for ni = 1:3
                            n     = hex.inBfsID[ni]
                            dot_vs = ρmi[1]*dyadG[ni,1] + ρmi[2]*dyadG[ni,2] + ρmi[3]*dyadG[ni,3]
                            z_vs   = temp * dot_vs * Jη₀divK * Vs
                            row_vs = n_surf + n
                            lock(row_locks[row_vs])
                            Z[row_vs, m] += z_vs
                            unlock(row_locks[row_vs])
                        end
                    end
                end
            end
        end
    end

    # ─── Phase 2: Z_VV + Z_SV (vol cols) ─────────────────────────────────────
    rank == 0 && println("  [SCFIE-PWCHex 2] Z_VV + Z_SV (vol cols)...")
    rank == 0 && flush(stdout)

    next_vv = Threads.Atomic{Int}(1)
    Threads.@threads for _ = 1:n_threads
        Z_ts_buf = zeros(CT, 3, 3)
        dyadG    = zeros(CT, 3, 3)

        while true
            tidx = Threads.atomic_add!(next_vv, 1)
            tidx > length(src_hexas_vol) && break

            js    = src_hexas_vol[tidx]
            hex_s = hexas[js]
            κₛ    = hex_s.κ
            rad_s = cbrt(hex_s.volume)

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
                    k, k², jk, Jη₀divK, div4π,
                )

                @inbounds for ni = 1:3
                    n   = hex_s.inBfsID[ni]
                    col = n_surf + n
                    col in local_vol_cols || continue
                    lock(row_locks[col])
                    for mi = 1:3
                        m   = hex_t.inBfsID[mi]
                        row = n_surf + m
                        if ti == js
                            val = Z_ts_buf[mi, ni] * κₛ
                            if mi == ni
                                selfImp = CT(1) / (im * omega) / hex_s.ε * hex_s.volume
                                val += selfImp
                            end
                            Z[row, col] += val
                        else
                            Z[row, col] += Z_ts_buf[mi, ni] * κₛ
                        end
                    end
                    unlock(row_locks[col])
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
                        R̂x = Rx * divR; R̂y = Ry * divR; R̂z = Rz * divR
                        GR   = exp(-jk * R) * div4π * divR * gq_hex_near.weight[gj]

                        RR11 = R̂x*R̂x; RR12 = R̂x*R̂y; RR13 = R̂x*R̂z
                        RR22 = R̂y*R̂y; RR23 = R̂y*R̂z; RR33 = R̂z*R̂z

                        dyadG[1,1] += GR * ((1 - RR11) * k² - (1 - 3RR11) * jkpR)
                        dyadG[2,2] += GR * ((1 - RR22) * k² - (1 - 3RR22) * jkpR)
                        dyadG[3,3] += GR * ((1 - RR33) * k² - (1 - 3RR33) * jkpR)
                        od12 = GR * (-RR12 * k² + 3RR12 * jkpR)
                        od13 = GR * (-RR13 * k² + 3RR13 * jkpR)
                        od23 = GR * (-RR23 * k² + 3RR23 * jkpR)
                        dyadG[1,2] += od12; dyadG[2,1] += od12
                        dyadG[1,3] += od13; dyadG[3,1] += od13
                        dyadG[2,3] += od23; dyadG[3,2] += od23
                    end

                    for mi = 1:3
                        m = tri.inBfsID[mi]
                        m == 0 && continue
                        lm    = tri.edgel[mi]
                        freeV = tri.vertices[:, mi]
                        ρmi   = SVector(rgi[1]-freeV[1], rgi[2]-freeV[2], rgi[3]-freeV[3])
                        temp  = gq_s.weight[gi] * lm / 2

                        for ni = 1:3
                            n   = hex_s.inBfsID[ni]
                            col = n_surf + n
                            col in local_vol_cols || continue
                            dot_sv = ρmi[1]*dyadG[1,ni] + ρmi[2]*dyadG[2,ni] + ρmi[3]*dyadG[3,ni]
                            z_sv   = temp * dot_sv * Jη₀divK * hex_s.volume * κₛ
                            lock(row_locks[m])
                            Z[m, col] += z_sv
                            unlock(row_locks[m])
                        end
                    end
                end
            end
        end
    end

    sync!(Z)

    rank == 0 && println("SCFIE-PWCHex MPI Assembly Completed.")
    rank == 0 && flush(stdout)

    return Z
end
