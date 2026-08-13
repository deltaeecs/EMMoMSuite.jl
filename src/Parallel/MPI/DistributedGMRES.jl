# DistributedGMRES.jl
# Phase 15.4 Method C: 分布式 Krylov GMRES，消除 I-3（Krylov 基向量全量复制）
#
# 核心思路：
#   - 保持 MPIMatrix 列分区不变（装配 API 不破坏）
#   - Krylov 基向量按行分区：每进程仅存 N/P 行 → O(N/P × restart × CT) 内存
#   - matvec: Allreduce-gather(N) + 现有列分区 mul!(内含 Allreduce) → 全量 y，取本地行
#   - dot:    分布式内积 → 1 次 Allreduce(scalar)
#   - axpy:   纯本地操作
#   - 收敛后 Allreduce-gather 全量 x，与当前 mpi_gmres! API 完全兼容

using LinearAlgebra, MPI
using Logging

# ── 行分区工具 ──────────────────────────────────────────────────────────────────

"""
    row_partition(N, rank, nproc) → UnitRange{Int}

计算 `rank` 进程在大小为 `N` 的全局向量上所拥有的行范围（1-indexed）。
负载均衡：各进程大小相差最多 1。
"""
function row_partition(N::Int, rank::Int, nproc::Int)
    base  = N ÷ nproc
    extra = N % nproc
    if rank < extra
        lo = rank * (base + 1) + 1
        hi = lo + base
    else
        lo = rank * base + extra + 1
        hi = lo + base - 1
    end
    return lo:hi
end

# ── 分布式基本运算 ──────────────────────────────────────────────────────────────

"""
    _dist_dot(u_local, v_local, comm) → CT

分布式内积：各进程贡献本地行的点积，Allreduce 求和。
"""
@inline function _dist_dot(u_local::AbstractVector, v_local::AbstractVector, comm::MPI.Comm)
    return MPI.Allreduce(dot(u_local, v_local), +, comm)
end

"""
    _dist_norm2(v_local, comm) → Real

分布式 ‖v‖² （实数）。
"""
@inline function _dist_norm2(v_local::AbstractVector, comm::MPI.Comm)
    return real(MPI.Allreduce(dot(v_local, v_local), +, comm))
end

"""
    _allreduce_gather!(buf, v_local, local_rows, comm)

通过 Allreduce(+) 将分布式向量拼接为全量向量 `buf`（长度 N）。
每个进程仅在 `local_rows` 处写入值，其余置零，Allreduce 求和即得完整向量。
"""
function _allreduce_gather!(buf::AbstractVector{CT},
                             v_local::AbstractVector{CT},
                             local_rows::UnitRange{Int},
                             comm::MPI.Comm) where {CT}
    fill!(buf, zero(CT))
    buf[local_rows] .= v_local
    MPI.Allreduce!(buf, +, comm)
    return buf
end

# ── Givens 旋转（支持复数） ────────────────────────────────────────────────────

"""
    _givens(a, b) → (c::Real, s::CT, r::CT)

计算 Givens 旋转使得：`[c s; -conj(s) c] * [a; b] = [r; 0]`
- c 为实余弦，s 为复正弦，r 为旋转后的对角元素
"""
@inline function _givens(a::CT, b::CT) where {CT}
    FT   = real(CT)
    absa = abs(a)
    absb = abs(b)
    n    = sqrt(absa^2 + absb^2)
    if n < eps(FT)
        return one(FT), zero(CT), a
    end
    c = FT(absa / n)
    if iszero(absa)
        # a = 0 的退化情形：旋转 90° 消去 b
        return zero(FT), one(CT), b
    end
    # s = sign(a) * conj(b) / n，r = n * sign(a)
    sa = CT(a / absa)          # sign(a)
    s  = CT(sa * conj(b) / n)
    r  = CT(n) * sa
    return c, s, r
end

"""
    _apply_givens(c, s, x, y) → (x', y')

将已计算的 Givens 旋转 (c, s) 作用于向量 [x; y]。
"""
@inline function _apply_givens(c::FT, s::CT, x::CT, y::CT) where {FT<:AbstractFloat, CT}
    return c * x + s * y, -conj(s) * x + c * y
end

# ── 上三角回代 ─────────────────────────────────────────────────────────────────

"""
    _upper_triangular_solve(H, g, j) → Vector

求解上三角系统 `H[1:j, 1:j] * y = g[1:j]`（原地回代）。
"""
function _upper_triangular_solve(H::Matrix{CT}, g::Vector{CT}, j::Int) where {CT}
    y = copy(g[1:j])
    for i in j:-1:1
        y[i] /= H[i, i]
        for k in 1:(i-1)
            y[k] -= H[k, i] * y[i]
        end
    end
    return y
end

# ── 核心实现 ───────────────────────────────────────────────────────────────────

"""
    distributed_gmres!(x, A::MPIMatrix, b; restart, maxiter, reltol, abstol,
                        verbose, log, initially_zero) → x [, history]

**分布式 Krylov GMRES**（Method C，消除 I-3）

## 与 mpi_gmres! 的区别
| 方面 | 旧（SPMD IterativeSolvers）| 新（DistributedGMRES）|
|------|---------------------------|----------------------|
| Krylov 基内存 | N × restart × CT × P（全量复制）| N × restart × CT（分布存储）|
| dot 操作 | 纯本地（O(N/P)，无通信）| 1 次 Allreduce(scalar)/步 |
| 额外通信 | 0 | 1 次 Allreduce(N) 每次 matvec（gather 基向量）|

## 通信分析（每次 Arnoldi 步）
- 1 × Allreduce(N CT)：gather 基向量 v_j
- 1 × Allreduce(N CT)：列分区 mul! 内部（与旧方案相同）
- j × Allreduce(scalar)：Modified Gram-Schmidt 内积（j ≤ restart）

## 内存（每进程）
- Krylov 基：`(restart+1) × n_local × CT`，`n_local ≈ N/P`
- 全量缓冲：2 × N × CT（gather 和 matvec 缓冲，不随 restart 增长）

## 参数
- `x`: 初始猜测（原地更新）；普通 `Vector{CT}`，长度 N
- `A`: 列分区 `MPIMatrix`
- `b`: 右端向量；普通 `Vector{CT}`，长度 N
- `restart`: Krylov 子空间大小
- `maxiter`: 最大迭代次数
- `reltol`: 相对收敛容差（基于 ‖b‖）
- `abstol`: 绝对收敛容差
- `verbose`: 打印收敛信息（仅 rank 0）
- `log`: 返回 `(x, history)` 而非仅 `x`
- `initially_zero`: 若为 `true`，跳过初始 r = b - A*x（节省 1 次 matvec）
- `Pl`: 左预条件（`nothing` = 不预条件；支持 DistributedBlockJacobiPreconditioner /
  DistributedDiagonalPreconditioner 及串行预条件回退）
"""
function distributed_gmres!(
    x           :: AbstractVector{CT},
    A           :: MPIMatrix,
    b           :: AbstractVector{CT};
    restart     :: Int   = min(30, size(A, 2)),
    maxiter     :: Int   = size(A, 2),
    reltol      :: Real  = 1e-6,
    abstol      :: Real  = 0.0,
    verbose     :: Bool  = false,
    log         :: Bool  = false,
    initially_zero :: Bool = false,
    Pl          :: Any   = nothing,
    kwargs...,       # 兼容旧 API 的额外 kwargs（忽略）
) where {CT}
    FT    = real(CT)
    comm  = A.comm
    rank  = MPI.Comm_rank(comm)
    nproc = MPI.Comm_size(comm)
    N     = length(b)

    restart = min(restart, maxiter, N)

    # ── 行分区 ────────────────────────────────────────────────────────────────
    local_rows = row_partition(N, rank, nproc)
    n_local    = length(local_rows)

    if rank == 0 && verbose
        @info "[DistGMRES] N=$N, restart=$restart, maxiter=$maxiter, reltol=$reltol, P=$nproc, n_local=$n_local"
    end

    # ── 预分配缓冲区 ──────────────────────────────────────────────────────────
    # 全量缓冲（不随 restart 增长，固定 2×N）
    v_gather  = zeros(CT, N)   # Allreduce-gather 缓冲（覆盖写）
    y_matvec  = zeros(CT, N)   # mul! 输出缓冲（覆盖写）
    # 预条件缓冲（仅 Pl 非空时分配）：z_pre 为预条件后的 Arnoldi 输入，r_full 为残差
    z_pre  = Pl === nothing ? y_matvec : zeros(CT, N)
    r_full = Pl === nothing ? y_matvec : zeros(CT, N)

    # 分布式 Krylov 基：各列长度 n_local ≈ N/P
    # 内存：(restart+1) * n_local * sizeof(CT) per process
    V = [zeros(CT, n_local) for _ in 1:(restart + 1)]

    # Hessenberg 矩阵（小矩阵，全量复制，很小）
    H  = zeros(CT, restart + 1, restart)

    # Givens 旋转系数
    cs = zeros(FT, restart)
    sn = zeros(CT, restart)

    # 最小二乘右端
    g  = zeros(CT, restart + 1)

    # 收敛历史
    residuals = FT[]
    sizehint!(residuals, min(maxiter + 3, 500))

    # ── 计算初始残差 r₀ = M⁻¹(b - A*x)（左预条件）────────────────────────────
    if initially_zero
        # x = 0, r = M⁻¹ b
        if Pl === nothing
            r_local = copy(b[local_rows])
        else
            apply_mpi_preconditioner!(r_full, Pl, b)
            r_local = r_full[local_rows]
        end
    else
        mul!(y_matvec, A, x)
        @. y_matvec = b - y_matvec
        if Pl === nothing
            r_local = y_matvec[local_rows]   # copy（UnitRange 索引）
        else
            apply_mpi_preconditioner!(r_full, Pl, y_matvec)
            r_local = r_full[local_rows]
        end
    end

    b_norm = Pl === nothing ? sqrt(_dist_norm2(b[local_rows], comm)) :
             sqrt(_dist_norm2(r_local, comm))
    r_norm = sqrt(_dist_norm2(r_local, comm))
    tol    = max(FT(reltol) * b_norm, FT(abstol))

    push!(residuals, r_norm)
    if rank == 0 && verbose
        @info "[DistGMRES]  iter    0  resnorm = $(r_norm)  tol = $(tol)"
    end

    # 初始已收敛
    if r_norm ≤ tol
        history = (isconverged = true, mvps = 0,
                   resnorm = r_norm, residuals = residuals)
        return log ? (x, history) : x
    end

    iter_total = 0
    converged  = false

    # ── 外层 restart 循环 ────────────────────────────────────────────────────
    while !converged && iter_total < maxiter

        # v₁ = r / ‖r‖（分布式归一化，仅本地操作）
        @. V[1] = r_local / r_norm

        # 重置 Hessenberg 和最小二乘右端
        fill!(H, zero(CT))
        fill!(cs, zero(FT))
        fill!(sn, zero(CT))
        fill!(g,  zero(CT))
        g[1] = CT(r_norm)

        j = 0   # 本次 restart 内迭代步数

        # ── 内层 Arnoldi 循环 ───────────────────────────────────────────────
        while !converged && j < restart && iter_total < maxiter
            j          += 1
            iter_total += 1

            # ─ Arnoldi: w = M⁻¹ * (A * V[j])（左预条件）──────────────────
            # Step 1: gather 分布式 V[j] → 全量 v_gather（Allreduce 方式）
            _allreduce_gather!(v_gather, V[j], local_rows, comm)

            # Step 2: 列分区 matvec（内含 Allreduce，结果全量复制到所有进程）
            mul!(y_matvec, A, v_gather)

            # Step 3: 左预条件 M⁻¹（分布式施加；无预条件时为恒等），取本地行
            if Pl === nothing
                w_local = view(y_matvec, local_rows)
            else
                apply_mpi_preconditioner!(z_pre, Pl, y_matvec)
                w_local = view(z_pre, local_rows)
            end

            # ─ Modified Gram-Schmidt 正交化 ────────────────────────────────
            for i in 1:j
                H[i, j] = _dist_dot(V[i], w_local, comm)
                @. w_local = w_local - H[i, j] * V[i]
            end

            # ‖w‖（分布式范数）
            H[j+1, j] = CT(sqrt(_dist_norm2(w_local, comm)))

            # 新基向量 V[j+1] = w / ‖w‖
            if abs(H[j+1, j]) > eps(FT)
                @. V[j+1] = w_local / H[j+1, j]
            else
                fill!(V[j+1], zero(CT))   # 退化情形（极罕见）
            end

            # ─ 应用之前的 Givens 旋转 ────────────────────────────────────
            for i in 1:(j - 1)
                H[i, j], H[i+1, j] = _apply_givens(cs[i], sn[i], H[i, j], H[i+1, j])
            end

            # ─ 计算本步 Givens 旋转并更新 ────────────────────────────────
            cs[j], sn[j], H[j, j] = _givens(H[j, j], H[j+1, j])
            H[j+1, j]  = zero(CT)
            g[j+1]     = -conj(sn[j]) * g[j]   # G = [c, s; -conj(s), c]，复数 Givens 旋转
            g[j]       =  cs[j] * g[j]           # cs[j] 为实数
            r_norm      = abs(g[j+1])

            push!(residuals, r_norm)
            if rank == 0 && verbose
                @info "[DistGMRES]  iter $(" "^(5 - ndigits(iter_total)))$(iter_total)  resnorm = $(r_norm)"
            end

            r_norm ≤ tol && (converged = true)
        end  # 内层 Arnoldi

        # ── 更新 x：x += V * y（上三角回代 + gather）─────────────────────
        y = _upper_triangular_solve(H, g, j)

        # dx_full[local_rows] = sum_i y[i] * V[i]，其余为 0，Allreduce 得全量
        fill!(v_gather, zero(CT))
        for i in 1:j
            v_gather[local_rows] .+= y[i] .* V[i]
        end
        MPI.Allreduce!(v_gather, +, comm)
        x .+= v_gather

        # ── restart：重新计算残差 ─────────────────────────────────────────
        if !converged && iter_total < maxiter
            mul!(y_matvec, A, x)
            @. y_matvec = b - y_matvec
            if Pl === nothing
                r_local = y_matvec[local_rows]
            else
                apply_mpi_preconditioner!(r_full, Pl, y_matvec)
                r_local = r_full[local_rows]
            end
            r_norm  = sqrt(_dist_norm2(r_local, comm))
            push!(residuals, r_norm)
            if rank == 0 && verbose
                @info "[DistGMRES] restart  resnorm = $(r_norm)"
            end
            r_norm ≤ tol && (converged = true)
        end

    end  # 外层 restart

    history = (isconverged = converged, mvps = iter_total,
               resnorm = r_norm, residuals = residuals)
    if rank == 0 && verbose
        @info "[DistGMRES] done: converged=$(converged), iters=$(iter_total), final_resnorm=$(r_norm)"
    end
    return log ? (x, history) : x
end
