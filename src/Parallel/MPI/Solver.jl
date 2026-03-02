# Solver.jl
# Phase 15.4 Method C: 分布式 Krylov GMRES（消除 I-3）
#
# 策略 (Distributed Krylov):
#   - Krylov 基向量按行分区：每进程仅存 N/P 行（见 DistributedGMRES.jl）
#   - mul!(y, A::MPIMatrix, x) 使用现有列分区实现（内含 Allreduce）
#   - dot/norm 通过 MPI.Allreduce(scalar) 计算（每步 j 次）
#   - 内存：O(N/P × restart × CT) per process（vs 旧方案 O(N × restart × CT × P)）
#   - 公共 API (mpi_gmres! / mpi_gmres) 签名不变

"""
    mpi_gmres!(x, A::MPIMatrix, b; restart, maxiter, reltol, abstol, verbose, log)

分布式 Krylov GMRES 求解 `A*x = b`，`A` 为列分区 `MPIMatrix`。

**工作原理** (Phase 15.4 Method C):
- Krylov 基向量 **行分区存储**：每进程持有 `N/P` 行，消除 I-3 全量复制内存
- `mul!(y, A, x)` 对列分区矩阵做局部乘积后 `MPI.Allreduce!`（与旧版相同）
- Modified Gram-Schmidt 内积通过 `MPI.Allreduce(scalar)` 完成
- 收敛后所有进程持有相同的完整解向量 `x`

**参数**:
- `x`: 初始猜测 (原地更新); 普通 `Vector{CT}`，长度 N
- `A`: 列分区 `MPIMatrix`
- `b`: 右端向量; 普通 `Vector{CT}`，长度 N
- `restart`: Krylov 子空间大小 (默认 min(30, N))
- `maxiter`: 最大迭代次数 (默认 N)
- `reltol`: 相对收敛容差 (默认 1e-6)
- `abstol`: 绝对收敛容差 (默认 0)
- `verbose`: rank 0 打印迭代信息 (默认 false)
- `log`: 返回 `(x, history)` 而非仅 `x` (默认 false)

**返回**: `log=false` 时返回 `x`; `log=true` 时返回 `(x, history)`
"""
function mpi_gmres!(
    x       :: AbstractVector,
    A       :: MPIMatrix,
    b       :: AbstractVector;
    restart :: Int  = min(30, size(A, 2)),
    maxiter :: Int  = size(A, 2),
    reltol  :: Real = 1e-6,
    abstol  :: Real = 0,
    verbose :: Bool = false,
    log     :: Bool = false,
    kwargs...,
)
    return distributed_gmres!(
        x, A, b;
        restart    = restart,
        maxiter    = maxiter,
        reltol     = reltol,
        abstol     = abstol,
        verbose    = verbose,
        log        = log,
        kwargs...,
    )
end

"""
    mpi_gmres(A::MPIMatrix, b; kwargs...)

返回 `x = A \\ b` 的 GMRES 近似解（从零初值出发）。
"""
function mpi_gmres(
    A::MPIMatrix,
    b::AbstractVector;
    kwargs...,
)
    x = zeros(eltype(b), length(b))
    return mpi_gmres!(x, A, b; initially_zero = true, kwargs...)
end
