# Solver.jl
# Phase 14.3: 分布式 GMRES via SPMD + MPIMatrix 列分区 mul!
#
# 策略 (SPMD — Single Program Multiple Data):
#   - 所有 MPI 进程同步执行相同的 GMRES 迭代树
#   - Krylov 向量 (x, b, Arnoldi 基向量) 在每个进程上完全复制 (长度 N)
#   - mul!(y, A::MPIMatrix, x) 内部用 Allreduce 使所有进程得到相同的 y
#   - dot / norm 操作对本地 Vector 直接计算 (各进程结果相同，无需 MPI 通信)
#   - 内存节省来自矩阵 A：每个进程只存 N×(N/P)，而不是 N×N
#
# 无需修改 IterativeSolvers.jl 内部; 只需正确的 LinearAlgebra dispatch 即可。

import IterativeSolvers

"""
    mpi_gmres!(x, A::MPIMatrix, b; restart, maxiter, reltol, abstol, verbose, log)

SPMD 分布式 GMRES 求解 `A*x = b`，`A` 为列分区 `MPIMatrix`。

**工作原理**:
- 所有 MPI 进程同步运行相同的 GMRES 迭代树
- `mul!(y, A, x)` 对列分区矩阵做局部乘积后 `MPI.Allreduce!`，所有进程获得相同 `y`
- Krylov 向量在所有进程上完全复制；输入 `x`, `b` 请为普通 `Vector`
- 收敛后所有进程持有相同的解向量 `x`

**参数**:
- `x`: 初始猜测 (将被原地更新); 普通 `Vector{CT}`
- `A`: 列分区 `MPIMatrix`
- `b`: 右端向量; 普通 `Vector{CT}`
- `restart`: Krylov 子空间大小 (默认 30)
- `maxiter`: 最大迭代次数 (默认 size(A,2))
- `reltol`: 相对收敛容差 (默认 1e-6)
- `abstol`: 绝对收敛容差 (默认 0)
- `verbose`: 是否在 rank 0 打印迭代信息 (默认 false)
- `log`: 是否返回收敛历史 (默认 false)

**返回**: `log=false` 时返回 `x`; `log=true` 时返回 `(x, history)`
"""
function mpi_gmres!(
    x::AbstractVector,
    A::MPIMatrix,
    b::AbstractVector;
    restart::Int   = min(30, size(A, 2)),
    maxiter::Int   = size(A, 2),
    reltol::Real   = 1e-6,
    abstol::Real   = 0,
    verbose::Bool  = false,
    log::Bool      = false,
    kwargs...,
)
    comm = A.comm
    rank = MPI.Comm_rank(comm)
    np   = MPI.Comm_size(comm)

    if rank == 0 && verbose
        N = size(A, 1)
        println("MPI GMRES: N=$N, restart=$restart, maxiter=$maxiter, reltol=$reltol, procs=$np, threads=$(Threads.nthreads())")
        flush(stdout)
    end

    # 所有进程同步执行 GMRES;  mul!(Ax, A, x) 内部走 Allreduce 分支
    result = IterativeSolvers.gmres!(
        x, A, b;
        restart        = restart,
        maxiter        = maxiter,
        reltol         = reltol,
        abstol         = abstol,
        log            = log,
        verbose        = false,   # 关闭 IterativeSolvers 内置日志（所有 rank 都会打印）
        kwargs...,
    )

    # 手动打印收敛结果（仅 rank 0）
    if rank == 0 && verbose
        if log
            history = result[2]
            println("GMRES converged: $(history.isconverged), iters=$(history.mvps), final_resnorm=$(last(history[:resnorm]))")
        end
        flush(stdout)
    end

    return result
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
