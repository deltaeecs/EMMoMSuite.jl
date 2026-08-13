# ScaLAPACKLU.jl — 通过本机 MinGW 的 ScaLAPACK（MSMPI 版）实现分布式稠密 LU
#
# 库：C:\\msys64\\mingw64\\bin\\libscalapack.dll（ScaLAPACK 2.2.2，链接 msmpi.dll）
# 接口：BLACS（blacs_gridinit_ 等）+ pzgesv_（复数部分主元 LU + 求解），2D 块循环布局。
# 可用环境变量 SCALAPACK_LIB_PATH 覆盖库路径。

using MPI

const SCALAPACK_LIB = get(
    ENV,
    "SCALAPACK_LIB_PATH",
    "C:\\msys64\\mingw64\\bin\\libscalapack.dll",
)

"""
    ScaLAPACKGrid

BLACS 进程网格（基于 MPI world 创建；`P*Q ≤ MPI 世界大小`，取最接近 sqrt 的因子）。
"""
struct ScaLAPACKGrid
    ictxt::Int32
    nprow::Int32
    npcol::Int32
    myrow::Int32
    mycol::Int32
    MB::Int32
    NB::Int32
    N::Int32
end

function _blacs_gridinit!(ictxt::Ref{Int32}, nprow::Int32, npcol::Int32)
    layout = Ref{UInt8}(UInt8('R'))
    ccall(
        (:blacs_gridinit_, SCALAPACK_LIB),
        Cvoid,
        (Ref{Int32}, Ref{UInt8}, Ref{Int32}, Ref{Int32}),
        ictxt, layout, Ref(nprow), Ref(npcol),
    )
    return nothing
end

"""
    _blacs_system_context() → Ref{Int32}

返回 BLACS 系统上下文（`BLACS_GET(-1, 0)`）。必须先取系统上下文再 `gridinit`，
否则 MSMPI 版 BLACS 无法正确登记进程网格。
"""
function _blacs_system_context()
    ictxt = Ref{Int32}(0)
    ccall(
        (:blacs_get_, SCALAPACK_LIB),
        Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Int32}),
        Ref{Int32}(-1), Ref{Int32}(0), ictxt,
    )
    return ictxt
end

function _blacs_gridinfo!(ictxt::Int32)
    nprow = Ref{Int32}(0); npcol = Ref{Int32}(0)
    myrow = Ref{Int32}(0); mycol = Ref{Int32}(0)
    ccall(
        (:blacs_gridinfo_, SCALAPACK_LIB),
        Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}),
        Ref(ictxt), nprow, npcol, myrow, mycol,
    )
    return nprow[], npcol[], myrow[], mycol[]
end

_blacs_gridexit!(ictxt::Int32) = ccall(
    (:blacs_gridexit_, SCALAPACK_LIB), Cvoid, (Ref{Int32},), Ref(ictxt)
)

"""
    init_grid(comm; MB=64, NB=64) → ScaLAPACKGrid

基于 MPI world 创建 BLACS 网格（行主序），块大小 MB×NB。
"""
function init_grid(comm; MB::Int = 64, NB::Int = 64)
    P = MPI.Comm_size(comm)
    nprow = Int32(1)
    best = typemax(Int)
    for r in 1:P
        P % r == 0 || continue
        c = P ÷ r
        d = abs(r - c)
        if d < best
            best = d
            nprow = Int32(r)
        end
    end
    npcol = Int32(P ÷ nprow)
    ictxt = _blacs_system_context()
    _blacs_gridinit!(ictxt, nprow, npcol)
    _, _, myrow, mycol = _blacs_gridinfo!(ictxt[])
    return ScaLAPACKGrid(ictxt[], nprow, npcol, myrow, mycol, Int32(MB), Int32(NB), Int32(0))
end

"""
    numroc(n, nb, iproc, isrcproc, nprocs) → Int

ScaLAPACK numroc 语义：进程 `iproc`（源进程 `isrcproc`，取 0）在 n 个元素、
块大小 nb、共 nprocs 个进程下的本地元素数。
"""
function numroc(n::Int, nb::Int, iproc::Int, isrcproc::Int, nprocs::Int)
    nblocks = cld(n, nb)
    count = 0
    for b in 1:nblocks
        ((b - 1 + isrcproc) % nprocs == iproc) || continue
        start = (b - 1) * nb + 1
        count += min(nb, n - start + 1)
    end
    return count
end

"""
    local_layout(N, MB, NB, grid) → (nlocal_rows, nlocal_cols, lld)

块循环布局下本进程的本地行数、列数与前导维度。
"""
function local_layout(N::Int, MB::Int, NB::Int, grid::ScaLAPACKGrid)
    nr = numroc(N, MB, Int(grid.myrow), 0, Int(grid.nprow))
    nc = numroc(N, NB, Int(grid.mycol), 0, Int(grid.npcol))
    lld = max(1, nr)
    return nr, nc, lld
end

"""
    distribute(A, grid; MB, NB) → A_local

把全量稠密矩阵 `A`（M×K，各秩持有完整副本）按块循环布局分发到本进程本地块
（RHS 传 `reshape(b, N, 1)`）。
"""
function distribute(A::Matrix{CT}, grid::ScaLAPACKGrid; MB::Int = 64, NB::Int = 64) where {CT<:Complex}
    M, K = size(A)
    myrow = Int(grid.myrow)
    mycol = Int(grid.mycol)
    nprow = Int(grid.nprow)
    npcol = Int(grid.npcol)
    nr = numroc(M, MB, myrow, 0, nprow)
    nc = numroc(K, NB, mycol, 0, npcol)
    Aloc = zeros(CT, max(1, nr), max(1, nc))
    lr0 = 1
    for bi in 1:cld(M, MB)
        (bi - 1) % nprow == myrow || continue
        i0 = (bi - 1) * MB + 1
        nrows = min(MB, M - i0 + 1)
        lc0 = 1
        for bj in 1:cld(K, NB)
            (bj - 1) % npcol == mycol || continue
            j0 = (bj - 1) * NB + 1
            ncols = min(NB, K - j0 + 1)
            Aloc[lr0:lr0 + nrows - 1, lc0:lc0 + ncols - 1] .=
                @view A[i0:i0 + nrows - 1, j0:j0 + ncols - 1]
            lc0 += ncols
        end
        lr0 += nrows
    end
    return Aloc
end

"""
    pzgesv!(Aloc, bloc, grid; MB, NB) → ipiv

ScaLAPACK `pzgesv`：复数分布式 LU + 求解（原地）。`bloc` 为块循环分布的右端/解。
"""
function pzgesv!(Aloc::Matrix{CT}, bloc::Vector{CT}, grid::ScaLAPACKGrid; MB::Int = 64, NB::Int = 64) where {CT<:Complex}
    N = Int(grid.N)
    nrhs = Int32(1)
    ia = Int32(1); ja = Int32(1); ib = Int32(1); jb = Int32(1)
    nr, nc, lld = local_layout(N, MB, NB, grid)
    desc = zeros(Int32, 9)
    descb = zeros(Int32, 9)
    info = Ref{Int32}(0)
    ccall(
        (:descinit_, SCALAPACK_LIB),
        Cvoid,
        (Ptr{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}),
        desc, Ref(Int32(N)), Ref(Int32(N)), Ref(Int32(MB)), Ref(Int32(NB)),
        Ref(Int32(0)), Ref(Int32(0)), Ref(grid.ictxt), Ref(Int32(lld)), info,
    )
    info[] == 0 || error("descinit failed: info=", info[])
    ccall(
        (:descinit_, SCALAPACK_LIB),
        Cvoid,
        (Ptr{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}),
        descb, Ref(Int32(N)), Ref(nrhs), Ref(Int32(MB)), Ref(Int32(NB)),
        Ref(Int32(0)), Ref(Int32(0)), Ref(grid.ictxt), Ref(Int32(length(bloc))), info,
    )
    info[] == 0 || error("descinit(b) failed: info=", info[])
    ipiv = zeros(Int32, N)
    ccall(
        (:pzgesv_, SCALAPACK_LIB),
        Cvoid,
        (Ref{Int32}, Ref{Int32}, Ptr{CT}, Ref{Int32}, Ref{Int32}, Ptr{Int32}, Ptr{Int32},
         Ptr{CT}, Ref{Int32}, Ref{Int32}, Ptr{Int32}, Ptr{Int32}),
        Ref(Int32(N)), Ref(nrhs), Aloc, ia, ja, desc, ipiv, bloc, ib, jb, descb, info,
    )
    info[] == 0 || error("pzgesv failed: info=", info[])
    return ipiv
end

"""
    gather_solution(bloc, grid; MB, NB, N) → x

从块循环分布的 `bloc`（N×1 右端/解）恢复全量解向量 `x`（各秩一致）。
"""
function gather_solution(bloc::Vector{CT}, grid::ScaLAPACKGrid; MB::Int = 64, NB::Int = 64) where {CT<:Complex}
    N = Int(grid.N)
    x = zeros(CT, N)
    mycol = Int(grid.mycol)
    myrow = Int(grid.myrow)
    nprow = Int(grid.nprow)
    if mycol == 0   # RHS 的单列块由进程列 0 持有
        for j in 1:N
            bi = (j - 1) ÷ MB + 1
            ((bi - 1) % nprow == myrow) || continue
            lr = ((bi - 1) ÷ nprow) * MB + (j - 1) % MB + 1
            x[j] = bloc[lr]
        end
    end
    MPI.Allreduce!(x, +, MPI.COMM_WORLD)
    return x
end

"""
    scalapack_lu_solve(A, b, comm; MB=64, NB=64) → x

高层接口：各秩持有完整 `A`/`b`，内部做块循环分发 → `pzgesv` → 恢复全量解。
"""
function scalapack_lu_solve(
    A::Matrix{CT},
    b::Vector{CT},
    comm;
    MB::Int = 64,
    NB::Int = 64,
) where {CT<:Complex}
    MB == NB || error("ScaLAPACK PZGESV requires square block decomposition (MB == NB); got MB=$MB, NB=$NB")
    N = size(A, 1)
    grid = init_grid(comm; MB = MB, NB = NB)
    grid = ScaLAPACKGrid(grid.ictxt, grid.nprow, grid.npcol, grid.myrow, grid.mycol, Int32(MB), Int32(NB), Int32(N))
    Aloc = distribute(A, grid; MB = MB, NB = NB)
    bloc_full = distribute(reshape(b, N, 1), grid; MB = MB, NB = NB)
    bloc = size(bloc_full, 2) >= 1 ? bloc_full[:, 1] : ComplexF64[]
    pzgesv!(Aloc, bloc, grid; MB = MB, NB = NB)
    x = gather_solution(bloc, grid; MB = MB, NB = NB)
    _blacs_gridexit!(grid.ictxt)
    return x
end
