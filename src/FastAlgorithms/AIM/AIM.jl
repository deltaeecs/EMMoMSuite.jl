module AIM

using LinearAlgebra
using SparseArrays
using StaticArrays
using FFTW
using MPI
using ...CoreModule
using ...Geometry
using ...BasisFunctions
using ...IntegralEquations
using ...IntegralEquations.Impedance: get_triangles_info
using ...IntegralEquations.EFIEModule: efie_interaction!, EFIE, efie_from_keta

export AIMOperator, AIMOperatorMPI

# Dunavant 7 点三角形求积（次数 5，重心坐标，权重和 = 1）
const _TRI7_L = (
    SVector{3,Float64}(1 / 3, 1 / 3, 1 / 3),
    SVector{3,Float64}(0.5, 0.5, 0.0),
    SVector{3,Float64}(0.5, 0.0, 0.5),
    SVector{3,Float64}(0.0, 0.5, 0.5),
    SVector{3,Float64}(1 / 6, 1 / 6, 2 / 3),
    SVector{3,Float64}(1 / 6, 2 / 3, 1 / 6),
    SVector{3,Float64}(2 / 3, 1 / 6, 1 / 6),
)
const _TRI7_W = (0.225, 0.13239415, 0.13239415, 0.13239415, 0.12593918, 0.12593918, 0.12593918)

"""
    AIMGrid{FT}

覆盖目标包围盒（外扩 2h）的均匀笛卡尔网格。节点间距 `h`，节点数 `n`。
"""
struct AIMGrid{FT<:Real}
    origin::SVector{3,FT}
    n::NTuple{3,Int}
    h::FT
end

@inline function node_index(g::AIMGrid, ix::Int, iy::Int, iz::Int)
    return ix + g.n[1] * (iy - 1) + g.n[1] * g.n[2] * (iz - 1)
end

@inline function node_coord(g::AIMGrid{FT}, idx::Int) where {FT}
    nx, ny, nz = g.n
    iz = div(idx - 1, nx * ny) + 1
    rem1 = idx - (iz - 1) * nx * ny
    iy = div(rem1 - 1, nx) + 1
    ix = rem1 - (iy - 1) * nx
    return SVector{3,FT}(
        g.origin[1] + (ix - 1) * g.h,
        g.origin[2] + (iy - 1) * g.h,
        g.origin[3] + (iz - 1) * g.h,
    )
end

function grid_from_mesh(mesh, h::FT) where {FT}
    nds = vertices(mesh)
    lo = minimum(nds, dims = 2) .- 2 * h
    hi = maximum(nds, dims = 2) .+ 2 * h
    n = Tuple(ceil(Int, (hi[i] - lo[i]) / h) + 1 for i in 1:3)
    return AIMGrid{FT}(SVector{3,FT}(lo), n, h)
end

@inline function greens(k::Real, R::Real)
    return R > 1e-14 ? exp(-im * k * R) / (4π * R) : zero(ComplexF64)
end

"""
    build_kernel(g::AIMGrid, k) -> Ghat

构造 Toeplitz→循环嵌入的格林函数核的 FFT（尺寸 2n）。
"""
function build_kernel(g::AIMGrid{FT}, k::Real) where {FT}
    P = (2g.n[1], 2g.n[2], 2g.n[3])
    K = zeros(Complex{FT}, P...)
    @inbounds for iz in 1:P[3], iy in 1:P[2], ix in 1:P[1]
        dx = (ix - 1 < g.n[1] ? ix - 1 : ix - 1 - 2g.n[1]) * g.h
        dy = (iy - 1 < g.n[2] ? iy - 1 : iy - 1 - 2g.n[2]) * g.h
        dz = (iz - 1 < g.n[3] ? iz - 1 : iz - 1 - 2g.n[3]) * g.h
        K[ix, iy, iz] = greens(k, sqrt(dx * dx + dy * dy + dz * dz))
    end
    return fft(K)
end

"""
    conv3!(out, U, Ghat, P, n)

均匀网格上的卷积 `V[i] = Σ_j G(r_i - r_j) U[j]`（FFT，循环嵌入），
`U`/`out` 为长度 `prod(n)` 的展平数组。
"""
function conv3!(out::AbstractVector{CT}, U::AbstractVector{CT}, Ghat::AbstractArray, P::NTuple{3,Int}, n::NTuple{3,Int}) where {CT}
    Up = zeros(CT, P...)
    Up[1:n[1], 1:n[2], 1:n[3]] .= reshape(U, n...)
    F = fft(Up)
    F .*= Ghat
    ifft!(F)
    out .= vec(view(F, 1:n[1], 1:n[2], 1:n[3]))
    return out
end

"""
    conv3!(out, U, Ghat, P, n, ws)

零分配版本：复用 `ws = (Up, pfft, pifft)` 的填充缓冲与 FFTW 计划。
"""
function conv3!(
    out::AbstractVector{CT},
    U::AbstractVector{CT},
    Ghat::AbstractArray,
    P::NTuple{3,Int},
    n::NTuple{3,Int},
    ws,
) where {CT}
    Up, pfft, pifft = ws
    fill!(Up, 0)
    Up[1:n[1], 1:n[2], 1:n[3]] .= reshape(U, n...)
    pfft * Up
    Up .*= Ghat
    pifft * Up
    out .= vec(view(Up, 1:n[1], 1:n[2], 1:n[3]))
    return out
end

"""
    build_projection(basis::RWGBasis, g::AIMGrid, tris)

对每个 RWG 基函数计算其在网格节点上的投影：
- `Pv[m][q]`：`∫ f_m Λ_i dS`（3 维向量，节点 `nodes[m][q]`）；
- `Pd[m][q]`：`∫ (∇·f_m) Λ_i dS`（标量）。
约定与 MLFMA `add_radiation_pattern_rwg!` 一致：`f·dS = sign·l/2·ρ·w_qp`，
`∇·f = sign·l/A`（常数）。
"""
function build_projection(basis::RWGBasis, g::AIMGrid{FT}, tris) where {FT}
    N = length(basis.functions)
    nodes = Vector{Vector{Int}}(undef, N)
    Pv = Vector{Vector{SVector{3,FT}}}(undef, N)
    Pd = Vector{Vector{FT}}(undef, N)
    for (m, bf) in enumerate(basis.functions)
        acc = Dict{Int,Tuple{SVector{3,FT},FT}}()
        for i_supp in 1:2
            tri_idx = bf.support[i_supp]
            tri_idx == 0 && continue
            verts = tris[tri_idx].vertices
            v_opp = verts[:, bf.local_edge_idx[i_supp]]
            sgn = FT(bf.signs[i_supp])
            l = FT(bf.edge_length)
            for q in 1:7
                L = _TRI7_L[q]
                w = _TRI7_W[q]
                r = verts * L
                rho = r - v_opp
                c = clamp.(floor.(Int, (r .- g.origin) ./ g.h) .+ 1, 1, g.n .- 1)
                u = (r .- g.origin) ./ g.h .- (c .- 1)
                for dz in 0:1, dy in 0:1, dx in 0:1
                    lam =
                        (dx == 0 ? 1 - u[1] : u[1]) *
                        (dy == 0 ? 1 - u[2] : u[2]) *
                        (dz == 0 ? 1 - u[3] : u[3])
                    idx = node_index(g, c[1] + dx, c[2] + dy, c[3] + dz)
                    pv = sgn * l / 2 * w * lam * rho
                    pd = sgn * l * w * lam
                    if haskey(acc, idx)
                        old = acc[idx]
                        acc[idx] = (old[1] + pv, old[2] + pd)
                    else
                        acc[idx] = (pv, pd)
                    end
                end
            end
        end
        ks = sort!(collect(keys(acc)))
        nodes[m] = ks
        Pv[m] = [acc[k][1] for k in ks]
        Pd[m] = [acc[k][2] for k in ks]
    end
    return nodes, Pv, Pd
end

function _box_overlap(b1, b2)
    return b1[1] <= b2[2] && b2[1] <= b1[2] &&
           b1[3] <= b2[4] && b2[3] <= b1[4] &&
           b1[5] <= b2[6] && b2[5] <= b1[6]
end

"""
    assemble_near_correction(op, basis, g, nodes, Pv, Pd, k, near_radius)

近场修正（空间哈希版）：对包围盒（外扩 `near_radius`）相交的基函数对，
`Z_near = Z_direct − Z_grid`（Z_grid 为网格近似，同一 G(0)=0 约定），
同时返回近场**直接**矩阵 `Z_direct`（供预条件使用）。
候选对由"基函数中心分箱（箱大小 2·near_radius）+ Chebyshev ≤ 2 邻箱"生成，
再经包围盒重叠过滤——结果与 O(N²) 全遍历逐项一致（测试验证）。
"""
function assemble_near_correction(
    op::EFIE,
    basis::RWGBasis,
    g::AIMGrid{FT},
    nodes,
    Pv,
    Pd,
    k::Real,
    near_radius::Real,
    ;
    rows = nothing,
) where {FT}
    N = length(basis.functions)
    rows === nothing && (rows = 1:N)
    tris = get_triangles_info(basis.mesh, basis)
    boxes = Vector{NTuple{6,FT}}(undef, N)
    cellof = Vector{NTuple{3,Int}}(undef, N)
    allc = reduce(hcat, [bf.center for bf in basis.functions])
    lo = minimum(allc, dims = 2)
    for m in 1:N
        bf = basis.functions[m]
        lo_b = fill(FT(Inf), 3)
        hi_b = fill(FT(-Inf), 3)
        for i_supp in 1:2
            tri_idx = bf.support[i_supp]
            tri_idx == 0 && continue
            for v in eachcol(tris[tri_idx].vertices)
                for d in 1:3
                    lo_b[d] = min(lo_b[d], v[d])
                    hi_b[d] = max(hi_b[d], v[d])
                end
            end
        end
        boxes[m] = (
            lo_b[1]-near_radius, hi_b[1]+near_radius,
            lo_b[2]-near_radius, hi_b[2]+near_radius,
            lo_b[3]-near_radius, hi_b[3]+near_radius,
        )
        cellof[m] = Tuple(floor.(Int, (allc[:, m] .- lo) ./ (2 * near_radius)))
    end

    cellmap = Dict{NTuple{3,Int},Vector{Int}}()
    for m in 1:N
        push!(get!(cellmap, cellof[m], Int[]), m)
    end

    # 近场校正的 EFIE 核须与算子参数一致（k/eta/factor），否则自定义 eta/factor
    # 的 EFIE 会出现近场校正因子与远场 jkη 因子不一致。
    efie_op = efie_from_keta(op.k, op.eta, op.factor)
    C = im * k * op.eta
    Is = Int[]
    Js = Int[]
    Vs = ComplexF64[]
    IsD = Int[]
    JsD = Int[]
    VsD = ComplexF64[]
    Z3 = zeros(ComplexF64, 3, 3)
    for m in rows
        bf_m = basis.functions[m]
        # 自对（m==n）：盒子必然重叠，直接计算
        Zdir = zero(ComplexF64)
        for a in 1:2, b in 1:2
            tm = bf_m.support[a]
            sn = bf_m.support[b]
            (tm == 0 || sn == 0) && continue
            fill!(Z3, 0)
            efie_interaction!(Z3, efie_op, tris[tm], tris[sn])
            Zdir +=
                Z3[bf_m.local_edge_idx[a], bf_m.local_edge_idx[b]] *
                bf_m.signs[a] * bf_m.signs[b]
        end
        Zgrid = zero(ComplexF64)
        for (qi, ni) in enumerate(nodes[m]), (qj, nj) in enumerate(nodes[m])
            G = greens(k, norm(node_coord(g, ni) - node_coord(g, nj)))
            Zgrid += dot(Pv[m][qi], Pv[m][qj]) * G
            Zgrid -= Pd[m][qi] * Pd[m][qj] * G / k^2
        end
        val = Zdir - C * Zgrid
        if abs(val) > 1e-300
            push!(Is, m)
            push!(Js, m)
            push!(Vs, val)
        end
        if abs(Zdir) > 1e-300
            push!(IsD, m)
            push!(JsD, m)
            push!(VsD, Zdir)
        end
        seen = Set{Int}()
        push!(seen, m)
        for dz in -2:2, dy in -2:2, dx in -2:2
            key = (cellof[m][1]+dx, cellof[m][2]+dy, cellof[m][3]+dz)
            for n in get(cellmap, key, Int[])
                n in seen && continue
                push!(seen, n)
                _box_overlap(boxes[m], boxes[n]) || continue
                bf_n = basis.functions[n]
                Zdir = zero(ComplexF64)
                for a in 1:2, b in 1:2
                    tm = bf_m.support[a]
                    sn = bf_n.support[b]
                    (tm == 0 || sn == 0) && continue
                    fill!(Z3, 0)
                    efie_interaction!(Z3, efie_op, tris[tm], tris[sn])
                    Zdir +=
                        Z3[bf_m.local_edge_idx[a], bf_n.local_edge_idx[b]] *
                        bf_m.signs[a] * bf_n.signs[b]
                end
                Zgrid = zero(ComplexF64)
                for (qi, ni) in enumerate(nodes[m]), (qj, nj) in enumerate(nodes[n])
                    G = greens(k, norm(node_coord(g, ni) - node_coord(g, nj)))
                    Zgrid += dot(Pv[m][qi], Pv[n][qj]) * G
                    Zgrid -= Pd[m][qi] * Pd[n][qj] * G / k^2
                end
                Zgrid *= C
                val = Zdir - Zgrid
                if abs(val) > 1e-300
                    push!(Is, m)
                    push!(Js, n)
                    push!(Vs, val)
                end
                if abs(Zdir) > 1e-300
                    push!(IsD, m)
                    push!(JsD, n)
                    push!(VsD, Zdir)
                end
            end
        end
    end
    return sparse(Is, Js, Vs, N, N), sparse(IsD, JsD, VsD, N, N)
end

"""
    AIMWS

AIM matvec 的按线程工作区（具体类型字段避免动态分配）。
"""
mutable struct AIMWS
    UA::Matrix{ComplexF64}
    UB::Vector{ComplexF64}
    VA::Matrix{ComplexF64}
    VB::Vector{ComplexF64}
    Up::Array{ComplexF64,3}
    pfft::Any
    pifft::Any
    UAs::Vector{Matrix{ComplexF64}}
    UBs::Vector{Vector{ComplexF64}}
end

"""
    AIMOperator(op::EFIE, basis::RWGBasis; h_ratio=0.1, near_radius=0.35)

IE-FFT/AIM 型算子：RWG 投影到均匀网格 + 标量格林函数 FFT 卷积 + 近场修正。
矩阵向量积 `Zx ≈ Z_near x + C·(A − B/k²)`，其中
`A_mn = Σ Πv_mi·G_ij·Πv_nj`、`B_mn = Σ Πd_mi·G_ij·Πd_nj`，`C = jkη`
（与 EFIE 直接求解器 `factor = jkη/16π` 的归一化等价，已用稠密矩阵实测校准）。
`h_ratio` 为网格间距（以 λ 计，默认 0.1）；`near_radius` 为近场修正的固定
物理距离（以 λ 计，默认 0.35 ≈ 3.5 个网格单元）：包围盒外扩该距离相交的基函数对使用
`Z_direct − Z_grid` 精确修正，避免网格近似承担近奇异对。
"""
struct AIMOperator{FT,CT} <: AbstractIntegralOperator
    freq::FT
    k::FT
    eta::FT
    grid::AIMGrid{FT}
    Ghat::Array{Complex{FT},3}
    nodes::Vector{Vector{Int}}
    Pv::Vector{Vector{SVector{3,FT}}}
    Pd::Vector{Vector{FT}}
    Z_near::SparseMatrixCSC{CT,Int}
    Z_near_direct::SparseMatrixCSC{CT,Int}
    work::Vector{AIMWS}
end

function AIMOperator(op::EFIE, basis::RWGBasis; h_ratio = 0.1, near_radius = 0.35)
    FT = eltype(basis.mesh.node)
    k = FT(op.k)
    λ = 2π / k
    g = grid_from_mesh(basis.mesh, FT(h_ratio * λ))
    tris = get_triangles_info(basis.mesh, basis)
    nodes, Pv, Pd = build_projection(basis, g, tris)
    Ghat = build_kernel(g, k)
    Z_near, Z_near_direct = assemble_near_correction(
        op, basis, g, nodes, Pv, Pd, k, FT(near_radius * λ)
    )
    CT = eltype(Z_near)
    NP = prod(g.n)
    P = size(Ghat)
    work = [
        AIMWS(
            zeros(ComplexF64, 3, NP),
            zeros(ComplexF64, NP),
            zeros(ComplexF64, 3, NP),
            zeros(ComplexF64, NP),
            zeros(ComplexF64, P...),
            nothing,
            nothing,
            [zeros(ComplexF64, 3, NP) for _ in 1:Threads.maxthreadid()],
            [zeros(ComplexF64, NP) for _ in 1:Threads.maxthreadid()],
        ) for _ in 1:Threads.nthreads()
    ]
    for w in work
        w.pfft = FFTW.plan_fft!(w.Up)
        w.pifft = FFTW.plan_ifft!(w.Up)
    end
    return AIMOperator{FT,CT}(
        FT(op.freq),
        k,
        FT(op.eta),
        g,
        Ghat,
        nodes,
        Pv,
        Pd,
        Z_near,
        Z_near_direct,
        work,
    )
end

Base.eltype(A::AIMOperator{FT,CT}) where {FT,CT} = CT
Base.size(A::AIMOperator) = size(A.Z_near)
Base.size(A::AIMOperator, i::Int) = size(A.Z_near, i)

function Base.:*(A::AIMOperator, x::AbstractVector)
    y = similar(x)
    mul!(y, A, x)
    return y
end

function _aim_ws(A)
    tid = Threads.threadid()
    ws = A.work
    if length(ws) < tid
        return ws[end]
    end
    return ws[tid]
end

function _aim_project!(UA, UB, A, x, rng)
    for m in rng
        xm = x[m]
        for q in eachindex(A.nodes[m])
            idx = A.nodes[m][q]
            pv = A.Pv[m][q]
            UA[1, idx] += xm * pv[1]
            UA[2, idx] += xm * pv[2]
            UA[3, idx] += xm * pv[3]
            UB[idx] += xm * A.Pd[m][q]
        end
    end
    return nothing
end

function _aim_test!(y, A, VA, VB, C, rng)
    for m in rng
        sA = zero(ComplexF64)
        sB = zero(ComplexF64)
        for q in eachindex(A.nodes[m])
            idx = A.nodes[m][q]
            pv = A.Pv[m][q]
            sA += pv[1] * VA[1, idx] + pv[2] * VA[2, idx] + pv[3] * VA[3, idx]
            sB += A.Pd[m][q] * VB[idx]
        end
        y[m] += C * (sA - sB / A.k^2)
    end
    return nothing
end

function LinearAlgebra.mul!(y::AbstractVector, A::AIMOperator, x::AbstractVector)
    mul!(y, A.Z_near, x)
    ws = _aim_ws(A)
    UA, UB, VA, VB, Up, pfft, pifft = ws.UA, ws.UB, ws.VA, ws.VB, ws.Up, ws.pfft, ws.pifft
    nx, ny, nz = A.grid.n
    NP = nx * ny * nz
    P = size(A.Ghat)
    fill!(UA, 0)
    fill!(UB, 0)
    nt = Threads.nthreads()
    if nt > 1
        UAs = ws.UAs
        UBs = ws.UBs
        for t in 1:length(UAs)
            fill!(UAs[t], 0)
            fill!(UBs[t], 0)
        end
        Threads.@threads :static for m in eachindex(x)
            tid = Threads.threadid()
            _aim_project!(UAs[tid], UBs[tid], A, x, m:m)
        end
        # 归约全部 maxthreadid 槽位（@threads 池实际可能使用 id > nthreads 的线程，
        # 未使用槽为零，无害）
        for t in 1:length(UAs)
            @views UA .+= UAs[t]
            @views UB .+= UBs[t]
        end
    else
        _aim_project!(UA, UB, A, x, eachindex(x))
    end
    fill!(VA, 0)
    fill!(VB, 0)
    for c in 1:3
        conv3!(
            reshape(view(VA, c, :), NP),
            reshape(view(UA, c, :), NP),
            A.Ghat,
            P,
            A.grid.n,
            (Up, pfft, pifft),
        )
    end
    conv3!(VB, UB, A.Ghat, P, A.grid.n, (Up, pfft, pifft))
    C = im * A.k * A.eta
    if nt > 1
        Threads.@threads :static for m in eachindex(x)
            _aim_test!(y, A, VA, VB, C, m:m)
        end
    else
        _aim_test!(y, A, VA, VB, C, eachindex(x))
    end
    return y
end

"""
    _row_partition(N, rank, nproc) → UnitRange{Int}

AIM 本地行分区（与 Parallel.row_partition 相同语义，避免跨模块依赖）。
"""
function _row_partition(N::Int, rank::Int, nproc::Int)
    base = N ÷ nproc
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

"""
    AIMOperatorMPI(op::EFIE, basis::RWGBasis; h_ratio=0.1, near_radius=0.35, comm=MPI.COMM_WORLD)

AIM 的 MPI+线程混合算子：每个秩构建完整投影/网格核（确定性一致），仅装配其**行分区**的
近场修正；matvec = 本地行近场 SpMV → Allreduce 完整 y_near；本地行投影 → Allreduce 网格
UA/UB（小数组）→ 各秩复制执行 FFT 卷积（网格小，代价低）→ 本地行测试 → Allreduce 完整 y_far。
秩内投影/测试/近场装配使用 `@threads`；结果与串行逐位一致（SPMD + 行分区）。
"""
struct AIMOperatorMPI{FT,CT} <: AbstractIntegralOperator
    freq::FT
    k::FT
    eta::FT
    grid::AIMGrid{FT}
    Ghat::Array{Complex{FT},3}
    rows::UnitRange{Int}
    nodes::Vector{Vector{Int}}
    Pv::Vector{Vector{SVector{3,FT}}}
    Pd::Vector{Vector{FT}}
    Z_near::SparseMatrixCSC{CT,Int}
    work::Vector{AIMWS}
    comm
end

function AIMOperatorMPI(
    op::EFIE,
    basis::RWGBasis;
    h_ratio = 0.1,
    near_radius = 0.35,
    comm = MPI.COMM_WORLD,
)
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    FT = eltype(basis.mesh.node)
    k = FT(op.k)
    λ = 2π / k
    g = grid_from_mesh(basis.mesh, FT(h_ratio * λ))
    tris = get_triangles_info(basis.mesh, basis)
    nodes, Pv, Pd = build_projection(basis, g, tris)
    Ghat = build_kernel(g, k)
    N = length(basis.functions)
    rows = _row_partition(N, rank, P)
    Z_near, _ = assemble_near_correction(
        op, basis, g, nodes, Pv, Pd, k, FT(near_radius * λ); rows = rows
    )
    Z_near = Z_near[rows, :]   # 裁剪到本地行（nloc × N）
    CT = eltype(Z_near)
    NP = prod(g.n)
    Pdim = size(Ghat)
    work = [
        AIMWS(
            zeros(ComplexF64, 3, NP),
            zeros(ComplexF64, NP),
            zeros(ComplexF64, 3, NP),
            zeros(ComplexF64, NP),
            zeros(ComplexF64, Pdim...),
            nothing,
            nothing,
            [zeros(ComplexF64, 3, NP) for _ in 1:Threads.maxthreadid()],
            [zeros(ComplexF64, NP) for _ in 1:Threads.maxthreadid()],
        ) for _ in 1:Threads.nthreads()
    ]
    for w in work
        w.pfft = FFTW.plan_fft!(w.Up)
        w.pifft = FFTW.plan_ifft!(w.Up)
    end
    return AIMOperatorMPI{FT,CT}(
        FT(op.freq), k, FT(op.eta), g, Ghat, rows, nodes, Pv, Pd, Z_near, work, comm
    )
end

Base.eltype(A::AIMOperatorMPI{FT,CT}) where {FT,CT} = CT
Base.size(A::AIMOperatorMPI) = (length(A.nodes), length(A.nodes))
Base.size(A::AIMOperatorMPI, i::Int) = size(A)[i]

function Base.:*(A::AIMOperatorMPI, x::AbstractVector)
    y = similar(x)
    mul!(y, A, x)
    return y
end

function LinearAlgebra.mul!(y::AbstractVector, A::AIMOperatorMPI, x::AbstractVector)
    N = length(A.nodes)
    nloc = length(A.rows)
    # 近场：本地行 SpMV
    yloc = zeros(ComplexF64, nloc)
    mul!(yloc, A.Z_near, x)

    # 网格：复制投影（O(N) 投影代价低，避免每次 matvec Allreduce 大网格数组）→
    # 复制卷积 → 本地行测试 → Allreduce y_far
    ws = _aim_ws(A)
    UA, UB, VA, VB, Up, pfft, pifft = ws.UA, ws.UB, ws.VA, ws.VB, ws.Up, ws.pfft, ws.pifft
    nx, ny, nz = A.grid.n
    NP = nx * ny * nz
    Pdim = size(A.Ghat)
    fill!(UA, 0)
    fill!(UB, 0)
    nt = Threads.nthreads()
    if nt > 1
        UAs = ws.UAs
        UBs = ws.UBs
        for t in 1:length(UAs)
            fill!(UAs[t], 0)
            fill!(UBs[t], 0)
        end
        Threads.@threads :static for m in eachindex(A.nodes)
            tid = Threads.threadid()
            _aim_project!(UAs[tid], UBs[tid], A, x, m:m)
        end
        for t in 1:length(UAs)
            @views UA .+= UAs[t]
            @views UB .+= UBs[t]
        end
    else
        _aim_project!(UA, UB, A, x, eachindex(A.nodes))
    end

    fill!(VA, 0)
    fill!(VB, 0)
    for c in 1:3
        conv3!(
            reshape(view(VA, c, :), NP),
            reshape(view(UA, c, :), NP),
            A.Ghat,
            Pdim,
            A.grid.n,
            (Up, pfft, pifft),
        )
    end
    conv3!(VB, UB, A.Ghat, Pdim, A.grid.n, (Up, pfft, pifft))
    C = im * A.k * A.eta
    y_far = zeros(ComplexF64, N)
    if nt > 1
        Threads.@threads :static for m in A.rows
            _aim_test!(y_far, A, VA, VB, C, m:m)
        end
    else
        _aim_test!(y_far, A, VA, VB, C, A.rows)
    end
    # 合并近场+远场本地行 → 一次 Allreduce 得到完整 y
    yloc .+= y_far[A.rows]
    fill!(y, 0)
    y[A.rows] .= yloc
    MPI.Allreduce!(y, +, A.comm)
    return y
end

end
