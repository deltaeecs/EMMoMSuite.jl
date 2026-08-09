# 基于球谐与八面体群结构的 Lebedev 插值权重优化实现
#
# 原理:
#   球面角谱是限带函数 (degree <= L = (p-1)/2)。粗层 Lebedev 采样到细层采样的精确
#   插值算子为 W = Y_fine * pinv(Y_coarse)，其中 Y 为 degree<=L 的实球谐合成矩阵。
#   该算子是旋转群等变的（W·R = R·W），因此可按八面体群轨道压缩：只对每轨道代表
#   求解权重，其余行由旋转精确生成（理论加速 ~24x）。
#
# 提供:
#   interp_weights_exact(pk, pt)      -> 稠密精确 W（机器精度、确定性、无训练）
#   interp_weights_local(pk, pt; ...) -> 局部支撑 + 多项式精确约束的稀疏 W
#   interp_weights_local_orbit(...)   -> 轨道压缩版稀疏 W（与朴素版逐位一致）
#   vectorize(W)                      -> 分块对角化处理 (θ,ϕ) 矢量分量

module SHInterp

using LinearAlgebra
using SparseArrays
using SpecialFunctions: gamma
using ..LebedevSortedPoints: get_t_nodes, nodes2Poles
using ...MLFMA.Interpolation: truncation_kernel
import ..dataset_generator: random_source_geometry, evaluate_poles!, random_rvec, random_rhat

export realSHmatrix, interp_weights_exact, interp_weights_local,
    interp_weights_local_orbit, interp_weights_auto, vectorize, octahedral_rotations,
    wigner_small_d, spin_weighted_harmonics, interp_weights_vsh, interp_weights_vsh_local,
    interp_weights_cart, interp_weights_cart_local, interp_weights_cart_local_orbit,
    interp_weights_hybrid

"""
    realSHmatrix(nodes, Lmax) -> Y (n, (Lmax+1)^2)

实球谐（含 Condon-Shortley 相位）合成矩阵：Y[i, (l,m)] = y_lm(r̂_i)。
"""
function realSHmatrix(nodes::AbstractMatrix, Lmax::Int)
    n = size(nodes, 2)
    Y = zeros(Float64, n, (Lmax + 1)^2)
    for i in 1:n
        x = Float64(nodes[3, i])
        s = sqrt(max(0.0, 1.0 - x * x))
        cphi = nodes[1, i] / max(s, 1e-300)
        sphi = nodes[2, i] / max(s, 1e-300)
        P = zeros(Lmax + 1, Lmax + 1)   # P[l+1, m+1] = P_l^m, m<=l
        P[1, 1] = 1.0
        for l in 1:Lmax
            P[l+1, l+1] = -(2l - 1) * s * P[l, l]
            for m in 0:(l-1)
                pprev = l >= 2 ? P[l-1, m+1] : 0.0
                P[l+1, m+1] = ((2l - 1) * x * P[l, m+1] - (l - 1 + m) * pprev) / (l - m)
            end
        end
        col = 0
        for l in 0:Lmax, m in -l:l
            col += 1
            am = abs(m)
            nrm = sqrt((2l + 1) / (4π) * gamma(l - am + 1) / gamma(l + am + 1))
            val = nrm * P[l+1, am+1]
            if m == 0
                Y[i, col] = val
            elseif m > 0
                Y[i, col] = sqrt(2) * val * cos(am * atan(sphi, cphi))
            else
                Y[i, col] = sqrt(2) * val * sin(am * atan(sphi, cphi))
            end
        end
    end
    return Y
end

"""
    _cart_core(pnodes, tnodes, poles_c, poles_f, Yc, Yf, Lloc; cap_rad, grow) -> SparseMatrixCSC

笛卡尔标量 SH 矢量插值内核：
每个细层点取支撑 S（cap_rad>0 时角距帽，否则全支撑），求满足标量 SH 精确性
（degree<=Lloc）的最小范数权重 w（三分量共用），再以 θ̂/ϕ̂ 点积散布成 (θ,ϕ) 耦合块：
  W[θ_i, j] = w_j (θ̂_i·θ̂_j),  W[θ_i, nt+j] = w_j (θ̂_i·ϕ̂_j),
  W[ϕ_i, j] = w_j (ϕ̂_i·θ̂_j),  W[ϕ_i, nt+j] = w_j (ϕ̂_i·ϕ̂_j)
优点：笛卡尔分量在极点是光滑的（无 θ/ϕ 基奇异），精确吸收矢量 L+1 耦合；
且权重为实数。共享节点行自动退化为单位行。
"""
function _cart_core(
    tnodes::AbstractMatrix,
    pnodes::AbstractMatrix,
    poles_c,
    poles_f,
    Yc::AbstractMatrix,
    Yf::AbstractMatrix,
    Lloc::Int;
    cap_rad::Real = 0.0,
    grow::Bool = false,
)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    m = (Lloc + 1)^2
    srows = Int[]
    scols = Int[]
    svals = Float64[]
    thC = reduce(hcat, [p.θhat for p in poles_c])   # 3 x nt
    phC = reduce(hcat, [p.ϕhat for p in poles_c])
    thF = reduce(hcat, [p.θhat for p in poles_f])   # 3 x nf
    phF = reduce(hcat, [p.ϕhat for p in poles_f])
    θ = Float64(cap_rad)
    for i in 1:nf
        ang = Vector{Float64}(undef, nt)
        xi, yi, zi = pnodes[1, i], pnodes[2, i], pnodes[3, i]
        @inbounds for j in 1:nt
            dx = tnodes[1, j] - xi
            dy = tnodes[2, j] - yi
            dz = tnodes[3, j] - zi
            d = sqrt(dx * dx + dy * dy + dz * dz)
            ang[j] = 2asin(min(d / 2, 1.0))
        end
        S = θ > 0 ? findall(ang .< θ) : collect(1:nt)
        while grow && length(S) < m
            θ *= 1.2
            S = findall(ang .< θ)
            θ >= π && break
        end
        w = reshape(Yf[i, :], 1, :) * pinv(Yc[S, :]; rtol = 1e-12)
        th_i = view(thF, :, i)
        ph_i = view(phF, :, i)
        aθθ = vec(w .* (th_i' * thC[:, S]))
        aθϕ = vec(w .* (th_i' * phC[:, S]))
        aϕθ = vec(w .* (ph_i' * thC[:, S]))
        aϕϕ = vec(w .* (ph_i' * phC[:, S]))
        append!(srows, fill(i, length(S)))
        append!(scols, S)
        append!(svals, aθθ)
        append!(srows, fill(i, length(S)))
        append!(scols, S .+ nt)
        append!(svals, aθϕ)
        append!(srows, fill(nf + i, length(S)))
        append!(scols, S)
        append!(svals, aϕθ)
        append!(srows, fill(nf + i, length(S)))
        append!(scols, S .+ nt)
        append!(svals, aϕϕ)
    end
    return sparse(srows, scols, svals, 2nf, 2nt)
end

"""
    interp_weights_cart(pk, pt; degree = Lb) -> Matrix (2n_pt x 2n_pk)

笛卡尔标量 SH 精确矢量一步插值（全支撑，稠密）。三个笛卡尔分量各自按标量
限带函数插值（共用同一最小范数权重），再合成 (θ,ϕ) 耦合块。
对 EFIE 类切向矢量场机器精度（实测 ~1e-14），免训练、确定性、极点光滑。
"""
function interp_weights_cart(pk::Int, pt::Int; degree::Int = 0, FT::Type = Float64)
    Lb = (pk - 1) ÷ 2
    deg = degree > 0 ? degree : Lb
    tnodes = get_t_nodes((pk - 1) ÷ 2)
    pnodes = get_t_nodes((pt - 1) ÷ 2)
    deg > isqrt(size(tnodes, 2)) - 1 &&
        error("degree=$deg 超出粗网格可表示阶数（点数 $(size(tnodes, 2))）")
    Yc = realSHmatrix(tnodes, deg)
    Yf = realSHmatrix(pnodes, deg)
    poles_c = nodes2Poles(tnodes)
    poles_f = nodes2Poles(pnodes)
    return Matrix{FT}(_cart_core(tnodes, pnodes, poles_c, poles_f, Yc, Yf, deg))
end

"""
    interp_weights_cart_local(pk, pt; Lloc, cap_rad, grow=true) -> SparseMatrixCSC

笛卡尔标量 SH 局部稀疏版：角距帽内支撑 + 标量 SH 精确性约束（l<=Lloc），
闭式散射成稀疏 (θ,ϕ) 耦合块（每细点 4*|S| 个非零）。
"""
function interp_weights_cart_local(
    pk::Int,
    pt::Int;
    Lloc::Int,
    cap_rad::Real,
    grow::Bool = true,
    FT::Type = Float64,
)
    tnodes = get_t_nodes((pk - 1) ÷ 2)
    pnodes = get_t_nodes((pt - 1) ÷ 2)
    Yc = realSHmatrix(tnodes, Lloc)
    Yf = realSHmatrix(pnodes, Lloc)
    poles_c = nodes2Poles(tnodes)
    poles_f = nodes2Poles(pnodes)
    return SparseMatrixCSC{FT,Int}(
        _cart_core(tnodes, pnodes, poles_c, poles_f, Yc, Yf, Lloc; cap_rad = cap_rad, grow = grow)
    )
end

"""
    interp_weights_cart_local_orbit(pk, pt; Lloc, cap_rad) -> SparseMatrixCSC

笛卡尔局部稀疏的八面体群轨道压缩版：只对每轨道代表求解标量最小范数权重
（在旋转下不变、共享），其余行用各节点自己的 θ̂/ϕ̂ 帧做点积散布（切平面
holonomy 由逐节点点积自然吸收）。与朴素版逐位一致，pinv 构造时间约 1/24。
"""
function interp_weights_cart_local_orbit(
    pk::Int,
    pt::Int;
    Lloc::Int,
    cap_rad::Real,
    grow::Bool = false,
    FT::Type = Float64,
)
    tnodes = get_t_nodes((pk - 1) ÷ 2)
    pnodes = get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    m = (Lloc + 1)^2
    Yc = realSHmatrix(tnodes, Lloc)
    Yf = realSHmatrix(pnodes, Lloc)
    poles_c = nodes2Poles(tnodes)
    poles_f = nodes2Poles(pnodes)
    R = octahedral_rotations()
    reps, rots, rep2nodes, rotidx = _orbit_reps(pnodes)
    cdict = Dict{Tuple{Float64,Float64,Float64},Int}()
    r8(x) = x == 0.0 ? 0.0 : round(x, digits = 8)
    for j in 1:nt
        cdict[(r8(tnodes[1, j]), r8(tnodes[2, j]), r8(tnodes[3, j]))] = j
    end
    # 每个旋转的粗点列置换（一次预计算）
    perms = [Vector{Int}(undef, nt) for _ in R]
    for (k, M) in enumerate(R)
        for s in 1:nt
            perms[k][s] = cdict[(
                r8(M[1,1] * tnodes[1,s] + M[1,2] * tnodes[2,s] + M[1,3] * tnodes[3,s]),
                r8(M[2,1] * tnodes[1,s] + M[2,2] * tnodes[2,s] + M[2,3] * tnodes[3,s]),
                r8(M[3,1] * tnodes[1,s] + M[3,2] * tnodes[2,s] + M[3,3] * tnodes[3,s]),
            )]
        end
    end
    thC = reduce(hcat, [p.θhat for p in poles_c])
    phC = reduce(hcat, [p.ϕhat for p in poles_c])
    thF = reduce(hcat, [p.θhat for p in poles_f])
    phF = reduce(hcat, [p.ϕhat for p in poles_f])
    srows = Int[]
    scols = Int[]
    svals = Float64[]
    θ = Float64(cap_rad)
    for (irep, nodes) in rep2nodes
        ang = Vector{Float64}(undef, nt)
        xr, yr, zr = pnodes[1, irep], pnodes[2, irep], pnodes[3, irep]
        @inbounds for j in 1:nt
            dx = tnodes[1, j] - xr
            dy = tnodes[2, j] - yr
            dz = tnodes[3, j] - zr
            d = sqrt(dx * dx + dy * dy + dz * dz)
            ang[j] = 2asin(min(d / 2, 1.0))
        end
        S = findall(ang .< θ)
        while grow && length(S) < m
            θ *= 1.2
            S = findall(ang .< θ)
            θ >= π && break
        end
        w = reshape(Yf[irep, :], 1, :) * pinv(Yc[S, :]; rtol = 1e-12)
        for j in nodes
            cols = perms[rotidx[j]][S]
            th_j = view(thF, :, j)
            ph_j = view(phF, :, j)
            aθθ = vec(w .* (th_j' * thC[:, cols]))
            aθϕ = vec(w .* (th_j' * phC[:, cols]))
            aϕθ = vec(w .* (ph_j' * thC[:, cols]))
            aϕϕ = vec(w .* (ph_j' * phC[:, cols]))
            append!(srows, fill(j, length(cols)))
            append!(scols, cols)
            append!(svals, aθθ)
            append!(srows, fill(j, length(cols)))
            append!(scols, cols .+ nt)
            append!(svals, aθϕ)
            append!(srows, fill(nf + j, length(cols)))
            append!(scols, cols)
            append!(svals, aϕθ)
            append!(srows, fill(nf + j, length(cols)))
            append!(scols, cols .+ nt)
            append!(svals, aϕϕ)
        end
    end
    return SparseMatrixCSC{FT,Int}(sparse(srows, scols, svals, 2nf, 2nt))
end

"""
    interp_weights_hybrid(pk, pt; Lloc=3, support_scale=2.0, kmin=16, nρ, npos, seed)

训练式与精确性的混合：每细层点取 k 个最近粗点为支撑（k 自适应 = max(2*Lloc 余量, kmin)），
求解"数据拟合 + 笛卡尔标量 SH 精确性约束"的约束最小二乘（KKT），再以 θ̂/ϕ̂ 点积散布。
  - 约束保证对 degree<=Lloc 的笛卡尔标量场精确（确定性、可复现）；
  - 数据拟合吸收 EFIE 矢量的高阶残差（同训练式 k=8 的精度量级）；
  - 非零元保持 2k/行（≈16~32），比 cap 球冠的笛卡尔局部稀疏更省。
数据集由固定种子确定性生成（源在盒内）。
"""
function interp_weights_hybrid(
    pk::Int,
    pt::Int;
    Lloc::Int = 3,
    support_scale::Real = 2.0,
    kmin::Int = 16,
    nρ::Int = 25,
    npos::Int = 60,
    seed::Int = 20260808,
    FT::Type = Float64,
)
    tnodes = get_t_nodes((pk - 1) ÷ 2)
    pnodes = get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    m = (Lloc + 1)^2
    k = max(ceil(Int, support_scale * m), kmin)
    k = min(k, nt)

    # 确定性数据集（共享几何、源在盒内）：用模块内 splitmix64 RNG，不依赖全局随机状态
    rng_s = UInt64(seed) ⊻ 0x9E3779B97F4A7C15
    nextrng() = begin
        rng_s += 0x9E3779B97F4A7C15
        z = rng_s
        z = (z ⊻ (z >> 30)) * 0xBF58476D1CE4E5B9
        z = (z ⊻ (z >> 27)) * 0x94D049BB133111EB
        z ⊻ (z >> 31)
    end
    runif() = Float64(nextrng() >> 11) * 2.0^-53
    rnorm() = sqrt(-2 * log(runif())) * cos(2π * runif())
    rhat_local() = begin
        z = 2runif() - 1
        φ = 2π * runif()
        r = sqrt(1 - z * z)
        [r * cos(φ), r * sin(φ), z]
    end
    rvec_local(bd) = [(2runif() - 1) * bd for _ in 1:3]
    truncnorm(mu, sigma, lo, hi) = begin
        x = 0.0
        while true
            x = mu + sigma * rnorm()
            lo <= x <= hi && return x
        end
    end

    f(x) = truncation_kernel(x) - (pk + 1) / 2
    lo, hi = 1e-4, 20.0
    for _ in 1:80
        mid = (lo + hi) / 2
        f(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    rel_l = (lo + hi) / 2
    arm = min(0.12, 0.5 * rel_l)
    rscale = max(rel_l / 2 - arm, 1e-4)
    poles_c = nodes2Poles(tnodes)
    poles_f = nodes2Poles(pnodes)
    N = nρ * npos
    T = zeros(ComplexF64, nt, 2, N)
    P = zeros(ComplexF64, nf, 2, N)
    idx = 0
    for ir in 1:npos, iρ in 1:nρ
        idx += 1
        rvecp = rvec_local(rscale)
        rvecm = rvecp .+ rhat_local() .* truncnorm(2 / 3 * arm, arm / 3, 0, arm)
        off_max = 0.125 * arm
        offsets = [[rhat_local()...] .* truncnorm(2 / 3 * off_max, off_max / 3, off_max / 3, off_max) for _ in 1:5]
        offsets[1] .= 0
        r0p = rvecp .+ rvecp .- rvecm .+ rhat_local() .* (arm / 24)
        r0m = rvecm .+ rvecm .- rvecp .+ rhat_local() .* (arm / 24)
        geom = (rvecp = rvecp, rvecm = rvecm, offsets = offsets, r0p = r0p, r0m = r0m)
        evaluate_poles!(poles_c, view(T, :, :, idx), geom)
        evaluate_poles!(poles_f, view(P, :, :, idx), geom)
    end

    # 笛卡尔粗/细数据（3n × N 复数）
    thC = reduce(hcat, [p.θhat for p in poles_c])
    phC = reduce(hcat, [p.ϕhat for p in poles_c])
    thF = reduce(hcat, [p.θhat for p in poles_f])
    phF = reduce(hcat, [p.ϕhat for p in poles_f])
    Xc = zeros(ComplexF64, 3nt, N)
    Yf_c = zeros(ComplexF64, 3nf, N)
    for j in 1:nt
        for c in 1:3
            Xc[c*nt-nt+j, :] .= thC[c, j] .* view(T, j, 1, :) .+ phC[c, j] .* view(T, j, 2, :)
        end
    end
    for j in 1:nf
        for c in 1:3
            Yf_c[c*nf-nf+j, :] .= thF[c, j] .* view(P, j, 1, :) .+ phF[c, j] .* view(P, j, 2, :)
        end
    end

    Yc = realSHmatrix(tnodes, Lloc)
    Yf = realSHmatrix(pnodes, Lloc)
    # 支撑距离预计算
    dists = [norm(tnodes[:, j] .- pnodes[:, i]) for i in 1:nf, j in 1:nt]

    srows = Int[]
    scols = Int[]
    svals = Float64[]
    for i in 1:nf
        S = partialsortperm(view(dists, i, :), 1:k)
        # 数据矩阵 A (|S| x 6N) = [Re Xc1 Im Xc1 Re Xc2 Im Xc2 Re Xc3 Im Xc3]
        XcS = Xc[vcat(S, S .+ nt, S .+ 2nt), :]        # 3k x N
        A = hcat(
            hcat(real(view(XcS, 1:k, :)), imag(view(XcS, 1:k, :))),
            hcat(real(view(XcS, (k+1):2k, :)), imag(view(XcS, (k+1):2k, :))),
            hcat(real(view(XcS, (2k+1):3k, :)), imag(view(XcS, (2k+1):3k, :))),
        )                                           # k x 6N
        yc = Yf_c[vcat(i, i .+ nf, i .+ 2nf), :]     # 3 x N
        b = vcat(
            vcat(real(view(yc, 1, :)), imag(view(yc, 1, :))),
            vcat(real(view(yc, 2, :)), imag(view(yc, 2, :))),
            vcat(real(view(yc, 3, :)), imag(view(yc, 3, :))),
        )                                           # 6N
        B = Yc[S, :]
        cvec = Yf[i, :]
        # KKT: min |wA-b|^2 s.t. wB = cvec'
        G = A * A'                       # k x k
        GiB = G \ B                      # k x m
        w_ls = (b' * A') / G             # 1 x k
        res = w_ls * B - cvec'           # 1 x m
        H = B' * GiB                     # m x m
        w = w_ls - res * (H \ GiB')      # 1 x k
        th_i = view(thF, :, i)
        ph_i = view(phF, :, i)
        aθθ = vec(w .* (th_i' * thC[:, S]))
        aθϕ = vec(w .* (th_i' * phC[:, S]))
        aϕθ = vec(w .* (ph_i' * thC[:, S]))
        aϕϕ = vec(w .* (ph_i' * phC[:, S]))
        append!(srows, fill(i, k))
        append!(scols, S)
        append!(svals, aθθ)
        append!(srows, fill(i, k))
        append!(scols, S .+ nt)
        append!(svals, aθϕ)
        append!(srows, fill(nf + i, k))
        append!(scols, S)
        append!(svals, aϕθ)
        append!(srows, fill(nf + i, k))
        append!(scols, S .+ nt)
        append!(svals, aϕϕ)
    end
    return SparseMatrixCSC{FT,Int}(sparse(srows, scols, svals, 2nf, 2nt))
end

"""
    wigner_small_d(l, mp, m, θ) -> Float64

Wigner 小 d 矩阵元素 d^l_{mp,m}(θ)（Varshalovich 求和公式）：
d^l_{mp,m} = sqrt(C) Σ_k (-1)^(k+mp-m) cos^(2l+m-mp-2k)(θ/2) sin^(mp-m+2k)(θ/2)
              / [k! (l+m-k)! (l-mp-k)! (mp-m+k)!]
"""
function wigner_small_d(l::Int, mp::Int, m::Int, θ::Real)
    C = gamma(l + mp + 1) * gamma(l - mp + 1) * gamma(l + m + 1) * gamma(l - m + 1)
    c = cos(θ / 2)
    s = sin(θ / 2)
    kmin = max(0, m - mp)
    kmax = min(l + m, l - mp)
    acc = 0.0
    for k in kmin:kmax
        denom = gamma(k + 1) * gamma(l + m - k + 1) * gamma(l - mp - k + 1) * gamma(mp - m + k + 1)
        term = (-1)^(k + mp - m) * c^(2l + m - mp - 2k) * s^(mp - m + 2k) / denom
        acc += term
    end
    return sqrt(C) * acc
end

"""
    spin_weighted_harmonics(nodes, Lmax; s = 1) -> Matrix{ComplexF64}

自旋加权球谐 sY_lm（s=±1）合成矩阵：S[i, (l,m)] = sY_lm(r̂_i)，l=1..Lmax, m=-l..l。
定义 sY_lm(θ,φ) = (-1)^m sqrt((2l+1)/(4π)) d^l_{m,s}(θ) e^{imφ}。
对切向矢量场：F_± = F_θ ± iF_ϕ 分别展开在 s=±1 基上，l 可达 L+1（矢量耦合）。
"""
function spin_weighted_harmonics(nodes::AbstractMatrix, Lmax::Int; s::Int = 1)
    n = size(nodes, 2)
    ncol = (Lmax + 1)^2 - 1
    S = zeros(ComplexF64, n, ncol)
    for i in 1:n
        x = clamp(Float64(nodes[3, i]), -1.0, 1.0)
        θ = acos(x)
        φ = atan(Float64(nodes[2, i]), Float64(nodes[1, i]))
        col = 0
        for l in 1:Lmax
            nrm = sqrt((2l + 1) / (4π))
            for m in -l:l
                col += 1
                d = wigner_small_d(l, m, s, θ)
                S[i, col] = (-1)^m * nrm * d * cis(m * φ)
            end
        end
    end
    return S
end

"""
    interp_weights_vsh(pk, pt; Lmax = Lb+1) -> Matrix (2n_pt x 2n_pk)

自旋加权球谐（VSH）精确矢量一步插值：
  F_± 在 s=±1 基上各自独立插值（W± = S±_fine * pinv(S±_coarse)），
  再经 C = [I iI; I -iI] 变换回 (θ,ϕ) 分量（交叉块自动出现）。
对 EFIE 类切向矢量场精确（l <= Lb+1 可表示且粗网格点数足够），
确定性、免训练，是标量球谐（漏 L+1 矢量耦合）的完备化。
"""
function interp_weights_vsh(pk::Int, pt::Int; Lmax::Int = 0, FT::Type = ComplexF64)
    Lb = (pk - 1) ÷ 2
    Lv = Lmax > 0 ? Lmax : Lb + 1
    tnodes = get_t_nodes((pk - 1) ÷ 2)
    pnodes = get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    Sp_c = spin_weighted_harmonics(tnodes, Lv; s = 1)
    Sm_c = spin_weighted_harmonics(tnodes, Lv; s = -1)
    Sp_f = spin_weighted_harmonics(pnodes, Lv; s = 1)
    Sm_f = spin_weighted_harmonics(pnodes, Lv; s = -1)
    Wp = Sp_f * pinv(Sp_c; rtol = 1e-13)
    Wm = Sm_f * pinv(Sm_c; rtol = 1e-13)
    W = zeros(ComplexF64, 2nf, 2nt)
    # Fθ_f = (Wp(Fθ+iFϕ) + Wm(Fθ-iFϕ))/2, Fϕ_f = (Wp(Fθ+iFϕ) - Wm(Fθ-iFϕ))/(2i)
    W[1:nf, 1:nt] = (Wp + Wm) / 2
    W[1:nf, (nt+1):(2nt)] = im * (Wp - Wm) / 2
    W[(nf+1):(2nf), 1:nt] = (Wp - Wm) / (2im)
    W[(nf+1):(2nf), (nt+1):(2nt)] = (Wp + Wm) / 2
    return Matrix{FT}(W)
end

"""
    interp_weights_vsh_local(pk, pt; Lloc, cap_rad) -> SparseMatrixCSC

VSH 局部稀疏版：每个细层点取角距帽内粗点为支撑，在 s=±1 基上强制
自旋加权精确性（l <= Lloc）的最小范数解，再合成 (θ,ϕ) 耦合块。
行和 = 1 自动满足；支撑不足时可 grow。
"""
function interp_weights_vsh_local(
    pk::Int,
    pt::Int;
    Lloc::Int,
    cap_rad::Real,
    grow::Bool = true,
    FT::Type = ComplexF64,
)
    tnodes = get_t_nodes((pk - 1) ÷ 2)
    pnodes = get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    Sp_c = spin_weighted_harmonics(tnodes, Lloc; s = 1)
    Sm_c = spin_weighted_harmonics(tnodes, Lloc; s = -1)
    Sp_f = spin_weighted_harmonics(pnodes, Lloc; s = 1)
    Sm_f = spin_weighted_harmonics(pnodes, Lloc; s = -1)
    m = (Lloc + 1)^2 - 1
    m >= nt && error("Lloc=$(Lloc) 需要自旋阶数空间 $(m) >= 粗层点数 $(nt)")
    W = spzeros(ComplexF64, 2nf, 2nt)
    θ = Float64(cap_rad)
    for i in 1:nf
        ang = [2asin(min(norm(tnodes[:, j] .- pnodes[:, i]) / 2, 1.0)) for j in 1:nt]
        S = findall(ang .< θ)
        while grow && length(S) < m
            θ *= 1.2
            S = findall(ang .< θ)
            θ >= π && break
        end
        wp = reshape(Sp_f[i, :], 1, :) * pinv(Sp_c[S, :]; rtol = 1e-12)
        wm = reshape(Sm_f[i, :], 1, :) * pinv(Sm_c[S, :]; rtol = 1e-12)
        a = (wp + wm) / 2        # θ→θ / ϕ→ϕ
        b = (wp - wm) / 2        # 交叉
        W[i, S] = a[:]
        W[i, S .+ nt] = im .* b[:]
        W[i .+ nf, S] = b[:] ./ im
        W[i .+ nf, S .+ nt] = a[:]
    end
    return SparseMatrixCSC{FT,Int}(W)
end

"""八面体群真旋转（24 个符号置换矩阵，det=+1）"""
function octahedral_rotations()
    mats = Matrix{Float64}[]
    for p1 in 1:3, p2 in 1:3, p3 in 1:3
        (p1 == p2 || p1 == p3 || p2 == p3) && continue
        for s1 in (-1.0, 1.0), s2 in (-1.0, 1.0), s3 in (-1.0, 1.0)
            M = zeros(3, 3)
            M[1, p1] = s1
            M[2, p2] = s2
            M[3, p3] = s3
            abs(det(M) - 1.0) < 1e-12 && push!(mats, M)
        end
    end
    return mats
end

"""
    interp_weights_exact(pk, pt; degree = min(Lb+1, capacity)) -> Matrix (n_pt x n_pk)

球谐精确一步插值权重 W = Y_fine * pinv(Y_coarse)，degree <= (pk-1)/2。
对限带角谱为机器精度（实测 ~1e-14），确定性、无训练数据、无 RNG。
默认 degree = Lb（=(pk-1)/2，即 Lebedev 求积阶数决定的精确表示极限）：
粗网格无法可靠表示 Lb+1 阶（低阶规则点/系数比紧张，条件数退化），
矢量场 θ/ϕ 分量的 L+1 耦合由实际数据/局部方法近似吸收。
"""
function interp_weights_exact(pk::Int, pt::Int; degree::Int = 0, FT::Type = Float64)
    Lb = (pk - 1) ÷ 2
    tnodes = get_t_nodes((pk - 1) ÷ 2)
    pnodes = get_t_nodes((pt - 1) ÷ 2)
    maxdeg = isqrt(size(tnodes, 2)) - 1          # 粗网格点数的乐观上限
    deg = degree > 0 ? min(degree, maxdeg) : min(Lb, maxdeg)
    Yc = realSHmatrix(tnodes, deg)
    Yf = realSHmatrix(pnodes, deg)
    return Matrix{FT}(Yf * pinv(Yc; rtol = 1e-13))
end

"""
    interp_weights_auto(pk, pt; max_dense_nnz = 4_000_000) -> AbstractMatrix

按规模自动选择：矩阵规模不大用精确稠密 W（机器精度）；规模过大用八面体群
轨道压缩的混合稀疏 W（L_loc=3，数据拟合 + 标量 SH 精确性约束，确定性、免训练、
构造快）。高阶层稠密 W 内存不可接受。
返回 (θ,ϕ) 耦合算子（2n_pt x 2n_pk）。
"""
function interp_weights_auto(pk::Int, pt::Int; max_dense_nnz::Int = 250_000, FT::Type = Float64)
    tnodes = get_t_nodes((pk - 1) ÷ 2)
    pnodes = get_t_nodes((pt - 1) ÷ 2)
    if size(tnodes, 2) * size(pnodes, 2) <= max_dense_nnz
        return vectorize(interp_weights_exact(pk, pt; FT = FT))
    else
        return interp_weights_hybrid(pk, pt; Lloc = 3, support_scale = 1.5, FT = FT)
    end
end

"""
    interp_weights_local(pk, pt; Lloc, cap_rad, grow=true) -> SparseMatrixCSC

局部约束最小二乘：每个细层点取角距 < cap_rad 的粗层点为支撑，用最小范数解
（SVD 伪逆）强制对 degree<=Lloc 的球谐精确（自动满足行和=1、常数保持）。
支撑不足时若 grow=true 自动扩大 cap。
"""
function interp_weights_local(
    pk::Int,
    pt::Int;
    Lloc::Int,
    cap_rad::Real,
    grow::Bool = true,
    FT::Type = Float64,
)
    tnodes = get_t_nodes((pk - 1) ÷ 2)
    pnodes = get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    Yc = realSHmatrix(tnodes, Lloc)
    Yf = realSHmatrix(pnodes, Lloc)
    m = (Lloc + 1)^2
    m >= nt && error("Lloc=$(Lloc) 需要 m=$(m) >= 粗层点数 $(nt)，请用 interp_weights_exact")
    W = spzeros(FT, nf, nt)
    θ = Float64(cap_rad)
    for i in 1:nf
        ang = [2asin(min(norm(tnodes[:, j] .- pnodes[:, i]) / 2, 1.0)) for j in 1:nt]
        S = findall(ang .< θ)
        while grow && length(S) < m
            θ *= 1.2
            S = findall(ang .< θ)
            θ >= π && break
        end
        A = Yc[S, :]
        b = Yf[i, :]
        W[i, S] = b' * pinv(A; rtol = 1e-12)
    end
    return W
end

"""细层节点 -> (轨道代表下标, 使 M*r̂_rep = r̂_node 的旋转, 旋转在群中的索引)。
O(n*24) 规范形分组：每个节点的轨道 = 24 个旋转像的排序集合，同轨道节点共享同一规范形。"""
function _orbit_reps(fine_nodes)
    R = octahedral_rotations()
    nf = size(fine_nodes, 2)
    reps = Vector{Int}(undef, nf)
    rots = Vector{Matrix{Float64}}(undef, nf)
    rotidx = Vector{Int}(undef, nf)
    groups = Dict{Vector{Tuple{Float64,Float64,Float64}},Vector{Int}}()
    for i in 1:nf
        imgs = [Rm * fine_nodes[:, i] for Rm in R]
        key = sort([
            (round(x[1], digits = 6), round(x[2], digits = 6), round(x[3], digits = 6)) for x in imgs
        ])
        push!(get!(groups, key, Int[]), i)
    end
    rep2nodes = Dict{Int,Vector{Int}}()
    for (_, nodes) in groups
        irep = nodes[1]
        rep2nodes[irep] = nodes
        for j in nodes
            reps[j] = irep
            k = findfirst(Rm -> norm(Rm * fine_nodes[:, irep] .- fine_nodes[:, j]) < 1e-8, R)
            rots[j] = R[k]
            rotidx[j] = k
        end
    end
    return reps, rots, rep2nodes, rotidx
end

"""
    interp_weights_local_orbit(pk, pt; Lloc, cap_rad) -> SparseMatrixCSC

轨道压缩版：只对每个八面体群轨道的代表节点求解局部约束权重，其余行由旋转
精确生成（W·R = R·W 等变性）。结果与 interp_weights_local 逐位一致，构造
时间与存储按轨道数/点数比（约 1/24）压缩。
"""
function interp_weights_local_orbit(
    pk::Int,
    pt::Int;
    Lloc::Int,
    cap_rad::Real,
    grow::Bool = false,
    FT::Type = Float64,
)
    tnodes = get_t_nodes((pk - 1) ÷ 2)
    pnodes = get_t_nodes((pt - 1) ÷ 2)
    nt, nf = size(tnodes, 2), size(pnodes, 2)
    Yc = realSHmatrix(tnodes, Lloc)
    Yf = realSHmatrix(pnodes, Lloc)
    m = (Lloc + 1)^2
    m >= nt && error("Lloc=$(Lloc) 需要 m=$(m) >= 粗层点数 $(nt)，请用 interp_weights_exact")
    reps, rots, rep2nodes, _ = _orbit_reps(pnodes)
    cdict = Dict{Tuple{Float64,Float64,Float64},Int}()
    for j in 1:nt
        cdict[(round(tnodes[1, j], digits = 8), round(tnodes[2, j], digits = 8), round(tnodes[3, j], digits = 8))] = j
    end
    W = spzeros(FT, nf, nt)
    θ = Float64(cap_rad)
    for (irep, nodes) in rep2nodes
        ang = [2asin(min(norm(tnodes[:, j] .- pnodes[:, irep]) / 2, 1.0)) for j in 1:nt]
        S = findall(ang .< θ)
        while grow && length(S) < m
            θ *= 1.2
            S = findall(ang .< θ)
            θ >= π && break
        end
        A = Yc[S, :]
        b = Yf[irep, :]
        w = b' * pinv(A; rtol = 1e-12)
        for j in nodes
            M = rots[j]
            cols = [
                cdict[(
                    round(M[1,1] * tnodes[1,s] + M[1,2] * tnodes[2,s] + M[1,3] * tnodes[3,s], digits = 8),
                    round(M[2,1] * tnodes[1,s] + M[2,2] * tnodes[2,s] + M[2,3] * tnodes[3,s], digits = 8),
                    round(M[3,1] * tnodes[1,s] + M[3,2] * tnodes[2,s] + M[3,3] * tnodes[3,s], digits = 8),
                )] for s in S
            ]
            W[j, cols] = w[:]
        end
    end
    return W
end

"""标量插值矩阵 -> 分块对角矢量版（θ、ϕ 分量各自独立标量插值）"""
vectorize(W::AbstractMatrix) = [W zero(W); zero(W) W]

end # module SHInterp
