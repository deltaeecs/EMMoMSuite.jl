module SphericalHarmonics

using LinearAlgebra
using SpecialFunctions: sphericalbesselj, sphericalbessely, loggamma
using ..Interpolation: GLPolesInfo

export ncoeff, lm_degree, spectral_bandwidth, sh_matrix, m2l_matrix

# ─────────────────────────────────────────────────────────────────────────────
# 归一化球谐（复谐、Condon–Shortley 相位）
#   Y_lm(θ,φ) = Θ_lm(cosθ) · exp(i m φ)，Θ 为实值归一化连带 Legendre
#   m ≥ 0: Θ = N_lm P_l^m；m < 0: Y_lm = (-1)^{|m|} conj(Y_{l,|m|})
#
# 谱域 M2L 调研产物（issue #22 问题 3）：原型阶段 2 的数学工具集。
# 集成阶段 3 结论为负结果——对角方案在 kR < L 放大区依赖 GL 采样求积的
# 广义求和/混叠相消，带限谱域往返不可复现（详见
# docs/dev/m2l_short_range_spectral.md §6）；本模块保留为分析/教学工具。
# GL(θ, L+1 点) × 均匀半格(φ, 2(L+1) 点) 极点网格（levelIntegralInfoCal）对
# degree ≤ L 的球谐正交归一精确（原型验证至 1e-14）。
# ─────────────────────────────────────────────────────────────────────────────

"""连带 Legendre P_l^m(x)（未归一化，Condon–Shortley 相位，m ≥ 0）。
Numerical Recipes 风格稳定递推，l ≤ ~40 精度充足。"""
function plgndr(l::Int, m::Int, x::Float64)
    (m < 0 || m > l) && return 0.0
    pmm = 1.0
    if m > 0
        somx2 = sqrt(max(0.0, (1.0 - x) * (1.0 + x)))
        fact = 1.0
        for i = 1:m
            pmm *= -fact * somx2
            fact += 2.0
        end
    end
    l == m && return pmm
    pmmp1 = x * (2m + 1) * pmm
    l == m + 1 && return pmmp1
    pll = 0.0
    for ll = m+2:l
        pll = ((2ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m)
        pmm = pmmp1
        pmmp1 = pll
    end
    return pll
end

_norm_factor(l::Int, m::Int) = sqrt((2l + 1) / (4π)) *
    exp(0.5 * (loggamma(l - abs(m) + 1) - loggamma(l + abs(m) + 1)))

"""Θ_lm(cosθ)：Y_lm = Θ_lm·e^{imφ} 的实值 θ 部分。"""
function theta_part(l::Int, m::Int, x::Float64)
    if m >= 0
        return _norm_factor(l, m) * plgndr(l, m, x)
    else
        # Y_{l,-|m|} = (-1)^{|m|} conj(Y_{l,|m|}) ⇒ Θ = (-1)^{|m|} N P_l^{|m|}
        return (-1.0)^(-m) * _norm_factor(l, m) * plgndr(l, -m, x)
    end
end

function ylm(l::Int, m::Int, θ::Float64, φ::Float64)
    return theta_part(l, m, cos(θ)) * exp(im * m * φ)
end

_spherical_h2(l::Int, x::Float64) =
    complex(sphericalbesselj(l, x), -sphericalbessely(l, x))

# 系数索引 a = (l, m)：l = 0:L, m = -l:l，线性序 l²+(m+l+1)
ncoeff(L::Int) = (L + 1)^2
lm_degree(i::Int) = isqrt(i - 1)

"""
    spectral_bandwidth(kw; digits = 2.0) -> Int

谱域 M2L 源系数列截止度 `L_eff`：盒超额带宽公式（Chew/Cha 型）
`τ = √3·kw + 2.16·d^{2/3}·(kw)^{1/3}` 在低 digits 下的取值。

谱域 M2L 的稠密矩阵 `M_ab` 对源系数 `b` 的放大随 `l_b` 增长
（`~h_{l_a+l_b}(kR)`）；物理（盒内支撑）聚合系数按 `(k·r_s)^l/(2l−1)!!`
衰减，`l > τ(digits)` 后残差/插值噪声占主导——超过 `L_eff` 的列置零
（原型验证：1e-3 相对白噪下无截止谱域 rel ~1e4–4e5，L_eff=5–6 → 2e-3–3e-3）。
"""
function spectral_bandwidth(kw::Real; digits::Real = 2.0)
    kw > 0 || return 0
    τ = sqrt(3.0) * kw + 2.16 * Float64(digits)^(2 / 3) * Float64(kw)^(1 / 3)
    return max(1, floor(Int, τ))
end

"""
    sh_matrix(poles::GLPolesInfo, L) -> (Y, YtW)

极点网格上的球谐矩阵 `Y[p, a]`（nPoles × (L+1)²）与带权分析矩阵
`YtW = Y'·Diagonal(W)`（(L+1)² × nPoles）。分析 `F̂ = YtW·F`、综合
`G = Y·Ĝ`——对限带（degree ≤ L）方向图精确（GL 网格性质）。
"""
function sh_matrix(poles::GLPolesInfo{FT}, L::Int) where {FT<:Real}
    nP = length(poles.r̂sθsϕs)
    M = ncoeff(L)
    Y = Matrix{Complex{FT}}(undef, nP, M)
    ia = 0
    for l = 0:L, m = -l:l
        ia += 1
        for ip = 1:nP
            info = poles.r̂sθsϕs[ip]
            r̂ = info.r̂
            θ = acos(clamp(r̂[3], -1.0, 1.0))
            φ = atan(r̂[2], r̂[1])
            Y[ip, ia] = Complex{FT}(ylm(l, m, θ, φ))
        end
    end
    YtW = adjoint(Y) * Diagonal(complex.(FT.(poles.Wθϕs), zeros(FT, nP)))
    return Y, YtW
end

# ─────────────────────────────────────────────────────────────────────────────
# Gaunt 系数（因子化 θ 求积；φ 由选择律解析处理）
#   G[a,b,c] = ∫ conj(Y_a) Y_b Y_c dΩ = 2π δ_{m_c, m_a−m_b} Σ_i w_i Θ_a Θ_b Θ_c
#   （细 GL 网格 Nθ = 2L+4 对 cosθ 度 ≤ 4L 精确）
# 按 (ai, bi) 稀疏存储，进程内按 L 缓存（构建耗时 O(L⁵) 但只发生一次）。
#
# **c 度上限为 L（不是 l_a+l_b）**：谱域 M2L 复刻的是对角方案逐点采样的
# 截断核 T_{≤L}(k̂) = Σ_{l≤L}(−j)^l(2l+1)h_l(kR)P_l 的配对 ⟨T_{≤L}·Y_b, Ȳ_a⟩，
# 故 c 求和与核同步截在 L。若按 Gaunt 完备性放行 c ∈ (L, l_a+l_b]，会在
# kR < L（h_l 发散区）注入 |h_{>L}(kR)|~1e3–1e4 量级的伪耦合（实测 GD2X
# 层3 (9,9) 元素虚部 −5.3e4 vs 正确值 ~1.5e3）——那是级数发散尾巴，不是
# 算子内容；对角方案的 L 截断恰好构成正则化（经验上其求积稳健）。
# ─────────────────────────────────────────────────────────────────────────────

struct _LM
    l::Int
    m::Int
end

const _GAUNT_CACHE = Dict{Int,Dict{Tuple{Int,Int},Vector{Pair{Int,Float64}}}}()
const _LM2_CACHE = Dict{Int,Vector{_LM}}()

function _lm2_list(L::Int)
    get!(() -> [_LM(l, m) for l = 0:(2 * L) for m = -l:l], _LM2_CACHE, L)
end

function _gaunt_table(L::Int)
    cached = get(_GAUNT_CACHE, L, nothing)
    cached !== nothing && return cached
    lms2 = _lm2_list(L)
    Nθ = 2 * L + 4
    # GL 节点（无 FastGaussQuadrature 依赖：node/weight 用对称性构造——
    # 这里直接调用 SpecialFunctions 的精确实现）
    xf, wf = _gauss_legendre(Nθ)
    Θf = Dict{Tuple{Int,Int},Vector{Float64}}()
    for lm in lms2
        Θf[(lm.l, lm.m)] = [theta_part(lm.l, lm.m, x) for x in xf]
    end
    cidx_of = Dict{_LM,Int}()
    for (ci, lm) in enumerate(lms2)
        cidx_of[lm] = ci
    end
    G = Dict{Tuple{Int,Int},Vector{Pair{Int,Float64}}}()
    src = [_LM(l, m) for l = 0:L for m = -l:l]
    Θ2 = zeros(Float64, Nθ)
    for (ai, A) in enumerate(src), (bi, B) in enumerate(src)
        mc = A.m - B.m
        entries = Pair{Int,Float64}[]
        Θa = Θf[(A.l, A.m)]
        Θb = Θf[(B.l, B.m)]
        for lc = max(abs(A.l - B.l), abs(mc)):min(A.l + B.l, L)
            abs(mc) > lc && continue
            Θc = Θf[(lc, mc)]
            @inbounds @simd for i = 1:Nθ
                Θ2[i] = Θb[i] * Θc[i]
            end
            s = 0.0
            @inbounds @simd for i = 1:Nθ
                s += wf[i] * Θa[i] * Θ2[i]
            end
            g = 2π * s
            abs(g) > 1e-14 && push!(entries, cidx_of[_LM(lc, mc)] => g)
        end
        isempty(entries) || (G[(ai, bi)] = entries)
    end
    _GAUNT_CACHE[L] = G
    return G
end

# Gauss–Legendre 节点/权重（Newton 迭代求根，对称配对；与
# FastGaussQuadrature 在 N ≤ 64 一致到 1e-15 量级，此处自包含实现）。
function _gauss_legendre(n::Int)
    x = zeros(n)
    w = zeros(n)
    m = (n + 1) >> 1
    for i = 1:m
        # 初值（Chebyshev 逼近）
        z = cos(π * (i - 0.25) / (n + 0.5))
        for _ = 1:100
            # P_n, P_n' 向下递推
            p0 = 1.0
            p1 = z
            for k = 2:n
                p0, p1 = p1, ((2k - 1) * z * p1 - (k - 1) * p0) / k
            end
            dp = n * (z * p1 - p0) / (z * z - 1.0)
            dz = p1 / dp
            z -= dz
            abs(dz) < 1e-16 && break
        end
        p0 = 1.0
        p1 = z
        for k = 2:n
            p0, p1 = p1, ((2k - 1) * z * p1 - (k - 1) * p0) / k
        end
        dp = n * (z * p1 - p0) / (z * z - 1.0)
        x[i] = -z
        x[n+1-i] = z
        w[i] = 2.0 / ((1.0 - z * z) * dp * dp)
        w[n+1-i] = w[i]
    end
    return x, w
end

"""
    m2l_matrix(L, Rvec, k; Leff = L) -> M（(L+1)² × (L+1)²，稠密）

谱域（多极子域）M2L 平移矩阵（行 = 接收系数 a，列 = 源系数 b）：

    M_ab(R) = 4π Σ_{c ≤ L} (-j)^{l_c} h_{l_c}⁽²⁾(kR) conj(Y_c(R̂)) G[a,b,c]

即 `⟨T_{≤L}·Y_b, Y_a⟩`（T 为 L 截断实谱核；不含 Green 常数 `-jk/(4π)`）。
c 求和与核同步截在 L（见 `_gaunt_table` 注释：Gaunt 完备范围 c ≤ l_a+l_b
中 (L, l_a+l_b] 的项在 kR < L 时是 h_l 发散区伪耦合，对角的 L 截断核不含
它们）。

`Leff` 为**行列对称**的谱截断（`spectral_bandwidth`）：`l_a, l_b > Leff`
的行与列均置零。

> ⚠ **分析工具，勿用于管线**（issue #22 问题 3 阶段 3 负结果）：
> 本矩阵是"精确 Gaunt 配对"对象，在理想化原型（带限 F、同网格综合）
> 中与对角方案一致；但真实管线中对角方案的逐点采样求积依赖 GL 采样
> 对截断级数尾部的**广义求和/混叠相消**（kR < L 放大区，实际叶层常态），
> 带限的分析→M→综合往返摧毁该相消（实测谱域 rel 1e2–1e3 vs 对角 3.6e-3）。
> 机制与数据见 `docs/dev/m2l_short_range_spectral.md` §6。
"""
function m2l_matrix(L::Int, Rvec, k::Real; Leff::Int = L)
    G = _gaunt_table(L)
    lms2 = _lm2_list(L)
    Rv = collect(float.(Rvec))
    R = norm(Rv)
    R̂ = Rv / R
    kR = Float64(k) * R
    θR = acos(clamp(R̂[3], -1.0, 1.0))
    φR = atan(R̂[2], R̂[1])
    src = [_LM(l, m) for l = 0:L for m = -l:l]
    # Y_c(R̂) 与 h_{l_c}(kR)，c 度 ≤ 2L
    nc2 = (2 * L + 1)^2
    Yc = Vector{ComplexF64}(undef, nc2)
    hcof = Dict{Int,ComplexF64}()
    ci = 0
    for l = 0:(2 * L), m = -l:l
        ci += 1
        Yc[ci] = conj(ylm(l, m, θR, φR))
        haskey(hcof, l) || (hcof[l] = _spherical_h2(l, kR))
    end
    Leff = clamp(Leff, 0, L)
    Mb = zeros(ComplexF64, ncoeff(L), ncoeff(L))
    for (bi, B) in enumerate(src)
        B.l > Leff && break   # 列按 l 分块连续排列，直接终止
        for (ai, A) in enumerate(src)
            A.l > Leff && continue   # 行对称截断
            ent = get(G, (ai, bi), nothing)
            ent === nothing && continue
            acc = 0.0 + 0.0im
            for (cidx, g) in ent
                lm = lms2[cidx]
                acc += ((-im)^lm.l) * hcof[lm.l] * Yc[cidx] * g
            end
            Mb[ai, bi] = 4π * acc
        end
    end
    return Mb
end

end
