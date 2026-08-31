# verify_m2l_spectral_shortrange.jl
#
# 阶段 2 独立原型验证（issue #22 问题 3）：谱域（多极子域）M2L 平移方案
# 在短距（nr=1 交互列表的最短距离区间 kR ∈ [2.5, 6]）的形状误差。
#
# ⚠ 阶段 3 集成结论（负结果，2025）：本原型的"通过"只在理想化设定
# （带限 F、同一 GL 网格上的合成=求积、无真实管线噪声）中成立，未迁移
# 到真实 MLFMA 管线。集成层实测（GD2X r=0.4/36×50/叶 0.1 m/nr=2，
# PMCHW 300 MHz εr=4）：对角 rel=3.6e-3 ✓，谱域 rel ≥ 3.6e2 ✗——
# 对角方案的逐点采样求积在 kR < L 放大区依赖 GL 采样对截断级数尾部的
# 广义求和/混叠相消才收敛于真值；带限的分析→平移→综合往返破坏该相消。
# 本脚本保留作为数学工具与理想化设定下结论的自洽性证据；机制分析见
# docs/dev/m2l_short_range_spectral.md §6。
#
# 设计依据：docs/dev/m2l_short_range_spectral.md
# 模式：单体脚本 —— 随机点源/接收点 + 精确 Green 对照，不加载 EMMoMSuite。
#
# 自检链（任一失败立即中止）：
#   1) 归一化球谐在 GL(θ)×均匀(φ) 极点网格上的正交归一性（degree ≤ L）
#   2) Gaunt 系数 vs 稠密 2D 数值积分（抽样）+ 选择律 + c=(0,0) 恒等
#   3) M 矩阵单谐解析锚点 M_{a,(0,0)} = -jk(-j)^{l_a} h_{l_a}(kR) conj(Y_a(R̂))/√(4π)
#   4) 大 kR 一致性：谱域 ≈ 对角 ≈ 精确（同 L）
# 主检验：kR ∈ {2.51, 3.55, 4.36}（偏移 (2,0,0)/(2,2,0)/(2,2,2)·w）、
#   L ∈ {9, 12, 13}：对角 T_L 复现失效，谱域形状误差 < 1e-2（争取 < 1e-3）。
#
# Usage: julia --project=scripts/verification verify_m2l_spectral_shortrange.jl
#        （或任意可解析 SpecialFunctions/FastGaussQuadrature 的环境）

using LinearAlgebra
using SpecialFunctions: sphericalbesselj, sphericalbessely, loggamma
using FastGaussQuadrature: gausslegendre
using Printf
using Random

const FAIL = Ref(false)
macro pass(name, cond)
    return quote
        ok = $(esc(cond))
        println((ok ? "  [PASS] " : "  [FAIL] ") * $(esc(name)))
        ok || (FAIL[] = true)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 归一化球谐（复谐、Condon–Shortley 相位）
#   Y_lm(θ,φ) = Θ_lm(cosθ) · exp(i m φ)，Θ 为实值归一化连带 Legendre
#   m ≥ 0: Θ = N_lm P_l^m；m < 0: Y_lm = (-1)^{|m|} conj(Y_{l,|m|})
# ─────────────────────────────────────────────────────────────────────────────

"""连带 Legendre P_l^m(x)（未归一化，Condon–Shortley 相位，m ≥ 0）。Numerical
Recipes 风格稳定递推，l ≤ ~40 精度充足。"""
function plgndr(l::Int, m::Int, x::Float64)
    if m < 0 || m > l
        return 0.0
    end
    pmm = 1.0
    if m > 0
        somx2 = sqrt(max(0.0, (1.0 - x) * (1.0 + x)))
        fact = 1.0
        for i = 1:m
            pmm *= -fact * somx2
            fact += 2.0
        end
    end
    if l == m
        return pmm
    end
    pmmp1 = x * (2m + 1) * pmm
    if l == m + 1
        return pmmp1
    end
    pll = 0.0
    for ll = m+2:l
        pll = ((2ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m)
        pmm = pmmp1
        pmmp1 = pll
    end
    return pll
end

norm_factor(l::Int, m::Int) = sqrt((2l + 1) / (4π)) *
    exp(0.5 * (loggamma(l - abs(m) + 1) - loggamma(l + abs(m) + 1)))

"""Θ_lm(cosθ)：Y_lm = Θ_lm·e^{imφ} 的实值 θ 部分。"""
function theta_part(l::Int, m::Int, x::Float64)
    if m >= 0
        return norm_factor(l, m) * plgndr(l, m, x)
    else
        # Y_{l,-|m|} = (-1)^{|m|} conj(Y_{l,|m|}) ⇒ Θ = (-1)^{|m|} N P_l^{|m|}
        return (-1.0)^(-m) * norm_factor(l, m) * plgndr(l, -m, x)
    end
end

ylm(l::Int, m::Int, θ::Float64, φ::Float64) =
    theta_part(l, m, cos(θ)) * exp(im * m * φ)

spherical_h2(l::Int, x::Float64) =
    complex(sphericalbesselj(l, x), -sphericalbessely(l, x))

# ─────────────────────────────────────────────────────────────────────────────
# 极点网格（复刻 Interpolation.levelIntegralInfoCal：GL(θ) × 均匀半格(φ)）
# ─────────────────────────────────────────────────────────────────────────────

struct Poles
    L::Int
    θs::Vector{Float64}
    φs::Vector{Float64}
    Wθ::Vector{Float64}
    Wφ::Vector{Float64}
    n::Int                     # 总极点数 2(L+1)^2
    # 展平顺序：外层 ϕ、内层 θ（与 levelIntegralInfoCal 一致）
    θp::Vector{Float64}
    φp::Vector{Float64}
    Wp::Vector{Float64}
    k̂::Matrix{Float64}         # 3 × n
end

function Poles(L::Int)
    x, w = gausslegendre(L + 1)          # cosθ 升序；权重和 = 2
    Xcosθ = -x                            # 复刻 lb=1, hb=-1 的倒序映射（顺序无关紧要）
    Wθ = w
    Xθs = acos.(Xcosθ)
    Mφ = 2 * (L + 1)
    Xϕs = [(j - 0.5) * 2π / Mφ for j = 1:Mφ]
    Wϕ = fill(2π / Mφ, Mφ)
    n = length(Xθs) * length(Xϕs)
    θp = Vector{Float64}(undef, n)
    φp = Vector{Float64}(undef, n)
    Wp = Vector{Float64}(undef, n)
    k̂ = Matrix{Float64}(undef, 3, n)
    i = 0
    for ϕ in Xϕs, θ in Xθs
        i += 1
        θp[i] = θ
        φp[i] = ϕ
        Wp[i] = Wθs(Xθs, θ, Wθ) * 2π / Mφ
        k̂[1, i] = sin(θ) * cos(ϕ)
        k̂[2, i] = sin(θ) * sin(ϕ)
        k̂[3, i] = cos(θ)
    end
    return Poles(L, Xθs, Xϕs, Wθ, Wϕ, n, θp, φp, Wp, k̂)
end

Wθs(Xθs, θ, Wθ) = Wθ[findmin(abs.(Xθs .- θ))[2]]  # 仅构造期使用

# 系数索引 a = (l, m)：l = 0:L, m = -l:l
ncoeff(L::Int) = (L + 1)^2
idx(l::Int, m::Int) = l^2 + (m + l + 1)
struct LM
    l::Int
    m::Int
end
lm_list(L::Int) = [LM(l, m) for l = 0:L for m = -l:l]

"""极点网格上的球谐矩阵 Ymat[p, a]（nPoles × (L+1)²）。"""
function ymat(p::Poles)
    L = p.L
    lms = lm_list(L)
    Y = Matrix{ComplexF64}(undef, p.n, ncoeff(L))
    for (a, lm) in enumerate(lms)
        for ip = 1:p.n
            Y[ip, a] = ylm(lm.l, lm.m, p.θp[ip], p.φp[ip])
        end
    end
    return Y
end

# ─────────────────────────────────────────────────────────────────────────────
# Gaunt 系数（因子化 θ 求积；φ 由选择律解析处理）
#   G[a,b,c] = ∫ conj(Y_a) Y_b Y_c dΩ = 2π δ_{m_c, m_a−m_b} Σ_i w_i Θ_a Θ_b Θ_c
#   （细 GL 网格 Nθ = 2L+4 对 cosθ 度 ≤ 4L 精确）
# 返回 Dict{(Int,Int)} — (a,b) => [(cidx, val), ...]
# ─────────────────────────────────────────────────────────────────────────────

function gaunt_table(L::Int)
    lms = lm_list(2 * L)
    Nθ = 2 * L + 4
    xf, wf = gausslegendre(Nθ)
    Θf = Dict{Tuple{Int,Int},Vector{Float64}}()
    for lm in lms
        Θf[(lm.l, lm.m)] = [theta_part(lm.l, lm.m, x) for x in xf]
    end
    G = Dict{Tuple{Int,Int},Vector{Pair{Int,Float64}}}()
    src = lm_list(L)
    cidx_of = Dict{LM,Int}()
    for (ci, lm) in enumerate(lm_list(2 * L))
        cidx_of[lm] = ci
    end
    Θ2 = zeros(Float64, Nθ)
    Θ3 = zeros(Float64, Nθ)
    for (ai, A) in enumerate(src), (bi, B) in enumerate(src)
        mc = A.m - B.m
        entries = Pair{Int,Float64}[]
        Θa = Θf[(A.l, A.m)]
        Θb = Θf[(B.l, B.m)]
        for lc = max(abs(A.l - B.l), abs(mc)):min(A.l + B.l, 2 * L)
            if abs(mc) > lc
                continue
            end
            C = LM(lc, mc)
            Θc = Θf[(lc, mc)]
            @inbounds @simd for i = 1:Nθ
                Θ2[i] = Θb[i] * Θc[i]
            end
            s = 0.0
            @inbounds @simd for i = 1:Nθ
                s += wf[i] * Θa[i] * Θ2[i]
            end
            g = 2π * s
            abs(g) > 1e-14 && push!(entries, cidx_of[C] => g)
        end
        isempty(entries) || (G[(ai, bi)] = entries)
    end
    return G, cidx_of, lm_list(2 * L), Θf, wf
end

# ─────────────────────────────────────────────────────────────────────────────
# M2L 谱域矩阵与两条远场路径
# ─────────────────────────────────────────────────────────────────────────────

const _GAUNT_CACHE = Dict{Int,Tuple}()
function _gaunt_cached(L::Int)
    get!(() -> gaunt_table(L), _GAUNT_CACHE, L)
end

"""
    m2l_matrix(L, Rvec, k; l2max) -> M[(L+1)² × (L+1)²]（稠密，行=a、列=b）

M_ab(R) = 4π Σ_c (-j)^{l_c} h_{l_c}⁽²⁾(kR) conj(Y_c(R̂)) G[a,b,c]

即**级数形式**（不含 Green 常数 -jk/16π²）平移算子与限带方向图的配对矩阵
⟨T_series·Y_b, Y_a⟩：对角方案逐点采样求积计算的是同一对象
（CG 常数在 farfield_case 中两条路径共享）。
输出行索引按 a ≤ (L+1)²（degree ≤ L 的接收侧投影）；c 的 degree 上限
l2max（默认 2L，即 M_ab 的完整 Gaunt 截断）。列 b 索引为 degree ≤ L 的源系数。
"""
function m2l_matrix(L::Int, Rvec::Vector{Float64}, k::Float64; l2max::Int = 2 * L)
    l2max = min(l2max, 2 * L)
    G, cidx_of, lms_all, Θf, wf = _gaunt_cached(L)
    R = norm(Rvec)
    R̂ = Rvec / R
    kR = k * R
    θR = acos(clamp(R̂[3], -1.0, 1.0))
    φR = atan(R̂[2], R̂[1])
    src = lm_list(L)
    Yc = Vector{ComplexF64}(undef, ncoeff(2 * L))
    hcof = Dict{Int,ComplexF64}()
    for lm in lm_list(2 * L)
        i = cidx_of[lm]
        Yc[i] = conj(ylm(lm.l, lm.m, θR, φR))
        hcof[lm.l] = spherical_h2(lm.l, kR)
    end
    Mb = zeros(ComplexF64, ncoeff(L), ncoeff(L))
    for (ai, A) in enumerate(src), (bi, B) in enumerate(src)
        ent = get(G, (ai, bi), nothing)
        ent === nothing && continue
        acc = 0.0 + 0.0im
        for (ci, g) in ent
            lm = lms_all[ci]
            acc += ((-im)^lm.l) * hcof[lm.l] * Yc[ci] * g
        end
        Mb[ai, bi] = 4π * acc
    end
    return Mb
end

"""对角平移函数 T_L(k̂_p, R)（复刻 Translation.jl 的级数，不含权重与 -jk/(4π)）。"""
function diagonal_T(L::Int, p::Poles, Rvec::Vector{Float64}, k::Float64)
    R = norm(Rvec)
    R̂ = Rvec / R
    kR = k * R
    h = [spherical_h2(l, kR) for l = 0:L]
    T = Vector{ComplexF64}(undef, p.n)
    for ip = 1:p.n
        cϕ = clamp(dot(view(p.k̂, :, ip), R̂), -1.0, 1.0)
        # Legendre 递推
        Plm1 = 1.0          # P_0
        Pl = cϕ             # P_1（L=0 时不使用）
        acc = 0.0 + 0.0im
        jt = im
        for l = 0:L
            jt *= -im
            P = l == 0 ? 1.0 : (l == 1 ? cϕ : Pl)
            acc += jt * (2l + 1) * h[l+1] * P
            # 递推推进
            if l >= 1
                Plm1, Pl = Pl, ((2l + 1) * cϕ * Pl - l * Plm1) / (l + 1)
            end
        end
        T[ip] = acc
    end
    return T
end

# ─────────────────────────────────────────────────────────────────────────────
# 自检
# ─────────────────────────────────────────────────────────────────────────────

function selfcheck(L::Int)
    println("— 自检（L = $L）—")
    p = Poles(L)
    Y = ymat(p)
    # 1) 正交归一性（网格精确性 ⇒ 机器精度）
    S = Y' * (p.Wp .* Y)
    ortho_err = maximum(abs, S - I)
    @pass "正交归一 Y^H W Y = I（max |·| = $(round(ortho_err, sigdigits=2)) < 1e-12）" ortho_err < 1e-12

    # 2) Gaunt：稠密 2D 数值积分抽样对照
    G, cidx_of, lms_all, Θf, wf = gaunt_table(L)
    # 稠密参照网格
    Nd = 3 * L + 8
    xd, wdd = gausslegendre(Nd)
    Mφd = 6 * L + 16
    φd = [(j - 0.5) * 2π / Mφd for j = 1:Mφd]
    wφd = 2π / Mφd
    rng = MersenneTwister(7)
    worst = 0.0
    for _ = 1:40
        A = LM(rand(rng, 0:L), 0); A = LM(A.l, rand(rng, -A.l:A.l))
        B = LM(rand(rng, 0:L), 0); B = LM(B.l, rand(rng, -B.l:B.l))
        mc = A.m - B.m
        lo = max(abs(A.l - B.l), abs(mc))
        hi = min(A.l + B.l, 2 * L)
        lo > hi && continue
        C = LM(rand(rng, lo:hi), mc)
        ref = 0.0 + 0.0im
        for (iφ, φ) in enumerate(φd)
            for i = 1:Nd
                ref +=
                    wdd[i] * wφd *
                    conj(ylm(A.l, A.m, acos(clamp(xd[i], -1, 1)), φ)) *
                    ylm(B.l, B.m, acos(clamp(xd[i], -1, 1)), φ) *
                    ylm(C.l, C.m, acos(clamp(xd[i], -1, 1)), φ)
            end
        end
        # 表值
        tbl = 0.0
        src = lm_list(L)
        ai = findfirst(x -> x == A, src)
        bi = findfirst(x -> x == B, src)
        if ai !== nothing && bi !== nothing
            for (ci, g) in get(G, (ai, bi), Pair{Int,Float64}[])
                lms_all[ci] == C && (tbl = g)
            end
        end
        # 宇称奇等选择律外的真零项：两侧都是噪声，跳过比值
        if abs(ref) < 1e-12
            abs(tbl) < 1e-12 || (worst = max(worst, 1.0))
            continue
        end
        worst = max(worst, abs(ref - tbl) / abs(ref))
    end
    @pass "Gaunt vs 稠密积分（抽样 20 组，worst rel = $(round(worst, sigdigits=2)) < 1e-8）" worst < 1e-8

    # 3) c=(0,0) 恒等：G[a,b,(0,0)] = δ_{ab}
    src = lm_list(L)
    c00 = cidx_of[LM(0, 0)]
    ok00 = true
    for (ai, A) in enumerate(src), (bi, B) in enumerate(src)
        v = 0.0
        for (ci, g) in get(G, (ai, bi), Pair{Int,Float64}[])
            ci == c00 && (v = g)
        end
        (ai == bi && abs(v - 1 / sqrt(4π)) > 1e-12) && (ok00 = false)
        (ai != bi && abs(v) > 1e-14) && (ok00 = false)
    end
    @pass "Gaunt c=(0,0) 恒等 = δ_{ab}/√(4π)" ok00

    # 4) 单谐解析锚点：M_{a,(0,0)}·√(4π) = 4π(-j)^{l_a} h_{l_a} conj(Y_a(R̂))
    k = 4π
    Rvec = [0.2, 0.35, 0.91]
    Mb = m2l_matrix(L, Rvec, k)
    R = norm(Rvec); R̂ = Rvec / R; kR = k * R
    θR = acos(clamp(R̂[3], -1, 1)); φR = atan(R̂[2], R̂[1])
    a00 = idx(0, 0)
    anchor_err = 0.0
    for (ai, A) in enumerate(src)
        lhs = Mb[ai, a00] * sqrt(4π)
        rhs = 4π * (-im)^A.l * spherical_h2(A.l, kR) * conj(ylm(A.l, A.m, θR, φR))
        anchor_err = max(anchor_err, abs(lhs - rhs) / max(abs(rhs), 1e-30))
    end
    @pass "M 单谐锚点（max rel = $(round(anchor_err, sigdigits=2)) < 1e-10）" anchor_err < 1e-10
    return p, Y
end

# ─────────────────────────────────────────────────────────────────────────────
# 物理检验：随机点源 → 随机接收点，精确 Green 对照
# ─────────────────────────────────────────────────────────────────────────────

function farfield_case(;
    L::Int,
    k::Float64,
    w::Float64,
    offset::Tuple{Int,Int,Int},
    nsrc::Int = 20,
    nrcv::Int = 30,
    seed::Int = 2024,
    CG::ComplexF64 = -im * k / (16π^2),   # 由大 kR 标定/验证
    F_noise::Float64 = 0.0,               # F 样本相对噪声（模拟层间插值高频残差）
    verbose::Bool = true,
)
    rng = MersenneTwister(seed)
    p = Poles(L)
    Y = ymat(p)

    # 源：随机点 + 随机复振幅，位于源盒（中心原点）
    src_pts = [(w / 2) .* (2 .* rand(rng, 3) .- 1) for _ = 1:nsrc]
    q = randn(rng, ComplexF64, nsrc)
    # 接收：目标盒（中心 offset·w）内随机点
    ct = w .* collect(offset)
    rcv_pts = [ct .+ (w / 2) .* (2 .* rand(rng, 3) .- 1) for _ = 1:nrcv]

    # 聚合方向图（与 Aggregation.jl 约定一致：e^{+jk k̂·(r'−c_s)}）
    F = zeros(ComplexF64, p.n)
    for ip = 1:p.n
        acc = 0.0 + 0.0im
        k̂ = view(p.k̂, :, ip)
        for n = 1:nsrc
            acc += q[n] * exp(im * k * dot(k̂, src_pts[n]))
        end
        F[ip] = acc
    end

    # 精确场
    y_ex = Vector{ComplexF64}(undef, nrcv)
    for (ir, r) in enumerate(rcv_pts)
        acc = 0.0 + 0.0im
        for n = 1:nsrc
            ρ = r .- src_pts[n]
            acc += q[n] * exp(-im * k * norm(ρ)) / (4π * norm(ρ))
        end
        y_ex[ir] = acc
    end

    # —— 路径 A：对角 T_L 逐点（复刻 Translation.jl）——
    Rvec = Float64.([offset[1], offset[2], offset[3]]) .* w
    R = norm(Rvec); kR = k * R
    T = diagonal_T(L, p, Rvec, k)
    # 非理想谱成分（模拟层间插值高频残差 / 极化坐标伪影）：只污染 F，
    # 精确场保持无噪 —— 误差全部来自平移算子对污染的放大
    Fn = F_noise > 0 ? F .+ F_noise * norm(F) .* (randn(rng, ComplexF64, p.n) ./ sqrt(p.n)) : F
    y_diag = zeros(ComplexF64, nrcv)
    for (ir, r) in enumerate(rcv_pts)
        acc = 0.0 + 0.0im
        for ip = 1:p.n
            phase = exp(-im * k * dot(view(p.k̂, :, ip), r .- ct))
            acc += p.Wp[ip] * T[ip] * Fn[ip] * phase
        end
        y_diag[ir] = CG * acc
    end

    # —— 路径 B：谱域 M2L（SHT → 稠密 M → 综合）——
    Mb = m2l_matrix(L, Rvec, k)
    F̂ = Y' * (p.Wp .* Fn)
    Ĝ = Mb * F̂
    G = Y * Ĝ
    y_spec = zeros(ComplexF64, nrcv)
    for (ir, r) in enumerate(rcv_pts)
        acc = 0.0 + 0.0im
        for ip = 1:p.n
            phase = exp(-im * k * dot(view(p.k̂, :, ip), r .- ct))
            acc += p.Wp[ip] * G[ip] * phase
        end
        y_spec[ir] = CG * acc
    end

    rel(v) = norm(v .- y_ex) / norm(y_ex)
    corr(v) = abs(adjoint(v) * y_ex) / (norm(v) * norm(y_ex) + 1e-300)
    verbose && @printf(
        "    L=%2d  kR=%5.2f  off=(%d,%d,%d)   rel_diag=%.3e (corr %.4f)   rel_spec=%.3e (corr %.6f)\n",
        L, kR, offset..., rel(y_diag), corr(y_diag), rel(y_spec), corr(y_spec)
    )
    return (kR = kR, rel_diag = rel(y_diag), rel_spec = rel(y_spec),
            corr_diag = corr(y_diag), corr_spec = corr(y_spec),
            maxT = maximum(abs.(T)))
end

function main()
    println("=" ^ 78)
    println(" 谱域（多极子域）短距 M2L 原型验证 — issue #22 问题 3 阶段 2")
    println("=" ^ 78)

    # —— 自检 ——
    for L in (9, 12)
        selfcheck(L)
    end
    FAIL[] && error("自检未通过，中止")

    # —— 常数标定：大 kR 处对角路径应精确，验证 CG = -jk/(16π²) ——
    println("— 常数标定（大 kR 一致性）—")
    k = 4π                     # GD2V fixture: λ_medium = 0.5, w = 0.1 → k·w = 1.2566
    w = 0.1
    for (off, L) in (((16, 0, 0), 12), ((8, 0, 0), 12))
        r = farfield_case(L = L, k = k, w = w, offset = off, verbose = false)
        @printf("    off=(%d,%d,%d) kR=%6.2f  rel_diag=%.2e  rel_spec=%.2e\n",
                off..., r.kR, r.rel_diag, r.rel_spec)
        @pass "大 kR 一致性 off=$(off)：rel_diag < 1e-2 且 rel_spec < 1e-2（谱域≈对角）" (
            r.rel_diag < 1e-2 && r.rel_spec < 1e-2 &&
            abs(r.rel_spec - r.rel_diag) < 1e-2
        )
    end

    # —— 主检验 1：nr=1 交互列表最短距离区间（理想 F）——
    println("— 主检验 1：kR ∈ [2.5, 6] 理想配对（GD2V fixture 叶盒 nr=1 偏移族）—")
    cases = [
        ((2, 0, 0), 9), ((2, 0, 0), 12), ((2, 0, 0), 13),
        ((2, 1, 0), 12),
        ((2, 2, 0), 9), ((2, 2, 0), 12), ((2, 2, 0), 13),
        ((2, 2, 2), 12),
        ((3, 0, 0), 12),
    ]
    worst_spec = 0.0
    worst_diag = 0.0
    for (off, L) in cases
        r = farfield_case(L = L, k = k, w = w, offset = off)
        worst_spec = max(worst_spec, r.rel_spec)
        worst_diag = max(worst_diag, r.rel_diag)
    end
    println()
    @pass "谱域形状误差 < 1e-2（worst rel_spec = $(round(worst_spec, sigdigits=3))）" worst_spec < 1e-2
    @pass "谱域形状误差 < 1e-3（争取目标）" worst_spec < 1e-3

    # —— 主检验 2：非理想谱成分放大（真实失效机理复现）——
    # 净室理想配对下对角级数并不立即失效（高 l 项与理想 F 的高频近零系数
    # 精确相消）；真实管线中 F 携带层间插值高频残差（实测 ~7e-4）与极化
    # 坐标伪影，被巨大的部分和样本（max|T_L| ~ 1e5–1e6）放大为 O(1)+ 误差。
    # 此处以相对白噪污染 F 复现该机理，谱域路径的算子有界、无放大。
    println("— 主检验 2：F 谱污染放大（对角复现失效，谱域免疫）—")
    println("    [噪 F: rel 1e-3 —— 模拟 GD2V 实测插值残差量级]")
    worst_diag_n = 0.0
    worst_spec_n = 0.0
    for (off, L, wcase, tag) in (
        ((2, 0, 0), 12, 0.1, "叶层 L=12"),
        ((2, 0, 0), 17, 0.2, "第3层 L=17"),
    )
        r = farfield_case(L = L, k = k, w = wcase, offset = off, F_noise = 1e-3)
        worst_diag_n = max(worst_diag_n, r.rel_diag)
        worst_spec_n = max(worst_spec_n, r.rel_spec)
        @printf("      %s  max|T_L|=%.2e\n", tag, r.maxT)
    end
    @pass "对角被噪 F 放大失效（worst rel_diag = $(round(worst_diag_n, sigdigits=3)) > 0.5）" worst_diag_n > 0.5
    @pass "谱域对噪 F 免疫（worst rel_spec = $(round(worst_spec_n, sigdigits=3)) < 1e-2）" worst_spec_n < 1e-2

    println("=" ^ 78)
    if FAIL[]
        println(" 结果：存在失败项 — 按约定停止，不合入主代码")
    else
        println(" 结果：全部通过 — 谱域方案可进入阶段 3（集成）")
    end
    println("=" ^ 78)
    return FAIL[]
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main() ? 1 : 0)
end
