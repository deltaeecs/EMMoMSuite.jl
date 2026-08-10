module dataset_generator

using LinearAlgebra
using ...MLFMA.Interpolation: truncation_kernel
using ..LebedevSortedPoints: get_t_nodes, nodes2Poles, p2nDict
using ....Utilities: Progress, next!, find_zero_bisection

export generate_dataset_on_pkpt, random_source_geometry, evaluate_poles!

const ws = [-4 / 5 9 / 20 9 / 20 9 / 20 9 / 20]

function rand_truncated_normal(mu, sigma, min_val, max_val)
    while true
        x = mu + sigma * randn()
        if min_val <= x <= max_val
            return x
        end
    end
end

"""
在球面生成随机向量
"""
function random_rvec(bd = 1; FT = Float64)
    [(2rand() - 1) * bd, (2rand() - 1) * bd, (2rand() - 1) * bd]
end

function random_rhat()
    z = 2rand() - 1
    phi = 2π * rand()
    r = sqrt(1 - z^2)
    x = r * cos(phi)
    y = r * sin(phi)
    return [x, y, z]
end

"""
    random_source_geometry(rvec; arm_max, off_max)

生成一个随机辐射源的几何（双端点 rvecp/rvecm、求积偏移、参考点 r0p/r0m）。
同一个几何必须同时用于粗层与细层采样；否则两层数据对应不同辐射函数，
插值矩阵将拟合"随机粗图案 -> 随机细图案"的噪声——这是原实现权重矩阵
全部失效的根本原因。
"""
function random_source_geometry(rvec; arm_max = 0.12, off_max = 0.03)
    rvecp = rvec
    rvecm = rvec .+ random_rhat() .* (rand_truncated_normal(2 / 3 * arm_max, arm_max / 3, 0, arm_max))
    offsets = [
        [random_rhat()...] .* (rand_truncated_normal(2 / 3 * off_max, off_max / 3, off_max / 3, off_max)) for
        _ in eachindex(ws)
    ]
    offsets[1] .= 0
    r0p = rvecp .+ rvecp .- rvecm .+ random_rhat() .* (arm_max / 24)
    r0m = rvecm .+ rvecm .- rvecp .+ random_rhat() .* (arm_max / 24)
    return (rvecp = rvecp, rvecm = rvecm, offsets = offsets, r0p = r0p, r0m = r0m)
end

"""
    evaluate_poles!(rHatsθsϕs, tArray, geom; k)

用给定几何（同一辐射源）在采样点集 rHatsθsϕs 上求 θ/ϕ 分量并写入 tArray。

对应论文式 (4-21) 的通用辐射函数：

```math
\\bm{\\mathcal{F}}(\\hat{k}) = C_f \\left(\\overline{I} - \\hat{k}\\hat{k}\\right) \\cdot
\\bm{\\hat{\\rho}}\\, e^{-{\\rm j}k \\hat{k} \\cdot (r_b - r')}
```

实现用 RWG 型源对（`rvecp`/`rvecm`）与 5 点求积权重 `ws` 离散：
对每个采样点 `ĥ` 累加 `θ̂/φ̂ · ρ̂ e^{jk ĥ·r}` 的加权差（实部虚部分别写入
`tArray[iPole, 1]` 与 `[iPole, 2]`）。**同一个几何必须同时用于粗层与细层**
采样，否则两层数据对应不同辐射函数，拟合出的插值矩阵将退化为噪声
（原实现权重矩阵全部失效的根因）。
"""
function evaluate_poles!(rHatsθsϕs, tArray, geom; k = 2π, FT = Float64)
    # 常数
    JK_0 = im * k
    rvecp, rvecm = geom.rvecp, geom.rvecm
    offsets, r0p, r0m = geom.offsets, geom.r0p, geom.r0m
    rp = copy(r0p)
    rm = copy(r0m)
    ρhatp_iw = copy(r0p)
    ρhatm_iw = copy(r0m)
    for iPole in eachindex(rHatsθsϕs)
        # 该多极子
        poler̂θϕ = rHatsθsϕs[iPole]
        for iw in eachindex(ws)
            rp .= rvecp .+ offsets[iw]
            rm .= rvecm .+ offsets[iw]
            ρhatp_iw .= rp .- r0p
            ρhatm_iw .= rm .- r0m
            # 公用的 指数项
            wpexptemp = ws[iw] * exp(JK_0 * dot(poler̂θϕ.r̂, rp))
            wmexptemp = ws[iw] * exp(JK_0 * dot(poler̂θϕ.r̂, rm))
            # 将结果写入目标数组
            tArray[iPole, 1] += dot(poler̂θϕ.θhat, ρhatp_iw) * wpexptemp
            tArray[iPole, 1] -= dot(poler̂θϕ.θhat, ρhatm_iw) * wmexptemp
            tArray[iPole, 2] += dot(poler̂θϕ.ϕhat, ρhatp_iw) * wpexptemp
            tArray[iPole, 2] -= dot(poler̂θϕ.ϕhat, ρhatm_iw) * wmexptemp
        end
    end # iPole
    return tArray
end

"""
    generate_dataset_on_poles(rHatsθsϕs, tArray; rvec = random_rvec(), FT = Precision.FT)

    simulate agg on basis functions.
"""
function generate_dataset_on_poles(
    rHatsθsϕs,
    tArray;
    rvec = random_rvec(),
    k = 2π,
    λ = 1.0,
    arm_max = 0.12 * λ,
    off_max = 0.03 * λ,
    geom = nothing,
    FT = Float64,
)
    # 传入 geom 时复用同一辐射源；否则生成新的随机源
    g = geom === nothing ? random_source_geometry(rvec; arm_max = arm_max, off_max = off_max) : geom
    evaluate_poles!(rHatsθsϕs, tArray, g; k = k)
    return tArray
end

"""
    generate_dataset_on_poles(rHatsθsϕs; rvec = random_rvec(), FT = Precision.FT)

TBW
"""
function generate_dataset_on_poles(rHatsθsϕs; rvec = random_rvec(), k = 2π, λ = 1.0, FT = Float64)
    # 目标数组
    tArray = zeros(Complex{FT}, length(rHatsθsϕs), 2)
    generate_dataset_on_poles(
        rHatsθsϕs,
        reshape(tArray, length(rHatsθsϕs), 2);
        rvec = rvec,
        k = k,
        λ = λ,
        arm_max = 0.12 * λ,
        off_max = 0.03 * λ,
        FT = FT,
    )
    return tArray
end

"""
    generate_dataset_on_pkpt(pk, pt, rel_l; FT = Precision.FT)

在多项式阶数为 `pk`（细层）与 `pt`（粗层）的两套 Lebedev（或 Fibonacci）
采样点上生成辐射函数数据集，用于求解层间插值矩阵（论文 4.3.1 节）。

对 `Nρ` 个随机 `ρ̂` 与 `Nr` 个随机盒内位置 `r_b − r'` 的组合，在细层/粗层
采样点上分别计算 `F(ĥ)` 的 θ/φ 分量，实部与虚部作为独立数据列（论文：
"实部虚部拆分后数据集翻倍"）。返回数组尺寸
`(N_p^{细}, 2, Nρ, Nr)`，后续由 `pinv2interpW` 按 80/20 划分训练/测试集。

# Arguments
- `pk`, `pt`: 细层、粗层多项式阶数（`p = 2τ+1`，奇数）。
- `rel_l`: 层盒子边长与波长之比（决定截断项，见
  [`Interpolation.truncation_kernel`](@ref)）。
- `k`: 波数（默认 `2π`，与 λ=1 匹配）。
- `λ`: 波长（默认 1）。
"""
function generate_dataset_on_pkpt(
    pk::T,
    pt::T,
    rel_l = find_zero_bisection(x -> truncation_kernel(x) - (pk + 1) ÷ 2, 0);
    k = 2π,
    λ = 1.0,
    FT = Float64,
) where {T<:Integer}
    # trunc
    τt = (pk - 1) ÷ 2
    τp = (pt - 1) ÷ 2

    # 多项式阶数
    pt = 2τt + 1
    # 高阶（>131）无 Lebedev 数据集：get_t_nodes 自动回退 Fibonacci 格点，无需报错

    # 生成基函数矢量
    ρhats = zeros(FT, 3, 50)
    for i in axes(ρhats, 2)
        ρhats[:, i] = random_rhat()
    end
    # 空间位置矢量：源对整体必须落在盒内（|rvec| + arm <= rel_l*λ/2 每坐标），
    # 否则数据带宽超出本层可表示阶数（arm 取物理值 0.12λ 与盒内余量的较小者）
    rbmrps = zeros(FT, 3, 500)
    arm_max = min(0.12 * λ, 0.5 * rel_l * λ)
    off_max = 0.125 * arm_max
    rscale = max(rel_l * λ / 2 - arm_max, 1e-4 * λ)
    @info "box size" (rel_l * λ / 2) "source radius" rscale
    for i in axes(rbmrps, 2)
        rbmrps[:, i] .= random_rvec() .* rscale
    end

    # nodes
    tnodes = get_t_nodes(τt; FT = FT)
    pnodes = get_t_nodes(τp; FT = FT)
    # rHatsθsϕs
    tr̂sθsϕs = nodes2Poles(tnodes)
    pr̂sθsϕs = nodes2Poles(pnodes)

    # 预分配内存
    tArray = zeros(Complex{FT}, length(tr̂sθsϕs), 2, size(ρhats, 2), size(rbmrps, 2))
    pArray = zeros(Complex{FT}, length(pr̂sθsϕs), 2, size(ρhats, 2), size(rbmrps, 2))

    # 开始计算
    # 关键：每个样本的源几何只生成一次，粗层与细层共用同一辐射函数；
    # 源几何按层盒子尺寸缩放，保证数据带宽 <= 本层可表示阶数。
    pmeter = Progress(size(rbmrps, 2), "计算数据集中…")
    for ir in axes(rbmrps, 2)#@threads 
        for iρ in axes(ρhats, 2)
            geom = random_source_geometry(
                rbmrps[:, ir];
                arm_max = arm_max,
                off_max = off_max,
            )
            @views generate_dataset_on_poles(
                tr̂sθsϕs,
                tArray[:, :, iρ, ir];
                k = k,
                FT = FT,
                geom = geom,
            )
            @views generate_dataset_on_poles(
                pr̂sθsϕs,
                pArray[:, :, iρ, ir];
                k = k,
                FT = FT,
                geom = geom,
            )
        end
        next!(pmeter)
    end

    return reshape(tArray, length(tr̂sθsϕs) * 2, :), reshape(pArray, length(pr̂sθsϕs) * 2, :)

end

end # module
