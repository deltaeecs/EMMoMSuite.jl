module dataset_generator

using LinearAlgebra
using ProgressMeter
using Roots
using ...MLFMA.Interpolation: truncation_kernel
using ..LebedevSortedPoints: get_t_nodes, nodes2Poles, p2nDict

export generate_dataset_on_pkpt

const ws = [-4/5 9/20 9/20 9/20 9/20]

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
function random_rvec(bd=1; FT = Float64)
    [(2rand()-1)*bd, (2rand()-1)*bd, (2rand()-1)*bd]
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
    generate_dataset_on_poles(rHatsθsϕs, tArray; rvec = random_rvec(), FT = Precision.FT)

    simulate agg on basis functions.
"""
function generate_dataset_on_poles(rHatsθsϕs, tArray; rvec = random_rvec(), k=2π, λ=1.0, FT = Float64)
    # 常数
    JK_0 = im * k
    # rvecp rmvec
    rvecp = rvec
    rvecm = rvec .+ random_rhat().*(rand_truncated_normal(0.08λ, 0.04λ, 0, 0.12λ))
    offsets = [[random_rhat()...].*(rand_truncated_normal(0.02λ, 0.01λ, 0.01λ, 0.03λ)) for _ in eachindex(ws)]
    offsets[1] .= 0

    r0p = rvecp .+ rvecp .- rvecm .+ random_rhat() .* 0.005λ
    r0m = rvecm .+ rvecm .- rvecp .+ random_rhat() .* 0.005λ

    rp  = copy(r0p)
    rm  = copy(r0m)

    ρhatp_iw = copy(r0p)
    ρhatm_iw = copy(r0m)
    
    for iPole in eachindex(rHatsθsϕs)
        # 该多极子
        poler̂θϕ =   rHatsθsϕs[iPole]
        for iw in eachindex(ws)

            rp .= rvecp .+ offsets[iw]
            rm .= rvecm .+ offsets[iw]

            ρhatp_iw .=  rp .- r0p
            ρhatm_iw .=  rm .- r0m
            
            # 公用的 指数项
            wpexptemp =   ws[iw]*exp(JK_0*dot(poler̂θϕ.r̂, rp))
            wmexptemp =   ws[iw]*exp(JK_0*dot(poler̂θϕ.r̂, rm))
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
    generate_dataset_on_poles(rHatsθsϕs; rvec = random_rvec(), FT = Precision.FT)
    
TBW
"""
function generate_dataset_on_poles(rHatsθsϕs; rvec = random_rvec(), k=2π, λ=1.0, FT = Float64)

    # 目标数组
    tArray = zeros(Complex{FT}, length(rHatsθsϕs), 2)
    generate_dataset_on_poles(rHatsθsϕs, reshape(tArray, length(rHatsθsϕs), 2); rvec = rvec, k=k, λ=λ, FT = FT)
    return tArray

end

"""
    generate_dataset_on_pkpt(pk, pt, rel_l; FT = Precision.FT)

TBW
"""
function generate_dataset_on_pkpt(pk::T, pt::T, rel_l = find_zero(x -> truncation_kernel(x) - (pk+1)÷2, 0); k=2π, λ=1.0, FT = Float64) where{T<:Integer}
    # trunc
    τt  =   (pk - 1) ÷ 2
    τp  =   (pt - 1) ÷ 2

    # 多项式阶数
    pt = 2τt+1
    # 若本层已超出Lebedev求积点取值范围则报错
    pt > maximum(keys(p2nDict)) && throw("多项式阶数已超出Lebedev求积点取值范围。")

    # 生成基函数矢量
    # ρhats, _ = getlbSortedData(13)
    ρhats = zeros(FT, 3, 50)
    for i in axes(ρhats, 2)
        ρhats[:, i] = random_rhat()
    end
    # 空间位置矢量
    # rbmrps = getlbSortedData(13)[1] .* (√3/2*rel_l*Params.λ_0)
    rbmrps = zeros(FT, 3, 500)
    @info "box size" (rel_l*λ/2)
    for i in axes(rbmrps, 2)
        rbmrps[:, i]   .= random_rvec() .* (rel_l*λ/2)
    end

    # nodes
    tnodes = get_t_nodes(τt; FT=FT)
    pnodes = get_t_nodes(τp; FT=FT)
    # rHatsθsϕs
    tr̂sθsϕs = nodes2Poles(tnodes)
    pr̂sθsϕs = nodes2Poles(pnodes)

    # 预分配内存
    tArray  = zeros(Complex{FT}, length(tr̂sθsϕs), 2, size(ρhats, 2), size(rbmrps, 2))
    pArray  = zeros(Complex{FT}, length(pr̂sθsϕs), 2, size(ρhats, 2), size(rbmrps, 2))

    # 开始计算
    pmeter =  Progress(size(rbmrps, 2), "计算数据集中…")
    for ir in axes(rbmrps, 2)#@threads 
        for iρ in axes(ρhats, 2)
            @views generate_dataset_on_poles(tr̂sθsϕs, tArray[:, :, iρ, ir]; rvec = rbmrps[:, ir], k=k, λ=λ, FT=FT)
            @views generate_dataset_on_poles(pr̂sθsϕs, pArray[:, :, iρ, ir]; rvec = rbmrps[:, ir], k=k, λ=λ, FT=FT)
        end
        next!(pmeter)
    end

    return reshape(tArray, length(tr̂sθsϕs)*2, :), reshape(pArray, length(pr̂sθsϕs)*2, :)

end

end # module
