module LVI

using SparseArrays
using LinearAlgebra
using HDF5
using ....Geometry
using ...MLFMA.Interpolation
using ...MLFMA.Interpolation: interp_type
using ..LebedevSortedPoints: getlbSortedData, nodes2Poles, p2nDict, high_order_nodes
using ..pinv2interpW: runpinvCal
using ..SHInterp:
    vectorize, interp_weights_exact, interp_weights_local, interp_weights_local_orbit,
    interp_weights_auto, interp_weights_hybrid
using ....Utilities: load_sparse_matrix

export LbPolesInfo, LbTrainedInterp1tepInfo

"""
多极子的极信息，即角谱空间采样信息，基于 Lebedev 采样点
Wθϕs    ::Vector{FT}， 权重向量
rHatsθsϕs  ::Vector{r̂θϕInfo{FT}}， 球面采样信息向量
"""
struct LbPolesInfo{FT<:Real} <: AbstractPolesInfo{FT}
    Wθϕs::Vector{FT}
    r̂sθsϕs::Vector{r̂θϕInfo{FT}}
    p::Int   # 多项式阶数（2τ+1），Lebedev 或高阶 Fibonacci 节点
end

"""
计算八叉树的积分相关信息，包括截断项、各层积分点和求积权重数据
输入:
levelCubeEdgel::FT,  层盒子边长, 一般叶层为0.25λ，其中 λ 为区域局部波长。
返回值
L           ::IT， 层 截断项
levelsPoles ::Vector{GLPolesInfo{FT}}，从叶层到第 “2” 层的角谱空间采样信息
"""
function Interpolation.levelIntegralInfoCal(
    levelCubeEdgel::FT,
    ::Val{:LbTrained1Step};
    λ = 1.0,
) where {FT<:Real}
    ## 计算截断项
    # truncation_kernel 的输入是盒子边长（以 λ 计）a/λ
    L = truncation_kernel(levelCubeEdgel / λ)
    truncL = ceil(Int, L)

    p = 2truncL + 1
    if p <= maximum(keys(p2nDict))
        # 读取 球 t 采样点信息并返回更新的 truncL
        nodes, weights = getlbSortedData(p)

        # 创建Poles实例保存
        r̂sθsϕs = nodes2Poles(nodes)
        poles = LbPolesInfo{FT}(weights, r̂sθsϕs, p)

        return truncL, poles
    else
        # 高阶无 Lebedev 数据集：Fibonacci 准均匀格点 + 等权重（无 GL 回退）
        nodes = high_order_nodes(p)
        n = size(nodes, 2)
        weights = fill(FT(4π / n), n)
        r̂sθsϕs = nodes2Poles(nodes)
        poles = LbPolesInfo{FT}(weights, r̂sθsϕs, p)
        @warn "本层 p=$p 超出 Lebedev 数据集上限（$(maximum(keys(p2nDict)))），改用 Fibonacci 格点（n=$n）。"
        return truncL, poles
    end
end

"""
保存整个方向的稀疏插值矩阵，存储形式定为稠密阵，因为稀疏矩阵元素密度较高时不如直接计算稠密阵乘积
θϕCSC   ::AbstractMatrix{FT} 插值矩阵，用于左乘本层多极子矩阵插值
θϕCSCT  ::AbstractMatrix{FT} 插值矩阵的转置，用于左乘本层多极子矩阵反插值
"""
mutable struct LbTrainedInterp1tepInfo{IT,FT<:Real} <: AbstractInterpInfo{IT,FT}
    θϕCSC::SparseMatrixCSC{FT,IT}
    θϕCSCT::SparseMatrixCSC{FT,IT}
    LbTrainedInterp1tepInfo{IT,FT}() where {IT,FT<:Real} = new{IT,FT}()
    LbTrainedInterp1tepInfo{IT,FT}(
        θϕCSC::AbstractArray,
        θϕCSCT::AbstractArray,
    ) where {IT,FT<:Real} = new{IT,FT}(θϕCSC, θϕCSCT)
end

LbTrainedInterp1tepInfo(θϕCSC::AbstractArray, θϕCSCT::AbstractArray) =
    LbTrainedInterp1tepInfo{Int,eltype(θϕCSC)}(θϕCSC, θϕCSCT)

"""
带参数的构造函数
"""
function LbTrainedInterp1tepInfo(
    pk::Int,
    pt::Int;
    FT = Float64,
    method::Symbol = :sh_auto,
    Lloc::Int = 0,
    cap_rad::Float64 = 1.0,
    depath = joinpath(@__DIR__, "../../../deps/InterpolationWeights/"),
)

    w = if method == :sh_exact
        # 球谐精确一步插值：确定性、机器精度、无需训练数据
        vectorize(interp_weights_exact(pk, pt; FT = FT))
    elseif method == :sh_auto
        # 默认修复路径：小规模用精确稠密 W，大规模用局部约束稀疏 W
        interp_weights_auto(pk, pt; FT = FT)
    elseif method == :sh_hybrid
        # 混合权重：数据拟合 + 笛卡尔标量 SH 精确性约束（确定性、稀疏、构造快）
        interp_weights_hybrid(pk, pt; Lloc = 3, support_scale = 1.5, FT = FT)
    elseif method == :sh_local
        Lloc > 0 || error(":sh_local 需要 Lloc > 0")
        vectorize(interp_weights_local(pk, pt; Lloc = Lloc, cap_rad = cap_rad, FT = FT))
    elseif method == :sh_local_orbit
        Lloc > 0 || error(":sh_local_orbit 需要 Lloc > 0")
        vectorize(interp_weights_local_orbit(pk, pt; Lloc = Lloc, cap_rad = cap_rad, FT = FT))
    else
        # 原训练式权重（IDW + 逐行 pinv）：h5 缓存不存在时现场训练
        fn = joinpath(depath, "$(pk)to$(pt).h5")
        if !isfile(fn)
            runpinvCal(pk, pt; FT = FT)
        end
        h5open(fn, "r") do file
            load_sparse_matrix(file, "data")
        end
    end
    θϕCSC = convert(SparseMatrixCSC{FT,Int}, w)
    θϕCSCT = sparse(transpose(θϕCSC))

    return LbTrainedInterp1tepInfo{Int,FT}(θϕCSC, θϕCSCT)

end

"""
Lebedev一步插值
"""
function Interpolation.interpolate(
    weights::LbTrainedInterp1tepInfo{IT,FT},
    data::AbstractArray,
) where {IT,FT}
    target = zeros(eltype(data), size(weights.θϕCSC, 1))
    interpolate!(target, weights, data)
    reshape(target, :, 2)
end

"""
Lebedev一步反插值
"""
function Interpolation.anterpolate(
    weights::LbTrainedInterp1tepInfo{IT,FT},
    data::AbstractArray,
) where {IT,FT}
    target = zeros(eltype(data), size(weights.θϕCSCT, 1))
    anterpolate!(target, weights, data)
    reshape(target, :, 2)
end

"""
Lebedev一步插值
"""
function Interpolation.interpolate!(
    target::AbstractArray,
    weights::LbTrainedInterp1tepInfo{IT,FT},
    data::AbstractArray,
) where {IT,FT}
    mul!(reshape(target, :), weights.θϕCSC, reshape(data, :))
    return target
end

"""
Lebedev一步反插值
"""
function Interpolation.anterpolate!(
    target::AbstractArray,
    weights::LbTrainedInterp1tepInfo{IT,FT},
    data::AbstractArray,
) where {IT,FT}
    mul!(reshape(target, :), weights.θϕCSCT, reshape(data, :))
    return target
end

function Interpolation.interpolationCSCMatCal(
    tLevelPoles::LbPolesInfo{FT},
    kLevelPoles::LbPolesInfo{FT},
    ::IT = 8,
) where {IT<:Integer,FT<:Real}
    # 多项式阶数（显式携带，支持高阶 Fibonacci 节点）
    pk = kLevelPoles.p
    pt = tLevelPoles.p
    # 插值矩阵
    LbTrainedInterp1tepInfo(pk, pt; FT = FT, method = :sh_auto)
end

Interpolation.interp_type(::LbPolesInfo{FT}) where {FT} = LbTrainedInterp1tepInfo{Int,FT}


function Interpolation.interpolationCSCMatCal(
    tLevelPoles::GLPolesInfo{FT},
    kLevelPoles::LbPolesInfo{FT},
    ::IT = 8,
) where {IT<:Integer,FT<:Real}
    # 多项式阶数（Lebedev/Fibonacci 层显式携带）
    pk = kLevelPoles.p
    pt = 2(length(tLevelPoles.Xθs) - 1) + 1
    # 插值矩阵
    LbTrainedInterp1tepInfo(pk, pt; FT = FT, method = :sh_auto)
end

end # module
