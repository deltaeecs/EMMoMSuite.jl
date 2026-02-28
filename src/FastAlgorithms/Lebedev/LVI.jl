module LVI

using SparseArrays
using LinearAlgebra
using HDF5
using ....Geometry
using ...MLFMA.Interpolation
using ..LebedevSortedPoints: getlbSortedData, nodes2Poles, p2nDict, n2pDict
using ..pinv2interpW: runpinvCal
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
    # Assuming truncation_kernel takes ka = k * a = 2pi/lambda * a
    L = truncation_kernel(levelCubeEdgel * 2π / λ)
    truncL = ceil(Int, L)

    if 2truncL + 1 < maximum(keys(p2nDict))
        # 读取 球 t 采样点信息并返回更新的 truncL
        nodes, weights = getlbSortedData(2truncL + 1)

        # 创建Poles实例保存
        r̂sθsϕs = nodes2Poles(nodes)
        poles = LbPolesInfo{FT}(weights, r̂sθsϕs)

        return truncL, poles
    else
        @warn "本层大小超出 Lebedev 求积极限，换回球面高斯求积。"
        return levelIntegralInfoCal(levelCubeEdgel, Val(:Lagrange2Step); λ = λ)
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
    depath = joinpath(@__DIR__, "../../../deps/InterpolationWeights/"),
)

    fn = joinpath(depath, "$(pk)to$(pt).h5")
    if !isfile(fn)
        runpinvCal(pk, pt; FT = FT)
    end
    w = h5open(fn, "r") do file
        load_sparse_matrix(file, "data")
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
    # 多极子数
    nt = length(tLevelPoles.Wθϕs)
    nk = length(kLevelPoles.Wθϕs)
    # 多项式阶数
    pk = n2pDict[nk]
    pt = n2pDict[nt]
    # 插值矩阵
    LbTrainedInterp1tepInfo(pk, pt; FT = FT)
end


function Interpolation.interpolationCSCMatCal(
    tLevelPoles::GLPolesInfo{FT},
    kLevelPoles::LbPolesInfo{FT},
    ::IT = 8,
) where {IT<:Integer,FT<:Real}
    # 多极子数
    nk = length(kLevelPoles.Wθϕs)
    # 多项式阶数
    pk = n2pDict[nk]
    pt = 2(length(tLevelPoles.Xθs) - 1) + 1
    # 插值矩阵
    LbTrainedInterp1tepInfo(pk, pt; FT = FT)
end

end # module
