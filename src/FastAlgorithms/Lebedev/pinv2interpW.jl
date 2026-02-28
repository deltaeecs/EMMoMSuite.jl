module pinv2interpW

using LinearAlgebra
using Statistics
using SparseArrays
using HDF5
using ProgressMeter
using NearestNeighbors
using Roots
using ...MLFMA.Interpolation: truncation_kernel
using ..LebedevSortedPoints: get_t_nodes
using ..dataset_generator: generate_dataset_on_pkpt
using ....Utilities: load_sparse_matrix, save_sparse_matrix

export runpinvCal, interpWeightsInitial

"""
在两层采样点之间初始化插值矩阵（采用反距离权重）
输入
tNodes::Matrix{T} 大小为 (3, nt) 的矩阵，表示 nt 个插值点集
pNodes::Matrix{T} 大小为 (3, np) 的矩阵，表示 np 个待插值点集
nInterp::Int, 插值点数
"""
function interpWeightsInitial(tNodes::Matrix{T}, pNodes::Matrix{T}; nInterp::Integer = 10) where {T}
    # 点数
    ptNodes = size(pNodes, 2)

    # 最近 nInterp 个结点计算
    kdtree = KDTree(tNodes)
    # 最近的 nInterp 个点的id, 笛卡尔距离
    idxs, dists = knn(kdtree, pNodes, nInterp, true)
    # 
    idxs = hcat(idxs...)
    dists = hcat(dists...)

    # 转换为球面距离
    sphdists = 2asin.(dists / 2)

    # 反距离插值权重
    interpWeits = 1 ./ sphdists
    for i = 1:size(interpWeits, 2)
        interpWeits[:, i] ./= sum(interpWeits[:, i])
    end
    interpWeitsDiagonal = deepcopy(interpWeits)
    # nan值由自插值引起，对同向 ( θ → θ, ϕ → ϕ ) 插值，变为1，对异向( θ → ϕ, ϕ → θ )插值变为 0 
    for ij in eachindex(interpWeits)
        isnan(interpWeits[ij]) && begin
            interpWeits[ij] = 1
            interpWeitsDiagonal[ij] = 0
        end
    end


    raws = repeat(1:ptNodes; inner = nInterp)

    # 同向 ( θ → θ, ϕ → ϕ ) 插值矩阵
    @views interpWeitsCSC = sparse(raws, idxs[:], interpWeits[:])
    dropzeros!(interpWeitsCSC)
    # 异向( θ → ϕ, ϕ → θ ) 插值矩阵
    @views interpWeitsCSCDiagonal = sparse(raws, idxs[:], interpWeitsDiagonal[:])
    dropzeros!(interpWeitsCSCDiagonal)
    # 总的插值矩阵
    interpWeits = [
        interpWeitsCSC interpWeitsCSCDiagonal
        interpWeitsCSCDiagonal interpWeitsCSC
    ]

    return interpWeits
end

function pinv2W!(w, nInterp, xx2D, yy2D)
    # 插值点满阵则全计算
    if nInterp >= (size(w, 2) ÷ 2)
        xxpinv = pinv(xx2D)
        wFinal = yy2D * xxpinv
        w.nzval .= 1
        for j in axes(w, 2)
            for i in axes(w, 1)
                (w[i, j] != 0) && begin
                    w[i, j] = wFinal[i, j]
                end
                continue
                w[i, j] = 0
            end
        end
    else# 非满阵计算
        # 进度条
        pmeter = Progress(size(w, 1), "Calculating Interpolation weights...")

        # 对行循环计算本行插值矩阵
        Threads.@threads for irow in axes(w, 1)
            # 更新进度条
            next!(pmeter)
            # 权重行
            wrow = w[irow, :]
            # 跳过不用计算的极点
            nnz(wrow) == 1 && continue
            # 非零元素的列索引
            nzind = wrow.nzind
            # 提取 本行 右侧项
            xx = xx2D[nzind, :]
            # 计算伪逆
            xxpinv = pinv(xx)
            # 计算权重
            wi = view(yy2D, irow:irow, :) * xxpinv
            # 写入结果
            w[irow, nzind] .= reshape(wi, :)
        end
    end

    return w

end

acc(w, x, y) = mean(abs, (w * x .- y), dims = 1) ./ maximum(abs, y, dims = 1)

function saveInterpW2file(
    τt,
    τp,
    ε,
    w,
    xx2Dte,
    yy2Dte;
    dpath = joinpath(@__DIR__, "../../../deps/InterpolationWeights/"),
)

    !isdir(dpath) && mkpath(dpath)

    fileε = try
        wfile = h5open(joinpath(dpath, "$(2τt+1)to$(2τp+1).h5"), "r") do file
            load_sparse_matrix(file, "data")
        end
        @info "elements" wfile = length(wfile.nzval) w = length(w.nzval)
        mean(acc(wfile, xx2Dte, yy2Dte))
    catch
        1.0
    end

    @info "旧的精度: $fileε, 新的精度: $ε"
    if fileε > ε
        h5open(joinpath(dpath, "$(2τt+1)to$(2τp+1).h5"), "w") do file
            @info "得到更精确结果！保存中…"
            save_sparse_matrix(file, "data", w)
            file["ε"] = ε
        end
        @info "已保存结果。"
    else
        @info "比上次结果差，不保存"
    end

    nothing

end

"""
采用伪逆计算插值矩阵，对稀疏矩阵要分行计算
f(k̂ₗ₋₁) = Wₗ₋₁ₗ f(k̂ₗ)
Wₗ₋₁ₗ = pinv(f(k̂ₗ)) f(k̂ₗ₋₁) 
"""
function calWFinal(τt, τp; nInterp, xx2D, yy2D, xx2Dte, yy2Dte, FT = Float64)
    # poles
    tnodes = get_t_nodes(τt; FT = FT)
    pnodes = get_t_nodes(τp; FT = FT)

    # 初始化权重
    w = interpWeightsInitial(tnodes, pnodes; nInterp = nInterp)

    # 计算权重
    pinv2W!(w, nInterp, xx2D, yy2D)

    # 计算误差
    ε = mean(acc(w, xx2Dte, yy2Dte))

    @info "$(2τt+1) → $(2τp+1)"
    @info "测试误差" nInterp = nInterp ε = ε

    saveInterpW2file(τt, τp, ε, w, xx2Dte, yy2Dte)

    return nothing
end


function runpinvCal(pk::T, pt::T; nInterp = pk < 20 ? 9 : 8, FT = Float64) where {T<:Integer}

    @info "Calculating interp weights $pk → $pt …"
    rel_l = find_zero(x -> truncation_kernel(x) - (pk + 1) / 2, 0)
    # 生成数据集
    tArray, pArray = generate_dataset_on_pkpt(pk, pt, rel_l; FT = FT)

    # 划分数据集
    flag = trunc(Int, 0.8 * size(tArray, 2))
    @views xx2D = vcat(real(tArray[:, 1:flag]), imag(tArray[:, 1:flag]))
    @views yy2D = vcat(real(pArray[:, 1:flag]), imag(pArray[:, 1:flag]))
    @views xx2Dte = tArray[:, (flag+1):end]
    @views yy2Dte = pArray[:, (flag+1):end]

    # trunc
    τt = (pk - 1) ÷ 2
    τp = (pt - 1) ÷ 2
    # 权重 计算
    calWFinal(
        τt,
        τp;
        nInterp = nInterp,
        xx2D = xx2D,
        yy2D = yy2D,
        xx2Dte = xx2Dte,
        yy2Dte = yy2Dte,
        FT = FT,
    )

    nothing

end

end # module
