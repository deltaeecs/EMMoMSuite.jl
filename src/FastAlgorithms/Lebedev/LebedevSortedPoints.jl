module LebedevSortedPoints

using LinearAlgebra
using StaticArrays
using ....Geometry
using ...MLFMA.Interpolation: octreeXWNCal

export getlbSortedData, get_t_nodes, nodes2Poles, p2nDict, n2pDict
export fibonacci_grid, high_order_nodes

const TargetDir = joinpath(@__DIR__, "../../../deps/sphere_lebedev/nodesSorted/")

function lbnt2fnDictConstruct(filedirs::String = TargetDir)
    if !isdir(filedirs)
        return Dict{Int,String}(), Dict{Int,Int}(), Dict{Int,Int}()
    end
    sphlebeFileNames = readdir(filedirs)
    t2fnDict = Dict{Int,String}()
    p2nDict = Dict{Int,Int}()
    n2pDict = Dict{Int,Int}()

    for filename in sphlebeFileNames
        try
            parts = split(filename, ".")
            p = parse(Int, parts[1])
            n = parse(Int, parts[2])
            t2fnDict[p] = filename
            p2nDict[p] = n
            n2pDict[n] = p
        catch
            nothing
        end
    end

    return t2fnDict, p2nDict, n2pDict
end

const lbnP2FILEDict, p2nDict, n2pDict = lbnt2fnDictConstruct()

"""
    fibonacci_grid(n) -> nodes (3, n)

Fibonacci 格点：准均匀、任意点数、无需数据集。黄金角 φ_i = i*π(3-√5)，
z_i = 1 - 2(i-0.5)/n。用于高阶（p > 131，无 Lebedev 数据集）的节点提供。
"""
function fibonacci_grid(n::Int)
    nodes = zeros(Float64, 3, n)
    ga = π * (3 - sqrt(5))
    for i in 1:n
        z = 1 - 2 * (i - 0.5) / n
        r = sqrt(1 - z * z)
        φ = ga * i
        nodes[1, i] = r * cos(φ)
        nodes[2, i] = r * sin(φ)
        nodes[3, i] = z
    end
    return nodes
end

"""
    high_order_nodes(p) -> nodes (3, n)

多项式阶数 p（奇数，>131，Lebedev 数据集不存在）的节点：
Fibonacci 格点，点数 n ≈ (4/3)*τ^2（τ=(p-1)/2，与 Lebedev 同效率）。
"""
function high_order_nodes(p::Int)
    τ = (p - 1) ÷ 2
    n = round(Int, 4 / 3 * τ^2)
    return fibonacci_grid(n)
end

function modiTgetFileName(p::Int, T2FILEDict::Dict)
    if isempty(T2FILEDict)
        error("Lebedev dictionary is empty. Check deps/sphere_lebedev/nodesSorted/")
    end
    p > maximum(keys(T2FILEDict)) &&
        throw(ArgumentError("p=$p is too large, no corresponding file."))

    filename = get(T2FILEDict, p) do
        modiTgetFileName(p + 1, T2FILEDict)
    end
    return filename
end

function getlbSortedData(p::Int; FT = Float64)
    filename = modiTgetFileName(p, lbnP2FILEDict)
    filepath = joinpath(TargetDir, filename)

    # Parse nNodes from filename or count lines
    nNodes = parse(Int, split(filename, ".")[2])

    nodes = zeros(FT, 3, nNodes)
    weights = zeros(FT, nNodes)

    open(filepath, "r") do file
        for ii = 1:nNodes
            line = readline(file)
            contents = split(line)
            if length(contents) >= 4
                nodes[1, ii] = parse(FT, contents[1])
                nodes[2, ii] = parse(FT, contents[2])
                nodes[3, ii] = parse(FT, contents[3])
                weights[ii] = parse(FT, contents[4])
            end
        end
    end

    return nodes, weights
end

function get_t_nodes(t; FT = Float64)
    p = 2t + 1
    nodes = if p <= maximum(keys(p2nDict))
        getlbSortedData(p; FT = FT)[1]
    else
        # 高阶（无 Lebedev 数据集）：Fibonacci 准均匀格点（任意点数）
        high_order_nodes(p)
    end
    return nodes
end

function nodes2Poles(nodes::Matrix{FT}) where {FT}
    n = size(nodes, 2)
    poles = Vector{r̂θϕInfo{FT}}(undef, n)
    for i = 1:n
        r̂ = SVector{3,FT}(nodes[:, i])
        poles[i] = r̂θϕInfo(r̂)
    end
    return poles
end

end # module
