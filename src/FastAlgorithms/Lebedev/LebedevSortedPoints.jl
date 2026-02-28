module LebedevSortedPoints

using LinearAlgebra
using StaticArrays
using ....Geometry
using ...MLFMA.Interpolation: octreeXWNCal

export getlbSortedData, get_t_nodes, nodes2Poles, p2nDict, n2pDict

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
        # θ direction
        Xcosθs, Wθs = octreeXWNCal(one(FT), -one(FT), t, :glq)
        Xθs = acos.(Xcosθs)
        # ϕ direction
        Xϕs, Wϕs = octreeXWNCal(zero(FT), convert(FT, 2π), t, :uni)

        reduce(hcat, [r̂θϕInfo(θ, ϕ).r̂ for ϕ in Xϕs for θ in Xθs])
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
