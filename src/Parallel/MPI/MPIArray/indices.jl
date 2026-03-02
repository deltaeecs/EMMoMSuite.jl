"""
    expandslice(idxs, bandwidth, bounds)

expand `idxs` with `bandwidth` in `bounds`.
"""
function expandslice(idxs, bandwidth, bounds)
    d = max(first(bounds), first(idxs) - bandwidth)
    u = min(last(bounds), last(idxs) + bandwidth)
    d:u
end

"""
    slicedim2bounds(sz::Int, nc::Int)

将区间 `1:sz` 划分为 `nc` 个区间并返回区间上下界。
"""
function slicedim2bounds(sz::Int, nc::Int)
    if sz >= nc
        chunk_size = div(sz, nc)
        remainder = rem(sz, nc)
        grid = zeros(Int64, nc + 1)
        for i = 1:(nc+1)
            grid[i] += (i - 1) * chunk_size + 1
            if i <= remainder
                grid[i] += i - 1
            else
                grid[i] += remainder
            end
        end
        return grid
    else
        return [[1:(sz+1);]; zeros(Int, nc - sz)]
    end
end

"""
    slicedim2bounds(dims, nc::Int)

将区间 `dims` 划分为 `nc` 个区间并返回区间上下界。
"""
function slicedim2partition(dims, nc::Int)
    dims = [dims...]
    chunks = ones(Int, length(dims))
    f = sort!(collect(keys(factor(nc))), rev = true)
    k = 1
    while nc > 1
        # repeatedly allocate largest factor to largest dim
        if nc % f[k] != 0
            k += 1
            if k > length(f)
                break
            end
        end
        fac = f[k]
        (d, dno) = findmax(dims)
        # resolve ties to highest dim
        dno = findlast(isequal(d), dims)
        if dims[dno] >= fac
            dims[dno] = div(dims[dno], fac)
            chunks[dno] *= fac
        end
        nc = div(nc, fac)
    end
    return Tuple(chunks)
end


@doc """
    sizeChunks2cuts(Asize, chunks)
    sizeChunks2cuts(Asize::Int, chunks)
    sizeChunks2cuts(Asize, chunks::Int)
    sizeChunks2cuts(Asize::Int, chunks::Int)

将数组大小 `Asize` 按 `chunks` 进行分块。
"""
function sizeChunks2cuts(Asize, chunks)
    map(slicedim2bounds, Asize, chunks)
end
function sizeChunks2cuts(Asize::Int, chunks)
    map(slicedim2bounds, (Asize,), chunks)
end
function sizeChunks2cuts(Asize, chunks::Int)
    map(slicedim2bounds, Asize, (chunks,))
end
function sizeChunks2cuts(Asize::Int, chunks::Int)
    map(slicedim2bounds, (Asize,), (chunks,))
end

@doc """
    sizeChunksCuts2indices(Asize, nchunk, cuts::Tuple)
    sizeChunksCuts2indices(Asize, nchunk, cuts::Vector{I}) where{I<:Integer}

根据数组大小 `Asize` 分块数量 `nchunk` 以及各块索引区间 `cuts` 计算各块的索引。
"""
function sizeChunksCuts2indices(Asize, nchunk, cuts::Tuple)
    n = length(Asize)
    idxs = Array{NTuple{n,UnitRange{Int}},n}(undef, nchunk...)
    for cidx in CartesianIndices(tuple(nchunk...))
        if n > 0
            idxs[cidx.I...] = ntuple(i -> (cuts[i][cidx[i]]:cuts[i][cidx[i]+1]-1), n)
        else
            throw("0 dim array not supported.")
        end
    end
    return idxs
end
function sizeChunksCuts2indices(Asize, nchunk, cuts::Vector{I}) where {I<:Integer}
    n = length(Asize)
    idxs = Array{NTuple{n,UnitRange{Int}},n}(undef, nchunk...)
    for cidx in CartesianIndices(tuple(nchunk...))
        idxs[cidx.I...] = (cuts[cidx[1]]:cuts[cidx[1]+1]-1,)
    end
    return idxs
end

"""
    sizeChunks2idxs(Asize, nchunk)

Borrowed form [DistributedArray.jl](https://github.com/JuliaParallel/DistributedArrays.jl), get the slice of matrix
size `Asize` on each dimension with `nchunk`.
"""
function sizeChunks2idxs(Asize, nchunk)
    cuts = sizeChunks2cuts(Asize, nchunk)
    return sizeChunksCuts2indices(Asize, nchunk, cuts)
end

@doc """
    indice2rank(indice::T, rank2indices::Dict{Integer, NTuple{1}}) where{T<:Integer}
    indice2rank(indice::NTuple{N, T}, rank2indices::Dict{Int, Tuple{Vararg{T2, N}}}) where{T<:Integer, N, T2}
    indice2ranks(indice::NTuple{N, Union{UnitRange{T}, Vector{T}}}, rank2indices) where{T, N}

Get the rank of `indice` form `rank2indices`
"""
function indice2rank(indice::T, rank2indices::Dict{Integer,NTuple{1}}) where {T<:Integer}

    # findall on a Dict returns keys whose values satisfy predicate.
    # NTuple{1, UnitRange}: check if indice falls within the UnitRange.
    re = findall(x -> any(r -> indice in r, x), rank2indices)

    isempty(re) && throw("No suitable rank found, please recheck!")
    length(re) > 1 && throw("Multi ranks found, please recheck!")

    return re[1]

end
function indice2rank(
    indice::NTuple{N,T},
    rank2indices::Dict{Int,Tuple{Vararg{T2,N}}},
) where {T<:Integer,N,T2}
    rks = map(i -> findall(x -> indice[i] in x[i], rank2indices), 1:N)
    re = intersect(rks...)
    # isempty(re) && throw("No suitable rank found, please recheck!")
    length(re) > 1 && throw("Multi ranks found, please recheck!")
    return re[1]

end
function indice2ranks(indice::NTuple{N,Union{UnitRange{T},Vector{T}}}, rank2indices) where {T,N}
    rks = map(i -> findall(x -> !isempty(intersect(indice[i], x[i])), rank2indices), 1:N)
    re = intersect(rks...)
    # isempty(re) && throw("No suitable rank found, please recheck!")
    return sort!(re)
end


"""
    Base.searchsortedfirst(a::Vector{T}, x) where {T<:UnitRange}

找出已排序、无重叠的 Vector{UnitRange} 中的某个 UnitRange 的开头元素，在整个区间中的位置。
"""
function Base.searchsortedfirst(a::Vector{T}, x) where {T<:UnitRange}
    gid = searchsortedfirst(a, x, by = first)
    return sum(length, view(a, 1:(gid-1)); init = 0)
end

@doc """
    grank2ghostindices(ghostranks, ghostindices::Tuple{Vararg{T1, N}}, rank2indices::Dict{Int, Tuple{Vararg{T2, N}}}; localrank = MPI.Comm_rank(MPI.COMM_WORLD)) where{N, T1, T2}
    grank2ghostindices(ghostranks, ghostindices::Tuple{Vararg{T1, N}}, rank2indices::Dict{Int, Tuple{Vararg{T2, N}}}; localrank = MPI.Comm_rank(MPI.COMM_WORLD)) where{N, T1<:Vector{Int}, T2}
    grank2ghostindices(ghostranks, ghostindices::Tuple{Vararg{T1, N}}, rank2indices::Dict{Int, Tuple{Vararg{T2, N}}}; localrank = MPI.Comm_rank(MPI.COMM_WORLD)) where{N, T1<:Vector{UnitRange}, T2}
    grank2ghostindices(ghostranks, ghostindices::NTuple{N, Union{UnitRange{Int}, Vector{Int}}}, rank2indices::Dict{Int, Tuple{Vararg{T2, N}}}; localrank = MPI.Comm_rank(MPI.COMM_WORLD)) where{N, T2}

get indices of ghost data in its hosting rank.
"""
function grank2ghostindices(
    ghostranks,
    ghostindices::Tuple{Vararg{T1,N}},
    rank2indices::Dict{Int,Tuple{Vararg{T2,N}}};
    localrank = MPI.Comm_rank(MPI.COMM_WORLD),
) where {N,T1,T2}

    grank2gindices = Dict{Int,Tuple{Vararg{T1,N}}}()
    for grank in ghostranks
        grank == localrank && continue
        grank2gindices[grank] = Tuple([
            intersect(rank2indices[grank][i], ghostindices[i]) .- (first(ghostindices[i]) - 1)
            for i = 1:N
        ])
    end

    return grank2gindices
end
function grank2ghostindices(
    ghostranks,
    ghostindices::Tuple{Vararg{T1,N}},
    rank2indices::Dict{Int,Tuple{Vararg{T2,N}}};
    localrank = MPI.Comm_rank(MPI.COMM_WORLD),
) where {N,T1<:Vector{Int},T2}

    grank2gindices = Dict{Int,Tuple{Vararg{T2,N}}}()
    for grank in ghostranks
        grank == localrank && continue
        intersectIndice = map(intersect, rank2indices[grank], ghostindices)
        grank2gindices[grank] = map(
            (gidc, intersidc) ->
                searchsortedfirst(gidc, intersidc[1]) .+ (0:(length(intersidc)-1)),
            ghostindices,
            intersectIndice,
        )
    end

    return grank2gindices
end
function grank2ghostindices(
    ghostranks,
    ghostindices::Tuple{Vararg{T1,N}},
    rank2indices::Dict{Int,Tuple{Vararg{T2,N}}};
    localrank = MPI.Comm_rank(MPI.COMM_WORLD),
) where {N,T1<:Vector{UnitRange},T2}

    grank2gindices = Dict{Int,Tuple{Vararg{T2,N}}}()
    for grank in ghostranks
        grank == localrank && continue
        intersectIndice = map(
            (rk2indice, gindice) -> reduce(vcat, map(g -> intersect(rk2indice, g), gindice)),
            rank2indices[grank],
            ghostindices,
        )
        grank2gindices[grank] = map(
            (gidc, intersidc) ->
                searchsortedfirst(gidc, intersidc[1]) .+ (0:(length(intersidc)-1)),
            ghostindices,
            intersectIndice,
        )
    end

    return grank2gindices
end
function grank2ghostindices(
    ghostranks,
    ghostindices::NTuple{N,Union{UnitRange{Int},Vector{Int}}},
    rank2indices::Dict{Int,Tuple{Vararg{T2,N}}};
    localrank = MPI.Comm_rank(MPI.COMM_WORLD),
) where {N,T2}

    grank2gindices = Dict{Int,Tuple{Vararg{T2,N}}}()
    for grank in ghostranks
        grank == localrank && continue
        intersectIndice = map(intersect, rank2indices[grank], ghostindices)
        grank2gindices[grank] = map(
            (gidc, intersidc) ->
                searchsortedfirst(gidc, intersidc[1]) .+ (0:(length(intersidc)-1)),
            ghostindices,
            intersectIndice,
        )
    end

    return grank2gindices
end


"""
    grank2gdataSize(ghostranks, ghostindices::NTuple{N, Union{UnitRange{Int}, Vector{Int}}}, rank2indices::Dict{Int, Tuple{Vararg{T2, N}}}; localrank = MPI.Comm_rank(MPI.COMM_WORLD)) where{N, T2}

get size ghost data.
"""
function grank2gdataSize(
    ghostranks,
    ghostindices::NTuple{N,Union{UnitRange{Int},Vector{Int}}},
    rank2indices::Dict{Int,Tuple{Vararg{T2,N}}};
    localrank = MPI.Comm_rank(MPI.COMM_WORLD),
) where {N,T2}

    grank2gsize = Dict{Int,Int}()
    for grank in ghostranks
        grank == localrank && continue
        grank2gdataSize = map(intersect, rank2indices[grank], ghostindices)
        grank2gsize[grank] = prod(length, grank2gdataSize)
    end

    return grank2gsize
end

"""
    grank2indices(ghostranks, ghostindices::NTuple{N, Union{UnitRange{T}, Vector{T}}}, rank2indices::Dict{Int, Tuple{Vararg{T2, N}}}; localrank = MPI.Comm_rank(MPI.COMM_WORLD)) where{N, T,  T2}

get global indices of ghost data in ghost ranks.
"""
function grank2indices(
    ghostranks,
    ghostindices::NTuple{N,Union{UnitRange{T},Vector{T}}},
    rank2indices::Dict{Int,Tuple{Vararg{T2,N}}};
    localrank = MPI.Comm_rank(MPI.COMM_WORLD),
) where {N,T,T2}

    grank2indices = Dict{Int,typeof(ghostindices)}()
    for grank in ghostranks
        grank == localrank && continue
        grank2indices[grank] = map(intersect, rank2indices[grank], ghostindices)
    end

    return grank2indices
end

@doc """
    intersectInIdc(idc::UnitRange, intersidc)
    intersectInIdc(idc::Vector, intersidc)

找出 `intersidc` 在 `idc` 中的位置。
"""
function intersectInIdc(idc::UnitRange, intersidc)
    intersidc .- (first(idc) - 1)
end
function intersectInIdc(idc::Vector, intersidc)
    map(i -> searchsortedfirst(idc, i), intersidc)
end


"""
    remoterank2indices(remoteranks, indices, rank2ghostindices::Dict{Int, T}; localrank = MPI.Comm_rank(MPI.COMM_WORLD)) where{T}

get remote rank and relative indices in data.
"""
function remoterank2indices(
    remoteranks,
    indices,
    rank2ghostindices::Dict{Int,T};
    localrank = MPI.Comm_rank(MPI.COMM_WORLD),
) where {T}

    rrank2indices = Dict{Int,T}()
    for rrank in remoteranks
        rrank == localrank && continue
        intersectIndice = map(intersect, rank2ghostindices[rrank], indices)
        rrank2indices[rrank] = map(intersectInIdc, indices, intersectIndice)
    end

    return rrank2indices
end

"""
    remoterank2indices(remoteranks, indices::Tuple{Vararg{T1, N}}, rank2ghostindices::Dict{Int, Tuple{Vararg{T2, N}}}; localrank = MPI.Comm_rank(MPI.COMM_WORLD)) where{N, T1<:UnitRange, T2<:Union{UnitRange, Vector}}

get indices of ghost data in its hosting rank.
"""
function remoterank2indices(
    remoteranks,
    indices::Tuple{Vararg{T1,N}},
    rank2ghostindices::Dict{Int,Tuple{Vararg{T2,N}}};
    localrank = MPI.Comm_rank(MPI.COMM_WORLD),
) where {N,T1<:UnitRange,T2<:Union{UnitRange,Vector}}

    rrank2indices = Dict{Int,Tuple{Vararg{T2,N}}}()
    for rrank in remoteranks
        rrank == localrank && continue
        intersectIndice = map(intersect, rank2ghostindices[rrank], indices)
        rrank2indices[rrank] = map(
            (gidc, intersidc) ->
                searchsortedfirst(gidc, intersidc[1]) .+ (intersidc .- intersidc[1]),
            indices,
            intersectIndice,
        )
    end

    return rrank2indices
end

"""
    get_rank2ghostindices(ghostranks, indices, ghostindices; comm = MPI.COMM_WORLD, rank = MPI.Comm_rank(comm), np = MPI.Comm_size(comm))

获取 ghostranks 到 ghostindices 的字典。
"""
function get_rank2ghostindices(
    ghostranks,
    indices,
    ghostindices;
    comm = MPI.COMM_WORLD,
    rank = MPI.Comm_rank(comm),
    np = MPI.Comm_size(comm),
)

    # 传输本地ghostdata的大小数据
    gsize = zeros(Int, length(ghostranks), length(indices))
    reqs = MPI.Request[]
    for (i, grk) in enumerate(ghostranks)
        push!(reqs, MPI.Isend([length(gi) for gi in ghostindices], grk, rank * np + grk, comm))
        push!(reqs, MPI.Irecv!(view(gsize, i, :), grk, grk * np + rank, comm))
    end
    MPI.Waitall(MPI.RequestSet(reqs), MPI.Status)

    # 分配用到的远程进程的 ghost indices
    rank2ghostindices = Dict{Int,typeof(ghostindices)}()
    for (i, grk) in enumerate(ghostranks)
        rank2ghostindices[grk] = map(sz -> zeros(Int, sz), Tuple(gsize[i, :]))
    end
    reqs = MPI.Request[]
    for rk in ghostranks
        rk == rank && continue
        map(gi -> push!(reqs, MPI.Isend(gi, rk, rk * np + rank, comm)), ghostindices)
        map(gi -> push!(reqs, MPI.Irecv!(gi, rk, rank * np + rk, comm)), rank2ghostindices[rk])
    end
    MPI.Waitall(MPI.RequestSet(reqs), MPI.Status)

    return rank2ghostindices

end
