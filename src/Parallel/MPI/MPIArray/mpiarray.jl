
SubOrArray{T,N} = Union{Array{T,N},SubArray{T,N,P,I,L}} where {T,N,P,I,L}

"""
	MPIArray{T, I, N, DT, IG}<:AbstractArray{T, N}

MPIArray used in MoM.

```
data::DT				the data stored in local rank;
indices::I				indice of `data` in global Array
dataOffset::OffsetArray{T, N, DT}	the OffsetArray of data stored in local rank
comm::MPI.Comm				`MPI.Comm`
myrank::Int				`MPI.Rank`
size::NTuple{N,Int}			global size of data
rank2indices::Dict{Int, I}		a dict record all MPI ranks to their datas global indices
ghostdata::Array{T, N}			the data stored or used in this rank
ghostindices::IG			the indices of ghostdata stored or used in this rank
grank2ghostindices::Dict{Int, I}	the data used but not stored in this rank, here provides the host rank and its indices in `ghostdata`
rrank2localindices::Dict{Int, IG}	the data host in this rank but used in other rank, here provides the remote rank and the indices used in `data`
```
"""
mutable struct MPIArray{T,I,N,DT,IG} <: AbstractArray{T,N}
    data::DT
    indices::I
    dataOffset::OffsetArray{T,N,DT}
    comm::MPI.Comm
    myrank::Int
    size::NTuple{N,Int}
    rank2indices::Dict{Int,I}
    ghostdata::Array{T,N}
    ghostindices::IG
    grank2ghostindices::Dict{Int,I}
    rrank2localindices::Dict{Int,IG}
end

const MPIVector{T,I,DT,IG} = MPIArray{T,I,1,DT,IG} where {T,I,DT,IG}
const MPIMatrix{T,I,DT,IG} = MPIArray{T,I,2,DT,IG} where {T,I,DT,IG}
const SubMPIVector{T,I,DT,IG,SI,L} =
    SubArray{T,1,MPIArray{T,I,MN,DT,IG},SI,L} where {T,I,MN,DT,IG,SI,L}
const SubMPIMatrix{T,I,DT,IG,SI,L} =
    SubArray{T,2,MPIArray{T,I,MN,DT,IG},SI,L} where {T,I,MN,DT,IG,SI,L}
const SubMPIArray{T,I,N,DT,IG,SI,L} =
    SubArray{T,N,MPIArray{T,I,MN,DT,IG},SI,L} where {T,I,N,MN,DT,IG,SI,L}
const SubOrMPIVector{T,I,DT,IG} =
    Union{MPIVector{T,I,DT,IG},SubMPIVector{T,I,DT,IG,SI,L}} where {T,I,DT,IG,SI,L}
const SubOrMPIMatrix{T,I,DT,IG} =
    Union{MPIMatrix{T,I,DT,IG},SubMPIMatrix{T,I,DT,IG,SI,L}} where {T,I,DT,IG,SI,L}
const SubOrMPIArray{T,I,N,DT,IG} =
    Union{MPIArray{T,I,N,DT,IG},SubMPIArray{T,I,N,DT,IG,SI,L}} where {T,I,N,DT,IG,SI,L}

Base.size(A::MPIArray) = A.size
Base.size(A::MPIArray, i::Integer) = A.size[i]
Base.length(A::MPIArray) = prod(A.size)
Base.eltype(::MPIArray{T,I,N}) where {T,I,N} = T
Base.getindex(A::MPIArray, I...) = getindex(A.dataOffset, I...)
Base.setindex!(A::MPIArray, X, I...) = setindex!(A.dataOffset, X, I...)
function Base.fill!(A::MPIArray, args...)
    fill!(A.data, args...)
    A
end
function Base.fill!(A::SubMPIArray, args...)
    fill!(getdata(A), args...)
    A
end
function Base.similar(A::MPIArray)
    B = deepcopy(A)
    fill!(B, 0)
    B
end

function Base.sum(x::T; root = -1) where {T<:SubOrMPIArray}

    if root == -1
        return MPI.Allreduce(sum(getdata(x)), +, getcomm(x))
    else
        return MPI.Reduce(sum(getdata(x)), +, root, getcomm(x))
    end

end

function getdata(A::MPIArray)
    A.data
end

function getdata(A::SubMPIArray)
    Alocalidcs = map(intersect, A.indices, A.parent.indices)
    view(A.parent.dataOffset, Alocalidcs...)
end

function getcomm(A::MPIArray)
    A.comm
end

function getcomm(A::SubMPIArray)
    A.parent.comm
end

function Base.copyto!(A::TA, B::TB) where {TA<:SubOrMPIArray,TB<:SubOrMPIArray}
    copyto!(getdata(A), getdata(B))
    A
end

Base.broadcast(f, A::TA, B::TB) where {TA<:SubOrMPIArray,TB<:SubOrMPIArray} =
    broadcast(f, getdata(A), getdata(B))
Base.broadcast!(
    f,
    dest::Td,
    A::TA,
    B::TB,
) where {Td<:SubOrMPIArray,TA<:SubOrMPIArray,TB<:SubOrMPIArray} =
    broadcast!(f, getdata(dest), getdata(A), getdata(B))


mpiarray(T::DataType, Asize::NTuple{N,Int}; args...) where {N} = mpiarray(T, Asize...; args...)

"""
	mpiarray(T::DataType, Asize::Vararg{Int, N}; buffersize = 0, comm = MPI.COMM_WORLD, partition = (1, MPI.Comm_size(comm))) where {N}

construct a mpi array with size `Asize` and distributed on MPI `comm` with partition.
"""
function mpiarray(
    T::DataType,
    Asize::Vararg{Int,N};
    buffersize = 0,
    comm = MPI.COMM_WORLD,
    partition = Tuple(map(i -> begin
        (i < length(Asize)) ? 1 : MPI.Comm_size(comm)
    end, 1:N)),
    rank = MPI.Comm_rank(comm),
    np = MPI.Comm_size(comm),
) where {N}

    allindices = sizeChunks2idxs(Asize, partition)
    rank2indices = Dict(zip(0:(np-1), allindices))
    indices = rank2indices[rank]

    rank2ghostindices = Dict{Int,eltype(values(rank2indices))}()
    for (k, v) in rank2indices
        rank2ghostindices[k] = map((indice, ub) -> expandslice(indice, buffersize, 1:ub), v, Asize)
    end

    ghostindices = map((indice, ub) -> expandslice(indice, buffersize, 1:ub), indices, Asize)
    ghostdata = zeros(T, map(length, ghostindices)...)

    ghostranks = indice2ranks(ghostindices, rank2indices)
    grank2gindices = grank2ghostindices(ghostranks, ghostindices, rank2indices; localrank = rank)

    remoteranks = indice2ranks(indices, rank2ghostindices)
    rrank2indices = remoterank2indices(remoteranks, indices, rank2ghostindices; localrank = rank)

    dataInGhostData = Tuple(map((i, gi) -> i .- (first(gi) - 1), indices, ghostindices))
    data = view(ghostdata, dataInGhostData...)

    A = MPIArray{T,typeof(indices),N,typeof(data),typeof(ghostindices)}(
        data,
        indices,
        OffsetArray(data, indices),
        comm,
        rank,
        Asize,
        rank2indices,
        ghostdata,
        ghostindices,
        grank2gindices,
        rrank2indices,
    )

    sync!(A)

    MPI.Barrier(comm)
    return A

end

"""
	sync!(A::MPIArray)

Synchronize data in `A` between MPI ranks.
"""
function sync!(A::MPIArray; comm = A.comm, rank = A.myrank, np = MPI.Comm_size(comm))

    # begin sync
    req_all = MPI.Request[]
    begin
        for (ghostrank, indices) in A.grank2ghostindices
            req =
                MPI.Irecv!(view(A.ghostdata, indices...), ghostrank, ghostrank * np + rank, A.comm)
            push!(req_all, req)
        end
        for (remoterank, indices) in A.rrank2localindices
            req = MPI.Isend(A.data[indices...], remoterank, rank * np + remoterank, A.comm)
            push!(req_all, req)
        end
    end
    MPI.Waitall(MPI.RequestSet(req_all), MPI.Status)

    MPI.Barrier(A.comm)

    nothing

end


"""
    gather(A::MPIArray; root = 0)

Gather data in `A` to `root`.
"""
function gather(A::MPIArray; root = 0)

    rank = MPI.Comm_rank(A.comm)

    Alc = rank == root ? zeros(eltype(A), A.size) : nothing

    reqs = [MPI.Isend(Array(A.data), A.comm; dest = root)]

    if rank == root
        append!(
            reqs,
            map(
                (rk, indices) -> MPI.Irecv!(view(Alc, indices...), A.comm; source = rk),
                keys(A.rank2indices),
                values(A.rank2indices),
            ),
        )
    else
        nothing
    end

    MPI.Waitall(MPI.RequestSet(reqs), MPI.Status)

    return Alc

end
