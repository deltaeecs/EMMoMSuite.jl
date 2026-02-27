"""
    ArrayChunk{T, N} <:AbstractArray{T, N}

Array chunk with custom indices.
```
data::Array{T, N}
indices::NTuple{N, Union{UnitRange{Int}, Vector{Int}}}
```
"""
struct ArrayChunk{T, N} <:AbstractArray{T, N}
    data::Array{T, N}
    indices::NTuple{N, Union{UnitRange{Int}, Vector{Int}}}
end

"""
    ArrayChunk(T::DataType, indices::Vararg{Union{UnitRange{Int}, Vector{Int}}, N}) where {N}

Initial ArrayChunk with datatype `T` and indics `indices`.
"""
function ArrayChunk(T::DataType, indices::Vararg{Union{UnitRange{Int}, Vector{Int}}, N}) where {N}
    size =  map(length, indices)
    data =  zeros(T, size...)
    ArrayChunk{T, N}(data, indices)
end

"""
    ArrayChunk(T::DataType, indices::Vararg{Union{UnitRange{Int}, Vector{Int}}, N}) where {N}

Initial ArrayChunk with data `data` and indics `indices`.
"""
function ArrayChunk(data::Array{T, N}, indices...) where {T,N}
    any(i -> size(data, i)!= length(indices[i]), 1:N) && throw("Demension miss match!")
    ArrayChunk{T, N}(data, indices)
end

import Base:size, show, display, eltype, length, fill!, getindex, setindex!, sum

size(A::T) where{T<:ArrayChunk}  = size(A.data)
size(A::T, d::Integer) where{T<:ArrayChunk}  = size(A.data, d)
eltype(::ArrayChunk{T, N}) where{T, N} = T
length(A::T) where{T<:ArrayChunk} = prod(size(A))
fill!(A::T, x) where{T<:ArrayChunk} = fill!(A.data, x)
sum(A::T) where{T<:ArrayChunk} = sum(A.data)

"""
    getindex(A::ArrayChunk{T, N}, idcs::Vararg{I, N}) where{T, N, I<:Integer}

overload getindex for `ArrayChunk`.
"""
function getindex(A::ArrayChunk{T, N}, idcs::Vararg{I, N}) where{T, N, I<:Integer}
    
    ptrs    =   map(searchsorted, A.indices, idcs)
    any(isempty, ptrs) && return zero(eltype(A.data))
    idcslc  =   map(first, ptrs)
    return A.data[idcslc...]

end

"""
    setindex!(A::ArrayChunk{T, N}, x, idcs::Vararg{I, N}) where{T, N, I<:Integer}

overload setindex! for `ArrayChunk`.
"""
function setindex!(A::ArrayChunk{T, N}, x, idcs::Vararg{I, N}) where{T, N, I<:Integer}

    ptrs    =   map(searchsorted, A.indices, idcs)
    any(isempty, ptrs) && begin @warn "$idcs not in indices of array!"; end
    idcslc  =   map(first, ptrs)
    setindex!(A.data, x, idcslc...)
    return nothing
end
