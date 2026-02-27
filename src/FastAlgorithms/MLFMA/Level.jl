module Level

using StaticArrays
using OffsetArrays
using ..Interpolation

export CubeInfo, AbstractLevel, LevelInfo

mutable struct CubeInfo{IT<:Integer, FT<:Real}
    kidsInterval    ::UnitRange{IT}
    bfInterval      ::UnitRange{IT}
    kidsIn8         ::Vector{IT}
    geoIDs          ::Vector{IT}
    neighbors       ::Vector{IT}
    farneighbors    ::Vector{IT}
    ID3D            ::MVector{3, IT}
    center          ::MVector{3, FT}
end

abstract type AbstractLevel end

mutable struct LevelInfo{IT<:Integer, FT<:Real, IPT} <: AbstractLevel
    ID          ::IT
    isleaf      ::Bool
    L           ::IT
    nCubes      ::IT
    cubes       ::Vector{CubeInfo{IT, FT}}
    cubeEdgel   ::FT
    poles       ::AbstractPolesInfo{FT}
    interpWθϕ   ::IPT
    aggS        ::Array{Complex{FT}, 3}
    disaggG     ::Array{Complex{FT}, 3}
    phaseShift2Kids     ::Array{Complex{FT}, 2}
    phaseShiftFromKids  ::Array{Complex{FT}, 2}
    αTrans      ::Array{Complex{FT}, 2}
    αTransIndex ::OffsetArray{IT, 3, Array{IT, 3}}

    function LevelInfo{IT, FT, IPT}() where {IT<:Integer, FT<:Real, IPT}
        new{IT, FT, IPT}()
    end
end

end
