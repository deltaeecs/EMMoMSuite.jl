module Octree

using StaticArrays
using ..Level

export OctreeInfo

struct OctreeInfo{FT<:Real, LT<:AbstractLevel}
    nLevels         ::Integer
    leafCubeEdgel   ::FT
    bigCubeLowerCoor::SVector{3, FT}
    levels          ::Dict{Int, LT}
end

end
