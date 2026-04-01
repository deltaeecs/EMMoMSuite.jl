module OctreeBuilder

using StaticArrays
using LinearAlgebra
using OffsetArrays
using Logging
using ....Geometry
using ..Level
using ..Octree
using ..Interpolation
using ..Precomputations
using ..Translation

export build_octree

"""
    build_octree(leafnodes, leafCubeEdgel; λ=1.0)

Construct the hierarchical Octree structure for the Multilevel Fast Multipole Algorithm (MLFMA).

# Algorithm Overview

1.  **Domain Definition**:
    - Computes the bounding box of all basis function centers (`leafnodes`).
    - Determines the size of the root cube ("Big Cube") to encompass the entire domain.
    - Calculates the number of levels \$L\$ such that the leaf level box size is approximately `leafCubeEdgel`.

2.  **Hierarchical Decomposition**:
    - **Leaf Level**: Assigns basis functions to leaf cubes based on their spatial coordinates.
    - **Upward Pass**: Recursively groups child cubes to form parent cubes at higher levels.

3.  **Connectivity & Neighbor Lists**:
    - **Parent-Child**: Establishes links between cubes at level \$l\$ and their children at level \$l+1\$.
    - **Near Neighbors**: Identifies cubes that share a boundary (vertex, edge, or face). These interactions are computed directly.
    - **Far Neighbors**: Identifies interaction list for MLFMA (children of parent's neighbors that are not self-neighbors).

4.  **Precomputations**:
    - **Interpolation/Anterpolation**: Computes matrices to shift between aggregation/disaggregation points and radiation patterns.
    - **Translation**: Computes the translation operators (transfer functions) for all possible translation vectors at each level.

# Arguments
- `leafnodes`: A `3 \\times N` matrix of basis function centers.
- `leafCubeEdgel`: Target edge length for leaf cubes (typically \$\\approx 0.25\\lambda\$).
- `λ`: Wavelength (used for wavenumber \$k\$ calculation).

# Returns
- `OctreeInfo`: The constructed octree data structure containing all levels and precomputed data.
- `leafsIDSorted`: Permutation vector sorting basis functions according to the octree structure (Morton order or similar).
"""
function build_octree(leafnodes::Matrix{FT}, leafCubeEdgel::FT; λ = 1.0, L_min::Int = 0, near_range::Int = 4) where {FT<:Real}
    println("Building Octree...")

    # 1. Set Big Cube
    nLevels, bigCubeLowerCoor, leafCubeEdgelUsed = setBigCube(leafnodes, leafCubeEdgel)
    if nLevels < 2
        # error("Too few levels ($nLevels). Check parameters!")
        # Allow 2 levels for testing
        @warn "nLevels = $nLevels. MLFMA might degenerate to near-field only."
        if nLevels < 1
            nLevels = 1 # Force at least 1 level?
        end
    end

    # 2. Create Leaf Level
    leafLevel, leafsIDSorted =
        setLevelInfo!(nLevels, leafnodes, leafCubeEdgelUsed, bigCubeLowerCoor; λ = λ, L_min = L_min, near_range = near_range)

    # Initialize levels dictionary
    levels = Dict{Int,LevelInfo{Int,FT,LagrangeInterpInfo{Int,FT}}}(nLevels => leafLevel)

    # Track sorted IDs
    levelsCubeIDSorted = Dict{Int,Vector{Int}}()
    levelsCubeIDSorted[nLevels+1] = leafsIDSorted

    # 3. Create Non-Leaf Levels
    for ilevel = (nLevels-1):-1:1
        ilevelCubeEdgel = leafCubeEdgelUsed * (2^(nLevels - ilevel))
        level, levelIDSorted =
            setLevelInfo!(ilevel, levels[ilevel+1], ilevelCubeEdgel, bigCubeLowerCoor; λ = λ, L_min = L_min, near_range = near_range)
        levels[ilevel] = level
        levelsCubeIDSorted[ilevel+1] = levelIDSorted
    end

    # 4. Reorder Cubes
    reOrderCubeID!(nLevels, levels, levelsCubeIDSorted)

    # 5. Set Basis Function Intervals
    setBFInterval!(nLevels, levels)

    # 6. Set Kids in 8
    setLevelsCubesKidsIn8!(nLevels, levels)

    # 7. Set Far Neighbors
    for ilevel = 1:(nLevels-1)
        setKidLevelFarNeighbors!(levels[ilevel], levels[ilevel+1])
    end

    # 8. Calculate Interpolation Matrices
    compute_interpolation_matrices!(nLevels, levels)

    # 9. Precompute Shift Factors
    k = 2π / λ
    compute_shift_factors!(nLevels, levels, k)

    # 10. Precompute Transfer Factors
    compute_translation_factors!(nLevels, levels, k; near_range = near_range)

    println("Octree built successfully.")
    return OctreeInfo(nLevels, leafCubeEdgel, bigCubeLowerCoor, levels), leafsIDSorted
end

function setBigCube(nodes::Matrix{FT}, leafCubeEdgel::FT) where {FT<:Real}
    xyzmin = minimum(nodes, dims = 2)
    xyzmax = maximum(nodes, dims = 2)
    Δxyz = xyzmax .- xyzmin .+ (sqrt(2.0) - 1.0) * leafCubeEdgel
    CubeEdgel = maximum(Δxyz)
    CubeCenter = SVector{3,FT}(xyzmin .+ xyzmax) ./ 2
    nLevels = ceil(Int, log2(CubeEdgel / leafCubeEdgel))
    if nLevels < 2
        nLevels = 2
    end

    bigCubeEdgel = leafCubeEdgel * 2^nLevels
    bigCubeLowerCoor = CubeCenter .- bigCubeEdgel / 2

    return nLevels, bigCubeLowerCoor, leafCubeEdgel
end

function setLevelInfo!(
    nLevels::Integer,
    leafnodes::Matrix{FT},
    cubeEdgel::FT,
    bigCubeLowerCoor::SVector{3,FT};
    λ = 1.0,
    L_min::Int = 0,
    near_range::Int = 4,
) where {FT<:Real}
    nleaves = size(leafnodes, 2)
    nodesInCubeID3D = zeros(Int, nleaves, 4)
    nodesInCubeID3D[:, 4] = 1:nleaves

    for ileaf = 1:nleaves
        nodesInCubeID3D[ileaf, 1:3] .=
            ceil.(Int, (leafnodes[:, ileaf] .- bigCubeLowerCoor) ./ cubeEdgel)
    end

    # Sort
    nodesInCubeID3D = sortslices(nodesInCubeID3D, dims = 1)
    kidsSorted = nodesInCubeID3D[:, 4]

    # Calculate intervals
    kidsIntervals = Int[]
    push!(kidsIntervals, 1)
    for ileaf = 1:(nleaves-1)
        if nodesInCubeID3D[ileaf, 1:3] != nodesInCubeID3D[ileaf+1, 1:3]
            push!(kidsIntervals, ileaf + 1)
        end
    end
    push!(kidsIntervals, nleaves + 1)

    nCubes = length(kidsIntervals) - 1
    kidsSlice = [kidsIntervals[i]:(kidsIntervals[i+1]-1) for i = 1:nCubes]
    cubesID3D = nodesInCubeID3D[kidsIntervals[1:(end-1)], 1:3]

    # Search neighbors
    cubesNeighbors = searchNearCubes(cubesID3D, nLevels; near_range = near_range)

    cubesInfo = Vector{CubeInfo{Int,FT}}(undef, nCubes)
    for icube = 1:nCubes
        center = bigCubeLowerCoor .+ (cubesID3D[icube, :] .- 0.5) .* cubeEdgel
        cubesInfo[icube] = CubeInfo{Int,FT}(
            kidsSlice[icube],
            0:0,
            Int[],
            Int[],
            cubesNeighbors[icube],
            Int[],
            MVector{3,Int}(cubesID3D[icube, :]),
            MVector{3,FT}(center),
        )
    end

    L, poles = levelIntegralInfoCal(cubeEdgel; λ = λ, L_min = L_min)

    level = LevelInfo{Int,FT,LagrangeInterpInfo{Int,FT}}()
    level.ID = nLevels
    level.isleaf = true
    level.L = L
    level.nCubes = nCubes
    level.cubes = cubesInfo
    level.cubeEdgel = cubeEdgel
    level.poles = poles

    return level, kidsSorted
end

function setLevelInfo!(
    levelID::Integer,
    kidLevel,
    cubeEdgel::FT,
    bigCubeLowerCoor::SVector{3,FT};
    λ = 1.0,
    L_min::Int = 0,
    near_range::Int = 4,
) where {FT<:Real}
    nkidCubes = length(kidLevel.cubes)
    kidCubesInCubeID3D = zeros(Int, nkidCubes, 4)
    kidCubesInCubeID3D[:, 4] = 1:nkidCubes

    for ikidCube = 1:nkidCubes
        kidCubesInCubeID3D[ikidCube, 1:3] .= ceil.(Int, kidLevel.cubes[ikidCube].ID3D ./ 2)
    end

    kidCubesInCubeID3D = sortslices(kidCubesInCubeID3D, dims = 1)
    kidCubesSorted = kidCubesInCubeID3D[:, 4]

    kidsIntervals = Int[]
    push!(kidsIntervals, 1)
    for ikidCube = 1:(nkidCubes-1)
        if kidCubesInCubeID3D[ikidCube, 1:3] != kidCubesInCubeID3D[ikidCube+1, 1:3]
            push!(kidsIntervals, ikidCube + 1)
        end
    end
    push!(kidsIntervals, nkidCubes + 1)

    nCubes = length(kidsIntervals) - 1
    kidsSlice = [kidsIntervals[i]:(kidsIntervals[i+1]-1) for i = 1:nCubes]
    kidsIn8 = [Vector{Int}(undef, length(kidSlice)) for kidSlice in kidsSlice]

    cubesID3D = kidCubesInCubeID3D[kidsIntervals[1:(end-1)], 1:3]

    cubesNeighbors = searchNearCubes(cubesID3D, levelID; near_range = near_range)

    cubesInfo = Vector{CubeInfo{Int,FT}}(undef, nCubes)
    for icube = 1:nCubes
        center = bigCubeLowerCoor .+ (cubesID3D[icube, :] .- 0.5) .* cubeEdgel
        cubesInfo[icube] = CubeInfo{Int,FT}(
            kidsSlice[icube],
            0:0,
            kidsIn8[icube],
            Int[],
            cubesNeighbors[icube],
            Int[],
            MVector{3,Int}(cubesID3D[icube, :]),
            MVector{3,FT}(center),
        )
    end

    L, poles = levelIntegralInfoCal(cubeEdgel; λ = λ, L_min = L_min)

    level = LevelInfo{Int,FT,LagrangeInterpInfo{Int,FT}}()
    level.ID = levelID
    level.isleaf = false
    level.L = L
    level.nCubes = nCubes
    level.cubes = cubesInfo
    level.cubeEdgel = cubeEdgel
    level.poles = poles

    return level, kidCubesSorted
end

function searchNearCubes(cubesID3D::Matrix{IT}, levelID::Integer; near_range::Int = 4) where {IT<:Integer}
    nCubes = size(cubesID3D, 1)
    maxCubes1D = 2^levelID
    cubesID1D = [
        ((row[1] - 1) * maxCubes1D^2 + (row[2] - 1) * maxCubes1D + row[3]) for
        row in eachrow(cubesID3D)
    ]

    neighbors = Vector{Vector{IT}}(undef, nCubes)

    for iCube = 1:nCubes
        cubeID3D = cubesID3D[iCube, :]
        neighborsOffsets = [(-near_range:near_range), (-near_range:near_range), (-near_range:near_range)]

        for ii = 1:3
            min_offset = max(-near_range, 1 - cubeID3D[ii])
            max_offset = min(near_range, maxCubes1D - cubeID3D[ii])
            neighborsOffsets[ii] = min_offset:max_offset
        end

        neighborsNonEmpty = IT[]
        for offsetx in neighborsOffsets[1]
            for offsety in neighborsOffsets[2]
                for offsetz in neighborsOffsets[3]
                    if offsetx == 0 && offsety == 0 && offsetz == 0
                        push!(neighborsNonEmpty, iCube)
                        continue
                    end

                    neighbor1DID =
                        (cubeID3D[1] - 1 + offsetx) * maxCubes1D^2 +
                        (cubeID3D[2] - 1 + offsety) * maxCubes1D +
                        cubeID3D[3] +
                        offsetz
                    idx = searchsortedfirst(cubesID1D, neighbor1DID)
                    if idx <= length(cubesID1D) && cubesID1D[idx] == neighbor1DID
                        push!(neighborsNonEmpty, idx)
                    end
                end
            end
        end
        sort!(neighborsNonEmpty)
        neighbors[iCube] = neighborsNonEmpty
    end
    return neighbors
end

function reOrderCubeID!(
    nLevels::Integer,
    levels::Dict{Int,LV},
    levelsCubeIDSorted::Dict{Int,Vector{Int}},
) where {LV<:AbstractLevel}
    for ilevel = 2:nLevels
        parentLevel = levels[ilevel-1]
        level = levels[ilevel]
        cubeIDSorted = levelsCubeIDSorted[ilevel]

        permute!(cubeIDSorted, vcat([cube.kidsInterval for cube in parentLevel.cubes]...))

        nCumKidCube = 0
        for pCube in parentLevel.cubes
            nkid = length(pCube.kidsInterval)
            pCube.kidsInterval = (nCumKidCube+1):(nCumKidCube+nkid)
            nCumKidCube += nkid
        end

        permute!(level.cubes, cubeIDSorted)

        nCubes = length(cubeIDSorted)
        oldNewIDpair = hcat(cubeIDSorted, 1:nCubes)
        oldNewIDpair = sortslices(oldNewIDpair, dims = 1)

        for cube in level.cubes
            cube.neighbors .= oldNewIDpair[cube.neighbors, 2]
            sort!(cube.neighbors)
        end
    end

    leafLevel = levels[nLevels]
    leafsIDSorted = levelsCubeIDSorted[nLevels+1]
    permute!(leafsIDSorted, vcat([cube.kidsInterval for cube in leafLevel.cubes]...))

    nCumKidCube = 0
    for lCube in leafLevel.cubes
        nkid = length(lCube.kidsInterval)
        lCube.kidsInterval = (nCumKidCube+1):(nCumKidCube+nkid)
        nCumKidCube += nkid
    end
end

function setBFInterval!(nLevels::Integer, levels::Dict{Int,LV}) where {LV<:AbstractLevel}
    leafCubes = levels[nLevels].cubes
    for cube in leafCubes
        cube.bfInterval = cube.kidsInterval
    end

    for ilevel = nLevels-1:-1:1
        tlevelCubes = levels[ilevel].cubes
        klevelCubes = levels[ilevel+1].cubes
        for tlevelCube in tlevelCubes
            tklevelCubes = view(klevelCubes, tlevelCube.kidsInterval)
            kidsbfstart = tklevelCubes[1].bfInterval[1]
            kidsbfend = tklevelCubes[end].bfInterval[end]
            tlevelCube.bfInterval = kidsbfstart:kidsbfend
        end
    end
end

function setLevelsCubesKidsIn8!(nLevels::Integer, levels::Dict{Int,LV}) where {LV<:AbstractLevel}
    for iLevel = (nLevels-1):-1:2
        tCubes = levels[iLevel].cubes
        kCubes = levels[iLevel+1].cubes
        for iCube in eachindex(tCubes)
            cube = tCubes[iCube]
            cubeID3D4 = cube.ID3D .* 4 .- 2
            for (ii, kCube) in enumerate(view(kCubes, cube.kidsInterval))
                relativeOffset = (kCube.ID3D .* 2 .- 1) .- cubeID3D4
                relativeOffsetTrunc = [relativeOffset[i] == -1 ? 0 : 1 for i = 1:3]
                relativeOffset1D =
                    1 +
                    relativeOffsetTrunc[1] +
                    relativeOffsetTrunc[2] * 2 +
                    relativeOffsetTrunc[3] * 4
                cube.kidsIn8[ii] = relativeOffset1D
            end
        end
    end
end

function setKidLevelFarNeighbors!(thisLevel, kidLevel)
    cubes = thisLevel.cubes
    kidCubes = kidLevel.cubes
    for icube in eachindex(cubes)
        cube = cubes[icube]
        neighbors = cube.neighbors
        kidCubesID = cube.kidsInterval

        neighborsKids = Int[]
        for iNeighbor in neighbors
            if iNeighbor != icube
                append!(neighborsKids, cubes[iNeighbor].kidsInterval)
            end
        end

        for iKidCube in kidCubesID
            kidCube = kidCubes[iKidCube]
            kidCubeNeighbors = kidCube.neighbors
            kidCube.farneighbors = sort!(setdiff(neighborsKids, kidCubeNeighbors))
        end
    end
end

end
