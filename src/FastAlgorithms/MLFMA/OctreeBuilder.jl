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
- `near_range`: Near-field radius in cube offsets, uniform for the whole tree.
  `nothing` (default) derives an **adaptive** radius from the leaf level that
  guarantees real-spectrum M2L convergence (`kR_min ≥ 0.55·L_leaf` — the leaf
  has the smallest `w` and is the binding level; calibrated on gate GD2V). An
  explicit `Int` pins the radius (backward compatible — a too-small value makes
  the diagonal M2L inaccurate and triggers a warning). Larger radii grow the
  direct near-field matrix as `(2nr+1)³` per cube — pair with
  `BlockJacobiPreconditioner`.

# Returns
- `OctreeInfo`: The constructed octree data structure containing all levels and precomputed data.
- `leafsIDSorted`: Permutation vector sorting basis functions according to the octree structure (Morton order or similar).
"""
function build_octree(
    leafnodes::Matrix{FT},
    leafCubeEdgel::FT;
    λ = 1.0,
    L_min::Int = 0,
    near_range::Union{Int,Nothing} = nothing,
    interp_method::Val = Val(:Lagrange2Step),
    m2l_stabilization::Symbol = :unscaled,
) where {FT<:Real}
    @info "Building Octree..."

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

    # 1.5 Tree-uniform near radius derived from the LEAF level (issue #22 #3).
    # The near/far list tiling requires ONE radius for the whole tree, and the
    # leaf level (smallest w) is the binding constraint for real-spectrum M2L
    # convergence: kR_min = (nr+1)·k·w_leaf ≥ 0.55·L_leaf.
    L_leaf = max(Interpolation.truncationLCal(leafCubeEdgelUsed; λ = λ), L_min)
    nr_tree = near_range === nothing ?
        adaptive_near_range(λ, leafCubeEdgelUsed, L_leaf) : near_range

    # 1.6 Efficiency guidance (issue #22, problem 3): the direct near-field
    # matrix grows as (2nr+1)^3 cubes per cube. Adaptive radii above ~2 mean
    # the leaf cubes are electrically small relative to the evanescent reach;
    # electrically larger leaf cubes keep nr small at equal M2L accuracy
    # (kR_min scales with w while L grows sublinearly).
    if near_range === nothing && nr_tree > 2
        @warn "Adaptive near_range = $nr_tree (> 2): the direct near-field matrix " *
              "covers (2·$nr_tree+1)^3 cubes per cube, which grows assembly time and " *
              "memory. This leaf cube (≈ $(round(leafCubeEdgelUsed / λ, digits = 2)) λ) is " *
              "electrically small; an electrically larger leaf cube (≈ 0.5λ–1λ keeps " *
              "near_range ≈ 2–3) trades pole count for a much smaller near-field " *
              "matrix (the adaptive radius follows kR_min ≥ 0.55·L). Pair the " *
              "near-field matrix with BlockJacobiPreconditioner. See issue #22 (problem 3)."
    end

    # 2. Create Leaf Level
    leafLevel, leafsIDSorted = setLevelInfo!(
        nLevels,
        leafnodes,
        leafCubeEdgelUsed,
        bigCubeLowerCoor;
        λ = λ,
        L_min = L_min,
        near_range = nr_tree,
        interp_method = interp_method,
    )

    # Initialize levels dictionary
    levels = Dict{Int,AbstractLevel}(nLevels => leafLevel)

    # Track sorted IDs
    levelsCubeIDSorted = Dict{Int,Vector{Int}}()
    levelsCubeIDSorted[nLevels+1] = leafsIDSorted

    # 3. Create Non-Leaf Levels
    for ilevel = (nLevels-1):-1:1
        ilevelCubeEdgel = leafCubeEdgelUsed * (2^(nLevels - ilevel))
        level, levelIDSorted = setLevelInfo!(
            ilevel,
            levels[ilevel+1],
            ilevelCubeEdgel,
            bigCubeLowerCoor;
            λ = λ,
            L_min = L_min,
            near_range = nr_tree,
            interp_method = interp_method,
        )
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
    compute_shift_factors!(nLevels, levels, k;
        near_range = nr_tree, m2l_stabilization = m2l_stabilization)

    # 10. Precompute Transfer Factors (tree-uniform near radius, leaf-derived)
    compute_translation_factors!(nLevels, levels, k;
        near_range = nr_tree, m2l_stabilization = m2l_stabilization)

    # M2L 有效距离校验（issue #22 问题 3）：对角实谱 M2L 的浮点收敛要求
    # 叶层最小 far 对距离 kR_min ≥ kr_factor·L（kr_factor = 0.55，经 GD2V 门校准）。
    # near_range = nothing（默认）时逐层自适应半径已按此构造，此处作为护栏复核；
    # 显式 near_range 过小时给出可行动警告，避免静默精度损失。
    leaf_level = levels[nLevels]
    if leaf_level.nCubes > 0
        kR_min = Inf
        n_far = 0
        for iCube in 1:leaf_level.nCubes
            cube = leaf_level.cubes[iCube]
            for fid in cube.farneighbors
                far = leaf_level.cubes[fid]
                off = (
                    cube.ID3D[1] - far.ID3D[1],
                    cube.ID3D[2] - far.ID3D[2],
                    cube.ID3D[3] - far.ID3D[3],
                )
                R = norm(off) * leaf_level.cubeEdgel
                kR_min = min(kR_min, k * R)
                n_far += 1
            end
        end
        if n_far > 0 && kR_min < 0.55 * leaf_level.L
            @warn "MLFMA 叶层远场对最小距离不足：kR_min = $(round(kR_min, digits = 2)) < " *
                  "0.55·L（L = $(leaf_level.L)）。对角实谱 M2L 在该距离下精度显著下降" *
                  "（近场倏逝波分量未被实角谱覆盖）。当前 near_range = $(repr(nr_tree))；" *
                  "建议增大 near_range 或减小叶层尺寸，使叶层最小 far 偏移满足 " *
                  "kR_min ≥ 0.55·L。"
        end
    end

    @info "Octree built successfully."
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
    near_range::Union{Int,Nothing} = nothing,
    interp_method::Val = Val(:Lagrange2Step),
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

    # Search neighbors (tree-uniform radius — see build_octree / issue #22 #3)
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

    L, poles = if interp_method == Val(:LbTrained1Step)
        levelIntegralInfoCal(cubeEdgel, Val(:LbTrained1Step); λ = λ)
    elseif interp_method == Val(:FFTSpectral)
        levelIntegralInfoCal(cubeEdgel, Val(:FFTSpectral); λ = λ, L_min = L_min)
    else
        levelIntegralInfoCal(cubeEdgel; λ = λ, L_min = L_min)
    end

    level = LevelInfo{Int,FT,interp_type(poles)}()
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
    near_range::Union{Int,Nothing} = nothing,
    interp_method::Val = Val(:Lagrange2Step),
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

    # Truncation number first: the tree-uniform near-field radius is derived
    # from it in build_octree (issue #22, problem 3)
    L, poles = if interp_method == Val(:LbTrained1Step)
        levelIntegralInfoCal(cubeEdgel, Val(:LbTrained1Step); λ = λ)
    elseif interp_method == Val(:FFTSpectral)
        levelIntegralInfoCal(cubeEdgel, Val(:FFTSpectral); λ = λ, L_min = L_min)
    else
        levelIntegralInfoCal(cubeEdgel; λ = λ, L_min = L_min)
    end

    # Search neighbors (tree-uniform radius — see build_octree / issue #22 #3)
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

    level = LevelInfo{Int,FT,interp_type(poles)}()
    level.ID = levelID
    level.isleaf = false
    level.L = L
    level.nCubes = nCubes
    level.cubes = cubesInfo
    level.cubeEdgel = cubeEdgel
    level.poles = poles

    return level, kidCubesSorted
end

function searchNearCubes(
    cubesID3D::Matrix{IT},
    levelID::Integer;
    near_range::Int = 4,
) where {IT<:Integer}
    nCubes = size(cubesID3D, 1)
    maxCubes1D = 2^levelID
    cubesID1D = [
        ((row[1] - 1) * maxCubes1D^2 + (row[2] - 1) * maxCubes1D + row[3]) for
        row in eachrow(cubesID3D)
    ]

    neighbors = Vector{Vector{IT}}(undef, nCubes)

    for iCube = 1:nCubes
        cubeID3D = cubesID3D[iCube, :]
        neighborsOffsets =
            [(-near_range:near_range), (-near_range:near_range), (-near_range:near_range)]

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
