module Aggregation

using LinearAlgebra
using StaticArrays
using ....CoreModule
using ....CoreModule: Constants
using ....Geometry
using ....BasisFunctions
using ....IntegralEquations
using ..Level
using ..Octree

using ....IntegralEquations.Impedance: get_triangle_info, get_triangles_info
using ....IntegralEquations.VEFIEModule: get_tetrahedra_info

export aggregate!, aggregate_leaf!, aggregate_upward!

function get_k(op::AbstractIntegralOperator)
    if hasfield(typeof(op), :k)
        return op.k
    elseif hasfield(typeof(op), :efie)
        return op.efie.k
    elseif hasfield(typeof(op), :freq)
        return 2π * op.freq / Constants.c0
    else
        error("Operator $(typeof(op)) does not have field k or efie or freq")
    end
end

"""
    aggregate!(octree::OctreeInfo, basis::AbstractBasisFunction, operator::AbstractIntegralOperator, x::AbstractVector, sorted_ids::Vector{Int})

Compute aggregation for all levels using coefficients x.
"""
function aggregate!(
    octree::OctreeInfo,
    basis::AbstractBasisFunction,
    operator::AbstractIntegralOperator,
    x::AbstractVector,
    sorted_ids::Vector{Int},
)
    aggregate!(octree, AbstractBasisFunction[basis], [num_basis(basis)], operator, x, sorted_ids)
end

function aggregate!(
    octree::OctreeInfo,
    bases::Vector{<:AbstractBasisFunction},
    offsets::Vector{Int},
    operator::AbstractIntegralOperator,
    x::AbstractVector,
    sorted_ids::Vector{Int},
)
    # 1. Leaf level aggregation (Radiation pattern of basis functions)
    leafLevel = octree.levels[octree.nLevels]
    aggregate_leaf!(leafLevel, bases, offsets, operator, x, sorted_ids)

    # 2. Upward pass (Child to Parent)
    for levelID = (octree.nLevels-1):-1:2
        parentLevel = octree.levels[levelID]
        childLevel = octree.levels[levelID+1]
        aggregate_upward!(parentLevel, childLevel)
    end
end

function get_basis_index(global_idx::Int, offsets::Vector{Int})
    idx = searchsortedfirst(offsets, global_idx)
    local_idx = idx == 1 ? global_idx : global_idx - offsets[idx-1]
    return idx, local_idx
end

"""
    aggregate_leaf!(level::LevelInfo, bases::Vector{<:AbstractBasisFunction}, offsets::Vector{Int}, operator::AbstractIntegralOperator, x::AbstractVector, sorted_ids::Vector{Int})

Compute radiation patterns of basis functions at the leaf level.
"""
function aggregate_leaf!(
    level::LevelInfo,
    bases::Vector{<:AbstractBasisFunction},
    offsets::Vector{Int},
    operator::AbstractIntegralOperator,
    x::AbstractVector,
    sorted_ids::Vector{Int};
    cube_filter = nothing,
)
    FT = eltype(level.cubeEdgel)
    CT = Complex{FT}

    k = get_k(operator)
    JK = im * k

    poles = level.poles
    nPoles = length(poles.r̂sθsϕs)
    nCubes = level.nCubes

    if !isdefined(level, :aggS)
        level.aggS = zeros(CT, nPoles, 2, nCubes)
    else
        fill!(level.aggS, zero(CT))
    end

    # Precompute element info and quadrature
    element_infos = Any[]
    gqs = Any[]

    for b in bases
        if b isa RWGBasis
            push!(element_infos, get_triangles_info(b.mesh, b))
            push!(gqs, GaussQuadratureInfo(:Triangle, 3, FT))
        elseif b isa SWGBasis
            if hasfield(typeof(operator), :permittivities)
                push!(element_infos, get_tetrahedra_info(b.mesh, b, operator.permittivities))
            else
                # Fallback
                push!(
                    element_infos,
                    get_tetrahedra_info(b.mesh, b, fill(1.0 + 0im, num_elements(b.mesh))),
                )
            end
            push!(gqs, GaussQuadratureInfo(:Tetrahedron, 5, FT))
        else
            error("Unsupported basis type")
        end
    end

    poles_r̂ = [p.r̂ for p in poles.r̂sθsϕs]
    poles_θhat = [p.θhat for p in poles.r̂sθsϕs]
    poles_ϕhat = [p.ϕhat for p in poles.r̂sθsϕs]

    # Loop over cubes
    Threads.@threads for iCube = 1:nCubes
        cube_filter !== nothing && !cube_filter(iCube) && continue
        cube = level.cubes[iCube]
        if isempty(cube.bfInterval)
            continue
        end

        cubeCenter = cube.center
        start_idx = first(cube.bfInterval)
        end_idx = last(cube.bfInterval)

        for i_sorted = start_idx:end_idx
            i_orig = sorted_ids[i_sorted]

            # Identify basis
            b_idx, i_local = get_basis_index(i_orig, offsets)
            basis = bases[b_idx]
            bf = basis.functions[i_local]

            coef = x[i_orig]
            if abs(coef) < 1e-12
                continue
            end

            elem_info = element_infos[b_idx]
            gq = gqs[b_idx]

            if basis isa RWGBasis
                add_radiation_pattern_rwg!(
                    level.aggS,
                    iCube,
                    basis,
                    bf,
                    elem_info,
                    gq,
                    k,
                    cubeCenter,
                    poles_r̂,
                    poles_θhat,
                    poles_ϕhat,
                    coef,
                )
            elseif basis isa SWGBasis
                add_radiation_pattern_swg!(
                    level.aggS,
                    iCube,
                    basis,
                    bf,
                    elem_info,
                    gq,
                    k,
                    cubeCenter,
                    poles_r̂,
                    poles_θhat,
                    poles_ϕhat,
                    coef,
                )
            end
        end
    end
end

function add_radiation_pattern_rwg!(
    aggS,
    iCube,
    basis,
    bf,
    elem_info,
    gq,
    k,
    cubeCenter,
    poles_r̂,
    poles_θhat,
    poles_ϕhat,
    coef,
)
    JK = im * k
    n_qp = length(gq.weight)
    nPoles = length(poles_r̂)

    for i_supp = 1:2
        tri_idx = bf.support[i_supp]
        if tri_idx == 0
            continue
        end

        tri = elem_info[tri_idx]
        tri_vertices = tri.vertices

        local_edge = bf.local_edge_idx[i_supp]
        opp_vertex_idx = local_edge
        v_opp = tri_vertices[:, opp_vertex_idx]

        sign = bf.signs[i_supp]

        for i_qp = 1:n_qp
            L = gq.coordinate[:, i_qp]
            r = tri_vertices * L
            rho = r - v_opp
            r_local = r - cubeCenter

            factor_vec = sign * bf.edge_length / 2 * gq.weight[i_qp] * coef

            for iPole = 1:nPoles
                r̂ = poles_r̂[iPole]
                phase = exp(JK * dot(r̂, r_local))
                vec = rho * factor_vec * phase

                aggS[iPole, 1, iCube] += dot(poles_θhat[iPole], vec)
                aggS[iPole, 2, iCube] += dot(poles_ϕhat[iPole], vec)
            end
        end
    end
end

function add_radiation_pattern_swg!(
    aggS,
    iCube,
    basis,
    bf,
    elem_info,
    gq,
    k,
    cubeCenter,
    poles_r̂,
    poles_θhat,
    poles_ϕhat,
    coef,
)
    JK = im * k
    n_qp = length(gq.weight)
    nPoles = length(poles_r̂)

    for i_supp = 1:2
        tet_idx = bf.support[i_supp]
        if tet_idx == 0
            continue
        end

        tet = elem_info[tet_idx]
        tet_vertices = tet.vertices

        # Kappa
        kappa = tet.κ

        local_face = bf.local_face_idx[i_supp]
        v_free = tet_vertices[:, local_face]

        sign = bf.signs[i_supp]

        for i_qp = 1:n_qp
            L = gq.coordinate[:, i_qp]
            r = tet_vertices * L
            rho = sign * (r - v_free)
            r_local = r - cubeCenter

            factor_vec = bf.area / 3 * gq.weight[i_qp] * coef * kappa

            for iPole = 1:nPoles
                r̂ = poles_r̂[iPole]
                phase = exp(JK * dot(r̂, r_local))
                vec = rho * factor_vec * phase

                aggS[iPole, 1, iCube] += dot(poles_θhat[iPole], vec)
                aggS[iPole, 2, iCube] += dot(poles_ϕhat[iPole], vec)
            end
        end
    end
end

"""
    aggregate_upward!(parentLevel::LevelInfo, childLevel::LevelInfo)

Aggregate radiation patterns from child level to parent level (Interpolation + Phase Shift).
"""
function aggregate_upward!(parentLevel::LevelInfo, childLevel::LevelInfo)
    FT = eltype(parentLevel.cubeEdgel)
    CT = Complex{FT}

    nPolesParent = length(parentLevel.poles.r̂sθsϕs)
    nCubesParent = parentLevel.nCubes

    if !isdefined(parentLevel, :aggS)
        parentLevel.aggS = zeros(CT, nPolesParent, 2, nCubesParent)
    else
        fill!(parentLevel.aggS, zero(CT))
    end

    parentAggS = parentLevel.aggS
    childAggS = childLevel.aggS

    interp = childLevel.interpWθϕ
    θCSC = interp.θCSC
    ϕCSC = interp.ϕCSC

    phaseShift = parentLevel.phaseShiftFromKids

    Threads.@threads for iCube = 1:nCubesParent
        parentCube = parentLevel.cubes[iCube]

        for iKid = 1:length(parentCube.kidsInterval)
            childID = parentCube.kidsInterval[iKid]
            childIn8 = parentCube.kidsIn8[iKid]

            aggChild = view(childAggS, :, :, childID)

            aggInterpPhi = ϕCSC * aggChild
            aggInterp = θCSC * aggInterpPhi

            shift = view(phaseShift, :, childIn8)

            for pol = 1:2
                @views parentAggS[:, pol, iCube] .+= shift .* aggInterp[:, pol]
            end
        end
    end
end

end
