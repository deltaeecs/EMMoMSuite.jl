module Interpolation

using StaticArrays
using SparseArrays
using FastGaussQuadrature
using LinearAlgebra
using ....Geometry

export AbstractPolesInfo, AbstractInterpInfo
export GLPolesInfo, LagrangeInterpInfo
export levelIntegralInfoCal, interpolationCSCMatCal, truncation_kernel
export interpolate, anterpolate, interpolate!, anterpolate!

abstract type AbstractPolesInfo{FT<:AbstractFloat} end
abstract type AbstractInterpInfo{IT<:Integer, FT<:Real} end

function interpolate end
function anterpolate end
function interpolate! end
function anterpolate! end

const NBDIGITS = 9.0 # Default value

"""
    GLPolesInfo{FT}

Gauss-Legendre Poles Information.
"""
struct GLPolesInfo{FT<:Real} <: AbstractPolesInfo{FT}
    Xθs::Vector{FT}
    Xϕs::Vector{FT}
    Wθϕs::Vector{FT}
    r̂sθsϕs::Vector{r̂θϕInfo{FT}}
end

"""
    LagrangeInterpInfo{IT, FT}

Lagrange Interpolation Information.
"""
mutable struct LagrangeInterpInfo{IT<:Integer, FT<:Real} <: AbstractInterpInfo{IT, FT}
    θCSC    ::SparseMatrixCSC{FT, IT}
    ϕCSC    ::SparseMatrixCSC{FT, IT}
    θCSCT   ::SparseMatrixCSC{FT, IT}
    ϕCSCT   ::SparseMatrixCSC{FT, IT}
    
    function LagrangeInterpInfo{IT, FT}() where {IT<:Integer, FT<:Real}
        new{IT, FT}()
    end
end

function integral1DXW(lb::FT, hb::FT, Nsample::IT, mod::Symbol) where{IT<:Integer, FT<:Real}
    X, W = zeros(FT, Nsample), zeros(FT, Nsample)

    if mod == :uni # Uniform integration
        dl = (hb-lb)/Nsample
        for j in 1:Nsample
            X[j] = lb + (j-1)*dl + dl/2
            W[j] = dl
        end
    elseif mod == :glq # Gauss-Legendre integration
        Dx      =   0.5 * (hb - lb)
        center  =   0.5 * (hb + lb)
        XGL, WGL    =   gausslegendre(Nsample)
        X   .=   center .+ Dx .* XGL
        W   .=   abs(Dx) .* WGL
    else
        error("Only :uni and :glq modes are supported")
    end
    return X, W
end

function octreeXWNCal(lb::FT, hb::FT, L::IT, mod::Symbol) where{IT<:Integer, FT<:Real}
    N = if mod == :uni
        2*(L + 1)
    elseif mod == :glq
        L + 1
    else
        error("Only :uni and :glq modes are supported")
    end
    Xs, Ws = integral1DXW(lb, hb, N, mod)
    return Xs, Ws
end

function truncation_kernel(rel_l)
    return 2π*rel_l*sqrt(3) + 2.16*NBDIGITS^(2.0/3.0)*(2π*rel_l)^(1/3)
end

function truncationLCal(cubel::FT; λ=1.0) where {FT<:Real}
    rel_l = cubel/λ
    L = floor(Int, truncation_kernel(rel_l))
    return L
end

function levelIntegralInfoCal(levelCubeEdgel::FT; λ=1.0) where{FT<:Real}
    ## Calculate truncation number
    L = truncationLCal(levelCubeEdgel; λ=λ)
    
    ## Integration points and weights
    # Theta direction (Gauss-Legendre)
    Xcosθs, Wθs = octreeXWNCal(one(FT), -one(FT), L, :glq)
    Xθs = acos.(Xcosθs)
    
    # Phi direction (Uniform)
    Xϕs, Wϕs = octreeXWNCal(zero(FT), convert(FT, 2π), L, :uni)
    
    # Calculate all poles info
    r̂sθsϕs = [r̂θϕInfo(θ, ϕ) for ϕ in Xϕs for θ in Xθs]
    Wθϕs = [Wθ * Wϕ for Wϕ in Wϕs for Wθ in Wθs]

    Poles = GLPolesInfo{FT}(Xθs, Xϕs, Wθϕs, r̂sθsϕs)
    
    return L, Poles
end

"""
    interpolationCSCMatCal(pLevelPoles::GLPolesInfo{FT}, tLevelPoles::GLPolesInfo{FT}, nlocalInterp::IT)

Calculate sparse interpolation matrices from parent level to current level (or vice versa).
"""
function interpolationCSCMatCal(pLevelPoles::GLPolesInfo{FT}, tLevelPoles::GLPolesInfo{FT}, nlocalInterp::IT) where {IT<:Integer, FT<:Real}
    
    # Parent and Current level theta coordinates
    pXθs    =   pLevelPoles.Xθs
    tXθs    =   tLevelPoles.Xθs
    npXθs   =   length(pXθs)
    ntXθs   =   length(tXθs)
    
    # Parent and Current level phi coordinates
    pXϕs    =   pLevelPoles.Xϕs
    tXϕs    =   tLevelPoles.Xϕs
    npXϕs   =   length(pXϕs)
    ntXϕs   =   length(tXϕs)

    ################################################################
    # Theta direction

    nlocalInterpTheta = nlocalInterp
    if nlocalInterpTheta > npXθs
        nlocalInterpTheta = npXθs
    end
    
    pθsIntθs    =   cooraInCoorb(pXθs, tXθs)

    interWθs    =   ones(FT, (nlocalInterpTheta, npXθs))
    interIDθs   =   zeros(IT, (nlocalInterpTheta, npXθs))
    RelativeOffsets =   (1:nlocalInterpTheta) .- nlocalInterpTheta ÷ 2
    
    @inbounds for ipXθs in 1:npXθs
        pθIntθ  =   pθsIntθs[ipXθs]
        iInterIDθs   =   pθIntθ .+  RelativeOffsets
        interIDθs[:,ipXθs] .=  iInterIDθs
    end

    @inbounds for i in 1:nlocalInterpTheta
        θsInterp    =   [pickθ(idx, tXθs) for idx in interIDθs[i, :]]
        sinHalfDiffθpLevel  =   sin.((pXθs .- θsInterp) ./ 2)
        for j in 1:nlocalInterpTheta
            if i != j
                θsInterpLocal       =   [pickθ(idx, tXθs) for idx in interIDθs[j, :]]
                sinHalfDiffθtLevel  =   sin.((θsInterpLocal .- θsInterp) ./ 2)
                interWθs[j, :]     .*=   sinHalfDiffθpLevel ./ sinHalfDiffθtLevel
            end
        end
    end

    @inbounds for i in 1:npXθs
        interWθs[:, i] ./=  sum(interWθs[:, i])
    end

    ################################################################
    # Phi direction

    nlocalInterpPhi = nlocalInterp
    if nlocalInterpPhi > npXϕs
        nlocalInterpPhi = npXϕs
    end
    
    pϕsIntϕs    =   cooraInCoorb(pXϕs, tXϕs)

    interWϕs    =   ones(FT, (nlocalInterpPhi, npXϕs))
    interIDϕs   =   zeros(IT, (nlocalInterpPhi, npXϕs))
    RelativeOffsets =   (1:nlocalInterpPhi) .- nlocalInterpPhi ÷ 2
    
    @inbounds for ipXϕs in 1:npXϕs
        pϕIntϕ  =   pϕsIntϕs[ipXϕs]
        iInterIDϕs   =   collect(pϕIntϕ .+  RelativeOffsets)
        interIDϕs[:,ipXϕs] .=  iInterIDϕs
    end

    @inbounds for i in 1:nlocalInterpPhi
        ϕsInterp    =   [pickϕ(idx, tXϕs) for idx in interIDϕs[i,:]]
        sinHalfDiffϕpLevel  =   sin.((pXϕs .- ϕsInterp) ./ 2)
        for j in 1:nlocalInterpPhi
            if i != j
                ϕsInterpLocal       =   [pickϕ(idx, tXϕs) for idx in interIDϕs[j,:]]
                sinHalfDiffϕtLevel  =   sin.((ϕsInterpLocal .- ϕsInterp) ./ 2)
                interWϕs[j, :]     .*=   sinHalfDiffϕpLevel ./ sinHalfDiffϕtLevel
            end
        end
    end

    @inbounds for i in 1:npXϕs
        interWϕs[:, i] ./=  sum(view(interWϕs, :, i))
    end

    ################################################################
    # Construct Sparse Matrices (2-step interpolation)
    
    npSample    =   npXθs*npXϕs
    ntSample    =   ntXθs*ntXϕs
    ntempSample =   ntXθs*npXϕs
    
    ### Step 1: Phi direction
    tSampleIndexes  =   reshape(collect(1:ntSample), ntXθs, ntXϕs)
    interIDGlobalϕs =   repeat(interIDϕs, inner = (1, ntXθs))
    
    for itθ in 1:ntXθs
        for ipϕ in 1:npXϕs
            interIDGlobalϕs[:, (ipϕ - 1)*ntXθs + itθ]    .=  [pickCycleVec(interIDϕ, tSampleIndexes[itθ,:])  for interIDϕ in interIDϕs[:,ipϕ]]
        end
    end
    
    rawIDϕs =   repeat(collect(IT, 1:ntempSample); inner = nlocalInterp)
    interWGlobalϕs  =   repeat(interWϕs, inner = (1, ntXθs))
    interpϕCSC  =   sparse(rawIDϕs, view(interIDGlobalϕs, :), view(interWGlobalϕs, :))


    #### Step 2: Theta direction
    for i in eachindex(interIDθs)
        ((interIDθs[i] < 1) | (interIDθs[i] > ntXθs)) && (interWθs[i] *= -1) 
    end
    
    tempSampleIndexes  =   reshape(collect(1:ntempSample), ntXθs, npXϕs)
    interIDGlobalθs =   repeat(interIDθs, outer = (1, npXϕs))
    
    halfnpϕ = npXϕs ÷ 2
    for ipϕ in 1:npXϕs
        for ipθ in 1:npXθs
            inGlobalIDs     =   zeros(IT, nlocalInterp)
            for jInter in 1:nlocalInterp
                interIDθ    =   interIDθs[jInter, ipθ]
                targetIdxInTempSampleIndexes  =  [interIDθ, ipϕ] 
                
                if (interIDθ < 1) | (interIDθ > ntXθs)
                    if interIDθ < 1
                        targetIdxInTempSampleIndexes[1] = -interIDθ + 1
                    elseif  interIDθ > ntXθs
                        targetIdxInTempSampleIndexes[1] = 2ntXθs + 1 -interIDθ
                    end
                    
                    if ipϕ <= halfnpϕ
                        targetIdxInTempSampleIndexes[2] =  ipϕ + halfnpϕ
                    else
                        targetIdxInTempSampleIndexes[2] =  ipϕ - halfnpϕ
                    end
                end
                inGlobalIDs[jInter] = tempSampleIndexes[targetIdxInTempSampleIndexes...]
            end
            interIDGlobalθs[:, (ipϕ - 1)*npXθs + ipθ]    .=  inGlobalIDs
        end
    end
    
    rawIDθs =   repeat(collect(IT, 1:npSample); inner = nlocalInterp)
    interWGlobalθs  =   repeat(interWθs; outer = (1, npXϕs))
    interpθCSC  =   sparse(rawIDθs, view(interIDGlobalθs, :), view(interWGlobalθs, :))

    dropzeros!(interpθCSC)
    dropzeros!(interpϕCSC)

    # Create struct
    info = LagrangeInterpInfo{IT, FT}()
    info.θCSC = interpθCSC
    info.ϕCSC = interpϕCSC
    info.θCSCT = sparse(transpose(interpθCSC))
    info.ϕCSCT = sparse(transpose(interpϕCSC))
    
    return info
end

function cooraInCoorb(coora::Vector{T}, coorb::Vector{T}) where {T<:Number}
    targetIDs = fill!(similar(coora, Int), 0)
    @inbounds for i in eachindex(coora)
        a = coora[i]
        for j in eachindex(coorb)
            b = coorb[j]
            if a >= b
                targetIDs[i] = j
                continue
            end
        end
    end
    return targetIDs
end

function pickϕ(index::Integer, ϕs::Vector{TT}) where {TT<:Real}
    re = zero(TT)
    ln = length(ϕs)
    if  index < 1
        re = (ϕs[index + ln] - 2π)
    elseif index > ln
        re = ϕs[index - ln] + 2π
    else
        re = ϕs[index]
    end
    return re
end

function pickCycleVec(index::Integer, cycleVec::Vector{T}) where {T<:Real}
    ln = length(cycleVec)
    if  index < 1
        re = cycleVec[index + ln]
    elseif index > ln
        re = cycleVec[index - ln]
    else
        re = cycleVec[index]
    end
    return re
end

function pickθ(index::Integer, θs::Vector{T}) where {T<:Real}
    re = zero(T)
    ln = length(θs)
    if  index < 1
        re = -θs[-index + 1]
    elseif index > ln
        re = 2pi - θs[2ln + 1 - index]
    else
        re = θs[index]
    end
    return re
end


end
