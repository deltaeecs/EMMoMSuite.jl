using StaticArrays
using LinearAlgebra

"""
    sphere2cart(r, θ, ϕ)

Convert spherical coordinates to Cartesian coordinates.
"""
function sphere2cart(r, θ, ϕ)
    sinθ, cosθ = sincos(θ)
    sinϕ, cosϕ = sincos(ϕ)
    x = r * sinθ * cosϕ
    y = r * sinθ * sinϕ
    z = r * cosθ
    return SVector(x, y, z)
end

"""
    r̂θϕInfo{FT}

Struct to hold spherical coordinate unit vectors.
"""
struct r̂θϕInfo{FT}
    r̂::SVector{3,FT}
    θhat::SVector{3,FT}
    ϕhat::SVector{3,FT}
end

"""
    r̂θϕInfo(θ::FT, ϕ::FT)

Constructor from angles.
"""
function r̂θϕInfo(θ::FT, ϕ::FT) where {FT}
    sinθ, cosθ = sincos(θ)
    sinϕ, cosϕ = sincos(ϕ)

    r̂ = SVector(sinθ * cosϕ, sinθ * sinϕ, cosθ)
    θhat = SVector(cosθ * cosϕ, cosθ * sinϕ, -sinθ)
    ϕhat = SVector(-sinϕ, cosϕ, zero(FT))

    return r̂θϕInfo{FT}(r̂, θhat, ϕhat)
end

function r̂θϕInfo{FT}(r̂::AbstractVector) where {FT}
    # Calculate theta, phi from r_hat and construct local basis
    r = norm(r̂)
    if r ≈ 0
        return r̂θϕInfo{FT}(SVector{3,FT}(0, 0, 0), SVector{3,FT}(0, 0, 0), SVector{3,FT}(0, 0, 0))
    end
    r̂_norm = SVector{3,FT}(r̂) / r
    θ = acos(r̂_norm[3])
    ϕ = atan(r̂_norm[2], r̂_norm[1])
    return r̂θϕInfo(θ, ϕ)
end

r̂θϕInfo(r̂::AbstractVector{FT}) where {FT} = r̂θϕInfo{FT}(r̂)

"""
    globalObs2LocalObs(r̂θϕs_obs, l2gRot)

Transform observation angles from global to local coordinates.
"""
function globalObs2LocalObs(r̂θϕs_obs::Matrix{r̂θϕInfo{FT}}, l2gRot::SMatrix{3,3,FT}) where {FT}
    r̂ObsLocal = [l2gRot' * r̂θϕ_obs.r̂ for r̂θϕ_obs in r̂θϕs_obs]
    return [r̂θϕInfo{FT}(r̂) for r̂ in r̂ObsLocal]
end

"""
    localObs2GlobalObs(r̂θϕs_obs, l2gRot)

Transform observation angles from local to global coordinates.
"""
function localObs2GlobalObs(r̂θϕs_obs::Matrix{r̂θϕInfo{FT}}, l2gRot::SMatrix{3,3,FT}) where {FT}
    r̂ObsGlobal = [l2gRot * r̂θϕ_obs.r̂ for r̂θϕ_obs in r̂θϕs_obs]
    return [r̂θϕInfo{FT}(r̂) for r̂ in r̂ObsGlobal]
end
