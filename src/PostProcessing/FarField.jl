module FarField

using StaticArrays
using LinearAlgebra
using Base.Threads

using ...Geometry
using ...BasisFunctions
using ...Utilities.Parameters
using ..CurrentOnGeos
using ..RadiationIntegral: r̂θϕInfo, raditionalIntegralNθϕCal, ∠Info
using ...CoreModule: num_elements

export farField

"""
    farField(θs_obs, ϕs_obs, ICoeff, basis, source)

Calculate Far Field pattern.
"""
function farField(θs_obs::Vector{FT}, ϕs_obs::Vector{FT}, 
                  ICoeff::Vector{CT}, basis::RWGBasis{IT, FT}, 
                  source) where {IT<:Integer, FT<:Real, CT<:Complex{FT}}
    
    # Calculate currents on geometry
    Jtris = geoElectricJCal(ICoeff, basis)
    
    k0 = get_k0()
    eta0 = get_eta0()
    jk0 = im * k0
    div4pi = 1.0 / (4 * pi)
    
    Nθ_obs = length(θs_obs)
    Nϕ_obs = length(ϕs_obs)
    nobs = Nθ_obs * Nϕ_obs
    
    # Precompute angle info
    θsobsInfo = [∠Info(θ) for θ in θs_obs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs_obs]
    
    # Create grid of observation directions
    r̂θsϕs = [r̂θϕInfo(θ, ϕ) for θ in θsobsInfo, ϕ in ϕsobsInfo]
    r̂θsϕs_flat = vec(r̂θsϕs)
    
    farEθsϕs = zeros(Complex{FT}, 2, Nθ_obs, Nϕ_obs)
    farE_flat = reshape(farEθsϕs, 2, nobs)
    
    # Construct TriangleInfo for radiation integral
    mesh = basis.mesh
    ntri = num_elements(mesh)
    trianglesInfo = Vector{TriangleInfo{IT, FT}}(undef, ntri)
    for t in 1:ntri
        trianglesInfo[t] = TriangleInfo(mesh, t)
    end
    
    @threads for ii in 1:nobs
        r_info = r̂θsϕs_flat[ii]
        
        # Calculate Radiation Integral
        Nθϕ = raditionalIntegralNθϕCal(r_info, trianglesInfo, Jtris)
        Nθϕ = raditionalIntegralNθϕCal(r_info, trianglesInfo, Jtris)
        
        # Calculate Far Field
        # E_far = (-j k0 eta0 / 4pi) * N
        # Legacy code uses div4pi.
        factor = -jk0 * eta0 * div4pi
        farEθϕ = factor .* Nθϕ
        
        # Add source field if applicable
        # farEθϕ .+= sourceFarEfield(source, r_info)
        
        farE_flat[:, ii] = farEθϕ
    end
    
    return farEθsϕs
end

end
