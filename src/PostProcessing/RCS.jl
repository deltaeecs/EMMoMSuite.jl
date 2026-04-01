module RCS

using StaticArrays
using LinearAlgebra
using Base.Threads

using ...Geometry
using ...BasisFunctions
using ...Utilities.Parameters
using ..CurrentOnGeos
using ..RadiationIntegral:
    r̂θϕInfo,
    raditionalIntegralNθϕCal,
    radiation_integral_rwg,
    radiation_integral_swg,
    radiation_integral_pwc,
    radiation_integral_pwc_hex,
    radiation_integral_rbf,
    ∠Info
using ...CoreModule: num_elements, Constants

export radarCrossSection

"""
    radarCrossSection(θs_obs, ϕs_obs, ICoeff, trianglesInfo, BFT)

Calculate the Radar Cross Section (RCS) of the target.

# Mathematical Formulation

The RCS (\$\\sigma\$) is defined as the limit of the area that intercepts the incident power to produce the scattered power density at a distance \$r\$:

```math
\\sigma(\\theta, \\phi) = \\lim_{r \\to \\infty} 4\\pi r^2 \\frac{|\\mathbf{E}^{scat}(\\theta, \\phi)|^2}{|\\mathbf{E}^{inc}|^2}
```

Using the radiation vector \$\\mathbf{N}\$:
```math
\\mathbf{N}(\\theta, \\phi) = \\int_S \\mathbf{J}(\\mathbf{r}') e^{jk \\hat{r} \\cdot \\mathbf{r}'} dS'
```

The scattered far-field is:
```math
\\mathbf{E}^{scat}_{far} \\approx -j k \\eta \\frac{e^{-jkr}}{4\\pi r} (N_\\theta \\hat{\\theta} + N_\\phi \\hat{\\phi})
```

Substituting this into the RCS definition (assuming unit incident field):
```math
\\sigma = \\frac{k^2 \\eta^2}{4\\pi} (|N_\\theta|^2 + |N_\\phi|^2)
```

# Arguments
- `θs_obs`: Vector of observation angles \$\\theta\$ (radians).
- `ϕs_obs`: Vector of observation angles \$\\phi\$ (radians).
- `ICoeff`: Vector of current coefficients (solution from MoM).
- `trianglesInfo`: Mesh geometry information.
- `BFT`: Basis function type (e.g., `RWG`).

# Returns
- `RCSθsϕs`: Array of component RCS (Theta/Phi).
- `RCS_total`: Total RCS (linear scale, \$\\text{m}^2\$).
- `RCS_dB`: Total RCS (dBsm, \$10 \\log_{10}(\\sigma)\$).
"""
function radarCrossSection(
    θs_obs::Vector{FT},
    ϕs_obs::Vector{FT},
    ICoeff::Vector{CT},
    basis::RWGBasis{IT,FT},
) where {IT<:Integer,FT<:Real,CT<:Complex{FT}}

    # Calculate currents on geometry
    Jtris = geoElectricJCal(ICoeff, basis)

    k0 = get_k0()
    eta0 = get_eta0()

    Nθ_obs = length(θs_obs)
    Nϕ_obs = length(ϕs_obs)
    nobs = Nθ_obs * Nϕ_obs

    # Precompute angle info
    θsobsInfo = [∠Info(θ) for θ in θs_obs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs_obs]

    # Create grid of observation directions
    r̂θsϕs = [r̂θϕInfo(θ, ϕ) for θ in θsobsInfo, ϕ in ϕsobsInfo]

    # Flatten for iteration
    r̂θsϕs_flat = vec(r̂θsϕs)

    # Use 2D array for linear indexing
    RCS_flat = zeros(FT, 2, nobs)

    # Construct TriangleInfo for radiation integral
    # (kept for backward compatibility, not used in radiation_integral_rwg path)
    mesh = basis.mesh
    ntri = num_elements(mesh)

    # We can use @threads for parallel execution
    @threads for ii = 1:nobs
        r_info = r̂θsϕs_flat[ii]

        # Calculate Radiation Integral
        # Nθϕ = raditionalIntegralNθϕCal(r_info, trianglesInfo, Jtris)
        Nθϕ = radiation_integral_rwg(r_info, basis, ICoeff)

        # Calculate RCS
        # Formula: (k0 * eta0)^2 / (4 * pi) * |N|^2
        factor = (k0 * eta0)^2 / (4 * pi)
        RCS_flat[1, ii] = factor * abs2(Nθϕ[1]) # Theta component
        RCS_flat[2, ii] = factor * abs2(Nθϕ[2]) # Phi component

        # Debug
        # println("ii: ", ii, " N: ", Nθϕ, " RCS: ", RCS_flat[:, ii])
    end

    # Reshape for output
    RCSθsϕs = reshape(RCS_flat, 2, Nθ_obs, Nϕ_obs)

    # Calculate total RCS
    RCS_total = RCSθsϕs[1, :, :] + RCSθsϕs[2, :, :]

    # Convert to dB
    RCS_dB = 10 .* log10.(RCS_total)

    return RCSθsϕs, RCS_total, RCS_dB
end

"""
    radarCrossSection(θs_obs, ϕs_obs, ICoeff, basis, permittivities)

Calculate the Radar Cross Section (RCS) for VEFIE using SWG basis.
"""
function radarCrossSection(
    θs_obs::Vector{FT},
    ϕs_obs::Vector{FT},
    ICoeff::Vector{CT},
    basis::SWGBasis{IT,FT},
    permittivities::Vector{CT},
) where {IT<:Integer,FT<:Real,CT<:Complex{FT}}

    Jtetras = geoElectricJCal(ICoeff, basis, permittivities)

    k0 = get_k0()
    eta0 = get_eta0()

    Nθ_obs = length(θs_obs)
    Nϕ_obs = length(ϕs_obs)
    nobs = Nθ_obs * Nϕ_obs

    # Precompute angle info
    θsobsInfo = [∠Info(θ) for θ in θs_obs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs_obs]

    # Create grid of observation directions
    r̂θsϕs = [r̂θϕInfo(θ, ϕ) for θ in θsobsInfo, ϕ in ϕsobsInfo]
    r̂θsϕs_flat = vec(r̂θsϕs)

    # Use 2D array for linear indexing
    RCS_flat = zeros(FT, 2, nobs)

    # We can use @threads for parallel execution
    @threads for ii = 1:nobs
        r_info = r̂θsϕs_flat[ii]

        # Legacy parity: precompute one weighted current vector per tetrahedron.
        Nθϕ = raditionalIntegralNθϕCal(r_info, basis, Jtetras)

        # Calculate RCS
        # Formula: (k0 * eta0)^2 / (4 * pi) * |N|^2
        # Note: For VEFIE, N is defined as Integral(J_eq * exp).
        # J_eq has units of A/m^2. Volume integral gives A*m.
        # N has units of A*m.
        # E_scat ~ -jk * eta * N / (4*pi*r)
        # RCS = 4*pi*r^2 * |E|^2 / |E_inc|^2
        # RCS = k^2 * eta^2 / (4*pi) * |N|^2
        # This matches the formula used for Surface.

        factor = (k0 * eta0)^2 / (4 * pi)
        RCS_flat[1, ii] = factor * abs2(Nθϕ[1]) # Theta component
        RCS_flat[2, ii] = factor * abs2(Nθϕ[2]) # Phi component
    end

    # Reshape for output
    RCSθsϕs = reshape(RCS_flat, 2, Nθ_obs, Nϕ_obs)

    # Calculate total RCS
    RCS_total = RCSθsϕs[1, :, :] + RCSθsϕs[2, :, :]

    # Convert to dB
    RCS_dB = 10 .* log10.(RCS_total)

    return RCSθsϕs, RCS_total, RCS_dB
end

"""
    radarCrossSection(θs_obs, ϕs_obs, ICoeff, basis::PWCBasis, permittivities)

Calculate the Radar Cross Section (RCS) for VEFIE using PWC basis.

# Legacy Parity
Matches `MoM_Kernels` RCS calculation for PWC tetrahedra (ConstBasisFunction).
"""
function radarCrossSection(
    θs_obs::Vector{FT},
    ϕs_obs::Vector{FT},
    ICoeff::Vector{CT},
    basis::PWCBasis{IT,FT},
    permittivities::Vector{CT},
) where {IT<:Integer,FT<:Real,CT<:Complex{FT}}

    k0 = get_k0()
    eta0 = get_eta0()

    Nθ_obs = length(θs_obs)
    Nϕ_obs = length(ϕs_obs)
    nobs = Nθ_obs * Nϕ_obs

    # Precompute angle info
    θsobsInfo = [∠Info(θ) for θ in θs_obs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs_obs]

    # Create grid of observation directions
    r̂θsϕs = [r̂θϕInfo(θ, ϕ) for θ in θsobsInfo, ϕ in ϕsobsInfo]
    r̂θsϕs_flat = vec(r̂θsϕs)

    RCS_flat = zeros(FT, 2, nobs)

    @threads for ii = 1:nobs
        r_info = r̂θsϕs_flat[ii]

        # Calculate Radiation Integral using PWC
        Nθϕ = radiation_integral_pwc(r_info, basis, ICoeff, permittivities)

        # RCS = k² η² / (4π) * |N|²
        factor = (k0 * eta0)^2 / (4 * pi)
        RCS_flat[1, ii] = factor * abs2(Nθϕ[1])
        RCS_flat[2, ii] = factor * abs2(Nθϕ[2])
    end

    RCSθsϕs = reshape(RCS_flat, 2, Nθ_obs, Nϕ_obs)
    RCS_total = RCSθsϕs[1, :, :] + RCSθsϕs[2, :, :]
    RCS_dB = 10 .* log10.(RCS_total)

    return RCSθsϕs, RCS_total, RCS_dB
end

"""
    radarCrossSection(θs_obs, ϕs_obs, ICoeff, basis::PWCHexBasis, permittivities)

Calculate the Radar Cross Section (RCS) for VEFIE using PWC hexahedra basis.
"""
function radarCrossSection(
    θs_obs::Vector{FT},
    ϕs_obs::Vector{FT},
    ICoeff::Vector{CT},
    basis::PWCHexBasis{IT,FT},
    permittivities::Vector{CT},
) where {IT<:Integer,FT<:Real,CT<:Complex{FT}}

    k0 = get_k0()
    eta0 = get_eta0()

    Nθ_obs = length(θs_obs)
    Nϕ_obs = length(ϕs_obs)
    nobs = Nθ_obs * Nϕ_obs

    θsobsInfo = [∠Info(θ) for θ in θs_obs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs_obs]
    r̂θsϕs = [r̂θϕInfo(θ, ϕ) for θ in θsobsInfo, ϕ in ϕsobsInfo]
    r̂θsϕs_flat = vec(r̂θsϕs)

    RCS_flat = zeros(FT, 2, nobs)

    @threads for ii = 1:nobs
        r_info = r̂θsϕs_flat[ii]
        Nθϕ = radiation_integral_pwc_hex(r_info, basis, ICoeff, permittivities)

        factor = (k0 * eta0)^2 / (4 * pi)
        RCS_flat[1, ii] = factor * abs2(Nθϕ[1])
        RCS_flat[2, ii] = factor * abs2(Nθϕ[2])
    end

    RCSθsϕs = reshape(RCS_flat, 2, Nθ_obs, Nϕ_obs)
    RCS_total = RCSθsϕs[1, :, :] + RCSθsϕs[2, :, :]
    RCS_dB = 10 .* log10.(RCS_total)

    return RCSθsϕs, RCS_total, RCS_dB
end

"""
    radarCrossSection(θs_obs, ϕs_obs, ICoeff, basis::RBFBasis, permittivities)

Calculate the Radar Cross Section (RCS) for VEFIE using RBF hexahedra basis.
"""
function radarCrossSection(
    θs_obs::Vector{FT},
    ϕs_obs::Vector{FT},
    ICoeff::Vector{CT},
    basis::RBFBasis{IT,FT},
    permittivities::Vector{CT},
) where {IT<:Integer,FT<:Real,CT<:Complex{FT}}

    k0 = get_k0()
    eta0 = get_eta0()

    Nθ_obs = length(θs_obs)
    Nϕ_obs = length(ϕs_obs)
    nobs = Nθ_obs * Nϕ_obs

    θsobsInfo = [∠Info(θ) for θ in θs_obs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs_obs]
    r̂θsϕs = [r̂θϕInfo(θ, ϕ) for θ in θsobsInfo, ϕ in ϕsobsInfo]
    r̂θsϕs_flat = vec(r̂θsϕs)

    RCS_flat = zeros(FT, 2, nobs)

    @threads for ii = 1:nobs
        r_info = r̂θsϕs_flat[ii]
        Nθϕ = radiation_integral_rbf(r_info, basis, ICoeff, permittivities)

        factor = (k0 * eta0)^2 / (4 * pi)
        RCS_flat[1, ii] = factor * abs2(Nθϕ[1])
        RCS_flat[2, ii] = factor * abs2(Nθϕ[2])
    end

    RCSθsϕs = reshape(RCS_flat, 2, Nθ_obs, Nϕ_obs)
    RCS_total = RCSθsϕs[1, :, :] + RCSθsϕs[2, :, :]
    RCS_dB = 10 .* log10.(RCS_total)

    return RCSθsϕs, RCS_total, RCS_dB
end

# ─── Phase 17.4: SCFIE unified RCS ────────────────────────────────────────────

"""
    radarCrossSection(θs_obs, ϕs_obs, ICoeff, basis_surf, basis_vol, permittivities)

SCFIE unified RCS: automatically splits the combined coefficient vector into
surface (RWG) and volume (SWG/PWC/PWCHex/RBF) parts, computes each radiation
integral, and returns the coherent superposition.

The coefficient vector layout must be:
  `ICoeff = [I_surf (1:n_surf); I_vol (n_surf+1:end)]`
where `n_surf = num_basis(basis_surf)`.
"""
function radarCrossSection(
    θs_obs::Vector{FT},
    ϕs_obs::Vector{FT},
    ICoeff::Vector{CT},
    basis_surf::RWGBasis{IT,FT},
    basis_vol::Union{SWGBasis{IT,FT},PWCBasis{IT,FT},PWCHexBasis{IT,FT},RBFBasis{IT,FT}},
    permittivities::Vector{CT},
) where {IT<:Integer,FT<:Real,CT<:Complex{FT}}

    n_surf = num_basis(basis_surf)
    n_vol = num_basis(basis_vol)
    @assert length(ICoeff) == n_surf + n_vol string(
        "SCFIE RCS: ICoeff length $(length(ICoeff)) ≠ n_surf($n_surf) + n_vol($n_vol)",
    )

    I_surf = ICoeff[1:n_surf]
    I_vol = ICoeff[n_surf+1:end]

    k0 = get_k0()
    eta0 = get_eta0()
    Nθ_obs = length(θs_obs)
    Nϕ_obs = length(ϕs_obs)
    nobs = Nθ_obs * Nϕ_obs

    θsobsInfo = [∠Info(θ) for θ in θs_obs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs_obs]
    r̂θsϕs = [r̂θϕInfo(θ, ϕ) for θ in θsobsInfo, ϕ in ϕsobsInfo]
    r̂θsϕs_flat = vec(r̂θsϕs)

    RCS_flat = zeros(FT, 2, nobs)

    if basis_vol isa SWGBasis
        _rad_vol = (r_info) -> radiation_integral_swg(r_info, basis_vol, I_vol, permittivities)
    elseif basis_vol isa PWCBasis
        _rad_vol = (r_info) -> radiation_integral_pwc(r_info, basis_vol, I_vol, permittivities)
    elseif basis_vol isa PWCHexBasis
        _rad_vol = (r_info) -> radiation_integral_pwc_hex(r_info, basis_vol, I_vol, permittivities)
    else  # RBFBasis
        _rad_vol = (r_info) -> radiation_integral_rbf(r_info, basis_vol, I_vol, permittivities)
    end

    @threads for ii = 1:nobs
        r_info = r̂θsϕs_flat[ii]

        Nθϕ_surf = radiation_integral_rwg(r_info, basis_surf, I_surf)
        Nθϕ_vol = _rad_vol(r_info)

        # Coherent superposition of surface and volume contributions
        Nθ = Nθϕ_surf[1] + Nθϕ_vol[1]
        Nϕ = Nθϕ_surf[2] + Nθϕ_vol[2]

        factor = (k0 * eta0)^2 / (4 * pi)
        RCS_flat[1, ii] = factor * abs2(Nθ)
        RCS_flat[2, ii] = factor * abs2(Nϕ)
    end

    RCSθsϕs = reshape(RCS_flat, 2, Nθ_obs, Nϕ_obs)
    RCS_total = RCSθsϕs[1, :, :] + RCSθsϕs[2, :, :]
    RCS_dB = 10 .* log10.(RCS_total)

    return RCSθsϕs, RCS_total, RCS_dB
end

"""
    radarCrossSection(θs_obs, ϕs_obs, ICoeff, basis, k0, eta0)

PMCHWT 双流远场 RCS。`ICoeff` 长度为 2N：前 N 个是等效电流 J 的系数，后 N 个是
等效磁流 M 的系数；两者均在同一 RWG basis 上展开。

散射远场（自由空间）公式（对应 e^{-jωt} 约定）：
```
E_θ = C × [η₀ N_θ + L_ϕ]
E_ϕ = C × [η₀ N_ϕ - L_θ]
```
其中
```
N = ∫ J e^{jk₀ r̂·r'} dS,   L = ∫ M e^{jk₀ r̂·r'} dS
```
RCS (m²) = k₀²/(4π) × (|E_θ|² + |E_ϕ|²)

# 参数
- `k0`:   自由空间波数 (rad/m)
- `eta0`: 自由空间波阻抗 (Ω)，通常 ≈ 377 Ω

# 返回
- `RCSθsϕs` : [2, Nθ, Nϕ]，分量 RCS (m²)
- `RCS_total`: [Nθ, Nϕ]，全极化 RCS (m²)
- `RCS_dB`  : [Nθ, Nϕ]，全极化 RCS (dBsm)
"""
function radarCrossSection(
    θs_obs::Vector{FT},
    ϕs_obs::Vector{FT},
    ICoeff::Vector{CT},
    basis::RWGBasis{IT,FT},
    k0::FT,
    eta0::FT,
) where {IT,FT,CT}
    N = num_basis(basis)
    @assert length(ICoeff) == 2N "PMCHWT: ICoeff 长度应为 2N=$(2N)，实际为 $(length(ICoeff))"

    # radiation_integral_rwg 内部通过全局 get_k0() 读取波数，必须在调用前同步
    freq_from_k0 = Float64(k0 * Constants.c0 / (2 * FT(π)))
    set_frequency!(freq_from_k0)

    I_J = ICoeff[1:N]
    I_M = ICoeff[N+1:2N]

    Nθ_obs = length(θs_obs)
    Nϕ_obs = length(ϕs_obs)
    nobs = Nθ_obs * Nϕ_obs

    θsobsInfo = [∠Info(θ) for θ in θs_obs]
    ϕsobsInfo = [∠Info(ϕ) for ϕ in ϕs_obs]

    r̂θsϕs = [r̂θϕInfo(θ, ϕ) for θ in θsobsInfo, ϕ in ϕsobsInfo]
    r̂θsϕs_flat = vec(r̂θsϕs)

    RCS_flat = zeros(FT, 2, nobs)

    @threads for ii = 1:nobs
        r_info = r̂θsϕs_flat[ii]

        # N = ∫ J e^{jk₀ r̂·r'} dS  (同现有 EFIE RCS 路径)
        Nθϕ = radiation_integral_rwg(r_info, basis, I_J)
        # L = ∫ M e^{jk₀ r̂·r'} dS  (M 展开在同一 RWG basis 上)
        Lθϕ = radiation_integral_rwg(r_info, basis, I_M)

        # 合并：E_θ = η₀ Nθ + Lϕ,  E_ϕ = η₀ Nϕ - Lθ
        E_theta = eta0 * Nθϕ[1] + Lθϕ[2]
        E_phi = eta0 * Nθϕ[2] - Lθϕ[1]

        factor = k0^2 / (4 * FT(π))
        RCS_flat[1, ii] = factor * abs2(E_theta)
        RCS_flat[2, ii] = factor * abs2(E_phi)
    end

    RCSθsϕs = reshape(RCS_flat, 2, Nθ_obs, Nϕ_obs)
    RCS_total = RCSθsϕs[1, :, :] + RCSθsϕs[2, :, :]
    RCS_dB = 10 .* log10.(RCS_total)

    return RCSθsϕs, RCS_total, RCS_dB
end

end
