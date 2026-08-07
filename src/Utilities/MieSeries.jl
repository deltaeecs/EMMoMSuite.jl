"""
    MieSeries — PEC 球 & 均匀介质球的双站 RCS Mie 级数
"""
module MieSeries

using SpecialFunctions
using LinearAlgebra
using ...CoreModule: Constants

export calculate_mie_rcs_pec_sphere,
    calculate_mie_rcs_dielectric_sphere, calculate_mie_rcs_pec_sphere_fullpol

_emsuite_to_bh_materials(eps_r, mu_r) = conj(Complex(eps_r)), conj(Complex(mu_r))

"""
    calculate_mie_rcs_pec_sphere(radius, freq, theta_range) -> rcs_E

计算 PEC 球双站 RCS (m²), E-plane (S₂) 分量。
入射波: +z 方向传播, x 极化。观测角 0=前向, π=后向.
"""
function calculate_mie_rcs_pec_sphere(radius, freq, theta_range)
    k = 2π * freq / Constants.c0
    x = k * radius
    n_max = ceil(Int, x + 4 * x^(1 / 3) + 2)

    a_n = zeros(ComplexF64, n_max)
    b_n = zeros(ComplexF64, n_max)

    # 时域约定 e^{jωt}: 出射波 ~ h_n^(2) = j_n - i y_n
    # PEC: a_n = -psi_n'/zeta_n', b_n = -psi_n/zeta_n
    for n = 1:n_max
        jn = sphericalbesselj(n, x)
        jnm1 = sphericalbesselj(n - 1, x)
        hn2 = jn - im * sphericalbessely(n, x)
        hn2m1 = jnm1 - im * sphericalbessely(n - 1, x)
        ψn = x * jn
        ψn′ = x * jnm1 - n * jn
        ζn = x * hn2
        ζn′ = x * hn2m1 - n * hn2
        a_n[n] = -ψn′ / ζn′
        b_n[n] = -ψn / ζn
    end

    rcs = zeros(Float64, length(theta_range))
    for (i, θ) in enumerate(theta_range)
        μ = cos(θ)
        S2 = zero(ComplexF64)
        p0 = 0.0
        p1 = 1.0
        for n = 1:n_max
            τn = n * μ * p1 - (n + 1) * p0
            c = (2n + 1) / (n * (n + 1))
            S2 += c * (a_n[n] * τn + b_n[n] * p1)
            pn = ((2n + 1) * μ * p1 - (n + 1) * p0) / n
            p0 = p1
            p1 = pn
        end
        rcs[i] = 4π / k^2 * abs2(S2)
    end
    return rcs
end

"""
    calculate_mie_rcs_pec_sphere_fullpol(radius, freq, theta_range)
    -> (rcs_S2, rcs_S1)

PEC sphere Mie RCS: both S₂-based (theta component, E-plane) and S₁-based
(phi component, H-plane), so the full bistatic sphere can be reconstructed as:

    σ_θθ(θ, φ) = rcs_S2[i] * cos²(φ)
    σ_φφ(θ, φ) = rcs_S1[i] * sin²(φ)
    σ_tot(θ, φ) = σ_θθ + σ_φφ

Convention: +z propagation, x-polarized incident wave.
"""
function calculate_mie_rcs_pec_sphere_fullpol(radius, freq, theta_range)
    k = 2π * freq / Constants.c0
    x = k * radius
    n_max = ceil(Int, x + 4 * x^(1 / 3) + 2)

    a_n = zeros(ComplexF64, n_max)
    b_n = zeros(ComplexF64, n_max)
    for n = 1:n_max
        jn = sphericalbesselj(n, x)
        jnm1 = sphericalbesselj(n - 1, x)
        hn2 = jn - im * sphericalbessely(n, x)
        hn2m1 = jnm1 - im * sphericalbessely(n - 1, x)
        ψn = x * jn
        ψn′ = x * jnm1 - n * jn
        ζn = x * hn2
        ζn′ = x * hn2m1 - n * hn2
        a_n[n] = -ψn′ / ζn′
        b_n[n] = -ψn / ζn
    end

    factor = 4π / k^2
    rcs_S2 = zeros(Float64, length(theta_range))
    rcs_S1 = zeros(Float64, length(theta_range))
    for (i, θ) in enumerate(theta_range)
        μ = cos(θ)
        S1 = zero(ComplexF64)
        S2 = zero(ComplexF64)
        p0 = 0.0
        p1 = 1.0
        for n = 1:n_max
            τn = n * μ * p1 - (n + 1) * p0
            c = (2n + 1) / (n * (n + 1))
            S1 += c * (a_n[n] * p1 + b_n[n] * τn)   # S1: a*π + b*τ
            S2 += c * (a_n[n] * τn + b_n[n] * p1)   # S2: a*τ + b*π
            pn = ((2n + 1) * μ * p1 - (n + 1) * p0) / n
            p0 = p1
            p1 = pn
        end
        rcs_S2[i] = factor * abs2(S2)
        rcs_S1[i] = factor * abs2(S1)
    end
    return rcs_S2, rcs_S1
end

function _calculate_mie_rcs_dielectric_sphere_bh(radius, freq, theta_range, eps_r, mu_r = 1.0)
    k0 = 2π * freq / Constants.c0
    x = k0 * radius
    m = sqrt(Complex(eps_r * mu_r))   # 相对折射率
    mx = m * x                          # 内部尺寸参数

    n_max = ceil(Int, abs(x) + 4.05 * abs(x)^(1 / 3) + 2) + 5

    # 对数导数下行递推; Dmx[n+1] = D_n(mx)
    function downward_D(z::Complex, nm::Int)
        D = zeros(ComplexF64, nm + 2)
        for n = nm+1:-1:1
            D[n] = n / z - 1.0 / (D[n+1] + n / z)
        end
        return D
    end
    Dmx = downward_D(Complex(mx), n_max)

    a_n = zeros(ComplexF64, n_max)
    b_n = zeros(ComplexF64, n_max)

    for n = 1:n_max
        jn = sphericalbesselj(n, x)
        jnm1 = sphericalbesselj(n - 1, x)
        yn = sphericalbessely(n, x)
        ynm1 = sphericalbessely(n - 1, x)

        ψn = x * jn
        ψnm1 = x * jnm1
        ξn = x * (jn + im * yn)     # h_n^(1), B&H 约定; |S|^2 不依赖时域符号
        ξnm1 = x * (jnm1 + im * ynm1)

        Dn = Dmx[n+1]
        An = Dn / m + n / x               # [Dn(mx)/m + n/x]
        Bn = m * Dn + n / x               # [m Dn(mx) + n/x]

        a_n[n] = (An * ψn - ψnm1) / (An * ξn - ξnm1)
        b_n[n] = (Bn * ψn - ψnm1) / (Bn * ξn - ξnm1)
    end

    nobs = length(theta_range)
    rcs_E = zeros(Float64, nobs)
    rcs_H = zeros(Float64, nobs)
    rcs_unpol = zeros(Float64, nobs)
    factor = 4π / k0^2

    for (i, θ) in enumerate(theta_range)
        μ = cos(θ)
        S1 = zero(ComplexF64)
        S2 = zero(ComplexF64)
        p0 = 0.0
        p1 = 1.0
        for n = 1:n_max
            τn = n * μ * p1 - (n + 1) * p0
            c = (2n + 1) / (n * (n + 1))
            S1 += c * (a_n[n] * p1 + b_n[n] * τn)
            S2 += c * (a_n[n] * τn + b_n[n] * p1)
            pn = ((2n + 1) * μ * p1 - (n + 1) * p0) / n
            p0 = p1
            p1 = pn
        end
        rcs_E[i] = factor * abs2(S2)
        rcs_H[i] = factor * abs2(S1)
        rcs_unpol[i] = 0.5 * (rcs_E[i] + rcs_H[i])
    end
    return rcs_E, rcs_H, rcs_unpol
end

"""
    calculate_mie_rcs_dielectric_sphere(radius, freq, theta_range, eps_r, mu_r=1.0)
    -> (rcs_E, rcs_H, rcs_unpol)

计算均匀介质球双站 RCS (m²), Bohren & Huffman (1983) Ch.4 公式.

EMMoMSuite 全库采用 `e^{-jωt}` 约定，因此被动有损介质满足 `imag(εr) < 0`。
本文件的 B&H 级数内核按其原始材料参数约定实现，因此在进入内核前将
EMMoMSuite 约定的 `eps_r`/`mu_r` 做复共轭映射。

约定: 入射波 +z 方向, x 极化 (PlaneWave(freq,0,0,[1,0,0]))
- rcs_E : E-plane φ=0, θ-component (S₂), 对应 RCSθsϕs[1,:,:] at φ=0
- rcs_H : H-plane φ=π/2, ϕ-component (S₁), 对应 RCSθsϕs[2,:,:] at φ=π/2
"""
function calculate_mie_rcs_dielectric_sphere(radius, freq, theta_range, eps_r, mu_r = 1.0)
    eps_bh, mu_bh = _emsuite_to_bh_materials(eps_r, mu_r)
    return _calculate_mie_rcs_dielectric_sphere_bh(radius, freq, theta_range, eps_bh, mu_bh)
end

end # module MieSeries
