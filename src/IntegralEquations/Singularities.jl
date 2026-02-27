module Singularities

using StaticArrays
using LinearAlgebra

export singularF1, singularF21, singularF22, faceSingularityIg, faceSingularityIgIvecg, compute_SSCg, greenfunc_star

const SglrOrder = 20

"""
    compute_SSCg(k::FT) where {FT}

Compute coefficients for Green's function expansion.
"""
function compute_SSCg(k::FT) where {FT}
    SSCg = zeros(Complex{FT}, SglrOrder)
    # exp(-jkr)/R = sum_{n=0} (-jk)^n / n! * R^{n-1}
    # SSCg[n+1] stores (-jk)^n / n!
    
    term = Complex{FT}(1.0)
    SSCg[1] = term
    
    for n in 1:SglrOrder-1
        term *= (-im * k) / n
        SSCg[n+1] = term
    end
    return SSCg
end

"""
    greenfunc_star(R::FT, k::FT) where {FT<:AbstractFloat}

Calculate the smooth part of Green's function (G - 1/R) using Taylor expansion.
"""
function greenfunc_star(R::FT, k::FT) where {FT<:AbstractFloat}
    # Taylor expansion of (exp(-jkr) - 1) / R
    # = (-jkr - (kr)^2/2 + ...) / R
    # = -jk - k^2 R / 2 + ...
    
    # We use the same order as Legacy (15)
    SglrOrder = 15
    
    minusJk = -im * k
    g_star = Complex{FT}(minusJk)
    temp0 = minusJk * R
    temp1 = Complex{FT}(minusJk)
    
    for i in 2:SglrOrder
        temp1 *= temp0 / i
        g_star += temp1
    end
    
    return g_star
end

"""
    singularF1(a::FT, b::FT, c::FT) where{FT<:AbstractFloat}

Calculate the analytical value of the singular integral F1 for a triangle with edge lengths a, b, c.
Integral of 1/R over the triangle.
"""
function singularF1(a::FT, b::FT, c::FT) where{FT<:AbstractFloat}
    s = (a + b + c) / 2
    return -4 * (
        1/a * log(1 - a/s) + 
        1/b * log(1 - b/s) + 
        1/c * log(1 - c/s)
    ) / 3
end

"""
    singularF21(a::FT, b::FT, c::FT, area2::FT) where{FT<:AbstractFloat}

Calculate the analytical value of the singular integral F2 for the self-term (m == n).
Integral of (rho_m . rho_n) / R over the triangle.
"""
function singularF21(a::FT, b::FT, c::FT, area2::FT) where{FT<:AbstractFloat}
    a2 = a^2; b2 = b^2; c2 = c^2
    s = (a + b + c) / 2
    return (
        (10 - 3*(a2-b2)/c2 - 3*(a2-c2)/b2)*a - 
        (5 - 3*(a2-b2)/c2 - 2*(b2-c2)/a2)*b -
        (5 - 3*(a2-c2)/b2 - 2*(c2-b2)/a2)*c +
        (a2 - 3*b2 - 3*c2 - 8*area2/a2)*2/a*log(1 - a/s) +
        (a2 - 2*b2 - 4*c2 + 6*area2/b2)*4/b*log(1 - b/s) + 
        (a2 - 4*b2 - 2*c2 + 6*area2/c2)*4/c*log(1 - c/s)
    ) / 30
end

"""
    singularF22(a::FT, b::FT, c::FT, area2::FT) where{FT<:AbstractFloat}

Calculate the analytical value of the singular integral F2 for the cross-term (m != n).
Integral of (rho_m . rho_n) / R over the triangle.
"""
function singularF22(a::FT, b::FT, c::FT, area2::FT) where{FT<:AbstractFloat}
    a2 = a^2; b2 = b^2; c2 = c^2
    s = (a + b + c) / 2
    return (
        (-10 - (a2-b2)/c2 - (a2-c2)/b2)*a +
        (5 + (a2-b2)/c2 - 6*(b2-c2)/a2)*b +
        (5 + (a2-c2)/b2 - 6*(c2-b2)/a2)*c + 
        (2*a2 - b2 - c2 + 4*area2/a2)*12/a*log(1 - a/s) +
        (9*a2 - 3*b2 - c2 + 4*area2/b2)*2/b*log(1 - b/s) +
        (9*a2 - b2 - 3*c2 + 4*area2/c2)*2/c*log(1 - c/s)
    ) / 60
end

"""
    faceSingularityIg(rgt, vertices, edgel, edgev, edgen, area, facen)

Calculate the singular integral of Green's function over a triangle.
"""
function faceSingularityIg(rgt::AbstractVector{FT}, vertices::AbstractMatrix{FT}, 
                          edgel::AbstractVector{FT}, edgev::AbstractMatrix{FT}, 
                          edgen::AbstractMatrix{FT}, area::FT, facen::AbstractVector{FT},
                          SSCg::AbstractVector{Complex{FT}}) where {FT<:Real}
    
    # Il_r stores integrals of R^n along edges
    # Indices: 1->n=-1, 2->n=0, ...
    Il_r = zeros(FT, SglrOrder)
    
    # IS_r stores integrals of R^n over surface
    IS_r = zeros(FT, SglrOrder)
    
    ISg = Complex{FT}(0.0)
    
    # Distance to plane
    dts = dot(facen, rgt - vertices[:, 1])
    r0gi = rgt - dts * facen
    dtsAbs = abs(dts)
    dts2 = dtsAbs^2
    
    # Loop over 3 edges
    # Edge 1: v2-v3. Vertices indices [2, 3]
    # Edge 2: v3-v1. Vertices indices [3, 1]
    # Edge 3: v1-v2. Vertices indices [1, 2]
    # Note: Legacy uses specific indexing.
    # edgev[:, j] is unit vector of edge j.
    # edgen[:, j] is normal to edge j (in plane).
    
    # Vertices mapping for edges
    # Edge 1 starts at v2?
    # Legacy: edgeNodei- = vertices[:, EDGEVmINTriVsID[edgej]]
    # EDGEVmINTriVsID = [2, 3, 1]
    
    edge_start_indices = [2, 3, 1]
    
    for edgej in 1:3
        edgeNodei_minus = vertices[:, edge_start_indices[edgej]]
        lj = abs(edgel[edgej])
        
        # lj_minus = (edgeNodei_minus - r0gi) . edgev[:, edgej]
        lj_minus = dot(edgeNodei_minus - r0gi, edgev[:, edgej])
        lj_plus = lj_minus + lj
        
        # p02jvec
        p02jvec = edgeNodei_minus - lj_minus * edgev[:, edgej] - r0gi
        p02jl = dot(p02jvec, edgen[:, edgej])
        
        R02 = p02jl^2 + dts2
        R_plus = sqrt(lj_plus^2 + R02)
        R_minus = sqrt(lj_minus^2 + R02)
        
        fj = 0.0
        betaj = 0.0
        epsilon_l = 1e-2 * lj
        
        if abs(p02jl) < epsilon_l
            if dtsAbs < epsilon_l
                # Point on edge
                continue
            else
                # Point on plane, projected on edge
                fj = log((lj_plus + R_plus) / (lj_minus + R_minus))
            end
        else
            fj = log((lj_plus + R_plus) / (lj_minus + R_minus))
            if dtsAbs < epsilon_l
                # Point on plane
                betaj = atan((p02jl * lj_plus) / R02) - atan((p02jl * lj_minus) / R02)
                IS_r[1] += p02jl * fj # n=-1 -> index 1
            else
                betaj = atan((p02jl * lj_plus) / (R02 + dtsAbs * R_plus)) - 
                        atan((p02jl * lj_minus) / (R02 + dtsAbs * R_minus))
                IS_r[1] += p02jl * fj - dtsAbs * betaj
            end
        end
        
        # Calculate Il_r
        Il_r[1] = fj # n=-1
        Il_r[2] = lj # n=0
        
        R_plus_n = 1.0
        R_minus_n = 1.0
        
        for n in 1:(SglrOrder-2)
            # n here corresponds to power n in R^n
            # We want Il_r[n+2]
            R_plus_n *= R_plus
            R_minus_n *= R_minus
            
            # Recursive formula:
            # int R^n dl = (l+ R+^n - l- R-^n + n R0^2 int R^{n-2}) / (n+1)
            # Il_r[n+2] corresponds to power n
            # Il_r[n] corresponds to power n-2
            
            Il_r[n+2] = (lj_plus * R_plus_n - lj_minus * R_minus_n + n * R02 * Il_r[n]) / (n+1)
        end
        
        # Add to IS_r
        for n in 1:(SglrOrder-2)
            # IS_r[n+2] corresponds to power n
            # IS_r[n+2] += p02jl * Il_r[n+2] / (n+2) ?
            # Wait, Legacy: ISr[n] += p02jl * Ilr[n]
            # Legacy n is power.
            # So IS_r[n+2] += p02jl * Il_r[n+2]
            IS_r[n+2] += p02jl * Il_r[n+2]
        end
    end
    
    # IS_r recursion
    IS_r[2] = area # n=0
    
    for n in 1:(SglrOrder-2)
        # IS_r[n+2] += n * dts^2 * IS_r[n]
        # IS_r[n+2] /= (n+2)
        IS_r[n+2] += n * dts2 * IS_r[n]
        IS_r[n+2] /= (n + 2)
    end
    
    # Sum Green's function
    for n in 0:(SglrOrder-1)
        # n is power of (-jk)
        # Corresponds to R^{n-1}
        # Index in IS_r is n+1
        ISg += SSCg[n+1] * IS_r[n+1]
    end
    
    return ISg
end

"""
    faceSingularityIgIvecg(rgt, vertices, edgel, edgev, edgen, area, facen, SSCg)

Calculate the singular integral of Green's function and its vector moment over a triangle.
Returns (ISg, IvecSg).
"""
function faceSingularityIgIvecg(rgt::AbstractVector{FT}, vertices::AbstractMatrix{FT}, 
                               edgel::AbstractVector{FT}, edgev::AbstractMatrix{FT}, 
                               edgen::AbstractMatrix{FT}, area::FT, facen::AbstractVector{FT},
                               SSCg::AbstractVector{Complex{FT}}) where {FT<:Real}
    
    # Il_r stores integrals of R^n along edges
    Il_r = zeros(FT, SglrOrder)
    # IS_r stores integrals of R^n over surface
    IS_r = zeros(FT, SglrOrder)
    
    ISg = Complex{FT}(0.0)
    IvecSg = zeros(Complex{FT}, 3)
    
    # Precompute SSCgdivnp1
    # SSCg[n+1] is coeff for R^{n-1}. (n=0 -> R^-1)
    # We need coeff(n)/(n+1) ?
    # Legacy loop uses SSCgdivnp1[n] for Il_r[n+1] (R^{n-1}).
    # So SSCgdivnp1[n] corresponds to power n-1.
    # So it should be coeff(n-1) / (n-1+2) ? No.
    # Let's assume Legacy SSCgdivnp1[n] = SSCg[n+1] / (n+1) ?
    # Actually, let's just use SSCg directly.
    # The term is coeff(n) * I_{R^n} ?
    # The formula for Ivec is sum u_j * sum coeff(n)/(n+2) * I_{l, n+1} ?
    # This is getting complicated to derive.
    # I will trust the Legacy loop structure:
    # Cdvnp1I = -Il_r[1] (R^-1 term)
    # Loop n=1..: subtract SSCg[n+1]/(n+1) * Il_r[n+1] (R^{n-1} term)
    # Wait, SSCg[n+1] is coeff for R^{n-1}.
    # So we are using coeff(n-1)/(n) ?
    # If n=1 (loop), Il_r[2] is R^0. coeff is SSCg[2] (-jk).
    # We divide by 2?
    
    dts = dot(facen, rgt - vertices[:, 1])
    r0gi = rgt - dts * facen
    dtsAbs = abs(dts)
    dts2 = dtsAbs^2
    
    edge_start_indices = [2, 3, 1]
    
    for edgej in 1:3
        edgeNodei_minus = vertices[:, edge_start_indices[edgej]]
        lj = abs(edgel[edgej])
        
        lj_minus = dot(edgeNodei_minus - r0gi, edgev[:, edgej])
        lj_plus = lj_minus + lj
        
        p02jvec = edgeNodei_minus - lj_minus * edgev[:, edgej] - r0gi
        p02jl = dot(p02jvec, edgen[:, edgej])
        
        R02 = p02jl^2 + dts2
        R_plus = sqrt(lj_plus^2 + R02)
        R_minus = sqrt(lj_minus^2 + R02)
        
        fj = 0.0
        betaj = 0.0
        epsilon_l = 1e-2 * lj
        
        if abs(p02jl) < epsilon_l
            if dtsAbs < epsilon_l
                continue
            else
                fj = log((lj_plus + R_plus) / (lj_minus + R_minus))
            end
        else
            fj = log((lj_plus + R_plus) / (lj_minus + R_minus))
            if dtsAbs < epsilon_l
                betaj = atan((p02jl * lj_plus) / R02) - atan((p02jl * lj_minus) / R02)
                IS_r[1] += p02jl * fj
            else
                betaj = atan((p02jl * lj_plus) / (R02 + dtsAbs * R_plus)) - 
                        atan((p02jl * lj_minus) / (R02 + dtsAbs * R_minus))
                IS_r[1] += p02jl * fj - dtsAbs * betaj
            end
        end
        
        Il_r[1] = fj
        Il_r[2] = lj
        
        R_plus_n = 1.0
        R_minus_n = 1.0
        
        for n in 1:(SglrOrder-2)
            R_plus_n *= R_plus
            R_minus_n *= R_minus
            Il_r[n+2] = (lj_plus * R_plus_n - lj_minus * R_minus_n + n * R02 * Il_r[n]) / (n+1)
        end
        
        for n in 1:(SglrOrder-2)
            IS_r[n+2] += p02jl * Il_r[n+2]
        end
        
        # IvecSg contribution
        # Legacy implementation:
        # Cdvnp1I = CT(-Ilᵣ[1])
        # for n in 1:(SglrOrder-3)
        #     Cdvnp1I -= SSCgdivnp1[n] * Ilᵣ[n + 1]
        # end
        # Mapping: Ilᵣ[k] -> Il_r[k+2]
        # SSCgdivnp1[n] -> SSCg[n+2]/(n+2)
        
        Cdvnp1I = Complex{FT}(-Il_r[3])
        for n in 1:(SglrOrder-3)
            Cdvnp1I -= (SSCg[n+1] / (n+1)) * Il_r[n+3]
        end
        
        # IvecSg += edgen * Cdvnp1I
        for k in 1:3
            IvecSg[k] += edgen[k, edgej] * Cdvnp1I
        end
    end
    
    IS_r[2] = area
    for n in 1:(SglrOrder-2)
        IS_r[n+2] += n * dts2 * IS_r[n]
        IS_r[n+2] /= (n + 2)
    end
    
    for n in 0:(SglrOrder-1)
        ISg += SSCg[n+1] * IS_r[n+1]
    end
    
    # IvecSg += facen * (dts * ISg)
    for k in 1:3
        IvecSg[k] += facen[k] * (dts * ISg)
    end
    
    return ISg, IvecSg
end

end
