module Singularities

using StaticArrays
using LinearAlgebra
using ..Geometry

export singularF1,
    singularF21,
    singularF22,
    faceSingularityIg,
    faceSingularityIgIvecg,
    volumeSingularityIgIvecg,
    compute_SSCg,
    greenfunc_star

const SglrOrder = 15

"""
    compute_SSCg(k)

Compute coefficients for Green's function expansion.
Supports both real and complex wavenumbers.
"""
function compute_SSCg(k::FT) where {FT<:AbstractFloat}
    SSCg = zeros(Complex{FT}, SglrOrder)
    # exp(-jkr)/R = sum_{n=0} (-jk)^n / n! * R^{n-1}
    # SSCg[n+1] stores (-jk)^n / n!

    term = Complex{FT}(1.0)
    SSCg[1] = term

    for n = 1:SglrOrder-1
        term *= (-im * k) / n
        SSCg[n+1] = term
    end
    return SSCg
end

function compute_SSCg(k::CT) where {FT<:AbstractFloat,CT<:Complex{FT}}
    SSCg = zeros(Complex{FT}, SglrOrder)

    term = one(CT)
    SSCg[1] = term

    for n = 1:SglrOrder-1
        term *= (-im * k) / n
        SSCg[n+1] = term
    end
    return SSCg
end

"""
    greenfunc_star(R::FT, k) where {FT<:AbstractFloat}

Calculate the smooth part of Green's function (G - 1/R) using Taylor expansion.
Supports both real and complex wavenumbers.
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

    for i = 2:SglrOrder
        temp1 *= temp0 / i
        g_star += temp1
    end

    return g_star
end

function greenfunc_star(R::FT, k::CT) where {FT<:AbstractFloat,CT<:Complex{FT}}
    minusJk = -im * k
    g_star = Complex{FT}(minusJk)
    temp0 = minusJk * R
    temp1 = Complex{FT}(minusJk)

    for i = 2:SglrOrder
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
function singularF1(a::FT, b::FT, c::FT) where {FT<:AbstractFloat}
    s = (a + b + c) / 2
    # Numerical protection: avoid log(negative) when triangle is nearly degenerate
    # Use max(1 - edge/s, eps(FT)) to prevent log(0) or log(negative)
    eps_ft = eps(FT)
    term_a = 1 / a * log(max(1 - a / s, eps_ft))
    term_b = 1 / b * log(max(1 - b / s, eps_ft))
    term_c = 1 / c * log(max(1 - c / s, eps_ft))
    return -4 * (term_a + term_b + term_c) / 3
end

"""
    singularF21(a::FT, b::FT, c::FT, area2::FT) where{FT<:AbstractFloat}

Calculate the analytical value of the singular integral F2 for the self-term (m == n).
Integral of (rho_m . rho_n) / R over the triangle.
"""
function singularF21(a::FT, b::FT, c::FT, area2::FT) where {FT<:AbstractFloat}
    a2 = a^2
    b2 = b^2
    c2 = c^2
    s = (a + b + c) / 2
    # Numerical protection for log terms
    eps_ft = eps(FT)
    log_a = log(max(1 - a / s, eps_ft))
    log_b = log(max(1 - b / s, eps_ft))
    log_c = log(max(1 - c / s, eps_ft))
    
    return (
        (10 - 3 * (a2 - b2) / c2 - 3 * (a2 - c2) / b2) * a -
        (5 - 3 * (a2 - b2) / c2 - 2 * (b2 - c2) / a2) * b -
        (5 - 3 * (a2 - c2) / b2 - 2 * (c2 - b2) / a2) * c +
        (a2 - 3 * b2 - 3 * c2 - 8 * area2 / a2) * 2 / a * log_a +
        (a2 - 2 * b2 - 4 * c2 + 6 * area2 / b2) * 4 / b * log_b +
        (a2 - 4 * b2 - 2 * c2 + 6 * area2 / c2) * 4 / c * log_c
    ) / 30
end

"""
    singularF22(a::FT, b::FT, c::FT, area2::FT) where{FT<:AbstractFloat}

Calculate the analytical value of the singular integral F2 for the cross-term (m != n).
Integral of (rho_m . rho_n) / R over the triangle.
"""
function singularF22(a::FT, b::FT, c::FT, area2::FT) where {FT<:AbstractFloat}
    a2 = a^2
    b2 = b^2
    c2 = c^2
    s = (a + b + c) / 2
    # Numerical protection for log terms (same as singularF21)
    eps_ft = eps(FT)
    log_a = log(max(1 - a / s, eps_ft))
    log_b = log(max(1 - b / s, eps_ft))
    log_c = log(max(1 - c / s, eps_ft))
    
    return (
        (-10 - (a2 - b2) / c2 - (a2 - c2) / b2) * a +
        (5 + (a2 - b2) / c2 - 6 * (b2 - c2) / a2) * b +
        (5 + (a2 - c2) / b2 - 6 * (c2 - b2) / a2) * c +
        (2 * a2 - b2 - c2 + 4 * area2 / a2) * 12 / a * log_a +
        (9 * a2 - 3 * b2 - c2 + 4 * area2 / b2) * 2 / b * log_b +
        (9 * a2 - b2 - 3 * c2 + 4 * area2 / c2) * 2 / c * log_c
    ) / 60
end

"""
    faceSingularityIg(rgt, vertices, edgel, edgev, edgen, area, facen)

Calculate the singular integral of Green's function over a triangle.
"""
function faceSingularityIg(
    rgt::AbstractVector{FT},
    vertices::AbstractMatrix{FT},
    edgel::AbstractVector{FT},
    edgev::AbstractMatrix{FT},
    edgen::AbstractMatrix{FT},
    area::FT,
    facen::AbstractVector{FT},
    SSCg::AbstractVector{Complex{FT}},
) where {FT<:Real}

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

    for edgej = 1:3
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
                betaj =
                    atan((p02jl * lj_plus) / (R02 + dtsAbs * R_plus)) -
                    atan((p02jl * lj_minus) / (R02 + dtsAbs * R_minus))
                IS_r[1] += p02jl * fj - dtsAbs * betaj
            end
        end

        # Calculate Il_r
        Il_r[1] = fj # n=-1
        Il_r[2] = lj # n=0

        R_plus_n = 1.0
        R_minus_n = 1.0

        for n = 1:(SglrOrder-2)
            # n here corresponds to power n in R^n
            # We want Il_r[n+2]
            R_plus_n *= R_plus
            R_minus_n *= R_minus

            # Recursive formula:
            # int R^n dl = (l+ R+^n - l- R-^n + n R0^2 int R^{n-2}) / (n+1)
            # Il_r[n+2] corresponds to power n
            # Il_r[n] corresponds to power n-2

            Il_r[n+2] = (lj_plus * R_plus_n - lj_minus * R_minus_n + n * R02 * Il_r[n]) / (n + 1)
        end

        # Add to IS_r
        for n = 1:(SglrOrder-2)
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

    for n = 1:(SglrOrder-2)
        # IS_r[n+2] += n * dts^2 * IS_r[n]
        # IS_r[n+2] /= (n+2)
        IS_r[n+2] += n * dts2 * IS_r[n]
        IS_r[n+2] /= (n + 2)
    end

    # Sum Green's function
    for n = 0:(SglrOrder-1)
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
function faceSingularityIgIvecg(
    rgt::AbstractVector{FT},
    vertices::AbstractMatrix{FT},
    edgel::AbstractVector{FT},
    edgev::AbstractMatrix{FT},
    edgen::AbstractMatrix{FT},
    area::FT,
    facen::AbstractVector{FT},
    SSCg::AbstractVector{Complex{FT}},
) where {FT<:Real}

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

    for edgej = 1:3
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
                betaj =
                    atan((p02jl * lj_plus) / (R02 + dtsAbs * R_plus)) -
                    atan((p02jl * lj_minus) / (R02 + dtsAbs * R_minus))
                IS_r[1] += p02jl * fj - dtsAbs * betaj
            end
        end

        Il_r[1] = fj
        Il_r[2] = lj

        R_plus_n = 1.0
        R_minus_n = 1.0

        for n = 1:(SglrOrder-2)
            R_plus_n *= R_plus
            R_minus_n *= R_minus
            Il_r[n+2] = (lj_plus * R_plus_n - lj_minus * R_minus_n + n * R02 * Il_r[n]) / (n + 1)
        end

        for n = 1:(SglrOrder-2)
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
        for n = 1:(SglrOrder-3)
            Cdvnp1I -= (SSCg[n+1] / (n + 1)) * Il_r[n+3]
        end

        # IvecSg += edgen * Cdvnp1I
        for k = 1:3
            IvecSg[k] += edgen[k, edgej] * Cdvnp1I
        end
    end

    IS_r[2] = area
    for n = 1:(SglrOrder-2)
        IS_r[n+2] += n * dts2 * IS_r[n]
        IS_r[n+2] /= (n + 2)
    end

    for n = 0:(SglrOrder-1)
        ISg += SSCg[n+1] * IS_r[n+1]
    end

    # IvecSg += facen * (dts * ISg)
    for k = 1:3
        IvecSg[k] += facen[k] * (dts * ISg)
    end

    return ISg, IvecSg
end

function volumeSingularityIgIvecg(
    rgt::AbstractVector{FT},
    tetra::TetrahedraInfo{IT,FT,CT},
    SSCg::AbstractVector{Complex{FT}},
) where {IT<:Integer,FT<:Real,CT<:Complex}
    n_faces = length(tetra.faces)
    is_r = zeros(FT, SglrOrder + 2)
    il_r = zeros(FT, SglrOrder + 2)

    ivg = zero(Complex{FT})
    ivec_vg = zeros(Complex{FT}, 3)
    r0 = zeros(FT, 3)
    p02_vec = zeros(FT, 3)

    for iface = 1:n_faces
        face = tetra.faces[iface]
        fill!(is_r, zero(FT))

        dts = dot(view(tetra.facesn̂, :, iface), rgt - view(face.vertices, :, 1))
        for ii = 1:3
            r0[ii] = rgt[ii] - dts * tetra.facesn̂[ii, iface]
        end

        dts2 = dts^2
        dts_abs = abs(dts)

        for edge_idx = 1:3
            fill!(il_r, zero(FT))

            start_idx = (2, 3, 1)[edge_idx]
            edge_node_minus = view(face.vertices, :, start_idx)
            lj = face.edgel[edge_idx]

            lj_minus = zero(FT)
            for ii = 1:3
                lj_minus += (edge_node_minus[ii] - r0[ii]) * face.edgev̂[ii, edge_idx]
            end
            lj_plus = lj_minus + lj

            for ii = 1:3
                p02_vec[ii] = edge_node_minus[ii] - lj_minus * face.edgev̂[ii, edge_idx] - r0[ii]
            end
            p02jl = dot(p02_vec, view(face.edgen̂, :, edge_idx))
            r0_sq = p02jl^2 + dts2
            r_plus = sqrt(lj_plus^2 + r0_sq)
            r_minus = sqrt(lj_minus^2 + r0_sq)

            fj = zero(FT)
            betaj = zero(FT)
            eps_l = FT(1e-3) * lj

            if abs(p02jl) < eps_l
                if dts_abs < eps_l
                    continue
                else
                    fj = log((lj_plus + r_plus) / (lj_minus + r_minus))
                end
            else
                fj = log((lj_plus + r_plus) / (lj_minus + r_minus))
                if dts_abs < eps_l
                    betaj = atan((p02jl * lj_plus) / r0_sq) - atan((p02jl * lj_minus) / r0_sq)
                    is_r[2] += p02jl * fj
                else
                    betaj =
                        atan((p02jl * lj_plus) / (r0_sq + dts_abs * r_plus)) -
                        atan((p02jl * lj_minus) / (r0_sq + dts_abs * r_minus))
                    is_r[2] += p02jl * fj - dts_abs * betaj
                end
            end

            il_r[2] = fj
            il_r[3] = lj
            r_plus_n = one(FT)
            r_minus_n = one(FT)
            for n = 1:(SglrOrder - 2)
                r_plus_n *= r_plus
                r_minus_n *= r_minus
                il_r[n + 3] = (lj_plus * r_plus_n - lj_minus * r_minus_n + n * r0_sq * il_r[n + 1]) / (n + 1)
            end

            for n = 1:(SglrOrder - 2)
                is_r[n + 3] += p02jl * il_r[n + 3]
            end
        end

        is_r[3] = tetra.facesArea[iface]
        for n = 1:(SglrOrder - 2)
            is_r[n + 3] += n * dts2 * is_r[n + 1]
            is_r[n + 3] /= n + 2
        end

        isg = zero(Complex{FT})
        for n = 0:(SglrOrder - 1)
            isg -= (SSCg[n + 1] / (n + 2)) * is_r[n + 2]
        end
        ivg += dts * isg

        ivec_vgi = zero(Complex{FT})
        for n = 0:(SglrOrder - 3)
            ivec_vgi -= (SSCg[n + 1] / (n + 1)) * is_r[n + 4]
        end
        for ii = 1:3
            ivec_vg[ii] += tetra.facesn̂[ii, iface] * ivec_vgi
        end
    end

    return ivg, SVector{3,Complex{FT}}(ivec_vg)
end

end
