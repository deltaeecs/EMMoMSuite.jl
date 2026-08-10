module MFIEModule

using ..CoreModule
using ..CoreModule: Constants
using ..Geometry
using ..BasisFunctions
using ..EFIEModule
using ..Kernels
using ..Impedance
using StaticArrays
using LinearAlgebra
using SparseArrays

import ..CoreModule: assemble_impedance_matrix

export MFIE, assemble_impedance_matrix, mfie_interaction!, assemble_K_offdiag

const MFIE_NEAR_SCALE = 1.5
const MFIE_MIN_GAP_SCALE = 0.5 * MFIE_NEAR_SCALE

"""
    MFIE{FT, CT} <: AbstractIntegralOperator

Magnetic Field Integral Equation (MFIE) operator.

Solves for the surface current density \$\\mathbf{J}\$ on a closed PEC object by enforcing the boundary condition \$\\hat{n} \\times \\mathbf{H}^{total} = \\mathbf{J}\$.

# Mathematical Formulation

The MFIE is given by:
```math
\\frac{1}{2}\\mathbf{J}(\\mathbf{r}) - \\hat{n}(\\mathbf{r}) \\times \\int_S \\mathbf{J}(\\mathbf{r}') \\times \\nabla' G(\\mathbf{r}, \\mathbf{r}') dS' = \\hat{n}(\\mathbf{r}) \\times \\mathbf{H}^{inc}(\\mathbf{r})
```

The impedance matrix element \$Z_{mn}\$ consists of a mass matrix term and the K-operator term.
To maintain consistency with EFIE for CFIE formulations, we scale the MFIE by the intrinsic impedance \$\\eta\$.

```math
Z_{mn} = \\eta \\left( \\frac{1}{2} \\int_S \\mathbf{f}_m(\\mathbf{r}) \\cdot \\mathbf{f}_n(\\mathbf{r}) dS - \\int_S \\mathbf{f}_m(\\mathbf{r}) \\cdot \\left( \\hat{n}(\\mathbf{r}) \\times \\int_{S'} \\mathbf{f}_n(\\mathbf{r}') \\times \\nabla' G(\\mathbf{r}, \\mathbf{r}') dS' \\right) dS \\right)
```

# Fields
- `freq`: Operating frequency (Hz).
- `k`: Wavenumber.
- `eta`: Intrinsic impedance.
- `gq_far`: Far-field quadrature information.
- `gq_near`: Near-field quadrature information for adjacent panels.
"""
struct MFIE{FT<:AbstractFloat,CT<:Complex} <: AbstractIntegralOperator
    freq::FT
    k::FT
    eta::FT
    gq_far::GaussQuadratureInfoStruct{FT,4,3}
    gq_near::GaussQuadratureInfoStruct{FT,7,3}
end

function MFIE(freq::FT) where {FT}
    k = FT(2π * freq / Constants.c0)
    eta = FT(Constants.eta0)

    # Legacy uses GQPNTri=4 for regular MFIE K-operator evaluations.
    # Adjacent panels are more sensitive, so we mirror EFIE's 7-point near rule.
    gq_far = GaussQuadratureInfo(:Triangle, 4, FT)
    gq_near = GaussQuadratureInfo(:Triangle, 7, FT)

    return MFIE{FT,Complex{FT}}(freq, k, eta, gq_far, gq_near)
end

"""
    assemble_impedance_matrix(mfie::MFIE, basis::RWGBasis)

Assemble the impedance matrix Z for the MFIE using RWG basis functions.
Z = eta * (0.5 * M + K)

Quad points are precomputed once to avoid per-pair heap allocations.
"""
function assemble_impedance_matrix(mfie::MFIE{FT,CT}, basis::RWGBasis{IT,FT}) where {IT,FT,CT}
    # Precompute quadrature points for all triangles (like EFIE does)
    mesh = basis.mesh
    nt = num_elements(mesh)
    gq = mfie.gq_far
    N_points = length(gq.weight)

    quad_points = Vector{SVector{N_points,SVector{3,FT}}}(undef, nt)
    Threads.@threads for t = 1:nt
        v_indices = mesh.triangles[:, t]
        v1 = SVector{3,FT}(mesh.node[:, v_indices[1]])
        v2 = SVector{3,FT}(mesh.node[:, v_indices[2]])
        v3 = SVector{3,FT}(mesh.node[:, v_indices[3]])
        quad_points[t] = SVector{N_points,SVector{3,FT}}(
            v1 * gq.coordinate[1, i] + v2 * gq.coordinate[2, i] + v3 * gq.coordinate[3, i] for
            i = 1:N_points
        )
    end

    # Wrapper with precomputed points
    wrapper = (Z, op, t1, t2) -> mfie_interaction!(Z, op, t1, t2, quad_points)
    return assemble_generic(mfie, basis, wrapper)
end

function mfie_interaction!(
    Z_local::AbstractMatrix{CT},
    mfie::MFIE{FT,CT},
    tri_test::TriangleInfo{IT,FT},
    tri_source::TriangleInfo{IT,FT},
) where {IT,FT,CT}
    # Check if self-interaction
    if tri_test.triID == tri_source.triID
        calc_self_term!(Z_local, mfie, tri_test)
    else
        calc_k_term!(Z_local, mfie, tri_test, tri_source)
    end
    return nothing
end

# Fast path with precomputed quad points (eliminates per-pair allocation)
function mfie_interaction!(
    Z_local::AbstractMatrix{CT},
    mfie::MFIE{FT,CT},
    tri_test::TriangleInfo{IT,FT},
    tri_source::TriangleInfo{IT,FT},
    quad_points::Vector,
) where {IT,FT,CT}
    if tri_test.triID == tri_source.triID
        calc_self_term!(Z_local, mfie, tri_test)
    elseif needs_near_quadrature(tri_test, tri_source)
        r_test = get_global_quad_points(tri_test, mfie.gq_near)
        r_src = get_global_quad_points(tri_source, mfie.gq_near)
        calc_k_term_fast!(Z_local, mfie, tri_test, tri_source, r_test, r_src, mfie.gq_near.weight)
    else
        r_test = quad_points[tri_test.triID]
        r_src = quad_points[tri_source.triID]
        calc_k_term_fast!(Z_local, mfie, tri_test, tri_source, r_test, r_src, mfie.gq_far.weight)
    end
    return nothing
end

function calc_self_term!(
    Z_local::AbstractMatrix{CT},
    mfie::MFIE{FT,CT},
    tri::TriangleInfo{IT,FT},
) where {IT,FT,CT}
    # Self-term is just the Mass Matrix term scaled by eta
    # Z_self = eta * 0.5 * <fm, fn>
    # Based on legacy code: Z = (eta / 8A) * lm * ln * sum(rho_m . rho_n * w)

    gq = mfie.gq_far
    r_quad = get_global_quad_points(tri, gq)
    w_quad = gq.weight
    n_pts = length(w_quad)

    # Precompute vertices as SVector for stack allocation
    v1 = SVector{3,FT}(tri.vertices[:, 1])
    v2 = SVector{3,FT}(tri.vertices[:, 2])
    v3 = SVector{3,FT}(tri.vertices[:, 3])

    eta_div_8A = mfie.eta / (8 * tri.area)

    for m = 1:3
        lm = tri.edgel[m]
        if tri.inBfsID[m] == 0
            continue
        end

        for n = 1:3
            ln = tri.edgel[n]
            if tri.inBfsID[n] == 0
                continue
            end

            sum_val = zero(CT)
            @inbounds for k = 1:n_pts
                rk = r_quad[k]
                # Compute rho vectors inline (no heap allocation)
                rho_m = rk - (m == 1 ? v1 : m == 2 ? v2 : v3)
                rho_n = rk - (n == 1 ? v1 : n == 2 ? v2 : v3)
                sum_val += dot(rho_m, rho_n) * w_quad[k]
            end

            Z_local[m, n] += sum_val * lm * ln * eta_div_8A
        end
    end

    return nothing
end

function calc_k_term!(
    Z_local::AbstractMatrix{CT},
    mfie::MFIE{FT,CT},
    tri_test::TriangleInfo{IT,FT},
    tri_source::TriangleInfo{IT,FT},
) where {IT,FT,CT}
    gq = needs_near_quadrature(tri_test, tri_source) ? mfie.gq_near : mfie.gq_far
    r_test = get_global_quad_points(tri_test, gq)
    r_src = get_global_quad_points(tri_source, gq)
    calc_k_term_fast!(Z_local, mfie, tri_test, tri_source, r_test, r_src, gq.weight)
end

function needs_near_quadrature(t1::TriangleInfo{IT,FT}, t2::TriangleInfo{IT,FT}) where {IT,FT}
    if EFIEModule.is_adjacent(t1, t2)
        return true
    end

    center_dist = norm(t1.center - t2.center)
    r1 = max(
        norm(t1.vertices[:, 1] - t1.center),
        norm(t1.vertices[:, 2] - t1.center),
        norm(t1.vertices[:, 3] - t1.center),
    )
    r2 = max(
        norm(t2.vertices[:, 1] - t2.center),
        norm(t2.vertices[:, 2] - t2.center),
        norm(t2.vertices[:, 3] - t2.center),
    )
    near_radius = r1 + r2
    center_dist <= FT(MFIE_NEAR_SCALE) * near_radius && return true

    min_gap = _triangle_triangle_distance(t1, t2)
    return min_gap <= FT(MFIE_MIN_GAP_SCALE) * near_radius
end

function _triangle_triangle_distance(t1::TriangleInfo{IT,FT}, t2::TriangleInfo{IT,FT}) where {IT,FT}
    verts1 = ntuple(i -> SVector{3,FT}(t1.vertices[:, i]), 3)
    verts2 = ntuple(i -> SVector{3,FT}(t2.vertices[:, i]), 3)

    min_dist_sq = typemax(FT)
    @inbounds for p in verts1
        min_dist_sq =
            min(min_dist_sq, _point_triangle_distance_sq(p, verts2[1], verts2[2], verts2[3]))
    end
    @inbounds for p in verts2
        min_dist_sq =
            min(min_dist_sq, _point_triangle_distance_sq(p, verts1[1], verts1[2], verts1[3]))
    end

    edges = ((1, 2), (2, 3), (3, 1))
    @inbounds for (i1, i2) in edges, (j1, j2) in edges
        min_dist_sq = min(
            min_dist_sq,
            _segment_segment_distance_sq(verts1[i1], verts1[i2], verts2[j1], verts2[j2]),
        )
    end

    return sqrt(max(min_dist_sq, zero(FT)))
end

function _point_triangle_distance_sq(
    p::SVector{3,FT},
    a::SVector{3,FT},
    b::SVector{3,FT},
    c::SVector{3,FT},
) where {FT}
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = dot(ab, ap)
    d2 = dot(ac, ap)
    if d1 <= zero(FT) && d2 <= zero(FT)
        return dot(ap, ap)
    end

    bp = p - b
    d3 = dot(ab, bp)
    d4 = dot(ac, bp)
    if d3 >= zero(FT) && d4 <= d3
        return dot(bp, bp)
    end

    vc = d1 * d4 - d3 * d2
    if vc <= zero(FT) && d1 >= zero(FT) && d3 <= zero(FT)
        v = d1 / (d1 - d3)
        proj = a + v * ab
        diff = p - proj
        return dot(diff, diff)
    end

    cp = p - c
    d5 = dot(ab, cp)
    d6 = dot(ac, cp)
    if d6 >= zero(FT) && d5 <= d6
        return dot(cp, cp)
    end

    vb = d5 * d2 - d1 * d6
    if vb <= zero(FT) && d2 >= zero(FT) && d6 <= zero(FT)
        w = d2 / (d2 - d6)
        proj = a + w * ac
        diff = p - proj
        return dot(diff, diff)
    end

    va = d3 * d6 - d5 * d4
    if va <= zero(FT) && (d4 - d3) >= zero(FT) && (d5 - d6) >= zero(FT)
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        proj = b + w * (c - b)
        diff = p - proj
        return dot(diff, diff)
    end

    n = cross(ab, ac)
    signed_dist = dot(ap, n)
    return signed_dist * signed_dist / dot(n, n)
end

function _segment_segment_distance_sq(
    p1::SVector{3,FT},
    q1::SVector{3,FT},
    p2::SVector{3,FT},
    q2::SVector{3,FT},
) where {FT}
    u = q1 - p1
    v = q2 - p2
    w = p1 - p2

    a = dot(u, u)
    b = dot(u, v)
    c = dot(v, v)
    d = dot(u, w)
    e = dot(v, w)
    D = a * c - b * b
    tol = sqrt(Base.eps(FT))

    sN = zero(FT)
    sD = D
    tN = zero(FT)
    tD = D

    if D <= tol
        sN = zero(FT)
        sD = one(FT)
        tN = e
        tD = c
    else
        sN = b * e - c * d
        tN = a * e - b * d
        if sN < zero(FT)
            sN = zero(FT)
            tN = e
            tD = c
        elseif sN > sD
            sN = sD
            tN = e + b
            tD = c
        end
    end

    if tN < zero(FT)
        tN = zero(FT)
        if -d < zero(FT)
            sN = zero(FT)
        elseif -d > a
            sN = sD
        else
            sN = -d
            sD = a
        end
    elseif tN > tD
        tN = tD
        if (-d + b) < zero(FT)
            sN = zero(FT)
        elseif (-d + b) > a
            sN = sD
        else
            sN = -d + b
            sD = a
        end
    end

    sc = abs(sN) <= tol ? zero(FT) : sN / sD
    tc = abs(tN) <= tol ? zero(FT) : tN / tD
    dP = w + sc * u - tc * v
    return dot(dP, dP)
end

"""
    calc_k_term_fast!(Z_local, mfie, tri_test, tri_source, r_test, r_src)

MFIE 的 K 算子核（论文式 (2-30) 与 (2-22)）：

```math
\\mathcal{K}[\\bm{X}](\\bm{r}) = \\int_S \\bm{X}(\\bm{r}') \\times \\nabla' G(R)\\, dS'
```

测试后展开为（`n̂_t` 为测试三角形外法向）：

```math
\\eta\\, \\frac{l_m l_n}{16\\pi} \\sum_{i,j} w_i w_j\\,
\\Big[(\\bm{\\rho}_m \\cdot \\bm{R})(\\hat{\\bm{n}}_t \\cdot \\bm{\\rho}_n)
- (\\hat{\\bm{n}}_t \\cdot \\bm{R})(\\bm{\\rho}_m \\cdot \\bm{\\rho}_n)\\Big]
\\left({\\rm j}k + \\frac{1}{R}\\right) \\frac{e^{-{\\rm j}kR}}{R^2}
```

实现优化：`(i,j)` 外层循环使 `R`、`e^{-jkR}/R²`、`(jk+1/R)` 每对只算一次，
`ρ_m/ρ_n` 提前提升到最内层循环之外；`η/(16π)` 与边长在最后统一乘。
"""
function calc_k_term_fast!(
    Z_local::AbstractMatrix{CT},
    mfie::MFIE{FT,CT},
    tri_test::TriangleInfo{IT,FT},
    tri_source::TriangleInfo{IT,FT},
    r_test,
    r_src,
    w,
) where {IT,FT,CT}
    n_pts = length(w)

    nt_hat = tri_test.facen̂
    JK_0 = im * mfie.k
    eta_div_16pi = mfie.eta / (16 * FT(π))

    # Precompute free vertices (SVector for stack allocation)
    v_test = SVector{3,SVector{3,FT}}(
        SVector{3,FT}(tri_test.vertices[:, 1]),
        SVector{3,FT}(tri_test.vertices[:, 2]),
        SVector{3,FT}(tri_test.vertices[:, 3]),
    )
    v_src = SVector{3,SVector{3,FT}}(
        SVector{3,FT}(tri_source.vertices[:, 1]),
        SVector{3,FT}(tri_source.vertices[:, 2]),
        SVector{3,FT}(tri_source.vertices[:, 3]),
    )

    # Loop order: (j, i) outer, (m, n) inner
    # This computes rvec/R/divr/temp/Green ONCE per (i,j) pair (7×7 = 49 times)
    # instead of 9× per (i,j) pair (7×7×9 = 441 times in old code)
    @inbounds for j = 1:n_pts
        rgj = r_src[j]
        wj = w[j]
        # Precompute source rho vectors for all 3 BFs
        rho_n1 = rgj - v_src[1]
        rho_n2 = rgj - v_src[2]
        rho_n3 = rgj - v_src[3]
        nt_dot_rho_n1 = dot(nt_hat, rho_n1)
        nt_dot_rho_n2 = dot(nt_hat, rho_n2)
        nt_dot_rho_n3 = dot(nt_hat, rho_n3)

        for i = 1:n_pts
            rgi = r_test[i]
            wi = w[i]

            rvec = rgi - rgj
            R = norm(rvec)
            R < 1e-12 && continue
            divr = one(FT) / R

            # Green's function: exp(-jkR)/R, gradient factor: (jk + 1/R)/R
            @fastmath G_over_R = exp(-JK_0 * R) * divr  # exp(-jkR)/R
            gw_ij = G_over_R * divr * wi * wj   # exp(-jkR)/R² * wi * wj
            temp = (JK_0 + divr) * gw_ij        # (jk + 1/R) * exp(-jkR)/R² * wi * wj

            # Precompute test rho and shared dot products
            rho_m1 = rgi - v_test[1]
            rho_m2 = rgi - v_test[2]
            rho_m3 = rgi - v_test[3]
            nt_dot_rvec = dot(nt_hat, rvec)

            # Accumulate all 3×3 (m,n) contributions at this (i,j)
            # Zmn += ((ρm·rvec)(n̂t·ρn) - (n̂t·rvec)(ρm·ρn)) * temp
            for (ni, rho_n, nt_rho_n) in
                ((1, rho_n1, nt_dot_rho_n1), (2, rho_n2, nt_dot_rho_n2), (3, rho_n3, nt_dot_rho_n3))
                for (mi, rho_m) in ((1, rho_m1), (2, rho_m2), (3, rho_m3))
                    term1 = dot(rho_m, rvec) * nt_rho_n
                    term2 = nt_dot_rvec * dot(rho_m, rho_n)
                    Z_local[mi, ni] += (term1 - term2) * temp
                end
            end
        end
    end

    # Apply edge lengths and constant factor
    @inbounds for n = 1:3
        ln = tri_source.edgel[n]
        for m = 1:3
            lm = tri_test.edgel[m]
            Z_local[m, n] *= lm * ln * eta_div_16pi
        end
    end

    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# PMCHWT support: principal-value K-operator (no mass-matrix self-term)
# ─────────────────────────────────────────────────────────────────────────────

"""
    assemble_K_offdiag(basis::RWGBasis, k) -> Matrix{Complex}

Assemble the **off-diagonal** (principal-value) K-operator matrix with
overall factor ``1 / (16\\pi)``.

The self-interaction (diagonal mass-matrix) term is deliberately omitted.
This is the building block for the PMCHWT ``Z^{EM}`` and ``Z^{HJ}`` blocks.

# Arguments
- `basis`: RWG basis on a closed triangular mesh
- `k`:     Wavenumber of the propagation medium (real for lossless)

# Returns
- ``N \\times N`` complex matrix ``K``, where the diagonal is identically zero.

# Note
`assemble_K_offdiag(basis, k0) + assemble_K_offdiag(basis, k1)` gives the
combined K-block needed for PMCHWT interior/exterior superposition.
"""
function assemble_K_offdiag(basis::RWGBasis{IT,FT}, k::FT) where {IT,FT}
    CT = Complex{FT}
    # eta = 1.0  →  eta_div_16pi = 1.0 / (16π), the pure K factor for PMCHWT
    gq_far = GaussQuadratureInfo(:Triangle, 4, FT)
    gq_near = GaussQuadratureInfo(:Triangle, 7, FT)
    mfie_k = MFIE{FT,CT}(zero(FT), k, one(FT), gq_far, gq_near)

    # Precompute quadrature points (same as the standard MFIE assembly)
    mesh = basis.mesh
    nt = num_elements(mesh)
    N_points = length(gq_far.weight)

    quad_points = Vector{SVector{N_points,SVector{3,FT}}}(undef, nt)
    Threads.@threads for t = 1:nt
        v_idx = mesh.triangles[:, t]
        v1 = SVector{3,FT}(mesh.node[:, v_idx[1]])
        v2 = SVector{3,FT}(mesh.node[:, v_idx[2]])
        v3 = SVector{3,FT}(mesh.node[:, v_idx[3]])
        quad_points[t] = SVector{N_points,SVector{3,FT}}(
            v1 * gq_far.coordinate[1, i] +
            v2 * gq_far.coordinate[2, i] +
            v3 * gq_far.coordinate[3, i] for i = 1:N_points
        )
    end

    # Interaction wrapper: skip self-term (do not add mass matrix)
    function k_offdiag_interaction!(Z_local, op, t_test, t_src, qpts)
        if t_test.triID == t_src.triID
            return nothing          # no self-term in PMCHWT K-block
        end
        r_test = qpts[t_test.triID]
        r_src = qpts[t_src.triID]
        calc_k_term_fast!(Z_local, op, t_test, t_src, r_test, r_src, op.gq_far.weight)
        return nothing
    end

    wrapper = (Z, op, t1, t2) -> k_offdiag_interaction!(Z, op, t1, t2, quad_points)
    # K is NOT symmetric (no symmetry exploitation)
    return assemble_generic(mfie_k, basis, wrapper, symmetric = false)
end

end
