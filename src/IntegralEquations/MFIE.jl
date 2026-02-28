module MFIEModule

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..Kernels
using ..Impedance
using StaticArrays
using LinearAlgebra
using SparseArrays

import ..CoreModule: assemble_impedance_matrix

export MFIE, assemble_impedance_matrix

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
- `gq_info`: Gauss quadrature information.
"""
struct MFIE{FT<:AbstractFloat, CT<:Complex} <: AbstractIntegralOperator
    freq::FT
    k::FT
    eta::FT
    gq_info::GaussQuadratureInfoStruct{FT, 4, 3}  # 4-point Gauss (matches Legacy GQPNTri=4)
end

function MFIE(freq::FT) where {FT}
    c0 = 299792458.0
    mu0 = 4π * 1e-7
    eps0 = 1.0 / (c0^2 * mu0)
    k = 2π * freq / c0
    eta = sqrt(mu0 / eps0)
    
    # Legacy uses GQPNTri=4 for MFIE (non-singular K-operator).
    # 4-point integrates degree ≤ 3 polynomials exactly on triangles,
    # sufficient for both K-operator and mass matrix self-term.
    gq_info = GaussQuadratureInfo(:Triangle, 4, FT)
    
    return MFIE{FT, Complex{FT}}(freq, k, eta, gq_info)
end

"""
    assemble_impedance_matrix(mfie::MFIE, basis::RWGBasis)

Assemble the impedance matrix Z for the MFIE using RWG basis functions.
Z = eta * (0.5 * M + K)

Quad points are precomputed once to avoid per-pair heap allocations.
"""
function assemble_impedance_matrix(mfie::MFIE{FT, CT}, basis::RWGBasis{IT, FT}) where {IT, FT, CT}
    # Precompute quadrature points for all triangles (like EFIE does)
    mesh = basis.mesh
    nt = num_elements(mesh)
    gq = mfie.gq_info
    N_points = length(gq.weight)

    quad_points = Vector{SVector{N_points, SVector{3, FT}}}(undef, nt)
    Threads.@threads for t in 1:nt
        v_indices = mesh.triangles[:, t]
        v1 = SVector{3, FT}(mesh.node[:, v_indices[1]])
        v2 = SVector{3, FT}(mesh.node[:, v_indices[2]])
        v3 = SVector{3, FT}(mesh.node[:, v_indices[3]])
        quad_points[t] = SVector{N_points, SVector{3, FT}}(
            v1 * gq.coordinate[1, i] + v2 * gq.coordinate[2, i] + v3 * gq.coordinate[3, i]
            for i in 1:N_points
        )
    end

    # Wrapper with precomputed points
    wrapper = (Z, op, t1, t2) -> mfie_interaction!(Z, op, t1, t2, quad_points)
    return assemble_generic(mfie, basis, wrapper)
end

function mfie_interaction!(Z_local::AbstractMatrix{CT}, mfie::MFIE{FT, CT}, tri_test::TriangleInfo{IT, FT}, tri_source::TriangleInfo{IT, FT}) where {IT, FT, CT}
    # Check if self-interaction
    if tri_test.triID == tri_source.triID
        calc_self_term!(Z_local, mfie, tri_test)
    else
        calc_k_term!(Z_local, mfie, tri_test, tri_source)
    end
    return nothing
end

# Fast path with precomputed quad points (eliminates per-pair allocation)
function mfie_interaction!(Z_local::AbstractMatrix{CT}, mfie::MFIE{FT, CT}, tri_test::TriangleInfo{IT, FT}, tri_source::TriangleInfo{IT, FT}, quad_points::Vector) where {IT, FT, CT}
    if tri_test.triID == tri_source.triID
        calc_self_term!(Z_local, mfie, tri_test)
    else
        r_test = quad_points[tri_test.triID]
        r_src = quad_points[tri_source.triID]
        calc_k_term_fast!(Z_local, mfie, tri_test, tri_source, r_test, r_src)
    end
    return nothing
end

function calc_self_term!(Z_local::AbstractMatrix{CT}, mfie::MFIE{FT, CT}, tri::TriangleInfo{IT, FT}) where {IT, FT, CT}
    # Self-term is just the Mass Matrix term scaled by eta
    # Z_self = eta * 0.5 * <fm, fn>
    # Based on legacy code: Z = (eta / 8A) * lm * ln * sum(rho_m . rho_n * w)
    
    gq = mfie.gq_info
    r_quad = get_global_quad_points(tri, gq)
    w_quad = gq.weight
    n_pts = length(w_quad)
    
    # Precompute vertices as SVector for stack allocation
    v1 = SVector{3, FT}(tri.vertices[:, 1])
    v2 = SVector{3, FT}(tri.vertices[:, 2])
    v3 = SVector{3, FT}(tri.vertices[:, 3])
    
    eta_div_8A = mfie.eta / (8 * tri.area)
    
    for m in 1:3
        lm = tri.edgel[m]
        if tri.inBfsID[m] == 0
            continue
        end
        
        for n in 1:3
            ln = tri.edgel[n]
            if tri.inBfsID[n] == 0
                continue
            end
            
            sum_val = zero(CT)
            @inbounds for k in 1:n_pts
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

function calc_k_term!(Z_local::AbstractMatrix{CT}, mfie::MFIE{FT, CT}, tri_test::TriangleInfo{IT, FT}, tri_source::TriangleInfo{IT, FT}) where {IT, FT, CT}
    # Slow path: compute quad points on the fly
    gq = mfie.gq_info
    r_test = get_global_quad_points(tri_test, gq)
    r_src = get_global_quad_points(tri_source, gq)
    calc_k_term_fast!(Z_local, mfie, tri_test, tri_source, r_test, r_src)
end

"""
    calc_k_term_fast!(Z_local, mfie, tri_test, tri_source, r_test, r_src)

Optimized K-term kernel. Key improvements over original:
1. Loop order: (i,j) outer, (m,n) inner → rvec/R/divr/temp computed once per (i,j)
2. Precomputed quad points → zero per-pair allocations
3. Pre-hoisted rho_m/rho_n vectors outside innermost loop
"""
function calc_k_term_fast!(Z_local::AbstractMatrix{CT}, mfie::MFIE{FT, CT}, tri_test::TriangleInfo{IT, FT}, tri_source::TriangleInfo{IT, FT}, r_test, r_src) where {IT, FT, CT}
    gq = mfie.gq_info
    w = gq.weight
    n_pts = length(w)

    nt_hat = tri_test.facen̂
    JK_0 = im * mfie.k
    eta_div_16pi = mfie.eta / (16 * FT(π))

    # Precompute free vertices (SVector for stack allocation)
    v_test = SVector{3, SVector{3, FT}}(
        SVector{3, FT}(tri_test.vertices[:, 1]),
        SVector{3, FT}(tri_test.vertices[:, 2]),
        SVector{3, FT}(tri_test.vertices[:, 3])
    )
    v_src = SVector{3, SVector{3, FT}}(
        SVector{3, FT}(tri_source.vertices[:, 1]),
        SVector{3, FT}(tri_source.vertices[:, 2]),
        SVector{3, FT}(tri_source.vertices[:, 3])
    )

    # Loop order: (j, i) outer, (m, n) inner
    # This computes rvec/R/divr/temp/Green ONCE per (i,j) pair (7×7 = 49 times)
    # instead of 9× per (i,j) pair (7×7×9 = 441 times in old code)
    @inbounds for j in 1:n_pts
        rgj = r_src[j]
        wj = w[j]
        # Precompute source rho vectors for all 3 BFs
        rho_n1 = rgj - v_src[1]
        rho_n2 = rgj - v_src[2]
        rho_n3 = rgj - v_src[3]
        nt_dot_rho_n1 = dot(nt_hat, rho_n1)
        nt_dot_rho_n2 = dot(nt_hat, rho_n2)
        nt_dot_rho_n3 = dot(nt_hat, rho_n3)

        for i in 1:n_pts
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
            for (ni, rho_n, nt_rho_n) in ((1, rho_n1, nt_dot_rho_n1),
                                           (2, rho_n2, nt_dot_rho_n2),
                                           (3, rho_n3, nt_dot_rho_n3))
                for (mi, rho_m) in ((1, rho_m1), (2, rho_m2), (3, rho_m3))
                    term1 = dot(rho_m, rvec) * nt_rho_n
                    term2 = nt_dot_rvec * dot(rho_m, rho_n)
                    Z_local[mi, ni] += (term1 - term2) * temp
                end
            end
        end
    end

    # Apply edge lengths and constant factor
    @inbounds for n in 1:3
        ln = tri_source.edgel[n]
        for m in 1:3
            lm = tri_test.edgel[m]
            Z_local[m, n] *= lm * ln * eta_div_16pi
        end
    end

    return nothing
end

end
