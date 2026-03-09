module EFIEModule

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..Kernels
using StaticArrays
using LinearAlgebra
using SparseArrays

using ..Impedance

# Include Singularities module
include("Singularities.jl")
using .Singularities

import ..CoreModule: assemble_impedance_matrix

export EFIE, assemble_impedance_matrix, assemble_impedance_matrix!, efie_interaction, efie_interaction!, efie_from_keta

"""
    EFIE{FT, CT} <: AbstractIntegralOperator

Electric Field Integral Equation (EFIE) operator.

Solves for the surface current density \$\\mathbf{J}\$ on a PEC object by enforcing the boundary condition \$\\hat{n} \\times \\mathbf{E}^{scat} = -\\hat{n} \\times \\mathbf{E}^{inc}\$.

# Mathematical Formulation

The impedance matrix element \$Z_{mn}\$ is derived using Galerkin's method with RWG basis functions:

```math
Z_{mn} = j\\omega\\mu \\int_{S_m} \\int_{S_n} \\mathbf{f}_m(\\mathbf{r}) \\cdot \\mathbf{f}_n(\\mathbf{r}') G(\\mathbf{r}, \\mathbf{r}') dS' dS + \\frac{1}{j\\omega\\epsilon} \\int_{S_m} \\int_{S_n} (\\nabla_s \\cdot \\mathbf{f}_m(\\mathbf{r})) (\\nabla_s \\cdot \\mathbf{f}_n(\\mathbf{r}')) G(\\mathbf{r}, \\mathbf{r}') dS' dS
```

where:
- \$\\mathbf{f}_m, \\mathbf{f}_n\$ are the testing and source RWG basis functions.
- \$G(\\mathbf{r}, \\mathbf{r}') = \\frac{e^{-jk|\\mathbf{r}-\\mathbf{r}'|}}{4\\pi|\\mathbf{r}-\\mathbf{r}'|}\$ is the free-space Green's function.
- \$S_m, S_n\$ are the supports of the testing and source functions.

# Fields
- `freq`: Operating frequency (Hz).
- `k`: Wavenumber (\$k = \\omega\\sqrt{\\mu\\epsilon}\$).
- `eta`: Intrinsic impedance (\$\\eta = \\sqrt{\\mu/\\epsilon}\$).
- `gq_info`: Gauss quadrature information for numerical integration.
- `C4divk2`: Precomputed constant \$4/k^2\$ (used in some formulations).
- `factor`: Precomputed scaling factor \$j k \\eta / (16\\pi)\$ (varies by implementation).
"""
struct EFIE{FT<:AbstractFloat,CT<:Complex,N_FAR,N_NEAR} <: AbstractIntegralOperator
    freq::FT
    k::FT
    eta::FT
    gq_far::GaussQuadratureInfoStruct{FT,N_FAR,3}
    gq_near::GaussQuadratureInfoStruct{FT,N_NEAR,3}
    C4divk2::FT
    factor::CT
    SSCg::Vector{CT}
end

function EFIE(freq::FT) where {FT}
    c0 = 299792458.0
    mu0 = 4π * 1e-7
    eps0 = 1.0 / (c0^2 * mu0)
    k = 2π * freq / c0
    eta = sqrt(mu0 / eps0)

    gq_far = GaussQuadratureInfo(:Triangle, 4, FT)
    gq_near = GaussQuadratureInfo(:Triangle, 7, FT)
    C4divk2 = 4 / k^2
    # Standard EFIE factor: jk * eta / 16pi (Legacy Parity)
    factor = im * k * eta / (16 * π)
    SSCg = compute_SSCg(k)

    return EFIE{FT,Complex{FT},4,7}(freq, k, eta, gq_far, gq_near, C4divk2, factor, SSCg)
end

"""
    efie_from_keta(k::FT, eta::FT, factor::Complex{FT}) where {FT}

Low-level constructor for EFIE-like operators with an **explicit** wavenumber,
intrinsic impedance, and overall scaling factor.

This is used internally by the PMCHWT formulation to construct L-operators for
interior/exterior regions and for the ``Z^{HM}`` block (which uses a different
factor from the standard EFIE).

# Arguments
- `k`:      Wavenumber of the propagation medium (must be real ≥ 0)
- `eta`:    Intrinsic impedance ``\\eta = \\sqrt{\\mu/\\varepsilon}`` (for record-keeping)
- `factor`: Scaling factor applied after numerical integration.
  - For ``Z^{EJ}`` block:  ``factor = jk\\eta/(16\\pi)``
  - For ``Z^{HM}`` block:  ``factor = jk/(\\eta \\cdot 16\\pi)``

# Returns
`EFIE{FT, Complex{FT}, 4, 7}` operator that can be passed to
`assemble_impedance_matrix`.
"""
function efie_from_keta(k::FT, eta::FT, factor::Complex{FT}) where {FT}
    gq_far  = GaussQuadratureInfo(:Triangle, 4, FT)
    gq_near = GaussQuadratureInfo(:Triangle, 7, FT)
    C4divk2 = FT(4) / k^2
    SSCg    = compute_SSCg(k)  # Singularities.compute_SSCg accessible here
    return EFIE{FT,Complex{FT},4,7}(FT(0), k, eta, gq_far, gq_near, C4divk2, factor, SSCg)
end
function _assemble_impedance_matrix!(
    Z::AbstractMatrix{CT},
    efie::EFIE{FT,CT},
    basis::RWGBasis{IT,FT};
    accumulate::Bool = false,
) where {IT,FT,CT}
    # Precompute quadrature points for Far Field
    mesh = basis.mesh
    nt = num_elements(mesh)
    gq = efie.gq_far
    N_points = length(gq.weight)

    # Precompute points
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

    # Wrapper
    interaction_wrapper = (Z, op, t1, t2) -> efie_interaction!(Z, op, t1, t2, quad_points)

    return assemble_generic!(Z, efie, basis, interaction_wrapper; symmetric = true, accumulate)
end

function assemble_impedance_matrix!(
    Z::AbstractMatrix{CT},
    efie::EFIE{FT,CT},
    basis::RWGBasis{IT,FT};
    accumulate::Bool = false,
) where {IT,FT,CT}
    return _assemble_impedance_matrix!(Z, efie, basis; accumulate)
end

function assemble_impedance_matrix(efie::EFIE{FT,CT}, basis::RWGBasis{IT,FT}) where {IT,FT,CT}
    Z = zeros(CT, num_basis(basis), num_basis(basis))
    return _assemble_impedance_matrix!(Z, efie, basis)
end

function efie_interaction!(
    Z_local::AbstractMatrix{CT},
    efie::EFIE{FT,CT},
    tri_test::TriangleInfo{IT,FT},
    tri_source::TriangleInfo{IT,FT},
) where {IT,FT,CT}
    # Slow path: compute points on the fly
    r_test = get_global_quad_points(tri_test, efie.gq_far)
    r_src = get_global_quad_points(tri_source, efie.gq_far)
    efie_interaction!(Z_local, efie, tri_test, tri_source, r_test, r_src)
end

function efie_interaction!(
    Z_local::Matrix{CT},
    efie::EFIE{FT,CT},
    tri_test::TriangleInfo{IT,FT},
    tri_source::TriangleInfo{IT,FT},
    quad_points::Vector,
) where {IT,FT,CT}
    # Fast path: use precomputed points
    r_test = quad_points[tri_test.triID]
    r_src = quad_points[tri_source.triID]
    efie_interaction!(Z_local, efie, tri_test, tri_source, r_test, r_src)
end

function efie_interaction!(
    Z_local::AbstractMatrix{CT},
    efie::EFIE{FT,CT},
    tri_test::TriangleInfo{IT,FT},
    tri_source::TriangleInfo{IT,FT},
    r_test,
    r_src,
) where {IT,FT,CT}
    if tri_test.triID == tri_source.triID
        # Self-term (Singular)
        calc_self_interaction!(Z_local, efie, tri_test)
    elseif is_adjacent(tri_test, tri_source)
        # Near-term (Singular/Near-Singular)
        # Enforce T_test.ID < T_src.ID to match Direct Solver and ensure symmetry
        # calc_near_interaction! is asymmetric due to semi-analytical integration.
        if tri_test.triID > tri_source.triID
            Z_temp = zeros(CT, 3, 3)
            calc_near_interaction!(Z_temp, efie, tri_source, tri_test)
            # Transpose back: Z_local[m, n] += Z_temp[n, m]
            for m = 1:3, n = 1:3
                Z_local[m, n] += Z_temp[n, m]
            end
        else
            calc_near_interaction!(Z_local, efie, tri_test, tri_source)
        end
    else
        # Far-term
        calc_interaction!(Z_local, efie, tri_test, tri_source, r_test, r_src)
    end

    # Apply factor
    Z_local .*= efie.factor

    if abs(Z_local[1, 1]) > 0
        # println("EFIE Z_local: $(Z_local[1,1])")
    end

    return nothing
end



function calc_self_interaction!(
    Z_local::AbstractMatrix{CT},
    efie::EFIE{FT,CT},
    tri::TriangleInfo{IT,FT},
) where {IT,FT,CT}
    gq = efie.gq_near
    C4divk2 = efie.C4divk2

    r_pts = get_global_quad_points(tri, gq)
    w_pts = gq.weight

    # Precompute F1 (Singular integral of 1/R)
    # Note: singularF1 takes edge lengths.
    # tri.edgel is [l1, l2, l3].
    sF1 = singularF1(tri.edgel[1], tri.edgel[2], tri.edgel[3])
    F1 = C4divk2 * sF1
    # println("DEBUG: sF1=", sF1, " C4divk2=", C4divk2, " F1=", F1)

    # Loop over quadrature points for Smooth part
    for j = 1:length(w_pts)
        rj = r_pts[j]
        wj = w_pts[j] * tri.area

        # Source basis function centers (vertices)
        v_src = tri.vertices

        for i = 1:length(w_pts)
            ri = r_pts[i]
            wi = w_pts[i] * tri.area

            R = norm(ri - rj)

            # Smooth Green's function: (e^-jkr - 1) / R
            # Use Taylor expansion for consistency with Legacy
            G_smooth = greenfunc_star(R, efie.k)

            # Basis function rho vectors
            # rho_n = r - v_n
            rho_src_1 = rj - v_src[:, 1]
            rho_src_2 = rj - v_src[:, 2]
            rho_src_3 = rj - v_src[:, 3]
            rhos_src = (rho_src_1, rho_src_2, rho_src_3)

            rho_test_1 = ri - v_src[:, 1]
            rho_test_2 = ri - v_src[:, 2]
            rho_test_3 = ri - v_src[:, 3]
            rhos_test = (rho_test_1, rho_test_2, rho_test_3)

            for n = 1:3
                for m = 1:3
                    # Term 1: Vector potential part (rho . rho')
                    term1 = dot(rhos_test[m], rhos_src[n])
                    # Term 2: Scalar potential part (-4/k^2)
                    term2 = -C4divk2

                    # For self-term, we subtract the singular part from the numerical integration
                    # The singular part of G is 1/R.
                    # So we integrate (G - 1/R) numerically.
                    # G - 1/R = G_smooth.

                    val = (term1 + term2) * G_smooth

                    Z_local[m, n] += val * wi * wj
                end
            end
        end
    end

    # Add Singular Part (Analytical)
    area2 = tri.area^2

    for n = 1:3
        for m = 1:3
            # We need to add \int \int (term1 + term2) * (1/R)
            # = \int \int (rho_m . rho_n) / R - C4divk2 * \int \int 1/R
            # = F2 - F1

            val_singular = zero(FT)

            if m == n
                # F21
                # Edges a, b, c starting from m
                # tri.edgel indices: [1, 2, 3]
                # If m=1: a=l1, b=l2, c=l3
                # If m=2: a=l2, b=l3, c=l1
                # If m=3: a=l3, b=l1, c=l2

                idx_a = m
                idx_b = mod1(m + 1, 3)
                idx_c = mod1(m + 2, 3)

                a = tri.edgel[idx_a]
                b = tri.edgel[idx_b]
                c = tri.edgel[idx_c]

                val_singular = singularF21(a, b, c, area2) - F1
            else
                # F22
                # Edges a, b, c. a is opposite to the third vertex.
                # Third vertex index k = 6 - m - n.
                k = 6 - m - n

                idx_a = k
                idx_b = mod1(k + 1, 3)
                idx_c = mod1(k + 2, 3)

                a = tri.edgel[idx_a]
                b = tri.edgel[idx_b]
                c = tri.edgel[idx_c]

                val_singular = singularF22(a, b, c, area2) - F1
            end

            # singularF1/F2 return values normalized by Area (I/A).
            # Z_local accumulates the full integral (I).
            # So we must multiply by Area^2.
            # WAIT: Legacy code does NOT multiply by Area^2 here?
            # Legacy: Ztemp += singularF21(...) - F1
            # Legacy: Ztemp *= lm*ln*JKηdiv16π
            # Legacy singularF21 returns value already scaled?
            # Let's check Legacy singularF21 implementation.
            # It returns ( ... ) / 30.
            # It seems it returns the integral value directly?
            # No, the formula involves area2.

            # Let's look at Legacy EFIEOnTris again.
            # Ztemp += singularF21(...) - F1
            # Ztri[mi, ni] = Ztemp * lm * ln * factor

            # In EMSuite, we are adding to Z_local which is accumulating the integral.
            # The smooth part is accumulating val * wi * wj.
            # wi = w_pts[i] * tri.area.
            # So smooth part is scaled by Area^2.

            # If singularF21 returns the integral value directly (not normalized), then we don't need to multiply by Area^2.
            # But if it returns I/Area^2 or something, we do.

            # Legacy singularF1: -4 * ( ... ) / 3.
            # This is \int \int 1/R dS dS'.
            # Wait, dimensionally:
            # 1/R ~ 1/L. dS dS' ~ L^4. Integral ~ L^3.
            # Formula: a * log(...) ~ L.
            # So singularF1 returns L.
            # But we expect L^3.
            # Ah, the formula for 1/R integral over two triangles is proportional to Area^2 / avg_R? No.

            # Let's check Eibert's paper or similar.
            # The analytical formulas usually give the integral value.
            # But wait, Legacy code passes `tri.area^2` to singularF21.

            # CRITICAL: In Legacy `EFIEOnTris`:
            # Ztemp += singularF21(...) - F1
            # It does NOT multiply by Area^2 after this.
            # It implies singularF21 and F1 return the FULL integral value.

            # However, in EMSuite `calc_self_interaction!`:
            # Z_local[m, n] += val_singular * area2
            # I added `* area2` because I thought they were normalized.
            # If Legacy doesn't multiply, maybe I shouldn't either.

            # Revert: Removing empirical factor 0.375.
            # We must find the root cause of the magnitude discrepancy.
            Z_local[m, n] += val_singular * area2
        end
    end

    # Multiply by edge lengths and divide by areas
    inv_areas = 1.0 / (tri.area * tri.area)
    for n = 1:3
        ln = tri.edgel[n]
        for m = 1:3
            lm = tri.edgel[m]
            Z_local[m, n] *= lm * ln * inv_areas
        end
    end

    return nothing
end

function is_adjacent(t1::TriangleInfo, t2::TriangleInfo)
    # Check if they share any vertices
    for i = 1:3
        v1 = t1.verticesID[i]
        for j = 1:3
            if v1 == t2.verticesID[j]
                return true
            end
        end
    end
    return false
end

function calc_near_interaction!(
    Z_local::AbstractMatrix{CT},
    efie::EFIE{FT,CT},
    tri_test::TriangleInfo{IT,FT},
    tri_source::TriangleInfo{IT,FT},
) where {IT,FT,CT}
    gq = efie.gq_near
    C4divk2 = efie.C4divk2

    # Quadrature points for Test Triangle
    r_test = get_global_quad_points(tri_test, gq)
    w_test = gq.weight

    for i = 1:length(w_test)
        ri = r_test[i]
        wi = w_test[i] * tri_test.area

        Ig, IvecSg = faceSingularityIgIvecg(
            ri,
            tri_source.vertices,
            tri_source.edgel,
            tri_source.edgev̂,
            tri_source.edgen̂,
            tri_source.area,
            tri_source.facen̂,
            efie.SSCg,
        )

        rho_test_1 = ri - tri_test.vertices[:, 1]
        rho_test_2 = ri - tri_test.vertices[:, 2]
        rho_test_3 = ri - tri_test.vertices[:, 3]
        rhos_test = (rho_test_1, rho_test_2, rho_test_3)

        rho_src_at_ri_1 = ri - tri_source.vertices[:, 1]
        rho_src_at_ri_2 = ri - tri_source.vertices[:, 2]
        rho_src_at_ri_3 = ri - tri_source.vertices[:, 3]
        rhos_src_at_ri = (rho_src_at_ri_1, rho_src_at_ri_2, rho_src_at_ri_3)

        for n = 1:3
            for m = 1:3
                # Legacy Parity: EFIEOnNearTris formula (EFIERWGTri.jl line 139):
                # Ztemp += ((ρm·ρn - C4divk²)*Ig - ρm·Ivec) * weight[gi] / tris.area
                # 
                # val = (ρm·ρn - C4/k²)*Ig - ρm·IvecSg
                # where IvecSg = ∫(r-r')G dS' (so -ρm·IvecSg gives the right sign)
                term1_int = -dot(rhos_test[m], IvecSg) + dot(rhos_test[m], rhos_src_at_ri[n]) * Ig
                term2_int = -C4divk2 * Ig

                val = term1_int + term2_int

                Z_local[m, n] += val * wi
            end
        end
    end

    # Multiply by edge lengths and divide by areas
    # Legacy Parity: EFIEOnNearTris does:
    #   Ztemp * weight[gi] / tris.area  (divides by A_src in loop)
    #   Ztemp *= lm * ln * JKηdiv16π   (no further area division)
    #
    # EMSuite: wi = w[i] * A_test, so accumulated = A_test * Σ w_i * val
    # Need to divide by A_test (cancel quadrature Jacobian) AND by A_src
    # (source basis normalization, matching Legacy's / tris.area)
    inv_areas = 1.0 / (tri_test.area * tri_source.area)
    for n = 1:3
        ln = tri_source.edgel[n]
        for m = 1:3
            lm = tri_test.edgel[m]
            Z_local[m, n] *= lm * ln * inv_areas
        end
    end

    return nothing
end

function calc_interaction!(
    Z_local::AbstractMatrix{CT},
    efie::EFIE{FT,CT},
    tri_test::TriangleInfo{IT,FT},
    tri_source::TriangleInfo{IT,FT},
    r_test,
    r_src,
) where {IT,FT,CT}
    gq = efie.gq_far
    C4divk2 = efie.C4divk2
    k = efie.k

    w_src = gq.weight
    w_test = gq.weight
    n_src = length(w_src)
    n_test = length(w_test)

    v_src_1 = tri_source.vertices[:, 1]
    v_src_2 = tri_source.vertices[:, 2]
    v_src_3 = tri_source.vertices[:, 3]
    v_test_1 = tri_test.vertices[:, 1]
    v_test_2 = tri_test.vertices[:, 2]
    v_test_3 = tri_test.vertices[:, 3]

    # Loop over quadrature points
    # Optimization: fully unrolled 3×3 dot products, no ternary branching,
    # tuple-indexed rho vectors for zero-cost access.
    @inbounds for j = 1:n_src
        rj = r_src[j]
        wj = w_src[j] * tri_source.area

        # Precompute source rho vectors as tuple (zero alloc, stack only)
        rho_src = (rj - v_src_1, rj - v_src_2, rj - v_src_3)

        for i = 1:n_test
            ri = r_test[i]
            wi = w_test[i] * tri_test.area

            # Green's function: exp(-jkR)/R
            diff = ri - rj
            R = sqrt(diff[1]^2 + diff[2]^2 + diff[3]^2)
            @fastmath G = exp(-im * k * R) / R

            # Precompute test rho vectors as tuple
            rho_test = (ri - v_test_1, ri - v_test_2, ri - v_test_3)

            val_common = G * wi * wj

            # Fully unrolled 3×3 accumulation — no branching, no @simd overhead
            for n = 1:3
                rho_n = rho_src[n]
                for m = 1:3
                    Z_local[m, n] += (dot(rho_test[m], rho_n) - C4divk2) * val_common
                end
            end
        end
    end

    # Multiply by edge lengths and divide by areas
    inv_areas = 1.0 / (tri_test.area * tri_source.area)
    @inbounds for n = 1:3
        ln = tri_source.edgel[n]
        for m = 1:3
            lm = tri_test.edgel[m]
            Z_local[m, n] *= lm * ln * inv_areas
        end
    end

    return nothing
end

end
