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
    gq_info::GaussQuadratureInfoStruct{FT, 7, 3}
end

function MFIE(freq::FT) where {FT}
    c0 = 299792458.0
    mu0 = 4π * 1e-7
    eps0 = 1.0 / (c0^2 * mu0)
    k = 2π * freq / c0
    eta = sqrt(mu0 / eps0)
    
    gq_info = GaussQuadratureInfo(:Triangle, 7, FT)
    
    return MFIE{FT, Complex{FT}}(freq, k, eta, gq_info)
end

"""
    assemble_impedance_matrix(mfie::MFIE, basis::RWGBasis)

Assemble the impedance matrix Z for the MFIE using RWG basis functions.
Z = eta * (0.5 * M + K)
"""
function assemble_impedance_matrix(mfie::MFIE{FT, CT}, basis::RWGBasis{IT, FT}) where {IT, FT, CT}
    return assemble_generic(mfie, basis, mfie_interaction!)
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

function calc_self_term!(Z_local::AbstractMatrix{CT}, mfie::MFIE{FT, CT}, tri::TriangleInfo{IT, FT}) where {IT, FT, CT}
    # Self-term is just the Mass Matrix term scaled by eta
    # Z_self = eta * 0.5 * <fm, fn>
    # Based on legacy code: Z = (eta / 8A) * lm * ln * sum(rho_m . rho_n * w)
    
    gq = mfie.gq_info
    r_quad = get_global_quad_points(tri, gq)
    w_quad = gq.weight
    
    # Precompute rho vectors at quad points
    # rho_i = r - v_i
    rhos = [Vector{SVector{3, FT}}(undef, length(w_quad)) for _ in 1:3]
    
    for k in 1:length(w_quad)
        rk = r_quad[k]
        rhos[1][k] = rk - tri.vertices[:, 1]
        rhos[2][k] = rk - tri.vertices[:, 2]
        rhos[3][k] = rk - tri.vertices[:, 3]
    end
    
    eta_div_8A = mfie.eta / (8 * tri.area)
    
    for m in 1:3
        lm = tri.edgel[m]
        # Check if basis function exists
        if tri.inBfsID[m] == 0
            continue
        end
        
        for n in 1:3
            ln = tri.edgel[n]
            if tri.inBfsID[n] == 0
                continue
            end
            
            sum_val = zero(CT)
            for k in 1:length(w_quad)
                sum_val += dot(rhos[m][k], rhos[n][k]) * w_quad[k]
            end
            
            Z_local[m, n] += sum_val * lm * ln * eta_div_8A
        end
    end
    
    return nothing
end

function calc_k_term!(Z_local::AbstractMatrix{CT}, mfie::MFIE{FT, CT}, tri_test::TriangleInfo{IT, FT}, tri_source::TriangleInfo{IT, FT}) where {IT, FT, CT}
    # K-term calculation based on legacy MFIERWGTri.jl
    # Zmn = eta * <fm, n x (fn x grad G)>
    
    gq = mfie.gq_info
    r_test = get_global_quad_points(tri_test, gq)
    w_test = gq.weight
    r_src = get_global_quad_points(tri_source, gq)
    w_src = gq.weight
    
    nt = tri_test.facen̂
    ns = tri_source.facen̂
    
    JK_0 = im * mfie.k
    eta_div_16pi = mfie.eta / (16 * π)
    
    # Precompute Green's function and weights
    gw = zeros(CT, length(w_test), length(w_src))
    for j in 1:length(w_src)
        for i in 1:length(w_test)
            R_vec = r_test[i] - r_src[j]
            R = norm(R_vec)
            if R > 1e-12
                gw[i, j] = (exp(-JK_0 * R) / R) * w_test[i] * w_src[j]
            end
        end
    end
    
    for n in 1:3
        ln = tri_source.edgel[n]
        if tri_source.inBfsID[n] == 0; continue; end
        vn = tri_source.vertices[:, n] # Free vertex for source
        
        for m in 1:3
            lm = tri_test.edgel[m]
            if tri_test.inBfsID[m] == 0; continue; end
            vm = tri_test.vertices[:, m] # Free vertex for test
            
            Zmn = zero(CT)
            
            for j in 1:length(w_src)
                rgj = r_src[j]
                rho_n = rgj - vn
                
                for i in 1:length(w_test)
                    rgi = r_test[i]
                    rho_m = rgi - vm
                    
                    rvec = rgi - rgj
                    R = norm(rvec)
                    if R < 1e-12; continue; end
                    
                    divr = 1.0 / R
                    temp = (JK_0 + divr) * divr * gw[i, j]
                    
                    # Legacy expansion:
                    # Zmn += ((ρmi ⋅ rvec) * (n̂t ⋅ ρnj) - (n̂t ⋅ rvec) * (ρmi ⋅ ρnj)) * temp
                    
                    term1 = dot(rho_m, rvec) * dot(nt, rho_n)
                    term2 = dot(nt, rvec) * dot(rho_m, rho_n)
                    
                    Zmn += (term1 - term2) * temp
                end
            end
            
            Z_local[m, n] += Zmn * lm * ln * eta_div_16pi
        end
    end
    
    return nothing
end

function get_global_quad_points(tri::TriangleInfo{IT, FT}, gq::GaussQuadratureInfoStruct{FT}) where {IT, FT}
    N = length(gq.weight)
    points = Vector{SVector{3, FT}}(undef, N)
    
    v1 = tri.vertices[:, 1]
    v2 = tri.vertices[:, 2]
    v3 = tri.vertices[:, 3]
    
    for i in 1:N
        l1 = gq.coordinate[1, i]
        l2 = gq.coordinate[2, i]
        l3 = gq.coordinate[3, i]
        
        points[i] = v1 * l1 + v2 * l2 + v3 * l3
    end
    
    return points
end

end
