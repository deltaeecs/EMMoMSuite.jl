using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.FastAlgorithms.MLFMA
using EMSuite.FastAlgorithms.Lebedev
using LinearAlgebra
using StaticArrays
using Printf
using SpecialFunctions

function debug_interaction()
    freq = 300e6
    k = 2 * pi * freq / 299792458.0
    eta = 376.73
    
    # Create two triangles far apart
    # Triangle 1 (Source) at origin
    v1 = [0.0 0.1 0.0; 0.0 0.0 0.1; 0.0 0.0 0.0]
    
    # Triangle 2 (Test) at (1.0, 0, 0)
    offset = [1.0, 0.0, 0.0]
    v2 = v1 .+ offset
    
    # Create Mesh
    vertices = hcat(v1, v2)
    elements = [1 4; 2 5; 3 6] 
    
    gq = GaussQuadratureInfo(:Triangle, 3)
    
    # Geometry Info
    area1 = 0.5 * 0.1 * 0.1
    area2 = area1
    l1 = 0.1 
    l2 = 0.1
    
    # Center of Tri 1
    c1 = [0.1/3, 0.1/3, 0.0]
    # Center of Tri 2
    c2 = c1 .+ offset
    
    # Direct Integration
    Z_direct = zero(ComplexF64)
    
    for i_qp in 1:3
        L_test = gq.coordinate[:, i_qp]
        r_test = v2 * L_test
        w_test = gq.weight[i_qp]
        
        rho_test = r_test - v2[:, 3] # v6
        f_test = (l2 / (2*area2)) * rho_test
        div_f_test = (l2 / area2) * 2 
        
        for j_qp in 1:3
            L_src = gq.coordinate[:, j_qp]
            r_src = v1 * L_src
            w_src = gq.weight[j_qp]
            
            rho_src = r_src - v1[:, 3] # v3
            f_src = (l1 / (2*area1)) * rho_src
            div_f_src = (l1 / area1) * 2
            
            R_vec = r_test - r_src
            R = norm(R_vec)
            G = exp(-im * k * R) / (4 * pi * R)
            
            # Vector part: j * omega * mu * (f_test . f_src) * G
            term_vec = im * k * eta * dot(f_test, f_src) * G
            
            # Scalar part: (1 / j * omega * eps) * (div f_test * div f_src) * G
            # 1 / (j w eps) = eta / (j k) = -j eta / k
            term_scalar = (eta / (im * k)) * (div_f_test * div_f_src) * G
            
            Z_direct += (term_vec + term_scalar) * w_test * area2 * w_src * area1
        end
    end
    
    println("Direct Z: ", Z_direct)
    println("Direct |Z|: ", abs(Z_direct))
    
    # MLFMA Calculation
    L = 3 
    p = 2*L + 1
    nodes, weights = getlbSortedData(p)
    nPoles = length(weights)
    
    r̂sθsϕs = [r̂θϕInfo(SVector{3, Float64}(nodes[:, i])) for i in 1:nPoles]
    poles = LbPolesInfo(weights, r̂sθsϕs)
    
    aggS = zeros(ComplexF64, nPoles, 3)
    
    # Aggregation Loop
    for j_qp in 1:3
        L_src = gq.coordinate[:, j_qp]
        r_src = v1 * L_src
        w_src = gq.weight[j_qp]
        
        rho_src = r_src - v1[:, 3]
        
        factor_vec = l1 / 2 * w_src
        factor_scalar = l1 * w_src # div f * A * w = (l/A) * A * w = l * w
        
        r_local = r_src - c1
        
        for iPole in 1:nPoles
            r̂ = poles.r̂sθsϕs[iPole].r̂
            θhat = poles.r̂sθsϕs[iPole].θhat
            ϕhat = poles.r̂sθsϕs[iPole].ϕhat
            
            phase = exp(im * k * dot(r̂, r_local))
            
            vec = rho_src * factor_vec * phase
            scalar = factor_scalar * phase
            
            aggS[iPole, 1] += dot(θhat, vec)
            aggS[iPole, 2] += dot(ϕhat, vec)
            aggS[iPole, 3] += scalar
        end
    end
    
    # 2. Translation (Box 1 -> Box 2)
    R_BA = c2 - c1
    Rab = norm(R_BA)
    R̂ab = R_BA / Rab
    
    trans_factor = zeros(ComplexF64, nPoles)
    # Use the standard factor
    const_factor_trans = -im * k / (16 * pi^2)
    
    x = k * Rab
    h2lxs = [sphericalbesselj(l, x) - im * sphericalbessely(l, x) for l in 0:L]
    
    for iPole in 1:nPoles
        r̂ = poles.r̂sθsϕs[iPole].r̂
        w_pole = poles.Wθϕs[iPole]
        cosϕ = dot(r̂, R̂ab)
        
        Pl = zeros(L+1)
        Pl[1] = 1.0
        Pl[2] = cosϕ
        for l in 1:L-1
            Pl[l+2] = ((2*l+1)*cosϕ*Pl[l+1] - l*Pl[l])/(l+1)
        end
        
        val = zero(ComplexF64)
        j_term = im
        for l in 0:L
            j_term *= -im
            val += j_term * (2*l + 1) * h2lxs[l+1] * Pl[l+1]
        end
        
        trans_factor[iPole] = val * const_factor_trans * w_pole
    end
    
    # Apply Translation
    disaggG = zeros(ComplexF64, nPoles, 3)
    for iPole in 1:nPoles
        disaggG[iPole, 1] = aggS[iPole, 1] * trans_factor[iPole]
        disaggG[iPole, 2] = aggS[iPole, 2] * trans_factor[iPole]
        disaggG[iPole, 3] = aggS[iPole, 3] * trans_factor[iPole]
    end
    
    # 3. Disaggregation (Box 2 -> Tri 2)
    Z_mlfma = zero(ComplexF64)
    const_factor_disagg = im * k * eta
    
    scalar_factor = -1.0 / (k^2)
    
    for i_qp in 1:3
        L_test = gq.coordinate[:, i_qp]
        r_test = v2 * L_test
        w_test = gq.weight[i_qp]
        
        rho_test = r_test - v2[:, 3]
        factor_vec = l2 / 2 * w_test
        factor_scalar = l2 * w_test
        
        r_local = r_test - c2
        
        for iPole in 1:nPoles
            r̂ = poles.r̂sθsϕs[iPole].r̂
            θhat = poles.r̂sθsϕs[iPole].θhat
            ϕhat = poles.r̂sθsϕs[iPole].ϕhat
            
            phase = exp(-im * k * dot(r̂, r_local))
            
            E_theta = disaggG[iPole, 1]
            E_phi = disaggG[iPole, 2]
            Phi_pole = disaggG[iPole, 3]
            
            E_inc = (E_theta * θhat + E_phi * ϕhat) * phase
            Phi_inc = Phi_pole * phase
            
            val_vec = dot(rho_test, E_inc) * factor_vec
            val_scalar = Phi_inc * factor_scalar * scalar_factor
            
            Z_mlfma += (val_vec + val_scalar) * const_factor_disagg
        end
    end
    
    println("MLFMA Z: ", Z_mlfma)
    println("MLFMA |Z|: ", abs(Z_mlfma))
    
    println("Ratio |MLFMA|/|Direct|: ", abs(Z_mlfma)/abs(Z_direct))
end

debug_interaction()