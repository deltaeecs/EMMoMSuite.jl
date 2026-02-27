module MieSeries

using SpecialFunctions
using LinearAlgebra

export calculate_mie_rcs_pec_sphere

"""
    calculate_mie_rcs_pec_sphere(radius, freq, theta_range)

Calculate the bistatic RCS of a PEC sphere using Mie series.
Returns RCS in m^2.

# Arguments
- `radius`: Sphere radius (m)
- `freq`: Frequency (Hz)
- `theta_range`: Vector of observation angles (radians). 0 is forward scattering, pi is backscattering.
"""
function calculate_mie_rcs_pec_sphere(radius, freq, theta_range)
    c0 = 299792458.0
    k = 2π * freq / c0
    x = k * radius
    
    n_max = ceil(Int, x + 4 * x^(1/3) + 2)
    
    # Precompute Bessel functions
    # We need j_n(x) and h_n(2)(x)
    # And their derivatives.
    # Riccati-Bessel functions:
    # psi_n(x) = x * j_n(x)
    # zeta_n(x) = x * h_n(2)(x)
    
    # Coefficients a_n and b_n
    a_n = zeros(ComplexF64, n_max)
    b_n = zeros(ComplexF64, n_max)
    
    for n in 1:n_max
        # Bessel functions
        jn = sphericalbesselj(n, x)
        jnm1 = sphericalbesselj(n-1, x)
        # Derivative of x*jn(x) = jn(x) + x*jn'(x)
        # Recurrence: jn'(x) = j(n-1)(x) - (n+1)/x * jn(x)
        # So (x*jn(x))' = jn(x) + x*(j(n-1)(x) - (n+1)/x * jn(x))
        #               = jn(x) + x*j(n-1)(x) - (n+1)*jn(x)
        #               = x*j(n-1)(x) - n*jn(x)
        
        psi_n = x * jn
        psi_n_prime = x * jnm1 - n * jn
        
        # Hankel functions (Second kind for outgoing e^jwt? No, e^-jkr is outgoing for e^jwt? 
        # Wait, EMSuite uses exp(-im * k * R).
        # Time convention: usually exp(j wt).
        # Green's function: exp(-j k R) / R.
        # This is consistent with exp(j wt).
        # For exp(j wt), outgoing wave is h_n(2).
        
        hn = sphericalbesselj(n, x) - im * sphericalbessely(n, x) # h_n(2) = j_n - i y_n
        hnm1 = sphericalbesselj(n-1, x) - im * sphericalbessely(n-1, x)
        
        zeta_n = x * hn
        zeta_n_prime = x * hnm1 - n * hn
        
        a_n[n] = -psi_n / zeta_n # TM mode?
        b_n[n] = -psi_n_prime / zeta_n_prime # TE mode?
        
        # Note: Definitions vary. 
        # Stratton: a_n -> Magnetic, b_n -> Electric
        # Balanis: a_n -> TM, b_n -> TE
        # For PEC sphere:
        # a_n = - j_n(x) / h_n(2)(x)  (This is wrong, needs Riccati)
        # Correct is:
        # a_n = - psi_n'(x) / zeta_n'(x)
        # b_n = - psi_n(x) / zeta_n(x)
        # Let's check Balanis "Advanced Engineering Electromagnetics"
        # For PEC Sphere:
        # b_n = - j_n(x) / h_n(2)(x) ? No.
        
        # Let's use the standard Mie coefficients for PEC (m -> infinity)
        # a_n (TM) = - psi_n'(x) / zeta_n'(x)
        # b_n (TE) = - psi_n(x) / zeta_n(x)
        
        # Wait, my code above has a_n = -psi/zeta.
        # Let's swap them to match standard notation if needed, or just use consistent S1/S2.
        # S1 = sum (2n+1)/(n(n+1)) (a_n pi_n + b_n tau_n)
        # S2 = sum (2n+1)/(n(n+1)) (a_n tau_n + b_n pi_n)
        
        # Let's stick to:
        # a_n = - psi_n'(x) / zeta_n'(x)
        # b_n = - psi_n(x) / zeta_n(x)
        
        a_n[n] = -psi_n_prime / zeta_n_prime
        b_n[n] = -psi_n / zeta_n
    end
    
    rcs = zeros(Float64, length(theta_range))
    
    for (i, theta) in enumerate(theta_range)
        mu = cos(theta)
        S1 = 0.0 + 0.0im
        S2 = 0.0 + 0.0im
        
        # Legendre polynomials and derivatives
        # pi_n = P_n^1(mu) / sin(theta)
        # tau_n = d/dtheta P_n^1(mu)
        
        # Recurrence for pi_n and tau_n
        # pi_0 = 0, pi_1 = 1
        # pi_n = (2n-1)/(n-1) * mu * pi_{n-1} - n/(n-1) * pi_{n-2}
        # tau_n = n * mu * pi_n - (n+1) * pi_{n-1}
        
        # Re-implement loop with proper recurrence
        p_1 = 1.0 # pi_1
        p_0 = 0.0 # pi_0
        
        for n in 1:n_max
            pi_n = p_1
            tau_n = n * mu * pi_n - (n+1) * p_0
            
            term = (2n+1) / (n*(n+1))
            S1 += term * (a_n[n] * pi_n + b_n[n] * tau_n)
            S2 += term * (a_n[n] * tau_n + b_n[n] * pi_n)
            
            # Update for next n
            # pi_{n+1} = (2n+1)/n * mu * pi_n - (n+1)/n * pi_{n-1}
            p_next = ((2n+1) * mu * p_1 - (n+1) * p_0) / n
            p_0 = p_1
            p_1 = p_next
        end
        
        # RCS = 4 pi |S|^2 / k^2
        # E-plane (phi=0): |S2|^2
        # H-plane (phi=pi/2): |S1|^2
        # For sphere, E-plane and H-plane are symmetric?
        # No, E-plane is parallel to polarization, H-plane is perpendicular.
        # Usually we plot both or average?
        # Let's return both.
        
        # But the function signature returns a single vector.
        # Let's return E-plane RCS (S2) as it's commonly plotted.
        
        rcs[i] = 4 * pi * abs2(S2) / k^2
    end
    
    return rcs
end

end
