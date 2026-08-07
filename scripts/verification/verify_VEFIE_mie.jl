using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.Solvers
using EMMoMSuite.PostProcessing
using EMMoMSuite.Utilities.Parameters
using EMMoMSuite.CoreModule.Sources
using StaticArrays
using LinearAlgebra
using SpecialFunctions
using Statistics
using Printf
using DelimitedFiles

function verify_vefie_mie()
    println("==================================================")
    println("   Verification: VEFIE Direct RCS vs Mie Series")
    println("==================================================")

    # 1. Parameters
    freq = 300e6
    c0 = 299792458.0
    lambda = c0 / freq
    k = 2 * π / lambda
    
    # Set global parameters
    EMMoMSuite.Utilities.Parameters.set_frequency!(freq)
    
    # 2. Mesh
    # Use the same Tetra.nas mesh as in verify_VEFIE_direct.jl
    mesh_file = joinpath(@__DIR__, "../../../deps/fixtures/AllinOne/meshfiles/Tetra.nas")
    if !isfile(mesh_file)
        println("Error: Mesh file not found: $mesh_file")
        return
    end
    
    println("Loading mesh from $mesh_file...")
    mesh = read_nas_mesh(mesh_file, scale=1e-3)
    
    # Center mesh
    min_c = minimum(mesh.node, dims=2)
    max_c = maximum(mesh.node, dims=2)
    center = (min_c + max_c) / 2
    mesh.node .-= center
    
    # Calculate equivalent radius for Mie comparison
    # Volume of sphere = 4/3 * pi * r^3
    # Sum of tetra volumes
    total_vol = 0.0
    nodes = mesh.node
    tetras = mesh.tetras
    for i in 1:size(tetras, 2)
        v_indices = tetras[:, i]
        r1 = nodes[:, v_indices[1]]
        r2 = nodes[:, v_indices[2]]
        r3 = nodes[:, v_indices[3]]
        r4 = nodes[:, v_indices[4]]
        mat = hcat(r2-r1, r3-r1, r4-r1)
        total_vol += abs(det(mat)) / 6.0
    end
    
    radius_equiv = (total_vol * 3 / (4 * π))^(1/3)
    println("Total Mesh Volume: $total_vol m^3")
    println("Equivalent Radius: $radius_equiv m")
    
    ka = k * radius_equiv
    println("Frequency: $(freq/1e6) MHz")
    println("ka: $ka")

    # 3. Basis
    println("Setting up SWG Basis...")
    basis = SWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")

    # 4. Permittivity
    eps_r = 2.0 + 0.0im
    permittivities = fill(eps_r, num_elements(mesh))

    # 5. Operator
    println("Setting up VEFIE...")
    vefie = VEFIE(freq, permittivities)
    
    # 6. Assembly
    println("Assembling Impedance Matrix...")
    Z = assemble_impedance_matrix(vefie, basis, permittivities)

    # 7. Excitation
    println("Setting up Excitation...")
    source = PlaneWave(freq, π, 0.0, [1.0, 0.0, 0.0]) # Incident from +z
    V = excitation_vector(vefie, source, basis, permittivities)

    # 8. Solve
    println("Solving...")
    println("Condition number of Z: ", cond(Z))
    println("Norm of V: ", norm(V))
    D = Z \ V
    println("Norm of Solution D: ", norm(D))
    
    # 9. Calculate RCS
    println("Calculating RCS...")
    theta_range = collect(0:1.0:180.0) .* (π/180.0)
    phi_range = [0.0] # E-plane
    
    _, _, rcs_db_matrix = radarCrossSection(theta_range, phi_range, D, basis, permittivities)
    rcs_db = vec(rcs_db_matrix)
    
    # 10. Mie Series (Dielectric Sphere)
    println("Calculating Mie Series Reference...")
    mie_rcs_db = calculate_mie_series_dielectric(radius_equiv, freq, 2.0, theta_range)
    
    # 11. Comparison
    println("\nTheta (deg) | VEFIE (dBsm) | Mie (dBsm) | Diff (dB)")
    println("---------------------------------------------------")
    
    rms_error = 0.0
    for i in 1:length(theta_range)
        theta_deg = rad2deg(theta_range[i])
        diff = abs(rcs_db[i] - mie_rcs_db[i])
        rms_error += diff^2
        
        if i % 10 == 1 || i == length(theta_range)
            @printf("%11.1f | %12.2f | %10.2f | %9.2f\n", theta_deg, rcs_db[i], mie_rcs_db[i], diff)
        end
    end
    rms_error = sqrt(rms_error / length(theta_range))
    println("---------------------------------------------------")
    println("RMS Error: $rms_error dB")
    
    # Save to file
    open("vefie_vs_mie.txt", "w") do io
        println(io, "Theta(deg) VEFIE(dBsm) Mie(dBsm)")
        for i in 1:length(theta_range)
            println(io, "$(rad2deg(theta_range[i])) $(rcs_db[i]) $(mie_rcs_db[i])")
        end
    end
    println("Results saved to vefie_vs_mie.txt")
end

function calculate_mie_series_dielectric(a, freq, eps_r, theta_range)
    # Simplified Mie Series for Dielectric Sphere
    # Based on standard formulas (e.g., Balanis, Harrington)
    
    k0 = 2π * freq / 299792458.0
    x = k0 * a
    m = sqrt(eps_r) # Refractive index
    
    n_max = ceil(Int, x + 4 * x^(1/3) + 2)
    
    rcs_db = zeros(length(theta_range))
    
    for (idx, theta) in enumerate(theta_range)
        # Scattering amplitude S1, S2
        S1 = 0.0 + 0.0im
        S2 = 0.0 + 0.0im
        
        for n in 1:n_max
            # Mie coefficients an, bn
            # Using Riccati-Bessel functions
            
            # Arguments
            z = x
            z1 = m * x
            
            # Spherical Bessel j_n(z) and y_n(z)
            # psi_n(z) = z * j_n(z)
            # xi_n(z) = z * h_n(1)(z) = z * (j_n(z) + i*y_n(z))
            
            jn_z = sphericalbesselj(n, z)
            yn_z = sphericalbessely(n, z)
            hn_z = jn_z + im * yn_z
            
            jn_z1 = sphericalbesselj(n, z1)
            
            # Derivatives
            # d/dz (z j_n(z)) = j_n(z) + z j'_n(z)
            # Recurrence: j'_n(z) = j_{n-1}(z) - (n+1)/z j_n(z)
            # So d/dz (z j_n(z)) = j_n(z) + z (j_{n-1}(z) - (n+1)/z j_n(z))
            #                    = j_n(z) + z j_{n-1}(z) - (n+1) j_n(z)
            #                    = z j_{n-1}(z) - n j_n(z)
            
            psi_n_z = z * jn_z
            psi_prime_n_z = z * sphericalbesselj(n-1, z) - n * jn_z
            
            xi_n_z = z * hn_z
            xi_prime_n_z = z * (sphericalbesselj(n-1, z) + im * sphericalbessely(n-1, z)) - n * hn_z
            
            psi_n_z1 = z1 * jn_z1
            psi_prime_n_z1 = z1 * sphericalbesselj(n-1, z1) - n * jn_z1
            
            # Coefficients
            # an = (m psi_n(mx) psi'_n(x) - psi_n(x) psi'_n(mx)) / (m psi_n(mx) xi'_n(x) - xi_n(x) psi'_n(mx))
            # bn = (psi_n(mx) psi'_n(x) - m psi_n(x) psi'_n(mx)) / (psi_n(mx) xi'_n(x) - m xi_n(x) psi'_n(mx))
            
            num_an = m * psi_n_z1 * psi_prime_n_z - psi_n_z * psi_prime_n_z1
            den_an = m * psi_n_z1 * xi_prime_n_z - xi_n_z * psi_prime_n_z1
            an = num_an / den_an
            
            num_bn = psi_n_z1 * psi_prime_n_z - m * psi_n_z * psi_prime_n_z1
            den_bn = psi_n_z1 * xi_prime_n_z - m * xi_n_z * psi_prime_n_z1
            bn = num_bn / den_bn
            
            # Angular functions pi_n and tau_n
            # pi_n = P_n^1(cos theta) / sin theta
            # tau_n = d/dtheta P_n^1(cos theta)
            
            # Using approximation or library?
            # Let's implement simple recurrence for Legendre Polynomials
            # P_n^1(x)
            
            ct = cos(theta)
            st = sin(theta)
            
            if abs(theta) < 1e-6
                pi_n = 0.5 * n * (n + 1)
                tau_n = 0.5 * n * (n + 1)
            elseif abs(theta - π) < 1e-6
                pi_n = 0.5 * (-1)^n * n * (n + 1)
                tau_n = -0.5 * (-1)^n * n * (n + 1)
            else
                # Recurrence for P_n^1
                # P_1^1 = -sin(theta)
                # P_2^1 = -3 sin(theta) cos(theta)
                # (n-m) P_n^m = (2n-1) x P_{n-1}^m - (n+m-1) P_{n-2}^m
                
                # We need pi_n = P_n^1 / sin(theta)
                # P_1^1 / sin = -1
                # P_2^1 / sin = -3 cos
                
                # Recurrence for pi_n:
                # pi_n(x) = (2n-1)/(n-1) * x * pi_{n-1} - n/(n-1) * pi_{n-2}
                
                if n == 1
                    pi_n = 1.0
                    tau_n = ct
                else
                    # Recompute for specific n (inefficient but safe)
                    p_prev2 = 0.0 # pi_0 doesn't exist really
                    p_prev1 = 1.0 # pi_1
                    
                    pi_val = 1.0
                    
                    for k_leg in 2:n
                        pi_val = ( (2*k_leg - 1) * ct * p_prev1 - k_leg * p_prev2 ) / (k_leg - 1)
                        p_prev2 = p_prev1
                        p_prev1 = pi_val
                    end
                    pi_n = pi_val
                    
                    # tau_n = n * ct * pi_n - (n+1) * pi_{n-1}
                    # Wait, relation:
                    # dP_n^1/dtheta = ...
                    # tau_n = n * cos(theta) * pi_n - (n+1) * pi_{n-1} ? No
                    # tau_n = n * ct * pi_n - (n+1) * p_prev2 ?
                    
                    # Standard relation:
                    # tau_n = n * ct * pi_n - (n+1) * pi_{n-1}
                    # where pi_{n-1} is the previous value
                    
                    tau_n = n * ct * pi_n - (n+1) * p_prev2
                end
            end
            
            factor = (2n + 1) / (n * (n + 1))
            S1 += factor * (an * pi_n + bn * tau_n)
            S2 += factor * (an * tau_n + bn * pi_n)
        end
        
        # RCS = 4 pi |S|^2 / k^2
        # But standard definition usually absorbs k?
        # sigma = 4 pi |S|^2 / k^2?
        # Wait, S here is dimensionless? No.
        # Standard Mie: E_s = E_0 exp(-ikr)/(-ikr) * S
        # sigma = 4 pi |S|^2 / k^2
        
        sigma = 4 * π * (abs(S1)^2 + abs(S2)^2)/2 / k0^2 # Average polarization?
        # For E-plane (phi=0), we use S2?
        # E_theta ~ S2, E_phi ~ S1
        # Incident x-pol.
        # E-plane is phi=0. E_theta is dominant.
        # H-plane is phi=90. E_phi is dominant.
        
        # Let's assume E-plane (phi=0) -> Theta component -> S2
        sigma_e = 4 * π * abs(S2)^2 / k0^2
        
        rcs_db[idx] = 10 * log10(sigma_e)
    end
    
    return rcs_db
end

verify_vefie_mie()
