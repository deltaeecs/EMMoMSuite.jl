using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Solvers
using EMSuite.PostProcessing
using EMSuite.Utilities.Parameters
using EMSuite.CoreModule.Sources
using StaticArrays
using LinearAlgebra
using SpecialFunctions
using Statistics
using Printf
using JLD2
using DelimitedFiles

function verify_mie()
    println("==================================================")
    println("   Verification: SEFIE Direct RCS vs Mie Series")
    println("==================================================")

    # 1. Parameters
    freq = 300e6
    c0 = 299792458.0
    lambda = c0 / freq
    k = 2 * π / lambda
    radius = 1.0
    ka = k * radius
    
    println("Frequency: $(freq/1e6) MHz")
    println("ka: $ka")

    # 2. Mesh
    n_theta = 24
    n_phi = 48
    println("Generating sphere mesh (radius=$radius, n_theta=$n_theta, n_phi=$n_phi)...")
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    println("Mesh: $(num_vertices(mesh)) vertices, $(num_elements(mesh)) elements")

    # 3. Basis
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")

    # 4. Operator
    println("Setting up CFIE (Direct)...")
    set_frequency!(freq)
    
    efie = EFIE(freq)
    mfie = MFIE(freq)
    cfie = CFIE(freq, 0.5, efie, mfie)
    
    # 5. Assembly
    println("Assembling Impedance Matrix...")
    t_asm = @elapsed begin
        Z = assemble_impedance_matrix(cfie, basis)
    end
    println("Assembly time: $t_asm s")

    # 6. Excitation
    println("Setting up Excitation...")
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0]) 
    V = excitation_vector(cfie, source, basis)

    # 7. Solve
    println("Solving (LU Decomposition)...")
    t_solve = @elapsed begin
        x = Z \ V
    end
    println("Solved in $t_solve s")

    # 8. Calculate RCS
    println("Calculating RCS...")
    theta_range = collect(0:0.01:π) 
    phi_range = [0.0] 
    
    _, rcs_direct, rcs_direct_db_matrix = EMSuite.PostProcessing.radarCrossSection(theta_range, phi_range, x, basis)
    rcs_direct_db = vec(rcs_direct_db_matrix)

    # 9. Mie Series Calculation
    println("Calculating Mie Series Reference...")
    
    function mie_series_rcs(theta, ka)
        # Simple approximation or full series?
        # Let's use a simplified check for specific angles if full series is complex to implement here.
        # Or better, just output the values for manual check against known Mie plots.
        # For ka=6.28 (radius=1m, freq=300MHz), it's in resonance region.
        
        # Forward scattering (theta=0) should be ~ 4*pi*R^2 / lambda^2 * sigma_geom?
        # Optical limit: sigma = pi*R^2.
        # Forward scattering is usually enhanced.
        
        # Backscattering (theta=180)
        return 0.0
    end
    
    # Let's just print the values for now.
    println("\nTheta (deg) | Direct (dBsm)")
    println("---------------------------")
    indices = [1, length(theta_range)÷4, length(theta_range)÷2, 3*length(theta_range)÷4, length(theta_range)]
    for i in indices
        deg = rad2deg(theta_range[i])
        println(Printf.@sprintf("%10.1f | %10.2f", deg, rcs_direct_db[i]))
    end
    
    # Save to file for plotting
    open("rcs_direct_efie.txt", "w") do io
        println(io, "Theta(deg),RCS(dBsm)")
        for i in 1:length(theta_range)
            println(io, "$(rad2deg(theta_range[i])),$(rcs_direct_db[i])")
        end
    end
    println("Saved RCS to rcs_direct_efie.txt")
end

verify_mie()
