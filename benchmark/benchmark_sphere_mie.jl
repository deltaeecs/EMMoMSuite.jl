using EMSuite
using LinearAlgebra
using StaticArrays
using Test
using Printf
using Statistics

# Benchmark: PEC Sphere Scattering (Mie Series)
# ---------------------------------------------
# This script benchmarks the accuracy and performance of the EFIE solver
# against the analytical Mie series solution for a PEC sphere.

function run_benchmark()
    println("==================================================")
    println("   Benchmark: PEC Sphere Scattering (Mie Series)  ")
    println("==================================================")

    # 1. Parameters
    freq = 300e6 # 300 MHz
    lambda = 299792458.0 / freq
    radius = 1.0 # 1 meter radius (ka = 2*pi*1/1 = 2*pi approx 6.28)
    
    println(@sprintf("Frequency: %.2f MHz", freq/1e6))
    println(@sprintf("Wavelength: %.4f m", lambda))
    println(@sprintf("Radius: %.4f m (ka = %.4f)", radius, 2*pi*radius/lambda))

    # 2. Mesh Generation
    println("\n[1/5] Generating Mesh...")
    # n_theta = 20, n_phi = 40 -> approx 800 faces?
    # Surface area = 4*pi*r^2 = 12.56
    # Target edge length = lambda / 10 = 0.1
    # Circumference = 2*pi*r = 6.28 -> 63 segments
    # Let's use a coarser mesh for speed in this test, or finer for accuracy.
    # Let's try n_theta=32, n_phi=64
    n_theta = 32
    n_phi = 64
    
    t_mesh = @elapsed begin
        mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    end
    println(@sprintf("      Mesh generated in %.4f s", t_mesh))
    println(@sprintf("      Vertices: %d", num_vertices(mesh)))
    println(@sprintf("      Elements: %d", num_elements(mesh)))

    # 3. Basis Setup
    println("\n[2/5] Setting up Basis...")
    t_basis = @elapsed begin
        basis = RWGBasis(mesh)
    end
    println(@sprintf("      Basis setup in %.4f s", t_basis))
    println(@sprintf("      Unknowns: %d", num_basis(basis)))

    # 4. Matrix Assembly
    println("\n[3/5] Assembling Impedance Matrix (CFIE)...")
    # Use CFIE to avoid internal resonances and singularity issues of EFIE
    efie = EFIE(freq)
    mfie = MFIE(freq)
    cfie = CFIE(freq, 0.5, efie, mfie)
    
    t_assembly = @elapsed begin
        Z = assemble_impedance_matrix(cfie, basis)
    end
    println(@sprintf("      Assembly completed in %.4f s", t_assembly))
    
    # 5. Solve
    println("\n[4/5] Solving Linear System...")
    
    # Source: From Bottom (-z), propagating +z (theta=0)
    # Note: PlaneWave(theta) defines propagation direction.
    # theta=0 -> k_dir = +z (Wave moves Up)
    # theta=pi -> k_dir = -z (Wave moves Down)
    # We use theta=0 to match "from bottom" description.
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    
    V = excitation_vector(cfie, source, basis)
    
    solver = GMRESSolver(tol=1e-4, maxiter=1000)
    
    t_solve = @elapsed begin
        I = solve!(solver, Z, V)
    end
    println(@sprintf("      Solved in %.4f s", t_solve))
    
    # 6. Post-Processing (RCS)
    println("\n[5/5] Calculating RCS and Comparing with Mie Series...")
    
    # Bistatic RCS in E-plane (phi=0)
    # Theta from 0 to pi (0 to 180 degrees)
    theta_range = collect(range(0, pi, length=181))
    phi_val = [0.0]
    
    set_frequency!(freq)
    
    t_rcs = @elapsed begin
        ff = farField(theta_range, phi_val, I, basis, nothing)
    end
    
    # Calculate RCS (dBsm)
    # farField returns 2 x N_theta x N_phi array
    # We want total RCS = 4*pi*(|E_theta|^2 + |E_phi|^2)
    
    rcs_num = zeros(Float64, length(theta_range))
    for i in 1:length(theta_range)
        E_theta = ff[1, i, 1]
        E_phi = ff[2, i, 1]
        rcs_num[i] = 4 * pi * (abs2(E_theta) + abs2(E_phi))
    end
    
    rcs_db_num = 10 .* log10.(rcs_num)
    
    # Analytical Mie Series
    # Note: Mie series usually assumes incident from -z?
    # If our source is from pi, then forward scatter is at pi?
    # Let's calculate Mie series.
    # calculate_mie_rcs_pec_sphere takes theta relative to forward direction?
    # Usually theta=0 is forward, theta=pi is back.
    # If incident is from pi (propagating +z), then forward is 0.
    # If incident is from 0 (propagating -z), then forward is pi.
    
    # Let's assume standard Mie: Incident +z (theta=0).
    # If we used source theta=pi, we might need to flip theta for comparison.
    # Let's use source theta=0 (propagating -z? or +z?)
    # PlaneWave(freq, theta, phi, pol) usually defines k vector.
    # If theta=0, k is +z.
    
    # Let's re-run solve with theta=0 if needed, or just interpret results.
    # If k is +z, then theta=0 is forward scatter.
    
    # Let's calculate Mie for the same angles.
    rcs_mie = calculate_mie_rcs_pec_sphere(radius, freq, theta_range)
    rcs_db_mie = 10 .* log10.(rcs_mie)
    
    # Error Calculation
    rmse = sqrt(mean((rcs_db_num .- rcs_db_mie).^2))
    println(@sprintf("      RCS Calculation in %.4f s", t_rcs))
    println(@sprintf("      RMSE (dB): %.4f", rmse))
    
    # Print some values
    println("\n      Theta (deg) | Numerical (dBsm) | Analytical (dBsm) | Diff (dB)")
    println("      -------------------------------------------------------------")
    indices = [1, 46, 91, 136, 181] # 0, 45, 90, 135, 180
    for i in indices
        deg = rad2deg(theta_range[i])
        diff = abs(rcs_db_num[i] - rcs_db_mie[i])
        println(@sprintf("      %11.1f | %16.4f | %17.4f | %9.4f", deg, rcs_db_num[i], rcs_db_mie[i], diff))
    end
    
    println("\n==================================================")
    println("   Benchmark Completed                            ")
    println("==================================================")
end

run_benchmark()
