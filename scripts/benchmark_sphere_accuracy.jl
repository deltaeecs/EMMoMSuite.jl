using EMSuite
using LinearAlgebra
using StaticArrays
using Test
using Printf
using Statistics

# Benchmark: Accuracy Comparison (EFIE, MFIE, CFIE)
# -------------------------------------------------

function run_accuracy_benchmark()
    println("==================================================")
    println("   Benchmark: Accuracy Comparison (Sphere)        ")
    println("==================================================")

    # 1. Parameters
    freq = 300e6 # 300 MHz
    lambda = 299792458.0 / freq
    radius = 0.5 # 0.5 meter radius (ka = pi approx 3.14)
    
    println(@sprintf("Frequency: %.2f MHz", freq/1e6))
    println(@sprintf("Wavelength: %.4f m", lambda))
    println(@sprintf("Radius: %.4f m (ka = %.4f)", radius, 2*pi*radius/lambda))

    # 2. Mesh Generation
    println("\n[1/4] Generating Mesh...")
    # Use a moderate mesh
    n_theta = 16
    n_phi = 32
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    println(@sprintf("      Vertices: %d", num_vertices(mesh)))
    println(@sprintf("      Elements: %d", num_elements(mesh)))

    # 3. Basis Setup
    println("\n[2/4] Setting up Basis...")
    basis = RWGBasis(mesh)
    println(@sprintf("      Unknowns: %d", num_basis(basis)))

    # 4. Solvers Loop
    solvers = [:EFIE, :MFIE, :CFIE]
    results = Dict()
    
    # Analytical Mie Series
    theta_range = collect(range(0, pi, length=181))
    phi_val = [0.0]
    rcs_mie = calculate_mie_rcs_pec_sphere(radius, freq, theta_range)
    rcs_db_mie = 10 .* log10.(rcs_mie)
    
    for solver_type in solvers
        println("\n--------------------------------------------------")
        println("   Running Solver: $solver_type")
        println("--------------------------------------------------")
        
        # Construct Operator
        if solver_type == :EFIE
            ie = EFIE(freq)
        elseif solver_type == :MFIE
            ie = MFIE(freq)
        elseif solver_type == :CFIE
            efie = EFIE(freq)
            mfie = MFIE(freq)
            ie = CFIE(freq, 0.5, efie, mfie)
        end
        
        # Assemble
        println("      Assembling Matrix...")
        t_asm = @elapsed Z = assemble_impedance_matrix(ie, basis)
        println(@sprintf("      Assembly: %.4f s", t_asm))
        
        # Excitation
        source = PlaneWave(freq, pi, 0.0, [1.0, 0.0, 0.0])
        V = excitation_vector(ie, source, basis)
        
        # Solve
        println("      Solving System...")
        # Use LU for small problems to be robust
        t_solve = @elapsed I = solve!(LUSolver(), Z, V)
        println(@sprintf("      Solve: %.4f s", t_solve))
        
        # RCS
        println("      Calculating RCS...")
        trianglesInfo = [get_triangle_info(mesh, basis, i) for i in 1:num_elements(mesh)]
        set_frequency!(freq)
        
        ff = farField(theta_range, phi_val, I, trianglesInfo, nothing, RWG)
        
        rcs_num = zeros(Float64, length(theta_range))
        for i in 1:length(theta_range)
            E_theta = ff[1, i, 1]
            E_phi = ff[2, i, 1]
            rcs_num[i] = 4 * pi * (abs2(E_theta) + abs2(E_phi))
        end
        rcs_db_num = 10 .* log10.(rcs_num)
        
        # Error
        rmse = sqrt(mean((rcs_db_num .- rcs_db_mie).^2))
        println(@sprintf("      RMSE (dB): %.4f", rmse))
        
        results[solver_type] = (rcs_db_num, rmse)
    end
    
    # 5. Summary
    println("\n==================================================")
    println("   Summary of Results (RMSE dB)                   ")
    println("==================================================")
    for solver_type in solvers
        rmse = results[solver_type][2]
        println(@sprintf("   %s: %.4f", solver_type, rmse))
    end
    
    # 6. MLFMA Check (Placeholder)
    println("\n[MLFMA Check]")
    println("   MLFMA integration is currently being verified separately.")
    # TODO: Add MLFMA benchmark once high-level API is confirmed.
    
end

run_accuracy_benchmark()
