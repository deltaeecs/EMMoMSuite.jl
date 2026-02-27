using EMSuite
using EMSuite.Geometry
using EMSuite.PostProcessing.RCS
using LinearAlgebra
using StaticArrays
using Printf

function run_benchmark()
    println("==================================================")
    println("   Benchmark: PEC Plate RCS (EFIE)                ")
    println("==================================================")
    
    freq = 300e6
    set_frequency!(freq)
    lambda = 299792458.0 / freq
    
    # Create mesh (1m x 1m, 30x30 segments)
    # Target edge length ~ lambda/10 = 0.1m
    # 1.0 / 0.1 = 10 segments.
    # Let's use 30x30 = 900 rectangles = 1800 triangles for a decent size.
    N = 30
    L = 1.0
    println("\n[1/5] Generating Mesh...")
    mesh = generate_rectangle_mesh(L, L, N, N)
    println("      Mesh: $(num_vertices(mesh)) vertices, $(num_elements(mesh)) elements")
    
    # Export mesh for legacy comparison
    nas_file = joinpath(@__DIR__, "plate_benchmark.nas")
    println("      Exporting mesh to $nas_file...")
    write_nas_mesh(nas_file, mesh)
    
    # Basis
    println("\n[2/5] Setting up Basis...")
    t_basis = @elapsed begin
        basis = RWGBasis(mesh)
    end
    println(@sprintf("      Basis setup in %.4f s", t_basis))
    n_unknowns = num_basis(basis)
    println("      Unknowns: $n_unknowns")
    
    # EFIE
    println("\n[3/5] Assembling Impedance Matrix (EFIE)...")
    ie = EFIE(freq)
    
    t_assembly = @elapsed begin
        Z = assemble_impedance_matrix(ie, basis)
    end
    println(@sprintf("      Assembly completed in %.4f s", t_assembly))
    
    # Excitation (Plane Wave, Normal Incidence)
    # Theta=0 (from +z?), Phi=0, Pol=x
    # Note: PlaneWave(freq, theta, phi, polarization_vector)
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(ie, source, basis)
    println("      Norm of V: $(norm(V))")
    println("      V[1]: ", V[1])
    println("      Z[1,1]: ", Z[1,1])
    println("      Sum Abs Z: $(sum(abs.(Z)))")
    
    # Solve
    println("\n[4/5] Solving Linear System (GMRES)...")
    solver = GMRESSolver(tol=1e-4, maxiter=1000)
    
    t_solve = @elapsed begin
        I = solve!(solver, Z, V)
    end
    println(@sprintf("      Solved in %.4f s", t_solve))
    println("      Norm of I: $(norm(I))")
    println("      I[1]: ", I[1])
    
    # Check Residual
    res = Z * I - V
    rel_res = norm(res) / norm(V)
    println("      Relative Residual norm(Z*I - V)/norm(V): $rel_res")
    println("      (Z*I)[1]: ", (Z*I)[1])
    
    # RCS Calculation
    println("\n[5/5] Calculating RCS...")
    
    # Observation angles
    # Legacy used 0 to 180 degrees
    theta_obs = collect(LinRange(0, pi, 181))
    phi_obs = [0.0]
    
    t_rcs = @elapsed begin
        # Returns: RCSθsϕs, RCS_total, RCS_dB
        _, _, rcs_db = radarCrossSection(theta_obs, phi_obs, I, basis)
    end
    println(@sprintf("      RCS calculated in %.4f s", t_rcs))
    
    # Save RCS
    rcs_file = joinpath(@__DIR__, "rcs_emsuite.txt")
    open(rcs_file, "w") do io
        for i in 1:length(theta_obs)
            deg = theta_obs[i] * 180 / pi
            val = rcs_db[i, 1]
            write(io, "$deg\t$val\n")
        end
    end
    println("      RCS saved to $rcs_file")
    
    println("\nBenchmark Completed.")
end

run_benchmark()
