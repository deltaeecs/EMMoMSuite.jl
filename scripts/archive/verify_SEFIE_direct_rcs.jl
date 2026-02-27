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

function verify_direct_rcs()
    println("==================================================")
    println("   Verification: SEFIE Direct RCS vs Legacy")
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
    println("Setting up EFIE (Direct)...")
    set_frequency!(freq)
    
    efie = EFIE(freq)
    # mfie = MFIE(freq)
    # cfie = CFIE(freq, 0.6, efie, mfie)
    
    # 5. Assembly
    println("Assembling Impedance Matrix...")
    t_asm = @elapsed begin
        Z = assemble_impedance_matrix(efie, basis)
    end
    println("Assembly time: $t_asm s")

    # 6. Excitation
    println("Setting up Excitation...")
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0]) 
    V = excitation_vector(efie, source, basis)

    # 7. Solve
    println("Solving (LU Decomposition)...")
    t_solve = @elapsed begin
        x = Z \ V
    end
    println("Solved in $t_solve s")

    # Save Matrix and Vector
    jldsave(joinpath(@__DIR__, "..", "emsuite_matrix.jld2"); Z=Z, V=V, I=x)

    # 8. Calculate RCS
    println("Calculating RCS...")
    theta_range = collect(0:0.01:π) 
    phi_range = [0.0] 
    
    _, rcs_direct, rcs_direct_db_matrix = EMSuite.PostProcessing.radarCrossSection(theta_range, phi_range, x, basis)
    rcs_direct_db = vec(rcs_direct_db_matrix)

    # 9. Load Legacy Reference
    println("Loading Legacy Reference...")
    legacy_file = joinpath(@__DIR__, "..", "legacy_rcs_reference.txt")
    if !isfile(legacy_file)
        error("Legacy reference file not found: $legacy_file. Run scripts/run_legacy_simulation.jl first.")
    end
    
    legacy_data = readdlm(legacy_file, ',', skipstart=1)
    rcs_legacy_db = legacy_data[:, 2]

    # 10. Compare
    if length(rcs_direct_db) != length(rcs_legacy_db)
        println("Warning: Length mismatch. Truncating to minimum.")
        n = min(length(rcs_direct_db), length(rcs_legacy_db))
        rcs_direct_db = rcs_direct_db[1:n]
        rcs_legacy_db = rcs_legacy_db[1:n]
        theta_range = theta_range[1:n]
    end

    rmse = sqrt(mean((rcs_direct_db .- rcs_legacy_db).^2))
    println("RMSE (dB) vs Legacy: $rmse")
    
    println("\nTheta (deg) | Direct (dBsm) | Legacy (dBsm) | Diff (dB)")
    println("------------------------------------------------------")
    indices = [1, length(theta_range)÷4, length(theta_range)÷2, 3*length(theta_range)÷4, length(theta_range)]
    for i in indices
        deg = rad2deg(theta_range[i])
        println(Printf.@sprintf("%10.1f | %10.2f | %10.2f | %10.2f", deg, rcs_direct_db[i], rcs_legacy_db[i], rcs_direct_db[i] - rcs_legacy_db[i]))
    end
    
    if rmse < 0.5
        println("\nSUCCESS: Direct RCS matches Legacy within tolerance.")
    else
        println("\nFAILURE: Direct RCS mismatch with Legacy.")
    end
end

verify_direct_rcs()
