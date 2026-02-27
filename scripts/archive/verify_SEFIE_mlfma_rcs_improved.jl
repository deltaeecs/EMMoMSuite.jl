using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Solvers
using EMSuite.FastAlgorithms.MLFMA
using EMSuite.PostProcessing
using EMSuite.Utilities.Parameters
using EMSuite.CoreModule.Sources
using StaticArrays
using LinearAlgebra
using SpecialFunctions
using IterativeSolvers
using Statistics
using Printf
using JLD2
using DelimitedFiles
using IncompleteLU

function verify_mlfma_rcs_improved()
    println("==================================================")
    println("   Verification: SEFIE MLFMA RCS vs Legacy (Improved)")
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

    # 4. MLFMA Operator
    println("Setting up MLFMA (CFIE)...")
    set_frequency!(freq)
    
    efie = EFIE(freq)
    mfie = MFIE(freq)
    cfie = CFIE(freq, 0.6, efie, mfie)
    
    l_min = lambda / 4.0
    mlfma_op = MLFMAOperator(cfie, basis, l_min)
    
    # 5. Excitation
    println("Setting up Excitation...")
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0]) 
    V = excitation_vector(cfie, source, basis)

    # 6. Solve
    println("Solving (GMRES with ILU)...")
    
    # Construct ILU Preconditioner from Near-Field Matrix
    println("Computing ILU factorization of Z_near...")
    # Use a small drop tolerance for better approximation
    P_ilu = ilu(mlfma_op.Z_near, τ=0.01) 
    
    t_start = time()
    # Tighter tolerance: 1e-6
    x = gmres(mlfma_op, V, Pl=P_ilu, restart=100, maxiter=200, reltol=1e-6, verbose=true)
    t_end = time()
    println("Solved in $(t_end - t_start) s")

    # 7. Calculate RCS
    println("Calculating RCS...")
    theta_range = collect(0:0.01:π) 
    phi_range = [0.0] 
    
    _, rcs_mlfma, rcs_mlfma_db_matrix = EMSuite.PostProcessing.radarCrossSection(theta_range, phi_range, x, basis)
    rcs_mlfma_db = vec(rcs_mlfma_db_matrix)

    # 8. Load Legacy Reference
    println("Loading Legacy Reference...")
    legacy_file = joinpath(@__DIR__, "..", "legacy_rcs_reference.txt")
    if !isfile(legacy_file)
        error("Legacy reference file not found: $legacy_file. Run scripts/run_legacy_simulation.jl first.")
    end
    
    legacy_data = readdlm(legacy_file, ',', skipstart=1)
    rcs_legacy_db = legacy_data[:, 2]

    # 9. Compare
    if length(rcs_mlfma_db) != length(rcs_legacy_db)
        println("Warning: Length mismatch. Truncating to minimum.")
        n = min(length(rcs_mlfma_db), length(rcs_legacy_db))
        rcs_mlfma_db = rcs_mlfma_db[1:n]
        rcs_legacy_db = rcs_legacy_db[1:n]
        theta_range = theta_range[1:n]
    end

    rmse = sqrt(mean((rcs_mlfma_db .- rcs_legacy_db).^2))
    println("RMSE (dB) vs Legacy: $rmse")
    
    println("\nTheta (deg) | MLFMA (dBsm) | Legacy (dBsm) | Diff (dB)")
    println("------------------------------------------------------")
    indices = [1, length(theta_range)÷4, length(theta_range)÷2, 3*length(theta_range)÷4, length(theta_range)]
    for i in indices
        deg = rad2deg(theta_range[i])
        println(Printf.@sprintf("%10.1f | %10.2f | %10.2f | %10.2f", deg, rcs_mlfma_db[i], rcs_legacy_db[i], rcs_mlfma_db[i] - rcs_legacy_db[i]))
    end
    
    if rmse < 1.0
        println("\nSUCCESS: MLFMA RCS matches Legacy within tolerance.")
    else
        println("\nFAILURE: MLFMA RCS mismatch with Legacy.")
    end
end

verify_mlfma_rcs_improved()
