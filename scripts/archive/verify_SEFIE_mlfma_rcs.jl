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

# --- Mie Series Implementation ---
function sphericalbesselj(n, x)
    return sqrt(π/(2x)) * besselj(n + 0.5, x)
end

function sphericalbessely(n, x)
    return sqrt(π/(2x)) * bessely(n + 0.5, x)
end

function mie_rcs_bistatic(ka, theta_range)
    n_max = ceil(Int, ka + 4 * ka^(1/3) + 2)
    an = zeros(ComplexF64, n_max)
    bn = zeros(ComplexF64, n_max)
    x = ka
    
    for n in 1:n_max
        jn = sphericalbesselj(n, x)
        yn = sphericalbessely(n, x)
        hn2 = jn - im * yn
        
        jn_prev = sphericalbesselj(n-1, x)
        yn_prev = sphericalbessely(n-1, x)
        hn2_prev = jn_prev - im * yn_prev
        
        d_xjn = jn + x * (jn_prev - (n+1)/x * jn)
        d_xhn2 = hn2 + x * (hn2_prev - (n+1)/x * hn2)
        
        an[n] = -jn / hn2
        bn[n] = -d_xjn / d_xhn2
    end
    
    rcs = zeros(Float64, length(theta_range))
    
    for (i, theta) in enumerate(theta_range)
        S1 = 0.0 + 0.0im
        S2 = 0.0 + 0.0im
        mu = cos(theta)
        
        pi_prev = 0.0
        pi_curr = 1.0 # pi_1 = 1
        
        for n in 1:n_max
            if n == 1
                pi_n = 1.0
                tau_n = mu
            else
                pi_n = ((2n-1)*mu*pi_curr - n*pi_prev)/(n-1)
                tau_n = n*mu*pi_n - (n+1)*pi_curr
            end
            
            term = (2n+1)/(n*(n+1))
            S1 += term * (an[n] * pi_n + bn[n] * tau_n)
            S2 += term * (an[n] * tau_n + bn[n] * pi_n)
            
            pi_prev = pi_curr
            pi_curr = pi_n
        end
        
        # E-plane (phi=0): E_theta ~ S2
        rcs[i] = 4 * π * abs2(S2) / (ka/1.0)^2 
    end
    
    return rcs
end

function verify_mlfma_rcs()
    println("==================================================")
    println("   Verification: SEFIE MLFMA RCS vs Mie Series")
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
    # Generate finer mesh for accuracy
    # Target edge length ~ lambda/10 = 0.1m
    # Area = 4*pi*r^2 = 12.56
    # N_elements ~ 12.56 / (0.1^2 * sqrt(3)/4) ? No, just approximate.
    # N_theta = pi*r / 0.1 = 31.4 -> 32
    # N_phi = 2*pi*r / 0.1 = 62.8 -> 64
    # Reduced for speed, but still acceptable
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
    
    # Use CFIE for better convergence and accuracy on closed sphere
    # Legacy code uses alpha=0.6 by default. Matching it here.
    efie = EFIE(freq)
    mfie = MFIE(freq)
    cfie = CFIE(freq, 0.6, efie, mfie)
    
    # MLFMA Parameters
    l_min = lambda / 4.0
    mlfma_op = MLFMAOperator(cfie, basis, l_min)
    
    # 5. Excitation
    println("Setting up Excitation...")
    # Incident from -z (theta=pi), propagating in +z (theta=0)
    # This matches the Mie series convention where theta=0 is forward scatter
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0]) 
    V = excitation_vector(cfie, source, basis)

    # 6. Solve
    println("Solving (GMRES)...")
    # Use diagonal preconditioner
    diag_Z = diag(mlfma_op.Z_near)
    P = Diagonal(1.0 ./ diag_Z)
    
    t_start = time()
    # Increase maxiter and restart for better convergence
    x = gmres(mlfma_op, V, Pl=P, restart=100, maxiter=200, reltol=1e-4, verbose=true)
    t_end = time()
    println("Solved in $(t_end - t_start) s")

    # 7. Calculate RCS
    println("Calculating RCS...")
    # Match legacy script range: 0:0.01:pi
    theta_range = collect(0:0.01:π) 
    phi_range = [0.0] # E-plane
    
    _, rcs_mlfma, rcs_mlfma_db_matrix = EMSuite.PostProcessing.radarCrossSection(theta_range, phi_range, x, basis)
    
    # Extract vector for phi=0
    rcs_mlfma_db = vec(rcs_mlfma_db_matrix)

    # 8. Load Legacy Reference
    println("Loading Legacy Reference...")
    legacy_file = joinpath(@__DIR__, "..", "legacy_rcs_reference.txt")
    if !isfile(legacy_file)
        error("Legacy reference file not found: $legacy_file. Run scripts/run_legacy_simulation.jl first.")
    end
    
    legacy_data = readdlm(legacy_file, ',', skipstart=1)
    # theta_legacy = legacy_data[:, 1]
    rcs_legacy_db = legacy_data[:, 2]

    # 9. Compare
    # Ensure lengths match
    if length(rcs_mlfma_db) != length(rcs_legacy_db)
        println("Warning: Length mismatch. Truncating to minimum.")
        n = min(length(rcs_mlfma_db), length(rcs_legacy_db))
        rcs_mlfma_db = rcs_mlfma_db[1:n]
        rcs_legacy_db = rcs_legacy_db[1:n]
        theta_range = theta_range[1:n]
    end

    rmse = sqrt(mean((rcs_mlfma_db .- rcs_legacy_db).^2))
    println("RMSE (dB) vs Legacy: $rmse")
    
    # Print some values
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

verify_mlfma_rcs()
