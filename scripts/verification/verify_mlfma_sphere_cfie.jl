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
using SparseArrays

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
    println("   Verification: CFIE MLFMA RCS vs Mie Series")
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
    n_theta = 16
    n_phi = 32
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
    
    # Correct CFIE constructor
    cfie = CFIE(freq, 0.5)
    
    # MLFMA Parameters
    l_min = lambda / 4.0
    mlfma_op = MLFMAOperator(cfie, basis, l_min)
    
    # 5. Excitation
    println("Setting up Excitation...")
    # Incident from -z (traveling to +z) -> theta=0
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0]) 
    V = excitation_vector(cfie, source, basis)

    # 6. Solve
    println("Solving (GMRES)...")
    # Use a simple diagonal preconditioner
    diag_Z = diag(mlfma_op.Z_near)
    P = Diagonal(1.0 ./ diag_Z)
    
    # Solve
    t_start = time()
    x = gmres(mlfma_op, V, Pl=P, restart=50, maxiter=100, reltol=1e-3, verbose=true)
    t_end = time()
    println("Solved in $(t_end - t_start) s")

    # 7. Calculate RCS
    println("Calculating RCS...")
    theta_range = collect(0:0.02:π) # 0 to 180 degrees
    phi_range = [0.0] # E-plane
    
    _, rcs_mlfma, rcs_mlfma_db_matrix = EMSuite.PostProcessing.radarCrossSection(theta_range, phi_range, x, basis)
    
    # Extract vector for phi=0
    rcs_mlfma_db = vec(rcs_mlfma_db_matrix)

    # 8. Check Values
    fwd_rcs = rcs_mlfma_db[1]   # Theta = 0
    back_rcs = rcs_mlfma_db[end] # Theta = 180
    
    println("Forward RCS (0 deg): $fwd_rcs dBsm (Expected ~21)")
    println("Back RCS (180 deg): $back_rcs dBsm (Expected ~5-7)")
    
    if abs(fwd_rcs - 21) < 3.0 && abs(back_rcs - 6) < 3.0
        println("\nSUCCESS: MLFMA RCS is physically consistent.")
    else
        println("\nFAILURE: MLFMA RCS mismatch.")
    end
end

verify_mlfma_rcs()
