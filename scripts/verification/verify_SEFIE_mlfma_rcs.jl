using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.Solvers
using EMMoMSuite.FastAlgorithms.MLFMA
using EMMoMSuite.PostProcessing
using EMMoMSuite.Utilities.Parameters
using EMMoMSuite.CoreModule.Sources
using StaticArrays
using LinearAlgebra
using SpecialFunctions
using IterativeSolvers
using Statistics
using Printf
using JLD2


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
        
        # Legendre Polynomials
        # P_n^1(cos theta) / sin theta = pi_n
        # d P_n^1 / d theta = tau_n
        
        pi_prev = 0.0
        pi_curr = 1.0 # pi_1 = 1
        
        for n in 1:n_max
            # Recurrence for pi_n
            # pi_n = P_n^1 / sin theta
            # pi_{n+1} = (2n+1)/(n+1) * mu * pi_n - (n+1)/n * pi_{n-1} -> This is for P_n
            # Correct recurrence for pi_n:
            # pi_n = (2n-1)/(n-1) * mu * pi_{n-1} - n/(n-1) * pi_{n-2}
            
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
        # H-plane (phi=90): E_phi ~ S1
        # RCS = 4 * pi * |S|^2 / k^2
        # Here we compute E-plane RCS
        rcs[i] = 4 * π * abs2(S2) / (ka/1.0)^2 # k = ka/a. If a=1, k=ka.
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
    mesh_file = joinpath(@__DIR__, "sphere_mesh_data.jld2")
    if !isfile(mesh_file)
        println("Generating mesh...")
        # Fallback or error
        error("Mesh file not found: $mesh_file")
    end
    
    # Load mesh (assuming JLD2 format from previous scripts)
    data = load(mesh_file)
    nodes = data["nodes"]
    triangles = data["triangles"]
    
    if size(nodes, 1) != 3; nodes = nodes'; end
    if size(triangles, 1) != 3; triangles = triangles'; end
    
    mesh = TriangleMesh(size(triangles, 2), Float64.(nodes), Int.(triangles))
    println("Mesh: $(size(mesh.node, 2)) vertices, $(size(mesh.triangles, 2)) elements")

    # 3. Basis
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")

    # 4. MLFMA Operator
    println("Setting up MLFMA...")
    set_frequency!(freq)
    efie = EFIE(freq)
    
    # MLFMA Parameters
    l_min = lambda / 4.0
    mlfma_op = MLFMAOperator(efie, basis, l_min)
    
    # 5. Excitation
    println("Setting up Excitation...")
    # Incident from +z (theta=0), x-polarized
    # k_dir = (0,0,-1) -> theta=180
    # Wait, standard RCS usually assumes incident from +z or -z.
    # Let's use incident from +z (propagating to -z).
    source = PlaneWave(freq, π, 0.0, [1.0, 0.0, 0.0]) 
    V = excitation_vector(efie, source, basis)

    # 6. Solve
    println("Solving (GMRES)...")
    # Use a simple diagonal preconditioner if possible, or just identity
    # Extract diagonal from Z_near for preconditioning
    # Z_near is sparse.
    diag_Z = diag(mlfma_op.Z_near)
    P = Diagonal(1.0 ./ diag_Z)
    
    # Solve
    t_start = time()
    x = gmres(mlfma_op, V, Pl=P, restart=50, maxiter=100, reltol=1e-4, verbose=true)
    t_end = time()
    println("Solved in $(t_end - t_start) s")

    # 7. Calculate RCS
    println("Calculating RCS...")
    theta_range = collect(0:0.02:π) # 0 to 180 degrees
    phi_range = [0.0] # E-plane
    
    # Use EMMoMSuite.PostProcessing.radarCrossSection
    # Returns: (RCS_components, RCS_total, RCS_dB)
    # Note: RCS_total is a Matrix [theta, phi]
    _, rcs_mlfma, rcs_mlfma_db_matrix = EMMoMSuite.PostProcessing.radarCrossSection(theta_range, phi_range, x, basis)
    
    # Extract vector for phi=0
    rcs_mlfma_db = vec(rcs_mlfma_db_matrix)

    # 8. Mie Series
    println("Calculating Mie Series...")
    rcs_mie = mie_rcs_bistatic(ka, theta_range)
    rcs_mie_db = 10 .* log10.(rcs_mie)

    # 9. Compare
    rmse = sqrt(mean((rcs_mlfma_db .- rcs_mie_db).^2))
    println("RMSE (dB): $rmse")
    
    # Print some values
    println("\nTheta (deg) | MLFMA (dBsm) | Mie (dBsm) | Diff (dB)")
    println("---------------------------------------------------")
    indices = [1, length(theta_range)÷4, length(theta_range)÷2, 3*length(theta_range)÷4, length(theta_range)]
    for i in indices
        deg = rad2deg(theta_range[i])
        println(Printf.@sprintf("%10.1f | %10.2f | %10.2f | %10.2f", deg, rcs_mlfma_db[i], rcs_mie_db[i], rcs_mlfma_db[i] - rcs_mie_db[i]))
    end
    
    if rmse < 3.0
        println("\nSUCCESS: MLFMA RCS matches Mie Series within tolerance.")
    else
        println("\nFAILURE: MLFMA RCS mismatch.")
    end
end

verify_mlfma_rcs()
