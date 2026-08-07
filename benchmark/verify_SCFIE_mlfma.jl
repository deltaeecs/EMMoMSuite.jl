using EMMoMSuite
using LinearAlgebra
using DataFrames, CSV
using Printf
using Statistics
using EMMoMSuite.Solvers

# Paths
const MOM_ALLINONE_DIR = joinpath(@__DIR__, "..", "deps", "fixtures", "AllinOne")
const OUTPUT_DIR = joinpath(@__DIR__, "../test_results/emsuite_verification")

if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
end

function verify_SCFIE_mlfma()
    println("=== Verifying SCFIE MLFMA (Jet) ===")
    
    # 1. Load Mesh
    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles/jet_100MHz.nas")
    println("Loading mesh: $mesh_file")
    
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    
    # 2. Parameters
    freq = 1e8
    set_frequency!(freq)
    lambda = 299792458.0 / freq
    k = 2π / lambda
    
    # 3. Basis
    println("Setting up RWG basis...")
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Number of unknowns: $N")
    
    # 4. Equation
    println("Setting up CFIE (alpha=0.5)...")
    cfie = CFIE(freq, 0.5)
    
    # 5. MLFMA Operator
    println("Setting up MLFMA Operator...")
    # Leaf size 0.25 lambda
    leaf_size = 0.25 * lambda
    Z_mlfma = MLFMAOperator(cfie, basis, leaf_size)
    
    # 6. Excitation
    # Match Legacy: Direction -x, Polarization +z
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(cfie, source, basis)
    
    # 7. Solve
    println("Solving with GMRES...")
    # Use GMRES with restart=50, tol=1e-3
    solver = GMRESSolver(restart=50, maxiter=200, tol=1e-3, verbose=true)
    
    # Preconditioner
    println("Setting up Near Field Preconditioner (LU)...")
    # Use LU factorization of Z_near
    P_near = lu(Z_mlfma.Z_near)
    
    # Wrap in a type that IterativeSolvers accepts
    # Note: We need to define this struct at top level if not already available
    # But since we are in a script, we can define it here or use the one from verify_SEFIE_mlfma.jl
    # Let's define it locally inside the function scope? No, must be global.
    
    I = solve!(solver, Z_mlfma, V, Pl=LUPreconditioner(P_near))
    
    # 8. RCS
    println("Calculating RCS...")
    θs_obs = collect(-180:0.5:180) .* (π/180)
    ϕs_obs = [0.0, π/2]
    
    RCS_res = radarCrossSection(θs_obs, ϕs_obs, I, basis)
    RCS_emsuite = RCS_res[2] # Total RCS linear scale (m²)
    
    # Save results
    RCS_Phi0 = abs.(RCS_emsuite[:, 1])
    RCS_Phi90 = abs.(RCS_emsuite[:, 2])
    
    df_out = DataFrame(
        Theta_Deg = θs_obs .* (180/π),
        RCS_Phi0_dB = 10log10.(RCS_Phi0),
        RCS_Phi90_dB = 10log10.(RCS_Phi90)
    )
    
    CSV.write(joinpath(OUTPUT_DIR, "SCFIE_MLFMA_Jet.csv"), df_out)
    println("Results saved to SCFIE_MLFMA_Jet.csv")
    
    # Compare with SCFIE Direct
    direct_file = joinpath(OUTPUT_DIR, "SCFIE_Direct_Jet.csv")
    if isfile(direct_file)
        df_direct = CSV.read(direct_file, DataFrame)
        println("Comparing with SCFIE Direct...")
        diff = df_out.RCS_Phi0_dB .- df_direct.RCS_Phi0_dB
        println("Mean Difference (MLFMA - Direct): $(mean(diff)) dB")
        println("RMSE: $(sqrt(mean(diff.^2))) dB")
    else
        println("Baseline file not found.")
    end
end

# Preconditioner Wrapper
struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

verify_SCFIE_mlfma()
