using EMSuite
using LinearAlgebra
using DataFrames, CSV
using Printf

# Paths
const MOM_ALLINONE_DIR = joinpath(@__DIR__, "../../MoM_AllinOne")
const LEGACY_BASELINE_DIR = joinpath(@__DIR__, "../test_results/legacy_baseline")
const OUTPUT_DIR = joinpath(@__DIR__, "../test_results/emsuite_verification")

if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
end

# Wrap in a type that IterativeSolvers accepts
struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

function verify_SEFIE_mlfma()
    println("=== Verifying SEFIE MLFMA (Jet) ===")
    
    # 1. Load Mesh
    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles/jet_100MHz.nas")
    println("Loading mesh: $mesh_file")
    
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    
    # 2. Parameters
    freq = 1e8
    set_frequency!(freq)
    lambda = 299792458.0 / freq
    
    # 3. Basis
    println("Setting up RWG basis...")
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Number of unknowns: $N")
    
    # 4. Equation
    println("Setting up EFIE...")
    efie = EFIE(freq)
    
    # 5. MLFMA Operator
    println("Setting up MLFMA Operator...")
    # Leaf size 0.35 lambda to improve near-field preconditioner quality
    leaf_size = 0.35 * lambda
    Z_mlfma = MLFMAOperator(efie, basis, leaf_size)
    
    # 6. Excitation
    # Match Legacy: Direction -x, Polarization +z
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(efie, source, basis)
    
    # 7. Solve
    println("Solving with GMRES...")
    # Use GMRES with restart=50, tol=1e-3 (Legacy uses 1e-3)
    solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
    
    # Preconditioner
    println("Setting up Near Field Preconditioner (LU)...")
    # Use LU factorization of Z_near
    P_near = lu(Z_mlfma.Z_near)
    
    P = LUPreconditioner(P_near)
    
    I = solve!(solver, Z_mlfma, V, Pl=P)
    
    # 8. RCS
    println("Calculating RCS...")
    θs_obs = collect(LinRange(-π, π, 721))
    ϕs_obs = [0.0, π/2]
    
    RCS_res = radarCrossSection(θs_obs, ϕs_obs, I, basis)
    RCS_emsuite = RCS_res[2] # Linear scale
    
    # 9. Compare with Legacy Direct
    baseline_file = joinpath(LEGACY_BASELINE_DIR, "SEFIE_Direct_Jet.csv")
    if isfile(baseline_file)
        df_base = CSV.read(baseline_file, DataFrame)
        RCS_legacy_phi0_dB = df_base.RCS_Phi0_dBsm
        
        RCS_emsuite_phi0_dB = 10log10.(RCS_emsuite[:, 1])
        
        # Calculate difference
        diff = RCS_emsuite_phi0_dB .- RCS_legacy_phi0_dB
        mean_diff = sum(diff) / length(diff)
        println("Mean Difference (EMSuite MLFMA - Legacy Direct): $(mean_diff) dB")
        
        rmse = sqrt(sum((diff .- mean_diff).^2) / length(diff))
        println("RMSE (after removing mean offset): $(rmse) dB")
        
        # Save results
        df_out = DataFrame(
            Theta_Rad = θs_obs,
            RCS_EMSuite_MLFMA_Phi0_dB = RCS_emsuite_phi0_dB,
            RCS_Legacy_Direct_Phi0_dB = RCS_legacy_phi0_dB,
            Diff_dB = diff
        )
        CSV.write(joinpath(OUTPUT_DIR, "SEFIE_MLFMA_Jet_Comparison.csv"), df_out)
        println("Comparison saved to SEFIE_MLFMA_Jet_Comparison.csv")
    else
        println("Baseline file not found.")
    end
end

verify_SEFIE_mlfma()
