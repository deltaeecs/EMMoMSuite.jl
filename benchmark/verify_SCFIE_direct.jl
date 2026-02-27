using EMSuite
using LinearAlgebra
using DataFrames, CSV
using Printf
using Statistics

# Paths
const MOM_ALLINONE_DIR = "f:/OneDrive/MoM/MoM_AllinOne"
const LEGACY_BASELINE_DIR = joinpath(@__DIR__, "../test_results/legacy_baseline")
const OUTPUT_DIR = joinpath(@__DIR__, "../test_results/emsuite_verification")

if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
end

function verify_SCFIE_direct()
    println("=== Verifying SCFIE Direct (Jet) ===")
    
    # 1. Load Mesh
    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles/jet_100MHz.nas")
    println("Loading mesh: $mesh_file")
    
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    n_elems = num_elements(mesh)
    
    # 2. Parameters
    freq = 1e8
    set_frequency!(freq)
    
    # 3. Basis
    println("Setting up RWG basis...")
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Number of unknowns: $N")
    
    # 4. Equation
    println("Setting up CFIE (alpha=0.5)...")
    cfie = CFIE(freq, 0.5)
    
    # 5. Assembly
    println("Assembling matrix...")
    Z = assemble_impedance_matrix(cfie, basis)
    
    # 6. Excitation
    # Match Legacy: Direction -x, Polarization +z
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(cfie, source, basis)
    
    # 7. Solve
    println("Solving...")
    I = Z \ V
    
    # 8. RCS
    println("Calculating RCS...")
    # Match Legacy resolution: 0.5 degree steps from -180 to 180?
    # Legacy file has 721 points (-180 to 180).
    # My previous run used 0:0.5:180 (361 points).
    θs_obs = collect(-180:0.5:180) .* (π/180)
    ϕs_obs = [0.0, π/2]
    
    RCS_res = radarCrossSection(θs_obs, ϕs_obs, I, basis)
    RCS_emsuite = RCS_res[2] # Total RCS linear scale (m²)
    
    # Save results
    # Ensure RCS is positive before log10
    RCS_Phi0 = abs.(RCS_emsuite[:, 1])
    RCS_Phi90 = abs.(RCS_emsuite[:, 2])
    
    df_out = DataFrame(
        Theta_Deg = θs_obs .* (180/π),
        RCS_Phi0_dB = 10log10.(RCS_Phi0),
        RCS_Phi90_dB = 10log10.(RCS_Phi90)
    )
    
    CSV.write(joinpath(OUTPUT_DIR, "SCFIE_Direct_Jet.csv"), df_out)
    println("Results saved to SCFIE_Direct_Jet.csv")
    
    # Compare with SEFIE Direct (just to see difference)
    sefie_file = joinpath(OUTPUT_DIR, "SEFIE_Direct_Jet_Comparison.csv")
    if isfile(sefie_file)
        df_sefie = CSV.read(sefie_file, DataFrame)
        # Interpolate or match indices?
        # Assuming same theta grid
        println("Comparing with SEFIE Direct...")
        diff = df_out.RCS_Phi0_dB .- df_sefie.RCS_EMSuite_Phi0_dB
        println("Mean Difference (CFIE - EFIE): $(mean(diff)) dB")
    end
end

verify_SCFIE_direct()
