using EMSuite
using LinearAlgebra
using DataFrames, CSV
using Printf

# Paths
const MOM_ALLINONE_DIR = "f:/OneDrive/MoM/MoM_AllinOne"
const LEGACY_BASELINE_DIR = joinpath(@__DIR__, "../test_results/legacy_baseline")
const OUTPUT_DIR = joinpath(@__DIR__, "../test_results/emsuite_verification")

if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
end

function verify_SEFIE_direct()
    println("=== Verifying SEFIE Direct (Jet) ===")
    println("Threads: $(Threads.nthreads())")
    
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
    println("Setting up EFIE...")
    efie = EFIE(freq)
    
    # 5. Assembly
    println("Assembling matrix...")
    t_assembly = @elapsed begin
        Z = assemble_impedance_matrix(efie, basis)
    end
    println("Assembly Time: $(t_assembly) seconds")
    
    # 6. Excitation
    # Legacy: PlaneWave(π/2, 0, 0f0, 1f0) -> Theta=90, Phi=0, Alpha=0, V=1.
    # Legacy PlaneWave:
    #   Direction: -r_hat. At (90, 0), r_hat is +x. So Direction is -x.
    #   Polarization: -theta_hat. At (90, 0), theta_hat is -z. So Polarization is +z.
    # EMSuite PlaneWave:
    #   Direction: k_dir. We want -x. So Theta=90, Phi=180 (π).
    #   Polarization: We want +z. So [0, 0, 1].
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(efie, source, basis)
    
    # 7. Solve
    println("Solving...")
    solver = LUSolver()
    t_solve = @elapsed begin
        I = solve!(solver, Z, V)
    end
    println("Solve Time: $(t_solve) seconds")
    
    # 8. RCS
    println("Calculating RCS...")
    θs_obs = collect(LinRange(-π, π, 721))
    ϕs_obs = [0.0, π/2]
    
    # Call RCS
    # radarCrossSection(θs_obs, ϕs_obs, ICoeff, basis)
    RCS_res = radarCrossSection(θs_obs, ϕs_obs, I, basis)
    # RCS_res returns (RCSθsϕs, RCSθsϕsdB, RCS, RCSdB)
    # RCS is [theta, phi]
    
    RCS_emsuite = RCS_res[2] # Total RCS linear scale (m²)
    
    # 9. Compare
    baseline_file = joinpath(LEGACY_BASELINE_DIR, "SEFIE_Direct_Jet.csv")
    if isfile(baseline_file)
        df_base = CSV.read(baseline_file, DataFrame)
        RCS_legacy_phi0_dB = df_base.RCS_Phi0_dBsm
        
        RCS_emsuite_phi0_dB = 10log10.(RCS_emsuite[:, 1])
        
        # Calculate difference
        diff = RCS_emsuite_phi0_dB .- RCS_legacy_phi0_dB
        mean_diff = sum(diff) / length(diff)
        println("Mean Difference (EMSuite - Legacy): $(mean_diff) dB")
        
        # RMSE after removing mean offset
        rmse = sqrt(sum((diff .- mean_diff).^2) / length(diff))
        println("RMSE (after removing mean offset): $(rmse) dB")
        
        # Save results
        df_out = DataFrame(
            Theta_Rad = θs_obs,
            RCS_EMSuite_Phi0_dB = RCS_emsuite_phi0_dB,
            RCS_Legacy_Phi0_dB = RCS_legacy_phi0_dB,
            Diff_dB = diff
        )
        CSV.write(joinpath(OUTPUT_DIR, "SEFIE_Direct_Jet_Comparison.csv"), df_out)
        println("Comparison saved to SEFIE_Direct_Jet_Comparison.csv")
    else
        println("Baseline file not found.")
    end
end

verify_SEFIE_direct()
