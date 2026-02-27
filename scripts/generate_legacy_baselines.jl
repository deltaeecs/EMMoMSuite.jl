using MoM_AllinOne
using DataFrames, CSV
using MoM_Visualizing

# Set up paths
const MOM_ALLINONE_DIR = "f:/OneDrive/MoM/MoM_AllinOne"
const OUTPUT_DIR = joinpath(@__DIR__, "../test_results/legacy_baseline")

if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
end

function run_case(case_name, mesh_file, freq, ie_type, solver_type, rtol=1e-3, restart=50)
    println("Running Legacy Case: $case_name")
    
    # Reset parameters
    setPrecision!(Float32)
    SimulationParams.SHOWIMAGE = false # Disable plotting during calculation

    # Parameters
    filename = joinpath(MOM_ALLINONE_DIR, "meshfiles", mesh_file)
    meshUnit = :m
    frequency = freq
    ieT = ie_type
    sbfT = :RWG
    vbfT = :nothing
    solverT = solver_type
    
    # Source
    source = PlaneWave(π/2, 0, 0f0, 1f0)

    # Observation
    θs_obs = LinRange{Precision.FT}(-π, π, 721)
    ϕs_obs = LinRange{Precision.FT}(0, π/2, 2)

    # Run Solver
    # We need to include the solver script from MoM_AllinOne
    # The solver scripts use variables defined in the global scope
    # So we need to make sure variables are available.
    # However, `include` works in the current scope.
    
    solver_script = solver_type == :direct ? "direct_solver.jl" : "fast_solver.jl"
    solver_path = joinpath(MOM_ALLINONE_DIR, "src", solver_script)
    
    # We wrap the include in a let block or similar if possible, but the solver scripts
    # likely expect global variables. Let's try to run them in the global scope of this function
    # or just include them. Since we are running cases sequentially, we can just overwrite variables.
    
    # Note: The solver scripts in MoM_AllinOne rely on `filename`, `meshUnit`, etc. being defined.
    # We need to make sure they are accessible.
    
    # To avoid scope issues with `include`, we might need to eval or just copy the logic.
    # But `include` should work if we are at top level. 
    # Since we are in a function, `include` will include into the function scope? 
    # No, `include` operates at global scope usually unless parsed.
    # Actually, `include` in a function includes into the global scope of the module.
    # This might be tricky if we want to run multiple cases.
    
    # Better approach: Write a separate script for each case or use a macro?
    # Or just set global variables and include.
    
    # Let's try setting global variables in the Main module and including.
    
    result = @eval Main begin
        filename = $filename
        meshUnit = $(QuoteNode(meshUnit))
        frequency = $frequency
        ieT = $(QuoteNode(ieT))
        sbfT = $(QuoteNode(sbfT))
        vbfT = $(QuoteNode(vbfT))
        solverT = $(QuoteNode(solverT))
        rtol = $rtol
        restart = $restart
        source = $source
        θs_obs = $θs_obs
        ϕs_obs = $ϕs_obs
        
        include($solver_path)
    end
    
    # Extract results from the return value of include
    # The solver scripts return (RCSθsϕs, RCSθsϕsdB, RCS, RCSdB)
    RCS = result[3]
    
    # Save results
    output_file = joinpath(OUTPUT_DIR, "$(case_name).csv")
    
    # RCS is likely a matrix [theta, phi]
    # We save theta and RCS for phi=0 and phi=90
    
    df = DataFrame(
        Theta_Rad = θs_obs,
        Theta_Deg = rad2deg.(θs_obs),
        RCS_Phi0_dBsm = 10log10.(RCS[:, 1]),
        RCS_Phi90_dBsm = 10log10.(RCS[:, 2])
    )
    
    CSV.write(output_file, df)
    println("Saved results to $output_file")
end

# 1. SEFIE Direct (Jet)
run_case("SEFIE_Direct_Jet", "jet_100MHz.nas", 1e8, :EFIE, :direct)

# 2. SEFIE MLFMA (Jet)
# run_case("SEFIE_MLFMA_Jet", "jet_100MHz.nas", 1e8, :EFIE, :gmres)

# 3. SCFIE Direct (Sphere)
run_case("SCFIE_Direct_Sphere", "sphere_600MHz.nas", 6e8, :CFIE, :direct)

# 4. SCFIE MLFMA (Sphere)
# run_case("SCFIE_MLFMA_Sphere", "sphere_600MHz.nas", 6e8, :CFIE, :gmres)

println("All legacy baselines generated.")
