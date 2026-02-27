
using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Solvers
using EMSuite.PostProcessing
using EMSuite.CoreModule.Sources
using LinearAlgebra
using Printf
# using CSV
# using DataFrames

function verify_vefie_plate()
    println("==================================================")
    println("   Verification: VEFIE Plate (Legacy Match)       ")
    println("==================================================")

    # 1. Parameters (Match Legacy examples_VEFIE_direct.jl)
    freq = 12e8 # 1.2 GHz
    
    # Set global parameters
    EMSuite.Utilities.Parameters.set_frequency!(freq)
    
    # 2. Mesh
    mesh_file = joinpath(@__DIR__, "../../../MoM_AllinOne/meshfiles/plate_1dot2GHz.nas")
    if !isfile(mesh_file)
        println("Error: Mesh file not found: $mesh_file")
        return
    end
    
    println("Loading mesh from $mesh_file...")
    mesh = read_nas_mesh(mesh_file)
    
    # Legacy uses meshUnit = :m. So we do NOT scale if it's already in m.
    # But we should check if the file is actually in m.
    # If Legacy says :m, it assumes the numbers in file are meters.
    # So we should NOT scale.
    
    println("Mesh: $(num_vertices(mesh)) vertices, $(num_elements(mesh)) tetrahedra")
    
    # 3. Basis
    println("Setting up SWG Basis...")
    basis = SWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")

    # 4. Permittivity
    # Legacy uses: setGeosPermittivity!(geosInfo, 2(1-0.0002im))
    eps_r = 2.0 * (1.0 - 0.0002im)
    permittivities = fill(eps_r, num_elements(mesh))

    # 5. Operator
    println("Setting up VEFIE...")
    vefie = VEFIE(freq, permittivities)
    
    # 6. Assembly
    println("Assembling Impedance Matrix...")
    Z = assemble_impedance_matrix(vefie, basis, permittivities)
    
    # 7. Excitation
    println("Setting up Excitation...")
    # Legacy: source = PlaneWave(π/4, 0, 0f0, 1f0) -> Theta=45deg, Phi=0.
    # Polarization: -Theta_hat = -[cos(45)cos(0), cos(45)sin(0), -sin(45)]
    # = -[0.70710678, 0, -0.70710678] = [-0.70710678, 0.0, 0.70710678]
    
    pol = [-0.70710678, 0.0, 0.70710678]
    source = PlaneWave(freq, π/4, 0.0, pol) 
    V = excitation_vector(vefie, source, basis, permittivities)
    
    # 8. Solve
    println("Solving...")
    solver = LUSolver()
    D = solve!(solver, Z, V)
    
    # 9. RCS
    println("Calculating RCS...")
    # theta_obs = collect(LinRange(-π, π, 721))
    theta_obs = [0.0] # Just one angle for debug
    phi_obs = zeros(length(theta_obs)) # Phi=0 cut
    
    # Let's calculate for Phi=0 first.
    rcs_comp, rcs_total, rcs_db = radarCrossSection(theta_obs, phi_obs, D, basis, permittivities)
    
    println("RCS at Theta=0: $(rcs_db[1]) dBsm")
end

verify_vefie_plate()
