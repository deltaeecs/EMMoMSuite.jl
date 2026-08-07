
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.Solvers
using EMMoMSuite.PostProcessing
using EMMoMSuite.CoreModule.Sources
using LinearAlgebra
using Printf

using EMMoMSuite.PostProcessing.RadiationIntegral: r̂θϕInfo, radiation_integral_rwg, radiation_integral_swg, ∠Info

function verify_vsie_plate()
    println("==================================================")
    println("   Verification: VSIE Plate+Metal (Legacy Match)  ")
    println("==================================================")

    # 1. Parameters
    freq = 12e8 # 1.2 GHz
    EMMoMSuite.Utilities.Parameters.set_frequency!(freq)
    
    # 2. Mesh
    mesh_file = joinpath(@__DIR__, "../../../deps/fixtures/AllinOne/meshfiles/plate_and_metal_1dot2GHz.nas")
    if !isfile(mesh_file)
        println("Error: Mesh file not found: $mesh_file")
        return
    end
    
    println("Loading mesh from $mesh_file...")
    # Legacy uses meshUnit = :m.
    surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file)
    
    println("Surface Mesh: $(num_elements(surf_mesh)) triangles")
    println("Volume Mesh: $(num_elements(vol_mesh)) tetrahedra")
    
    # 3. Basis
    println("Setting up Bases...")
    basis_surf = RWGBasis(surf_mesh)
    basis_vol = SWGBasis(vol_mesh)
    
    println("Surface Unknowns: $(num_basis(basis_surf))")
    println("Volume Unknowns: $(num_basis(basis_vol))")
    
    # 4. Permittivity
    # Legacy uses: setGeosPermittivity!(geosInfo, 2(1-0.0002im))
    eps_r = 2.0 * (1.0 - 0.0002im)
    permittivities = fill(eps_r, num_elements(vol_mesh))
    
    # 5. Operator
    println("Setting up VSIE (SCFIE with alpha=0.0)...")
    # SCFIE(freq, permittivities; alpha=0.0) -> EFIE formulation
    vsie = SCFIE(freq, permittivities; alpha=0.0)
    
    # 6. Assembly
    println("Assembling Impedance Matrix...")
    # We need a function that takes both bases.
    # assemble_impedance_matrix(vsie, basis_surf, basis_vol) ?
    # Let's check SCFIE.jl
    
    Z = assemble_impedance_matrix(vsie, basis_surf, basis_vol)
    
    # 7. Excitation
    println("Setting up Excitation...")
    pol = [-0.70710678, 0.0, 0.70710678]
    source = PlaneWave(freq, π/4, 0.0, pol) 
    
    # We need excitation vector for both
    # Construct CFIE and VEFIE operators to call excitation_vector
    cfie = CFIE(freq, vsie.alpha)
    vefie = VEFIE(freq, permittivities)
    
    V_surf = excitation_vector(cfie, source, basis_surf)
    V_vol = excitation_vector(vefie, source, basis_vol, permittivities)
    # Check VEFIE excitation_vector signature.
    # It seems excitation_vector(op, source, basis) is available.
    # But wait, VEFIE excitation vector is Integral(f . E_inc).
    # It does NOT depend on permittivity.
    # So excitation_vector(vefie, source, basis_vol) should work if defined.
    # Or excitation_vector(source, basis_vol).
    
    # Let's try excitation_vector(vefie, source, basis_vol)
    # If not defined, try excitation_vector(source, basis_vol)
    
    V = [V_surf; V_vol]
    
    # 8. Solve
    println("Solving...")
    solver = LUSolver()
    I = solve!(solver, Z, V)
    
    # Split solution
    N_surf = num_basis(basis_surf)
    I_surf = I[1:N_surf]
    I_vol = I[N_surf+1:end]
    
    # 9. RCS
    println("Calculating RCS...")
    theta_obs = [0.0] # Just one angle for debug
    phi_obs = zeros(length(theta_obs))
    
    # radarCrossSection needs to handle both?
    # Or we calculate N_surf and N_vol and sum them?
    # RCS module might not have a combined function yet.
    # We can calculate N separately and sum.
    
    # N_total = N_surf + N_vol
    # But we need N vector (theta, phi components).
    
    # Let's check RCS.jl
    
    # If not available, we can implement a simple check here.
    
    # Calculate N_surf
    # We need to call radiation_integral_rwg directly?
    # Or use radarCrossSection for each and combine?
    # radarCrossSection returns RCS, not N.
    
    # We can use radiation_integral_rwg and radiation_integral_swg from RadiationIntegral module.
    
    k0 = EMMoMSuite.Utilities.Parameters.get_k0()
    eta0 = EMMoMSuite.Utilities.Parameters.get_eta0()
    factor = (k0 * eta0)^2 / (4 * pi)
    
    r_info = r̂θϕInfo(∠Info(0.0), ∠Info(0.0)) # Theta=0, Phi=0
    
    N_surf = radiation_integral_rwg(r_info, basis_surf, I_surf)
    N_vol = radiation_integral_swg(r_info, basis_vol, I_vol, permittivities)
    
    N_total = N_surf + N_vol
    
    rcs_val = factor * (abs2(N_total[1]) + abs2(N_total[2]))
    rcs_db = 10 * log10(rcs_val)
    
    println("RCS at Theta=0: $rcs_db dBsm")

    
end

verify_vsie_plate()
