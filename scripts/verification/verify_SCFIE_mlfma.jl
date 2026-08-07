
using Pkg
Pkg.activate(joinpath(@__DIR__, "../../"))
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.Solvers
using EMMoMSuite.PostProcessing
using EMMoMSuite.CoreModule.Sources
using EMMoMSuite.FastAlgorithms.MLFMA
using LinearAlgebra
using SparseArrays
using StaticArrays
using Printf

using EMMoMSuite.PostProcessing.RadiationIntegral: r̂θϕInfo, radiation_integral_rwg, radiation_integral_swg, ∠Info

function verify_scfie_mlfma()
    println("==================================================")
    println("   Verification: VSIE MLFMA (Plate+Metal)         ")
    println("==================================================")

    # 1. Parameters
    freq = 12e8 # 1.2 GHz
    EMMoMSuite.Utilities.Parameters.set_frequency!(freq)
    lambda = 299792458.0 / freq
    
    # 2. Mesh
    mesh_file = joinpath(@__DIR__, "../../../deps/fixtures/AllinOne/meshfiles/plate_and_metal_1dot2GHz.nas")
    if !isfile(mesh_file)
        println("Error: Mesh file not found: $mesh_file")
        return
    end
    
    println("Loading mesh from $mesh_file...")
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
    eps_r = 2.0 * (1.0 - 0.0002im)
    permittivities = fill(eps_r, num_elements(vol_mesh))
    
    # 5. Operator
    println("Setting up VSIE (SCFIE with alpha=1.0)...")
    vsie = SCFIE(freq, permittivities; alpha=1.0)
    
    # 6. MLFMA Operator
    println("Setting up MLFMA Operator...")
    leaf_size = 0.23 * lambda # Match legacy
    mlfma_op = MLFMAOperator(vsie, [basis_surf, basis_vol], leaf_size)
    
    println("MLFMA Octree Levels: $(mlfma_op.octree.nLevels)")
    println("Near Field Matrix: $(size(mlfma_op.Z_near)) with $(nnz(mlfma_op.Z_near)) non-zeros")
    
    # 7. Excitation
    println("Setting up Excitation...")
    pol = [-0.70710678, 0.0, 0.70710678]
    source = PlaneWave(freq, π/4, 0.0, pol) 
    
    cfie = CFIE(freq, vsie.alpha)
    vefie = VEFIE(freq, permittivities)
    
    V_surf = excitation_vector(cfie, source, basis_surf)
    V_vol = excitation_vector(vefie, source, basis_vol, permittivities)
    V = [V_surf; V_vol]
    
    println("RHS Norm: $(norm(V))")
    if norm(V) == 0
        println("Error: RHS is zero!")
        return
    end
    
    # 8. Solve
    println("Solving with GMRES...")
    # Use a simple identity preconditioner or diagonal if available
    # For now, no preconditioner
    solver_conf = GMRESSolver(tol=1e-3, maxiter=200, restart=50)
    I = solve!(solver_conf, mlfma_op, V)
    
    # Split solution
    N_surf = num_basis(basis_surf)
    I_surf = I[1:N_surf]
    I_vol = I[N_surf+1:end]
    
    # 9. RCS
    println("Calculating RCS... (Skipped for debugging)")
    theta_obs = [0.0]
    # ...
    phi_obs = zeros(length(theta_obs))
    
    # Manual RCS calculation for mixed bases
    k = 2π * freq / 299792458.0
    rcs_db = Float64[]
    
    for i in 1:length(theta_obs)
        θ = theta_obs[i]
        ϕ = phi_obs[i]
        
        angle_theta = ∠Info(θ)
        angle_phi = ∠Info(ϕ)
        angle_info = r̂θϕInfo(angle_theta, angle_phi)
        
        # Surface Contribution
        N_s = radiation_integral_rwg(angle_info, basis_surf, I_surf)
        
        # Volume Contribution
        N_v = radiation_integral_swg(angle_info, basis_vol, I_vol, permittivities)
        
        # Total N
        N_total = N_s + N_v
        
        # E_far = -j k eta / (4 pi r) * exp(-jkr) * (I - rr) * N
        # RCS = 4 pi r^2 |E_far|^2 / |E_inc|^2
        # |E_far| propto k eta / (4 pi) * |N_transverse|
        # RCS = 4 pi * (k eta / 4 pi)^2 * |N_transverse|^2
        # RCS = k^2 eta^2 / (4 pi) * |N_transverse|^2
        
        # Wait, radiation_integral returns N vector.
        # We need transverse component.
        # N_total is already [N_theta, N_phi]
        N_theta = N_total[1]
        N_phi = N_total[2]
        
        N_mag2 = abs2(N_theta) + abs2(N_phi)
        
        eta = 376.730313668
        sigma = (k * eta)^2 / (4π) * N_mag2
        
        sigma_db = 10 * log10(sigma)
        push!(rcs_db, sigma_db)
        
        @printf("Theta: %.2f, RCS: %.4f dBsm\n", rad2deg(θ), sigma_db)
    end
    
    # Compare with Direct Solver result (-15.35 dBsm)
    expected = -15.35
    diff = abs(rcs_db[1] - expected)
    println("Difference from Direct Solver: $(diff) dB")
    
    if diff < 1.0
        println("VERIFICATION PASSED")
    else
        println("VERIFICATION FAILED")
    end
end

verify_scfie_mlfma()
