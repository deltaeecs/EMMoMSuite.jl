using EMMoMSuite
using EMMoMSuite.PostProcessing.RadiationIntegral
using EMMoMSuite.Utilities.Parameters: get_k0, get_eta0, set_frequency!
using LinearAlgebra
using StaticArrays
using Printf
using CSV
using DataFrames

# Verification: VSEFIE (Surface-Volume) against MoM_AllinOne/FEKO
# ---------------------------------------------------------------

function verify_vsefie_allinone()
    println("==================================================")
    println("   Verification: VSEFIE (Plate + Dielectric)      ")
    println("==================================================")
    
    # 1. Parameters
    freq = 1.2e9
    epsilon_r = 2.0 - 0.0004im
    
    # Set global parameters
    set_frequency!(freq)
    
    # 2. Mesh
    # EMMoMSuite is in MoM/EMMoMSuite
    # MoM_AllinOne is in MoM/MoM_AllinOne
    # @__DIR__ is MoM/EMMoMSuite/scripts/verification
    # dirname(@__DIR__) is MoM/EMMoMSuite/scripts
    # dirname(dirname(@__DIR__)) is MoM/EMMoMSuite
    # dirname(dirname(dirname(@__DIR__))) is MoM
    
    root_dir = dirname(dirname(dirname(@__DIR__)))
    mesh_file = joinpath(root_dir, "deps", "fixtures", "AllinOne", "meshfiles", "plate_and_metal_1dot2GHz.nas")
    println("Loading mesh from: $mesh_file")
    
    surface_mesh, volume_mesh = read_mixed_nas_mesh(mesh_file, scale=1.0)
    
    # Check mesh bounds
    min_bound = minimum(surface_mesh.node, dims=2)
    max_bound = maximum(surface_mesh.node, dims=2)
    println("Mesh Bounds: $min_bound to $max_bound")
    println("Mesh Size: $(max_bound - min_bound)")
    
    println("Surface Mesh: $(num_vertices(surface_mesh)) vertices, $(num_elements(surface_mesh)) triangles")
    println("Volume Mesh: $(num_vertices(volume_mesh)) vertices, $(num_elements(volume_mesh)) tetrahedra")
    
    # 3. Basis
    surface_basis = RWGBasis(surface_mesh)
    volume_basis = SWGBasis(volume_mesh)
    
    println("Surface Unknowns: $(num_basis(surface_basis))")
    println("Volume Unknowns: $(num_basis(volume_basis))")
    
    # 4. Operator
    nt_vol = num_elements(volume_mesh)
    permittivities = fill(Complex(epsilon_r), nt_vol)
    
    # alpha=0.0 for EFIE on surface
    scfie = SCFIE(freq, permittivities; alpha=0.0)
    
    # 5. Assembly
    println("Assembling Impedance Matrix...")
    t_start = time()
    Z = assemble_impedance_matrix(scfie, surface_basis, volume_basis)
    t_end = time()
    println(@sprintf("Assembly Time: %.4f s", t_end - t_start))
    println("Matrix size: $(size(Z))")
    
    # Debug: Check for singularity
    diags = diag(Z)
    min_diag = minimum(abs.(diags))
    println("Min diagonal magnitude: $min_diag")
    
    zero_diags = findall(x -> abs(x) < 1e-12, diags)
    if !isempty(zero_diags)
        println("Found $(length(zero_diags)) zero diagonals at indices: $(zero_diags[1:min(10, end)])...")
    end
    
    # 6. Excitation
    theta_inc = pi/2
    phi_inc = 0.0
    # theta_hat at pi/2, 0 is [0, 0, -1]
    pol = [0.0, 0.0, -1.0]
    
    source = PlaneWave(freq, theta_inc, phi_inc, pol)
    
    println("Excitation Vector...")
    V = excitation_vector(source, surface_basis, volume_basis)
    
    # Scale Volume part of V by 1/(im * omega) to match VEFIE scaling
    omega = 2π * freq
    N_surf = num_basis(surface_basis)
    V[N_surf+1:end] ./= (im * omega)
    
    # 7. Solve
    println("Solving...")
    I = Z \ V
    println("Solved. Current magnitude: $(norm(I))")
    
    # 8. RCS Calculation
    println("Calculating RCS...")
    
    # Load Reference Data
    ref_file = joinpath(root_dir, "deps", "fixtures", "AllinOne", "deps", "compare_feko", "plate_metal_1dot2GHzRCS.csv")
    df = CSV.read(ref_file, DataFrame; delim=' ', ignorerepeated=true)
    
    # Extract FEKO data (First 721 points are phi=0 cut)
    feko_rcs_m2 = df[1:721, "in"]
    feko_theta_deg = df[1:721, "THETA"]
    
    # Calculate RCS at these angles
    k0 = get_k0()
    eta0 = get_eta0()
    factor = (k0 * eta0)^2 / (4 * pi)
    
    N_surf = num_basis(surface_basis)
    I_surf = I[1:N_surf]
    I_vol = I[N_surf+1:end]
    
    println("I_surf norm: ", norm(I_surf))
    println("I_vol norm: ", norm(I_vol))
    println("Factor: ", factor)

    rcs_emsuite = Float64[]
    
    # Loop over angles
    for i in 1:length(feko_theta_deg)
        theta_deg = feko_theta_deg[i]
        theta = deg2rad(theta_deg)
        phi = 0.0 # phi=0 cut
        
        r_info = r̂θϕInfo(∠Info(theta), ∠Info(phi))
        
        N_s = radiation_integral_rwg(r_info, surface_basis, I_surf)
        N_v = radiation_integral_swg(r_info, volume_basis, I_vol, permittivities)
        
        # Scale N_v by j*omega because radiation_integral_swg returns N / (j*omega)
        # J_eq = j * omega * kappa * D
        # radiation_integral_swg computes Integral(kappa * D * exp)
        omega = k0 * 299792458.0
        N_v = (im * omega) * N_v
        
        if i == 1
            println("Angle 1: Theta=$theta_deg")
            println("N_s: ", N_s)
            println("N_v: ", N_v)
        end

        N_total = N_s + N_v
        
        sigma = factor * (abs2(N_total[1]) + abs2(N_total[2]))
        push!(rcs_emsuite, sigma)
    end
    
    # Compare
    diff = norm(rcs_emsuite - feko_rcs_m2) / norm(feko_rcs_m2)
    println("RCS Relative Difference (Norm): $diff")
    
    # Print some values
    println("Theta (deg) | FEKO (m^2) | EMMoMSuite (m^2)")
    for i in 1:50:length(feko_theta_deg)
        @printf("%10.2f | %10.4e | %10.4e\n", feko_theta_deg[i], feko_rcs_m2[i], rcs_emsuite[i])
    end
    
    if diff < 0.1
        println("Verification PASSED!")
    else
        println("Verification FAILED!")
    end
end

verify_vsefie_allinone()

