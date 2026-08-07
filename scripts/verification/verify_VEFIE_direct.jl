using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.Solvers
using EMMoMSuite.PostProcessing
using EMMoMSuite.CoreModule.Sources
using LinearAlgebra
using Printf

function verify_vefie_direct()
    println("==================================================")
    println("   Verification: VEFIE Direct Solver              ")
    println("==================================================")

    # 1. Parameters
    freq = 300e6
    c0 = 299792458.0
    lambda = c0 / freq
    
    # Set global parameters
    EMMoMSuite.Utilities.Parameters.set_frequency!(freq)
    println("Global parameters set: k0=$(EMMoMSuite.Utilities.Parameters.get_k0())")
    
    # 2. Mesh
    mesh_file = joinpath(@__DIR__, "../../../deps/fixtures/AllinOne/meshfiles/Tetra.nas")
    if !isfile(mesh_file)
        println("Error: Mesh file not found: $mesh_file")
        return
    end
    
    println("Loading mesh from $mesh_file...")
    mesh = read_nas_mesh(mesh_file)
    
    # Scale mesh (assuming mm to m)
    # Check bounding box first
    nodes = mesh.node
    max_coord = maximum(nodes)
    if max_coord > 10.0
        println("Scaling mesh by 0.001 (mm -> m)...")
        mesh.node .*= 0.001
    end
    
    # Center the mesh
    min_c = minimum(mesh.node, dims=2)
    max_c = maximum(mesh.node, dims=2)
    center = (min_c + max_c) / 2
    println("Centering mesh... Old Center: $center")
    mesh.node .-= center
    println("New Center: $((minimum(mesh.node, dims=2) + maximum(mesh.node, dims=2)) / 2)")
    
    println("Mesh: $(num_vertices(mesh)) vertices, $(num_elements(mesh)) tetrahedra")
    
    if !(mesh isa TetrahedraMesh)
        println("Error: Mesh is not a TetrahedraMesh.")
        return
    end

    # 3. Basis
    println("Setting up SWG Basis...")
    basis = SWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")

    # 4. Permittivity
    # Assume homogeneous dielectric sphere with eps_r = 2.0
    eps_r = 2.0 + 0.0im
    permittivities = fill(eps_r, num_elements(mesh))

    # 5. Operator
    println("Setting up VEFIE...")
    vefie = VEFIE(freq, permittivities)
    
    # 6. Assembly
    println("Assembling Impedance Matrix...")
    t_asm = @elapsed begin
        Z = assemble_impedance_matrix(vefie, basis, permittivities)
    end
    println("        Done in $(t_asm) s")
    println("        Matrix Size: $(size(Z))")
    
    # 7. Excitation
    println("Setting up Excitation...")
    # Plane wave incident from +z, x-polarized
    source = PlaneWave(freq, π, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(vefie, source, basis, permittivities)
    println("        RHS Norm: $(norm(V))")
    
    # 8. Solve
    println("Solving...")
    solver = LUSolver()
    t_solve = @elapsed begin
        D = solve!(solver, Z, V)
    end
    println("        Done in $(t_solve) s")
    println("        Solution Norm: $(norm(D))")
    
    # 9. Check Symmetry?
    # VEFIE matrix should be symmetric if formulated correctly.
    is_sym = isapprox(Z, transpose(Z), rtol=1e-5)
    println("        Symmetric: $is_sym")
    
    if is_sym
        println("SUCCESS: VEFIE Matrix is symmetric.")
    else
        println("WARNING: VEFIE Matrix is NOT symmetric.")
    end
    
    # 10. RCS Calculation
    println("Calculating RCS...")
    theta_range = collect(0:1.0:180.0) .* (π/180.0)
    phi_range = [0.0, π/2.0] # E-plane (0) and H-plane (90)
    
    # Need to pass permittivities to radarCrossSection for VEFIE
    # But radarCrossSection signature might need update or we call radiation_integral directly?
    # I updated RCS.jl to accept permittivities.
    
    RCS_components, RCS_total, RCS_dB = radarCrossSection(theta_range, phi_range, D, basis, permittivities)
    
    # Save RCS to file
    open("rcs_vefie_direct.txt", "w") do io
        println(io, "Theta(deg) Phi(deg) RCS(dBsm)")
        for p in 1:length(phi_range)
            phi_deg = rad2deg(phi_range[p])
            for t in 1:length(theta_range)
                theta_deg = rad2deg(theta_range[t])
                val = RCS_dB[t, p]
                println(io, "$theta_deg $phi_deg $val")
            end
        end
    end
    println("        RCS saved to rcs_vefie_direct.txt")

end

verify_vefie_direct()
