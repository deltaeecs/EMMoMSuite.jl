using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using LinearAlgebra
using StaticArrays
using Printf
using CSV
using DataFrames

# Import necessary modules/types that might not be exported
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.CoreModule
using EMSuite.PostProcessing.RCS

function run_reproduction()
    # 1. Load Mesh
    # Assuming the script is run from EMSuite/scripts/
    mesh_file = joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles", "Tri.nas")
    if !isfile(mesh_file)
        # Try absolute path if relative fails (e.g. if running from VS Code root)
        mesh_file = "f:\\OneDrive\\MoM\\MoM_AllinOne\\meshfiles\\Tri.nas"
    end
    
    println("Loading mesh from $mesh_file")
    mesh = read_nas_mesh(mesh_file, scale=0.001)
    println("Mesh loaded: $(num_vertices(mesh)) vertices, $(num_elements(mesh)) elements")

    # Check mesh bounds
    verts = mesh.node
    min_v = minimum(verts, dims=2)
    max_v = maximum(verts, dims=2)
    println("Mesh Bounds:")
    println("  X: $(min_v[1]) to $(max_v[1])")
    println("  Y: $(min_v[2]) to $(max_v[2])")
    println("  Z: $(min_v[3]) to $(max_v[3])")

    # 2. Setup Parameters
    freq = 1.0e9 # 1 GHz
    c0 = 299792458.0
    lambda = c0 / freq
    println("Frequency: $freq Hz, Lambda: $lambda m")
    set_frequency!(freq)

    # 3. Basis Functions
    println("Initializing RWG Basis...")
    basis = RWGBasis(mesh)
    n_unknowns = num_basis(basis)
    println("Number of unknowns: $n_unknowns")

    # 4. Integral Equation
    println("Initializing EFIE...")
    ie = EFIE(freq)

    # 5. Assembly
    println("Assembling Impedance Matrix...")
    @time Z = assemble_impedance_matrix(ie, basis)
    println("Z matrix size: $(size(Z))")
    println("Z matrix norm: $(norm(Z))")
    println("Z matrix sample (1,1): $(Z[1,1])")

    # 6. Excitation
    # Normal incidence (from +z direction, propagating in -z direction)
    # theta = pi, phi = 0. Polarization along x.
    theta_inc = pi
    phi_inc = 0.0
    pol = [1.0, 0.0, 0.0]
    source = PlaneWave(freq, theta_inc, phi_inc, pol)

    println("Assembling Excitation Vector...")
    @time V = excitation_vector(ie, source, basis)
    println("Norm of V: $(norm(V))")
    if norm(V) == 0
        println("WARNING: Excitation vector is zero!")
    end

    # 7. Solve
    println("Solving System...")
    @time I = solve!(LUSolver(), Z, V)
    println("Norm of I: $(norm(I))")
    if norm(I) == 0
        println("WARNING: Current vector is zero!")
    end

    # 8. RCS Calculation
    println("Calculating RCS...")
    
    # Prepare TriangleInfo with Basis IDs
    nt = num_elements(mesh)
    # We need to use the constructor from Geometry module
    trianglesInfo = [TriangleInfo(mesh, i) for i in 1:nt]
    
    # Populate inBfsID
    for bf in basis.functions
        # Plus triangle
        t_plus = bf.support[1]
        if t_plus > 0
            local_edge = bf.local_edge_idx[1]
            ids = Vector(trianglesInfo[t_plus].inBfsID)
            ids[local_edge] = bf.id # Positive ID
            trianglesInfo[t_plus].inBfsID = SVector{3, Int}(ids)
        end
        
        # Minus triangle
        t_minus = bf.support[2]
        if t_minus > 0
            local_edge = bf.local_edge_idx[2]
            ids = Vector(trianglesInfo[t_minus].inBfsID)
            ids[local_edge] = -bf.id # Negative ID
            trianglesInfo[t_minus].inBfsID = SVector{3, Int}(ids)
        end
    end
    
    # Observation angles
    theta_obs = collect(range(-pi, pi, length=361))
    phi_obs = [0.0] # Observation in xz plane (phi=0)

    # Calculate RCS
    # radarCrossSection(θs_obs, ϕs_obs, ICoeff, trianglesInfo, BFT)
    # BFT is RWG
    RCS_components, RCS_total, RCS_dB = radarCrossSection(theta_obs, phi_obs, I, trianglesInfo, RWG)
    
    # 9. Save Results
    output_file = joinpath(@__DIR__, "rcs_plate_emsuite.csv")
    println("Saving results to $output_file")
    
    df = DataFrame(Theta_deg = rad2deg.(theta_obs), RCS_dB = vec(RCS_dB))
    CSV.write(output_file, df)
    
    println("Reproduction script completed.")
end

run_reproduction()
