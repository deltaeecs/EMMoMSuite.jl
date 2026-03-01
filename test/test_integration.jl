using Test
using EMSuite
using LinearAlgebra
using StaticArrays

@testset "Full Integration Test" begin
    # 1. Setup Parameters
    freq = 300e6 # 300 MHz
    set_frequency!(freq)
    k0 = get_k0()

    # 2. Create Geometry (Simple Plate: 2 triangles)
    # Square plate in xy plane, z=0, side length 1.0
    # (-0.5, -0.5) to (0.5, 0.5)
    nodes = Matrix([
        -0.5 -0.5 0.0
        0.5 -0.5 0.0
        0.5 0.5 0.0
        -0.5 0.5 0.0
    ]') # 3x4 matrix

    # Two triangles: (1,2,3) and (1,3,4)
    elements = Matrix([
        1 2 3
        1 3 4
    ]') # 3x2 matrix

    mesh = TriangleMesh(size(elements, 2), nodes, elements)
    @test num_vertices(mesh) == 4
    @test num_elements(mesh) == 2

    # 3. Basis Functions
    basis = RWGBasis(mesh)
    n_unknowns = num_basis(basis)
    @test n_unknowns > 0

    # 4. Integral Equation (EFIE)
    ie = EFIE(freq)

    # 5. Matrix Assembly
    Z = assemble_impedance_matrix(ie, basis)
    @test size(Z) == (n_unknowns, n_unknowns)

    # 6. Excitation (Plane Wave)
    # Incident from +z direction, polarized in x
    theta_inc = 0.0
    phi_inc = 0.0
    # E_inc = exp(-jk z) x_hat
    # At z=0, E_inc = x_hat

    V = zeros(ComplexF64, n_unknowns)
    # Simple excitation vector calculation (approximation for test)
    # In a real scenario, we would use a proper excitation function
    # For now, let's just fill V with ones to test the solver
    fill!(V, 1.0 + 0.0im)

    # 7. Solve
    I = solve!(LUSolver(), Z, V)
    @test length(I) == n_unknowns
    @test !any(isnan.(I))

    # 8. Post-Processing
    # Calculate RCS at a few angles
    theta_obs = [0.0, pi / 4, pi / 2]
    phi_obs = [0.0]

    rcs_comp, rcs_total, rcs_db = radarCrossSection(theta_obs, phi_obs, I, basis)
    @test length(rcs_total) == length(theta_obs) * length(phi_obs)
    @test all(rcs_total .>= 0)

    # 9. I/O
    # Export VTK
    vtk_file = "test_integration_output"

    # Calculate J for visualization
    Jtris = geoElectricJCal(I, basis)
    Jmag = vec(sqrt.(sum(abs2.(Jtris), dims = 1)))

    save_vtk(vtk_file, mesh, Jmag; data_name = "CurrentMagnitude")
    @test isfile(vtk_file * ".vtu")
    rm(vtk_file * ".vtu")

    # Export Results
    h5_file = "test_integration_results.h5"
    save_results_hdf5(h5_file; freq = freq, theta = theta_obs, phi = phi_obs, rcs = rcs_total)
    @test isfile(h5_file)
    rm(h5_file)

    txt_file = "test_integration_rcs.txt"
    # Save Total RCS
    save_RCS_txt(txt_file, theta_obs, phi_obs, rcs_total)
    @test isfile(txt_file)
    rm(txt_file)

    println("Integration test completed successfully.")
end
