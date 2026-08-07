using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.Solvers
using EMMoMSuite.PostProcessing
using EMMoMSuite.Utilities.Parameters
using EMMoMSuite.CoreModule.Sources
using StaticArrays
using LinearAlgebra
using Printf
using DelimitedFiles

function verify_vefie_legacy()
    println("==================================================")
    println("   Verification: VEFIE Direct RCS vs Legacy")
    println("==================================================")

    # 1. Parameters
    freq = 300e6
    EMMoMSuite.Utilities.Parameters.set_frequency!(freq)
    
    # 2. Mesh
    mesh_file = joinpath(@__DIR__, "../../../deps/fixtures/AllinOne/meshfiles/Tetra.nas")
    # Legacy uses meshUnit=:mm, which scales by 1e-3.
    mesh = read_nas_mesh(mesh_file, scale=1e-3)
    
    # 3. Basis
    basis = SWGBasis(mesh)
    N = num_basis(basis)
    println("Unknowns: $N")

    # 4. Permittivity
    eps_r = 2.0 + 0.0im
    permittivities = fill(eps_r, num_elements(mesh))

    # 5. Operator
    vefie = VEFIE(freq, permittivities)
    
    # 6. Assembly
    println("Assembling Matrix...")
    Z = assemble_impedance_matrix(vefie, basis, permittivities)
    
    # Apply Legacy Scaling Factor (1/16pi)
    # Legacy code seems to scale the impedance matrix by 1/16pi relative to standard units.
    # This results in currents that are 16pi times larger, and RCS that is (16pi)^2 ~ 34dB higher.
    # To match Legacy, we apply this factor.
    # Z .*= (1.0 / (16π))
    
    println("EMMoMSuite Z norm: ", norm(Z))
    max_val, max_idx = findmax(abs.(Z))
    println("Max Z element: ", max_val, " at ", max_idx)
    println("Z[max_idx]: ", Z[max_idx])

    # 7. Excitation
    println("Excitation...")
    # Legacy: PlaneWave(pi, 0, 0, 1) -> Incident from +z, x-pol.
    # EMMoMSuite: PlaneWave(freq, pi, 0, [1,0,0]) -> Incident from +z, x-pol.
    source = PlaneWave(freq, π, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(vefie, source, basis, permittivities)
    println("EMMoMSuite V norm: ", norm(V))

    # 8. Solve
    println("Solving...")
    D = Z \ V
    
    # 9. Calculate RCS
    println("Calculating RCS...")
    theta_range = collect(0:1.0:180.0) .* (π/180.0)
    phi_range = [0.0]
    
    _, _, rcs_db_matrix = radarCrossSection(theta_range, phi_range, D, basis, permittivities)
    
    # Load Legacy Data
    legacy_file = joinpath(@__DIR__, "legacy_vefie_rcs.txt")
    if isfile(legacy_file)
        legacy_data = readdlm(legacy_file, skipstart=1)
        legacy_theta = legacy_data[:, 1]
        legacy_rcs = legacy_data[:, 2]
        
        println("\nComparison:")
        println("Theta | EMMoMSuite | Legacy | Diff")
        for i in 1:10:181
            t = legacy_theta[i]
            ems = rcs_db_matrix[i, 1]
            leg = legacy_rcs[i]
            diff = abs(ems - leg)
            @printf("%5.1f | %7.2f | %7.2f | %5.2f\n", t, ems, leg, diff)
        end
    else
        println("Legacy data not found at $legacy_file")
    end
end

verify_vefie_legacy()
