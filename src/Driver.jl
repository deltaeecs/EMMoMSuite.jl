module Driver

using ..EMSuite
using ..EMSuite.CoreModule: PlaneWave, SimulationResult
using Logging

export run_simulation

function run_simulation(config_path::String)
    # 1. Load Configuration
    config = load_config(config_path)
    
    # 2. Setup Logging
    log_file = joinpath(config.output.directory, config.simulation.name * ".log")
    init_logging(log_file)
    
    @info "Starting simulation: $(config.simulation.name)"
    @info "Frequency: $(config.simulation.frequency) Hz"
    
    # 3. Load Geometry
    @info "Loading mesh from: $(config.geometry.mesh_file)"
    
    if !isfile(config.geometry.mesh_file)
        error("Mesh file not found: $(config.geometry.mesh_file)")
    end

    if endswith(config.geometry.mesh_file, ".nas")
        mesh = read_nas_mesh(config.geometry.mesh_file)
    elseif endswith(config.geometry.mesh_file, ".msh")
        mesh = read_msh_mesh(config.geometry.mesh_file)
    else
        error("Unsupported mesh format: $(config.geometry.mesh_file)")
    end
    
    # Apply scaling if needed (not implemented in Mesh types yet, so we skip or assume 1.0)
    if config.geometry.unit_scale != 1.0
        @warn "Mesh scaling not yet implemented. Ignoring unit_scale."
    end
    
    # 4. Setup Basis
    @info "Setting up basis functions: $(config.basis.type)"
    if config.basis.type == "RWG"
        basis = RWGBasis(mesh)
    elseif config.basis.type == "SWG"
        basis = SWGBasis(mesh)
    else
        error("Unsupported basis type: $(config.basis.type)")
    end
    
    # 5. Define Operator & Solver
    @info "Assembling system matrix..."
    # For now, assume EFIE
    efie = EFIE(config.simulation.frequency)
    
    # Matrix Assembly
    Z = assemble_impedance_matrix(efie, basis)
    
    # 6. Excitation
    @info "Setting up excitation..."
    if config.excitation.type == "PlaneWave"
        source = PlaneWave(config.simulation.frequency, config.excitation.theta, config.excitation.phi, config.excitation.polarization)
    else
        error("Unsupported excitation type: $(config.excitation.type)")
    end
    
    V = excitation_vector(source, basis)
    
    # 7. Solve
    @info "Solving linear system using $(config.simulation.solver_type)..."
    if config.simulation.solver_type == "GMRES"
        solver = GMRESSolver(tol=config.solver.tolerance, maxiter=config.solver.max_iter, restart=config.solver.restart)
    elseif config.simulation.solver_type == "BiCGSTAB"
        solver = BiCGSTABSolver(tol=config.solver.tolerance, maxiter=config.solver.max_iter)
    elseif config.simulation.solver_type == "LU"
        solver = LUSolver()
    else
        error("Unsupported solver type: $(config.simulation.solver_type)")
    end
    
    I = solve!(solver, Z, V)
    
    # 8. Save Results
    @info "Saving results..."
    result = SimulationResult(config, I)
    saved_path = save_result(result)
    
    @info "Simulation completed successfully. Results saved to: $saved_path"
    return result
end

end
