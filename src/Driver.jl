module Driver

using ..EMMoMSuite
using ..EMMoMSuite.CoreModule: PlaneWave, SimulationResult
using Logging

export run_simulation

"""
    run_simulation(config_path::String)

从配置文件启动完整 MoM 仿真流程：
加载配置 → 初始化日志 → 读取几何网格 → 构造基函数 → 装配阻抗矩阵 →
施加激励 → 求解 → 后处理（RCS/方向图/近场）→ 保存结果。

返回 `SimulationResult`，其中包含仿真配置与求解结果。
"""
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

    # Apply scaling if needed
    if config.geometry.unit_scale != 1.0
        @warn "Mesh scaling not yet implemented. Ignoring unit_scale."
    end

    # 4. Setup Basis
    @info "Setting up basis functions: $(config.basis.type)"
    if config.basis.type == "RWG"
        basis = RWGBasis(mesh)
    elseif config.basis.type == "SWG"
        basis = SWGBasis(mesh)
    elseif config.basis.type == "PWC"
        basis = PWCBasis(mesh)
    else
        error("Unsupported basis type: $(config.basis.type)")
    end

    # 5. Define Operator & Solve
    ie_type = config.simulation.ie_type
    @info "Integral equation type: $ie_type"

    # Setup operator
    freq = config.simulation.frequency

    if ie_type == "EFIE"
        @info "Assembling EFIE system matrix..."
        op = EFIE(freq)
        Z = assemble_impedance_matrix(op, basis)
    elseif ie_type == "MFIE"
        @info "Assembling MFIE system matrix..."
        op = MFIE(freq)
        Z = assemble_impedance_matrix(op, basis)
    elseif ie_type == "CFIE"
        alpha = config.simulation.cfie_alpha
        @info "Assembling CFIE system matrix (α=$alpha)..."
        op = CFIE(freq, alpha)
        Z = assemble_impedance_matrix(op, basis)
    elseif ie_type == "VEFIE"
        perms = config.simulation.permittivities
        @info "Assembling VEFIE system matrix..."
        op = VEFIE(freq, perms)
        Z = assemble_impedance_matrix(op, basis)
    elseif ie_type == "SCFIE"
        # Surface-Volume coupled IE: requires both surface and volume meshes
        error("SCFIE requires explicit surface+volume mesh setup via scripting API.")
    else
        error("Unsupported integral equation type: $ie_type")
    end

    # 6. Excitation
    @info "Setting up excitation..."
    if config.excitation.type == "PlaneWave"
        source = PlaneWave(
            freq,
            config.excitation.theta,
            config.excitation.phi,
            config.excitation.polarization,
        )
    else
        error("Unsupported excitation type: $(config.excitation.type)")
    end

    if ie_type == "VEFIE" && basis isa SWGBasis
        V = excitation_vector(op, source, basis, config.simulation.permittivities)
    elseif ie_type == "VEFIE" && basis isa PWCBasis
        V = excitation_vector(op, source, basis)
    else
        V = excitation_vector(source, basis)
    end

    # 7. Solve
    @info "Solving linear system using $(config.simulation.solver_type)..."
    if config.simulation.solver_type == "GMRES"
        solver = GMRESSolver(
            tol = config.solver.tolerance,
            maxiter = config.solver.max_iter,
            restart = config.solver.restart,
        )
    elseif config.simulation.solver_type == "BiCGSTAB"
        solver = BiCGSTABSolver(tol = config.solver.tolerance, maxiter = config.solver.max_iter)
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
