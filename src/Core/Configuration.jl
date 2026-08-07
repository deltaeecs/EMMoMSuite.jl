module Configuration

using TOML

export EMMoMSuiteConfig, load_config
export SimulationConfig, GeometryConfig, BasisConfig, ExcitationConfig, SolverConfig, OutputConfig

# --- Configuration Structs ---

"""
    SimulationConfig

General simulation settings.
"""
Base.@kwdef struct SimulationConfig
    name::String = "Simulation"
    frequency::Float64
    solver_type::String = "GMRES"
    ie_type::String = "EFIE"
    cfie_alpha::Float64 = 0.5
    permittivities::Vector{ComplexF64} = ComplexF64[]
end

"""
    GeometryConfig

Geometry and mesh settings.
"""
Base.@kwdef struct GeometryConfig
    mesh_file::String
    unit_scale::Float64 = 1.0
end

"""
    BasisConfig

Basis function settings.
"""
Base.@kwdef struct BasisConfig
    type::String = "RWG"
end

"""
    ExcitationConfig

Excitation source settings.
"""
Base.@kwdef struct ExcitationConfig
    type::String = "PlaneWave"
    theta::Float64 = 0.0
    phi::Float64 = 0.0
    polarization::Vector{Float64} = [1.0, 0.0, 0.0]
end

"""
    SolverConfig

Linear solver settings.
"""
Base.@kwdef struct SolverConfig
    tolerance::Float64 = 1e-4
    max_iter::Int = 1000
    restart::Int = 50
end

"""
    OutputConfig

Output and result saving settings.
"""
Base.@kwdef struct OutputConfig
    directory::String = "results"
    save_fields::Bool = true
    save_rcs::Bool = true
    format::String = "h5"
end

"""
    EMMoMSuiteConfig

Top-level configuration struct.
"""
Base.@kwdef struct EMMoMSuiteConfig
    simulation::SimulationConfig
    geometry::GeometryConfig
    basis::BasisConfig
    excitation::ExcitationConfig
    solver::SolverConfig
    output::OutputConfig
end

# --- Loading Logic ---

"""
    load_config(path::String)

Load configuration from a TOML file.
"""
function load_config(path::String)
    if !isfile(path)
        error("Configuration file not found: $path")
    end

    data = TOML.parsefile(path)

    # Helper to safely extract sections
    function get_section(dict, key, default = Dict())
        return get(dict, key, default)
    end

    # Parse Simulation
    sim_data = get_section(data, "simulation")
    simulation = SimulationConfig(;
        name = get(sim_data, "name", "Simulation"),
        frequency = Float64(sim_data["frequency"]),
        solver_type = get(sim_data, "solver_type", "GMRES"),
        ie_type = get(sim_data, "ie_type", "EFIE"),
        cfie_alpha = Float64(get(sim_data, "cfie_alpha", 0.5)),
        permittivities = Vector{ComplexF64}(get(sim_data, "permittivities", ComplexF64[])),
    )

    # Parse Geometry
    geo_data = get_section(data, "geometry")
    geometry = GeometryConfig(;
        mesh_file = geo_data["mesh_file"],
        unit_scale = get(geo_data, "unit_scale", 1.0),
    )

    # Parse Basis
    basis_data = get_section(data, "basis")
    basis = BasisConfig(; type = get(basis_data, "type", "RWG"))

    # Parse Excitation
    exc_data = get_section(data, "excitation")
    excitation = ExcitationConfig(;
        type = get(exc_data, "type", "PlaneWave"),
        theta = get(exc_data, "theta", 0.0),
        phi = get(exc_data, "phi", 0.0),
        polarization = Vector{Float64}(get(exc_data, "polarization", [1.0, 0.0, 0.0])),
    )

    # Parse Solver
    sol_data = get_section(data, "solver")
    solver = SolverConfig(;
        tolerance = get(sol_data, "tolerance", 1e-4),
        max_iter = get(sol_data, "max_iter", 1000),
        restart = get(sol_data, "restart", 50),
    )

    # Parse Output
    out_data = get_section(data, "output")
    output = OutputConfig(;
        directory = get(out_data, "directory", "results"),
        save_fields = get(out_data, "save_fields", true),
        save_rcs = get(out_data, "save_rcs", true),
        format = get(out_data, "format", "h5"),
    )

    return EMMoMSuiteConfig(simulation, geometry, basis, excitation, solver, output)
end

end
