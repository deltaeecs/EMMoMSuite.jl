module Results

using ..Configuration

export SimulationResult

"""
    SimulationResult

Container for simulation results.
"""
struct SimulationResult
    config::EMMoMSuiteConfig
    currents::Vector{Complex{Float64}}
    metrics::Dict{String,Any}

    function SimulationResult(config::EMMoMSuiteConfig, currents::Vector{Complex{Float64}})
        new(config, currents, Dict{String,Any}())
    end
end

end
