module ResultIO

using DelimitedFiles
using HDF5
using TOML
using ...CoreModule: SimulationResult

export save_RCS_txt, save_results_hdf5, save_result

"""
    save_RCS_txt(filename, theta, phi, rcs_data)

Save RCS data to a text file.
Format: theta (deg) rows, phi (deg) columns (if 2D) or similar.
"""
function save_RCS_txt(filename::String, theta::AbstractVector, phi::AbstractVector, rcs_data::AbstractMatrix)
    open(filename, "w") do io
        # Write header
        write(io, "# RCS Data\n")
        write(io, "# Theta (deg): $(rad2deg.(theta))\n")
        write(io, "# Phi (deg): $(rad2deg.(phi))\n")
        
        # Write data
        writedlm(io, rcs_data)
    end
end

"""
    save_results_hdf5(filename, data)

Save results to an HDF5 file.
"""
function save_results_hdf5(filename::String, data::Dict)
    h5open(filename, "w") do file
        for (k, v) in data
            file[string(k)] = v
        end
    end
end

function save_results_hdf5(filename::String; kwargs...)
    h5open(filename, "w") do file
        for (k, v) in kwargs
            file[string(k)] = v
        end
    end
end

"""
    save_result(result::SimulationResult)

Save the full simulation result to disk using the path specified in the configuration.
"""
function save_result(result::SimulationResult)
    # Determine output path
    out_dir = result.config.output.directory
    if !isdir(out_dir)
        mkpath(out_dir)
    end
    
    filename = joinpath(out_dir, result.config.simulation.name * "_result.h5")
    
    h5open(filename, "w") do file
        # Save currents
        file["currents"] = result.currents
        
        # Save metrics
        g_metrics = create_group(file, "metrics")
        for (k, v) in result.metrics
            g_metrics[string(k)] = v
        end
        
        # Save config as TOML string
        io = IOBuffer()
        # We need a way to serialize config to TOML. 
        # Since config is a struct, we might need to convert it to Dict first or use a serializer.
        # For now, let's just save the raw config file content if we had it, but we have the struct.
        # Let's try to serialize the struct fields.
        # Or just save the important parts.
        
        # Simple serialization of config fields
        g_config = create_group(file, "config")
        g_config["simulation_name"] = result.config.simulation.name
        g_config["frequency"] = result.config.simulation.frequency
        # Add more fields as needed
    end
    
    return filename
end

end
