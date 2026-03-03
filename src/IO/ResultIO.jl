module ResultIO

using DelimitedFiles
using HDF5
using TOML
using ...CoreModule: SimulationResult

export save_RCS_txt, save_results_hdf5, save_result, save_RCS_csv

"""
    save_RCS_txt(filename, theta, phi, rcs_data)

Save RCS data to a text file.
Format: theta (deg) rows, phi (deg) columns (if 2D) or similar.
"""
function save_RCS_txt(
    filename::String,
    theta::AbstractVector,
    phi::AbstractVector,
    rcs_data::AbstractMatrix,
)
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
    save_RCS_csv(filename, θs, ϕs, rcs_comp, rcs_total, rcs_dB; linear=false)

Save RCS results to a CSV file suitable for comparison with CST, FEKO, etc.

# Arguments
- `filename`: Output file path (`.csv` extension appended if missing)
- `θs`: Elevation angle vector (radians)
- `ϕs`: Azimuth angle vector (radians)
- `rcs_comp`: Component RCS array of size `(2, Nθ, Nϕ)` (linear m²)
- `rcs_total`: Total RCS matrix `(Nθ, Nϕ)` (linear m²)
- `rcs_dB`: Total RCS in dBsm `(Nθ, Nϕ)`
- `linear=false`: When `false` (default) saves dBsm values; when `true` saves linear m²

# Format (1D, Nϕ==1)
`theta_deg, rcs_theta_dBsm, rcs_phi_dBsm, rcs_total_dBsm`

# Format (2D, Nϕ>1)
`theta_deg, phi_deg, rcs_theta_dBsm, rcs_phi_dBsm, rcs_total_dBsm`
"""
function save_RCS_csv(
    filename::String,
    θs::AbstractVector{<:Real},
    ϕs::AbstractVector{<:Real},
    rcs_comp::AbstractArray{<:Real,3},
    rcs_total::AbstractMatrix{<:Real},
    rcs_dB::AbstractMatrix{<:Real};
    linear::Bool = false,
)
    fn = endswith(filename, ".csv") ? filename : filename * ".csv"

    Nθ = length(θs)
    Nϕ = length(ϕs)

    open(fn, "w") do io
        if Nϕ == 1
            # 1D sweep: columns = theta_deg, θ-comp, φ-comp, total
            if linear
                write(io, "theta_deg,rcs_theta_m2,rcs_phi_m2,rcs_total_m2\n")
            else
                write(io, "theta_deg,rcs_theta_dBsm,rcs_phi_dBsm,rcs_total_dBsm\n")
            end
            for iθ in 1:Nθ
                θ_deg = rad2deg(θs[iθ])
                if linear
                    c1 = rcs_comp[1, iθ, 1]
                    c2 = rcs_comp[2, iθ, 1]
                    ct = rcs_total[iθ, 1]
                else
                    c1 = 10 * log10(max(rcs_comp[1, iθ, 1], 1e-30))
                    c2 = 10 * log10(max(rcs_comp[2, iθ, 1], 1e-30))
                    ct = rcs_dB[iθ, 1]
                end
                write(io, string(round(θ_deg,digits=6), ",", round(c1,digits=6), ",", round(c2,digits=6), ",", round(ct,digits=6), "\n"))
            end
        else
            # 2D grid: columns = theta_deg, phi_deg, θ-comp, φ-comp, total
            if linear
                write(io, "theta_deg,phi_deg,rcs_theta_m2,rcs_phi_m2,rcs_total_m2\n")
            else
                write(io, "theta_deg,phi_deg,rcs_theta_dBsm,rcs_phi_dBsm,rcs_total_dBsm\n")
            end
            for iθ in 1:Nθ, iϕ in 1:Nϕ
                θ_deg = rad2deg(θs[iθ])
                ϕ_deg = rad2deg(ϕs[iϕ])
                if linear
                    c1 = rcs_comp[1, iθ, iϕ]
                    c2 = rcs_comp[2, iθ, iϕ]
                    ct = rcs_total[iθ, iϕ]
                else
                    c1 = 10 * log10(max(rcs_comp[1, iθ, iϕ], 1e-30))
                    c2 = 10 * log10(max(rcs_comp[2, iθ, iϕ], 1e-30))
                    ct = rcs_dB[iθ, iϕ]
                end
                write(io, string(round(θ_deg,digits=6), ",", round(ϕ_deg,digits=6), ",", round(c1,digits=6), ",", round(c2,digits=6), ",", round(ct,digits=6), "\n"))
            end
        end
    end
    return fn
end


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
