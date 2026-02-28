module Parameters

export SimulationParameters, set_frequency!, get_k0, get_eta0, get_omega

"""
    SimulationParameters

Global container for simulation parameters.

# Fields
- `frequency`: Operating frequency \$f\$ (Hz).
- `omega`: Angular frequency \$\\omega = 2\\pi f\$ (rad/s).
- `k0`: Free-space wavenumber \$k_0 = \\omega/c_0\$ (rad/m).
- `eta0`: Free-space intrinsic impedance \$\\eta_0 = \\sqrt{\\mu_0/\\epsilon_0}\$ (\$\\Omega\$).
"""
mutable struct SimulationParameters
    frequency::Float64
    omega::Float64
    k0::Float64
    eta0::Float64

    function SimulationParameters()
        new(0.0, 0.0, 0.0, 376.73031346177)
    end
end

const GLOBAL_PARAMS = SimulationParameters()

function set_frequency!(freq::Float64)
    GLOBAL_PARAMS.frequency = freq
    GLOBAL_PARAMS.omega = 2 * pi * freq
    c0 = 299792458.0
    GLOBAL_PARAMS.k0 = GLOBAL_PARAMS.omega / c0
    GLOBAL_PARAMS.eta0 = 376.73031346177
end

get_k0() = GLOBAL_PARAMS.k0
get_eta0() = GLOBAL_PARAMS.eta0
get_omega() = GLOBAL_PARAMS.omega

end
