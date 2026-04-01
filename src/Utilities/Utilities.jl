module Utilities

include("LightweightSupport.jl")
using .LightweightSupport
export Progress, next!, @showprogress, find_zero_bisection, knn_bruteforce

include("Logging.jl")
using .LoggingUtils
export init_logging

include("Parameters.jl")
using .Parameters
export SimulationParameters, set_frequency!, get_k0, get_eta0, get_omega

include("HDF5Utils.jl")
using .HDF5Utils
export save_sparse_matrix, load_sparse_matrix

include("MieSeries.jl")
using .MieSeries
export calculate_mie_rcs_pec_sphere, calculate_mie_rcs_dielectric_sphere,
       calculate_mie_rcs_pec_sphere_fullpol

end
