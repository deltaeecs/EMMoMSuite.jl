module CoreModule

using LinearAlgebra

include("Interfaces.jl")
include("Types.jl")
include("Constants.jl")
include("Materials.jl")

export AbstractMesh, vertices, elements, dimension, num_vertices, num_elements
export AbstractBasisFunction, num_basis, support, evaluate
export AbstractIntegralOperator, kernel, impedance_element, assemble_impedance_matrix
export AbstractSolver, solve!
export AbstractSource, incident_field, excitation_vector

export Constants

# Re-export Materials
using .Materials
export AbstractMaterial, PEC, PMC, Dielectric, ImpedanceSurface
export permittivity, permeability, impedance

# Re-export Configuration
include("Configuration.jl")
using .Configuration
export EMMoMSuiteConfig, load_config
export SimulationConfig, GeometryConfig, BasisConfig, ExcitationConfig, SolverConfig, OutputConfig

# Re-export Results
include("Results.jl")
using .Results
export SimulationResult

# Re-export Sources
include("Sources.jl")
using .Sources
export PlaneWave, DeltaGapSource, incident_field

end
