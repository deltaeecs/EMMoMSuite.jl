"""
    EMSuite

A comprehensive Computational Electromagnetics (CEM) library in Julia.

# Overview
EMSuite provides a modular framework for solving electromagnetic scattering and radiation problems using the Method of Moments (MoM). It supports surface and volume integral equations, fast algorithms (MLFMA), and parallel computing.

# Key Modules
- **Geometry**: Mesh handling (Triangle, Tetrahedra) and geometric queries.
- **BasisFunctions**: RWG (surface), SWG (volume), and other basis functions.
- **IntegralEquations**: EFIE, MFIE, CFIE, and VIE formulations.
- **Solvers**: Direct (LU) and Iterative (GMRES, BiCGSTAB) solvers.
- **FastAlgorithms**: Multilevel Fast Multipole Algorithm (MLFMA).
- **PostProcessing**: RCS, Near-field, and Far-field calculations.
- **Parallel**: MPI and Threading support.

# Quick Start
```julia
using EMSuite

# 1. Load Mesh
mesh = read_nas_mesh("plate.nas")

# 2. Setup Basis
basis = RWGBasis(mesh)

# 3. Define Operator
freq = 300e6
efie = EFIE(freq)

# 4. Assemble Matrix
Z = assemble_impedance_matrix(efie, basis)

# 5. Solve
V = excitation_vector(PlaneWave(freq), basis)
I = solve!(GMRESSolver(), Z, V)

# 6. Post-Process
rcs = radarCrossSection(theta, phi, I, mesh, RWG)
```
"""
module EMSuite

using LinearAlgebra
using SparseArrays
using StaticArrays
using Logging

# Export Core interfaces
export AbstractMesh, AbstractBasisFunction, AbstractIntegralOperator, AbstractSolver, AbstractSource

# Include submodules
include("Core/Core.jl")
include("Utilities/Utilities.jl")
include("Geometry/Geometry.jl")
include("BasisFunctions/BasisFunctions.jl")
include("IntegralEquations/IntegralEquations.jl")
include("Solvers/Solvers.jl")
include("PostProcessing/PostProcessing.jl")
include("IO/IO.jl")

# Re-export Core symbols
using .CoreModule
export AbstractMesh, vertices, elements, dimension, num_vertices, num_elements
export AbstractBasisFunction, num_basis, support, evaluate
export AbstractIntegralOperator, kernel, impedance_element
export AbstractSolver, solve!
export AbstractSource, incident_field, excitation_vector, PlaneWave, DeltaGapSource
export EMSuiteConfig, load_config
export SimulationConfig, GeometryConfig, BasisConfig, ExcitationConfig, SolverConfig, OutputConfig

# Re-export Utilities symbols
using .Utilities
export init_logging,
    @showprogress, SimulationParameters, set_frequency!, get_k0, get_eta0, get_omega
export calculate_mie_rcs_pec_sphere

# Re-export Geometry symbols
using .Geometry
export TriangleMesh,
    TetrahedraMesh,
    HexahedraMesh,
    TriangleInfo,
    read_nas_mesh,
    read_mixed_nas_mesh,
    read_msh_mesh,
    generate_rectangle_mesh,
    generate_cylinder_mesh,
    generate_sphere_mesh

# Re-export PostProcessing symbols
using .PostProcessing
export radarCrossSection, farField, geoElectricJCal, calculate_near_field

# Re-export IO symbols
using .IO
export save_RCS_txt, save_results_hdf5, save_result, save_vtk


# Re-export BasisFunctions symbols
using .BasisFunctions
export RWGBasis,
    RWG,
    SWGBasis,
    SWG,
    evaluate_swg,
    PWCBasis,
    PWCHexBasis,
    PWC,
    RBFBasis,
    RBF,
    get_triangle_info,
    get_tetrahedra_info,
    get_hexahedra_info

# Re-export IntegralEquations symbols
using .IntegralEquations
export EFIE, MFIE, CFIE, VEFIE, SCFIE, assemble_impedance_matrix

# Re-export Solvers symbols
using .Solvers
export LUSolver, GMRESSolver, BiCGSTABSolver, solve!
export DiagonalPreconditioner,
    IdentityPreconditioner, ILUPreconditioner, SPAIPreconditioner, BlockJacobiPreconditioner

# Re-export FastAlgorithms symbols
include("FastAlgorithms/FastAlgorithms.jl")
using .FastAlgorithms
export FastAlgorithms, MLFMA, MLFMAOperator, get_leaf_intervals, MLFMAOperatorMPI

# Re-export Parallel symbols
include("Parallel/Parallel.jl")
using .Parallel
export init_parallel!, mpi_rank, nprocs, barrier, @root, num_threads, thread_id
export assemble_impedance_matrix_parallel

# Re-export Driver
include("Driver.jl")
using .Driver
export run_simulation

end # module
