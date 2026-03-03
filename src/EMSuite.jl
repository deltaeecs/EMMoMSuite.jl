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

# 1. 设置频率（更新全局 k0/η0）
freq = 300e6
set_frequency!(freq)

# 2. 加载网格与基函数
mesh  = read_nas_mesh("plate.nas")
basis = RWGBasis(mesh)

# 3. 定义算子 + 装配矩阵
efie = EFIE(freq)
Z    = assemble_impedance_matrix(efie, basis)

# 4. 激励（+z 方向入射，x 极化）
V = excitation_vector(PlaneWave(freq, π/2, π, [1.0, 0.0, 0.0]), basis)

# 5. 求解
I = solve!(LUSolver(), Z, V)

# 6. RCS
θ = collect(range(0.0, π, length=181))
_, _, rcs_db = radarCrossSection(θ, [0.0], I, basis)
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
include("Materials/Materials.jl")
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
export SimulationResult

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
    write_nas_mesh,
    read_msh_mesh,
    generate_rectangle_mesh,
    generate_cylinder_mesh,
    generate_sphere_mesh,
    generate_ellipsoid_mesh,
    generate_cone_mesh,
    generate_torus_mesh,
    generate_box_volume_mesh,
    generate_box_tet_mesh,
    extract_surface,
    translate_mesh, scale_mesh, rotate_mesh, transform_mesh, merge_meshes,
    MeshQualityReport, mesh_quality,
    remove_duplicate_nodes, fix_element_orientation, detect_degenerates,
    read_stl_mesh, write_stl_mesh,
    tet_volume, area,
    BRepFace, BRepSolid, CSGNode,
    box_solid, solid_volume, solid_surface_area, check_manifold, convert_to_triangle_mesh,
    intersect_solids, union_solids, subtract_solid, csg_volume,
    generate_gmsh_sphere, generate_gmsh_box, generate_gmsh_from_file,
    surface_mesh_gmsh, surface_mesh,
    mesh_face_labels, label_mesh_tags, propagate_labels,
    tet_mesh_gmsh, tet_mesh,
    read_hex_mesh, validate_mesh,
    BoundMesh, bind_materials, validate_bindings, element_material

# Re-export Materials symbols
using .MaterialsModule
export MaterialModel, Isotropic, Anisotropic
export MaterialEntry, MaterialLibrary
export add_material!, get_material
export save_library, load_library, load_builtin_library
export DebyePole, LorentzPole
export DebyeModel, DrudeModel, LorentzModel
export eval_permittivity

# Re-export PostProcessing symbols
using .PostProcessing
export radarCrossSection, farField, geoElectricJCal, geoVolumeCurrentCal, calculate_near_field
export field_cut_line, field_cut_plane
export antenna_directivity, input_impedance, beam_metrics
export absorbed_power, sar

# Re-export IO symbols
using .IO
export save_RCS_txt, save_RCS_csv, save_results_hdf5, save_result, save_vtk, save_vtk_multi


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
export mpi_gmres!, mpi_gmres
export MPIMatrix

# Re-export Driver
include("Driver.jl")
using .Driver
export run_simulation

end # module
