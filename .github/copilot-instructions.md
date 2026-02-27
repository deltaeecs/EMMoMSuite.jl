# Refactoring Prompts for EMSuite.jl

This file contains a checklist of tasks for the refactoring process, based on `REFACTORING_PLAN.md`.
Use this as a guide for the AI assistant to track progress.

## Phase 1: Infrastructure Setup (Completed)

### 1.1 Project Initialization
- [x] Create new project directory `EMSuite`
- [x] Initialize `Project.toml` with basic dependencies
- [x] Set up Git repository (init)
- [x] Create `.gitignore`
- [x] Create `README.md` with project overview
- [x] Create `LICENSE` file
- [x] Create `.github/workflows/CI.yml` for CI/CD

### 1.2 Core Module Design
- [x] Create `src/EMSuite.jl` (Main entry point)
- [x] Create `src/Core/` directory
- [x] Implement `src/Core/Interfaces.jl` (Abstract types: `AbstractMesh`, `AbstractBasisFunction`, etc.)
- [x] Implement `src/Core/Types.jl` (Common types)
- [x] Implement `src/Core/Constants.jl` (Physical constants)
- [x] Create `src/Utilities/` directory
- [x] Implement `src/Utilities/Logging.jl` (Logging system)
- [x] Implement `src/Utilities/Parameters.jl` (Parameter management)

### 1.3 Documentation Framework
- [x] Create `docs/` directory structure
- [x] Set up `docs/make.jl` using Documenter.jl
- [x] Create `docs/src/index.md`
- [x] Create `docs/Project.toml`

## Phase 2: Geometry and Basis Functions (3-4 Weeks)

### 2.1 Geometry Module Refactoring
- [x] Create `src/Geometry/` directory
- [x] Migrate mesh types from `MoM_Basics` to `src/Geometry/MeshTypes.jl`
- [x] Migrate mesh I/O (Nastran, etc.) to `src/Geometry/MeshIO.jl`
- [x] Implement `src/Geometry/CoordinateTransforms.jl`
- [x] Implement `src/Geometry/GaussQuadrature.jl`
- [x] Add tests for Geometry module in `test/test_geometry.jl`

### 2.2 Basis Functions Module Refactoring
- [x] Create `src/BasisFunctions/` directory
- [x] Migrate RWG basis functions to `src/BasisFunctions/RWG.jl`
- [x] Migrate SWG basis functions to `src/BasisFunctions/SWG.jl`
- [x] Migrate other basis functions (RBF, PWC)
- [x] Implement `src/BasisFunctions/BasisUtilities.jl`
- [x] Add tests for BasisFunctions module in `test/test_basis_functions.jl`

## Phase 3: Integral Equations and Matrix Assembly (4-5 Weeks)

### 3.1 Integral Equations
- [x] Create `src/IntegralEquations/` directory
- [x] Implement `src/IntegralEquations/EFIE.jl`
- [x] Implement `src/IntegralEquations/MFIE.jl`
- [x] Implement `src/IntegralEquations/CFIE.jl`
- [x] Implement `src/IntegralEquations/Impedance.jl` (Matrix assembly)
- [x] Add tests in `test/test_integral_equations.jl`

### 3.2 Direct Solvers
- [x] Create `src/Solvers/` directory
- [x] Implement `src/Solvers/DirectSolvers.jl`
- [x] Add tests in `test/test_solvers.jl`

## Phase 4: MLFMA Fast Algorithm (5-6 Weeks)

### 4.1 MLFMA Core
- [x] Create `src/FastAlgorithms/MLFMA/` directory
- [x] Implement Octree structure in `src/FastAlgorithms/MLFMA/Octree.jl`
- [x] Implement Aggregation in `src/FastAlgorithms/MLFMA/Aggregation.jl`
- [x] Implement Translation in `src/FastAlgorithms/MLFMA/Translation.jl`
- [x] Implement Disaggregation in `src/FastAlgorithms/MLFMA/Disaggregation.jl`
- [x] Add tests in `test/test_mlfma.jl`

### 4.2 Lebedev Integration
- [x] Migrate Lebedev code to `src/FastAlgorithms/Lebedev/`
- [x] Integrate with MLFMA

## Phase 5: Solvers and Parallel Computing (4-5 Weeks)

### 5.1 Iterative Solvers
- [x] Integrate `IterativeSolvers.jl` functionality
- [x] Implement GMRES, BiCGSTAB wrappers in `src/Solvers/IterativeSolvers/`
- [x] Implement Preconditioners

### 5.2 Parallel Computing
- [x] Create `src/Parallel/` directory
- [x] Implement MPI support in `src/Parallel/MPI/`
- [x] Implement Threading support in `src/Parallel/Threading.jl`
- [x] Add tests in `test/test_parallel.jl`

## Phase 6: Post-Processing and Visualization (3-4 Weeks)

### 6.1 Post-Processing
- [x] Create `src/PostProcessing/` directory
- [x] Implement RCS calculation in `src/PostProcessing/RCS.jl`
- [x] Implement Near/Far field calculation
- [x] Add tests in `test/test_postprocessing.jl`

### 6.2 Visualization (Replaced with VTK Export)
- [x] Remove `src/Visualization/` (Decoupled visualization)
- [x] Implement VTK export in `src/IO/VTKExport.jl`

## Phase 7: I/O and Utilities (Completed)

### 7.1 I/O
- [x] Create `src/IO/` directory
- [x] Implement Result I/O (HDF5, CSV, VTK)

### 7.2 Utilities
- [x] Finalize Utility modules

## Phase 10: Workflow and Configuration (Restructuring) (New)

### 10.1 Configuration Management
- [x] Design `Configuration` structs in `src/Core/Configuration.jl`
- [x] Implement `load_config` with `TOML`
- [x] Integrate `Configuration` into `EMSuite` module

### 10.2 Logging & Diagnostics
- [x] Refactor `src/Utilities/Logging.jl` to support file output and structured logging
- [x] Integrate `ProgressMeter.jl` for long-running tasks

### 10.3 Data Management
- [x] Create `SimulationResult` struct in `src/Core/Types.jl` or `src/IO/Results.jl`
- [x] Implement `save_result` function (Migrated to HDF5)

### 10.4 Workflow Orchestration
- [x] Create `src/Driver.jl` or `src/App.jl`
- [x] Implement `run_simulation(config_path::String)`

## Phase 11: Dependency Optimization (Completed)

- [x] Replace `JLD2` with `HDF5` for language-agnostic data storage
- [x] Remove heavy dependencies (`Distributions`, `LaTeXStrings`)
- [x] Migrate interpolation weights to HDF5 format

## Phase 8: Integration and Testing (3-4 Weeks)

- [x] Run full integration tests
- [ ] Perform benchmarks
- [ ] Complete documentation

### 8.1 Comprehensive Verification
- [ ] Verify RCS with finer mesh (Sphere)
- [ ] Verify SWG basis functions (Dielectric)
- [ ] Verify PWC basis functions
- [ ] Verify RBF basis functions
- [ ] Verify Mixed Basis scenarios

## Phase 9: Release Preparation (2 Weeks)

- [ ] Finalize CHANGELOG and README
- [ ] Prepare for registration
