# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added
- **Core**: Initial project structure and core interfaces (`AbstractMesh`, `AbstractBasisFunction`, etc.).
- **Geometry**: Support for Triangle, Tetrahedra, and Hexahedra meshes. Mesh I/O for Nastran (.nas) and Gmsh (.msh).
- **Basis Functions**: Implementation of RWG, SWG, RBF, and PWC basis functions.
- **Integral Equations**: EFIE, MFIE, and CFIE formulations.
- **Solvers**: Direct (LU) and Iterative (GMRES, BiCGSTAB) solvers.
- **Fast Algorithms**: Multilevel Fast Multipole Algorithm (MLFMA) with Lebedev integration.
- **Parallel**: MPI support for distributed matrix assembly and solving.
- **Post-Processing**: RCS, Near Field, and Far Field calculations. VTK export.
- **Documentation**: Comprehensive documentation using Documenter.jl.

### Changed
- Refactored legacy code (`MoM_Basics`, `MoM_Kernels`, `MoM_MPI`) into a unified `EMSuite` package.
- Improved MLFMA accuracy and calibrated factors to match legacy behavior.
- Standardized I/O and logging.

### Fixed
- Fixed singular term scaling in EFIE.
- Fixed edge sorting logic in RWG basis functions.
- Verified parallel consistency with serial execution.
