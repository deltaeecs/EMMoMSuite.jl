# EMSuite.jl

Welcome to the documentation for **EMSuite.jl**, a comprehensive Computational Electromagnetics (CEM) library written in Julia.

## Overview

EMSuite provides a modular and high-performance framework for solving electromagnetic scattering and radiation problems using the **Method of Moments (MoM)**. It is designed for researchers and engineers who need a flexible platform for algorithm development and large-scale simulation.

## Key Features

*   **Versatile Geometry Support**: Handles triangular (surface) and tetrahedral (volume) meshes.
*   **Advanced Integral Equations**: Supports EFIE, MFIE, CFIE, and VIE formulations.
*   **Fast Algorithms**: Implements Multilevel Fast Multipole Algorithm (MLFMA) for accelerating large-scale problems.
*   **Parallel Computing**: Built-in support for MPI and multi-threading to leverage modern clusters.
*   **Rich Basis Functions**: Includes RWG, SWG, PWC, and Rooftop basis functions.
*   **Robust Solvers**: Integrated direct (LU) and iterative (GMRES, BiCGSTAB) solvers with preconditioning.

## Documentation Structure

*   **[User Guide](guide/quick_start.md)**: Step-by-step instructions for installation, quick start, and running examples.
*   **[Algorithm Theory](theory/electromagnetics.md)**: Detailed mathematical derivations of the underlying physics and numerical methods.
*   **[API Reference](api/public_api.md)**: Comprehensive documentation of types and functions.

## Getting Started

To install the package, simply run:

```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/EMSuite.jl")
```

See the [Installation Guide](guide/installation.md) for more details.
