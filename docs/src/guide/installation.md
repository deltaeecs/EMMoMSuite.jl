# Installation Guide

## Prerequisites

*   **Julia**: Version 1.9 or higher is recommended. Download from [julialang.org](https://julialang.org/downloads/).
*   **MPI** (Optional): For parallel computing features. Ensure an MPI implementation (e.g., MPICH, OpenMPI, or Microsoft MPI on Windows) is installed and available in your system path.

## Installation Steps

### Method 1: From GitHub (Recommended for Users)

Open the Julia REPL and enter the package manager by pressing `]`. Then run:

```julia
pkg> add https://github.com/yourusername/EMSuite.jl
```

### Method 2: Development Mode (Recommended for Contributors)

If you plan to modify the source code:

1.  Clone the repository:
    ```bash
    git clone https://github.com/yourusername/EMSuite.jl.git
    cd EMSuite.jl
    ```

2.  Start Julia in the project directory:
    ```bash
    julia --project=.
    ```

3.  Instantiate dependencies:
    ```julia
    using Pkg
    Pkg.instantiate()
    ```

## Verifying Installation

To verify that the package is installed correctly, run the test suite:

```julia
pkg> test EMSuite
```

