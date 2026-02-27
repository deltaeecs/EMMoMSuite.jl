# Quick Start

This guide will walk you through a simple simulation: calculating the Radar Cross Section (RCS) of a Perfect Electric Conductor (PEC) plate.

## 1. Prepare the Mesh

First, you need a mesh file. EMSuite supports Nastran (`.nas`) and Gmsh (`.msh`) formats.
For this example, assume we have a file named `plate.nas`.

## 2. Write the Simulation Script

Create a file named `run_plate.jl`:

```julia
using EMSuite
using StaticArrays

# 1. Initialize Simulation Parameters
freq = 300e6 # 300 MHz
set_frequency!(freq)

# 2. Load Geometry
# Ensure 'plate.nas' exists in your working directory
mesh = read_nas_mesh("plate.nas")
println("Mesh loaded: $(num_elements(mesh)) triangles")

# 3. Setup Basis Functions
# Use RWG basis functions for surface currents
basis = RWGBasis(mesh)
println("Basis functions: $(num_basis(basis)) unknowns")

# 4. Define Integral Equation Operator
# Electric Field Integral Equation (EFIE)
operator = EFIE(freq)

# 5. Assemble Impedance Matrix
println("Assembling matrix...")
Z = assemble_impedance_matrix(operator, basis)

# 6. Define Excitation
# Plane wave incident from theta=0, phi=0 (z-direction)
# Polarization along x-axis
inc_wave = PlaneWave(theta=0.0, phi=0.0, pol=SVector(1.0, 0.0, 0.0))
V = excitation_vector(inc_wave, basis)

# 7. Solve Linear System
println("Solving system...")
solver = GMRESSolver(tol=1e-4, maxiter=100)
I_coeff = solve!(solver, Z, V)

# 8. Post-Processing (RCS)
println("Calculating RCS...")
theta_obs = collect(0:0.01:pi) # Observation angles (0 to 180 degrees)
phi_obs = [0.0]                # Phi = 0 plane

rcs_data, rcs_total, rcs_db = radarCrossSection(theta_obs, phi_obs, I_coeff, mesh, RWG)

# 9. Save Results
save_results_hdf5("plate_results.h5", I_coeff, rcs_db)
println("Done!")
```

## 3. Run the Simulation

Execute the script from the terminal:

```bash
julia --project=. run_plate.jl
```

## 4. Visualize Results

You can load the saved `plate_results.h5` file for plotting, or use the built-in VTK export to visualize currents:

```julia
# Add this to your script to export currents for ParaView
save_vtk("plate_currents", mesh, abs.(I_coeff))
```

