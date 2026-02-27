using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Solvers
using LinearAlgebra

# Load Mesh
mesh_file = "f:/OneDrive/MoM/MoM_AllinOne/meshfiles/jet_100MHz.nas"
if !isfile(mesh_file)
    error("Mesh file not found: $mesh_file")
end

println("Loading mesh from: $mesh_file")
mesh = read_nas_mesh(mesh_file)
println("Mesh loaded: $(length(mesh.triangles)) triangles")

# Basis
freq = 3.0e8
basis = RWGBasis(mesh)
println("Basis functions: $(length(basis.functions))")

# Assembly
println("Starting Assembly with $(Threads.nthreads()) threads...")
t_start = time()
Z = assemble_impedance_matrix(EFIE(freq), basis)
t_end = time()
println("Assembly Time: $(t_end - t_start) s")
