# Diagnose Z matrix element differences between EMSuite and Legacy formula
using EMSuite
using LinearAlgebra
using Printf

const MOM_ALLINONE_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne")

println("=" ^ 60)
println("Z Matrix Element Diagnosis")
println("=" ^ 60)

# Load mesh and build basis
mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
mesh = read_nas_mesh(mesh_file, scale=1.0)
freq = 1e8
set_frequency!(freq)
basis = RWGBasis(mesh)
N = num_basis(basis)
println("Unknowns: $N")

# Build EFIE and Z matrix
efie = EFIE(freq)
Z = assemble_impedance_matrix(efie, basis)

# Print diagonal and a few off-diagonal elements
println("\n--- First 10 diagonal elements ---")
for i in 1:min(10, N)
    @printf("Z[%d,%d] = %+.6e %+.6ej  (|Z|=%.4e)\n", i, i, real(Z[i,i]), imag(Z[i,i]), abs(Z[i,i]))
end

# Find statistics
diag_vals = [abs(Z[i,i]) for i in 1:N]
println("\nDiag |Z| stats: min=$(minimum(diag_vals)), max=$(maximum(diag_vals)), mean=$(sum(diag_vals)/N)")

# Print some off-diagonal elements (near the diagonal)
println("\n--- Near-diagonal elements (row 1) ---")
for j in 1:min(20, N)
    if j != 1 && abs(Z[1, j]) > 0
        @printf("Z[1,%d] = %+.6e %+.6ej  (|Z|=%.4e)\n", j, real(Z[1,j]), imag(Z[1,j]), abs(Z[1,j]))
    end
end

# Excitation vector
source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
V = excitation_vector(efie, source, basis)

println("\n--- First 10 excitation vector elements ---")
for i in 1:min(10, N)
    @printf("V[%d] = %+.6e %+.6ej  (|V|=%.4e)\n", i, real(V[i]), imag(V[i]), abs(V[i]))
end

# Solution coefficients
I_coeff = Z \ V
println("\n--- First 10 solution coefficients ---")
for i in 1:min(10, N)
    @printf("I[%d] = %+.6e %+.6ej  (|I|=%.4e)\n", i, real(I_coeff[i]), imag(I_coeff[i]), abs(I_coeff[i]))
end

# Check near-interaction vs far-interaction contribution for bf=1
bf1 = basis.functions[1]
println("\n--- Basis function 1 info ---")
println("  support: $(bf1.support)")
println("  edge_length: $(bf1.edge_length)")
println("  signs: $(bf1.signs)")

# Print the support triangle info
for k in 1:2
    t = bf1.support[k]
    tri = EMSuite.IntegralEquations.EFIEModule.Impedance.Geometry.TriangleInfo(mesh, t)
    println("  T$(k) (id=$t): area=$(tri.area), edgel=$(tri.edgel)")
end

# Print the triangle areas statistics
println("\n--- Triangle area statistics ---")
ntri = EMSuite.CoreModule.num_elements(mesh)
areas = Float64[]
for t in 1:ntri
    tri = EMSuite.IntegralEquations.EFIEModule.Impedance.Geometry.TriangleInfo(mesh, t)
    push!(areas, tri.area)
end
println("  ntri: $ntri")
println("  min area: $(minimum(areas))")
println("  max area: $(maximum(areas))")
println("  mean area: $(sum(areas)/length(areas))")
println("  1/mean_area: $(1/(sum(areas)/length(areas)))")
println("  1.25 * mean_area: $(1.25 * sum(areas)/length(areas))")
