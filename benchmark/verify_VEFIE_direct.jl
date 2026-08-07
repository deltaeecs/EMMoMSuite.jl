using EMMoMSuite
using LinearAlgebra
using StaticArrays
using CSV, DataFrames
using Statistics

# Parameters
freq = 1.2e9 # 1.2 GHz
mesh_file = joinpath(@__DIR__, "..", "deps", "fixtures", "AllinOne", "meshfiles", "plate_1dot2GHz.nas")

# 1. Load Mesh
println("Loading mesh: $mesh_file")
mesh = read_nas_mesh(mesh_file)
println("Mesh loaded. Type: $(typeof(mesh))")
println("Vertices: $(size(mesh.node, 2))")
if mesh isa TetrahedraMesh
    println("Tetrahedra: $(size(mesh.tetras, 2))")
else
    error("Expected TetrahedraMesh")
end

# 2. Setup Basis
println("Setting up SWG basis...")
basis = SWGBasis(mesh)
n_unknowns = num_basis(basis)
println("Number of unknowns: $n_unknowns")

# 3. Setup VEFIE
println("Setting up VEFIE...")
vefie = VEFIE(freq)

# 4. Material Properties
# Example uses 2(1-0.0002im)
eps_r_val = 2.0 * (1 - 0.0002im)
permittivities = fill(eps_r_val, size(mesh.tetras, 2))

# 5. Assembly
println("Assembling Impedance Matrix...")
Z = assemble_impedance_matrix(vefie, basis, permittivities)
println("Z size: $(size(Z))")

# 6. Excitation (Plane Wave)
println("Computing Excitation...")
theta_inc = π/4
phi_inc = 0.0
k = vefie.k
k_hat = [sin(theta_inc)*cos(phi_inc), sin(theta_inc)*sin(phi_inc), cos(theta_inc)]
E0 = [0.0, 1.0, 0.0] # E_y polarization (phi=90) or E_theta?
# If theta=45, phi=0. k_hat = (0.707, 0, 0.707).
# E0 = (0, 1, 0) is perpendicular to k_hat. So it's valid.
# This corresponds to TE (E perpendicular to plane of incidence xz).

V = zeros(ComplexF64, n_unknowns)
gq = vefie.gq_info
for m in 1:n_unknowns
    bf = basis.functions[m]
    val = 0.0 + 0.0im
    for i_supp in 1:2
        tet_idx = bf.support[i_supp]
        if tet_idx == 0; continue; end
        
        verts = mesh.node[:, mesh.tetras[:, tet_idx]]
        v1=verts[:,1]; v2=verts[:,2]; v3=verts[:,3]; v4=verts[:,4]
        vol = abs(dot(v2-v1, cross(v3-v1, v4-v1))) / 6.0
        
        for i_qp in 1:length(gq.weight)
            L = gq.coordinate[:, i_qp]
            w = gq.weight[i_qp]
            r = verts * L
            
            f_m, _ = evaluate_swg(bf, i_supp, r, verts, vol)
            
            phase = exp(-im * k * dot(k_hat, r))
            E_inc = E0 * phase
            
            val += dot(f_m, E_inc) * w * vol
        end
    end
    V[m] = val
end

# 7. Solve
println("Solving system...")
I = Z \ V
println("Solved.")

# 8. Save Results
println("Saving coefficients...")
df = DataFrame(Real = real.(I), Imag = imag.(I))
CSV.write("VEFIE_Direct_Coeffs.csv", df)
println("Done.")
