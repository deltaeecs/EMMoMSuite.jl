# Quick MatVec accuracy test for all MLFMA operators
using EMSuite
using LinearAlgebra

println("=" ^ 60)
println("MLFMA MatVec Accuracy Tests")
println("=" ^ 60)

mesh_file = joinpath(@__DIR__, "..", "..", "MoM_Basics", "meshfiles", "TriTetra.nas")
tetra_file = joinpath(@__DIR__, "..", "..", "MoM_Basics", "meshfiles", "Tetra.nas")

surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file, scale=0.001)

# --- Test SCFIE MLFMA ---
for (label, freq) in [("2GHz", 2e9), ("4GHz", 4e9)]
    set_frequency!(freq)
    λ = 299792458.0 / freq
    eps_r = fill(2.0*(1-0.001im), num_elements(vol_mesh))
    sb = RWGBasis(surf_mesh)
    vb = SWGBasis(vol_mesh)
    n = num_basis(sb) + num_basis(vb)

    scfie = SCFIE(freq, eps_r; alpha=0.5)
    Z_d = assemble_impedance_matrix(scfie, sb, vb)
    op = MLFMAOperator(scfie, [sb, vb], 0.25λ)

    x = randn(ComplexF64, n)
    err = norm(Z_d*x - op*x) / norm(Z_d*x)
    println("SCFIE MLFMA [$label] N=$n, Rel Err: $(round(err*100, digits=4))%")
end

# --- Test EFIE MLFMA ---
begin
    freq = 4e9
    set_frequency!(freq)
    λ = 299792458.0 / freq
    basis = RWGBasis(surf_mesh)
    n = num_basis(basis)
    efie = EFIE(freq)
    Z_d = assemble_impedance_matrix(efie, basis)
    op = MLFMAOperator(efie, basis, 0.25λ)
    x = randn(ComplexF64, n)
    err = norm(Z_d*x - op*x) / norm(Z_d*x)
    println("EFIE MLFMA [4GHz] N=$n, Rel Err: $(round(err*100, digits=4))%")
end

# --- Test VEFIE MLFMA (TriTetra vol) ---
begin
    freq = 4e9
    set_frequency!(freq)
    λ = 299792458.0 / freq
    basis = SWGBasis(vol_mesh)
    n = num_basis(basis)
    eps_r = fill(complex(2.0, -0.002), num_elements(vol_mesh))
    vefie = VEFIE(freq, eps_r)
    Z_d = assemble_impedance_matrix(vefie, basis)
    op = MLFMAOperator(vefie, basis, 0.25λ)
    x = randn(ComplexF64, n)
    err = norm(Z_d*x - op*x) / norm(Z_d*x)
    println("VEFIE MLFMA [TriTetra 4GHz] N=$n, Rel Err: $(round(err*100, digits=4))%")
end

# --- Test VEFIE MLFMA (pure Tetra) ---
if isfile(tetra_file)
    tet_mesh = read_nas_mesh(tetra_file, scale=0.001)
    freq = 1.5e9
    set_frequency!(freq)
    λ = 299792458.0 / freq
    basis = SWGBasis(tet_mesh)
    n = num_basis(basis)
    eps_r = fill(complex(2.0, -0.002), num_elements(tet_mesh))
    vefie = VEFIE(freq, eps_r)
    Z_d = assemble_impedance_matrix(vefie, basis)
    op = MLFMAOperator(vefie, basis, 0.25λ)
    x = randn(ComplexF64, n)
    err = norm(Z_d*x - op*x) / norm(Z_d*x)
    println("VEFIE MLFMA [Tetra 1.5GHz] N=$n, Rel Err: $(round(err*100, digits=4))%")
end

println("\nDone.")
