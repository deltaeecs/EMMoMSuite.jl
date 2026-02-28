using EMSuite, LinearAlgebra, Printf, Statistics

# ===== Case 1: Plate EFIE (small, N≈2640) =====
set_frequency!(3e8)
mesh_plate = read_nas_mesh(joinpath(@__DIR__, "plate_benchmark.nas"), scale=1.0)
basis_plate = RWGBasis(mesh_plate)
println("Plate N = ", num_basis(basis_plate))

# Warmup
assemble_impedance_matrix(EFIE(3e8), basis_plate)

times_plate = [(@elapsed assemble_impedance_matrix(EFIE(3e8), basis_plate)) for _ in 1:3]
baseline_plate = 1.02  # Pre-optimization EMSuite baseline (s)
@printf("Plate EFIE: %.3fs median (baseline %.3fs)  ratio=%.2fx vs baseline\n",
        median(times_plate), baseline_plate, baseline_plate / median(times_plate))

# ===== Case 2: Sphere CFIE (medium) to test merged assembly =====
# Use sphere mesh if available; else skip
sphere_mesh_path = joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles", "sphere_600MHz.nas")
if isfile(sphere_mesh_path)
    set_frequency!(6e8)
    mesh_sph = read_nas_mesh(sphere_mesh_path, scale=1.0)
    basis_sph = RWGBasis(mesh_sph)
    println("\nSphere N = ", num_basis(basis_sph))
    
    # Warmup
    assemble_impedance_matrix(CFIE(6e8, 0.5), basis_sph)
    
    times_cfie = [(@elapsed assemble_impedance_matrix(CFIE(6e8, 0.5), basis_sph)) for _ in 1:3]
    times_efie = [(@elapsed assemble_impedance_matrix(EFIE(6e8), basis_sph)) for _ in 1:3]
    times_mfie = [(@elapsed assemble_impedance_matrix(MFIE(6e8), basis_sph)) for _ in 1:3]
    
    @printf("Sphere EFIE:  %.3fs median\n", median(times_efie))
    @printf("Sphere MFIE:  %.3fs median\n", median(times_mfie))
    @printf("Sphere CFIE:  %.3fs median (separate EFIE+MFIE + in-place @.)\n", median(times_cfie))
    separate_time = median(times_efie) + median(times_mfie)
    @printf("Merged vs Separate EFIE+MFIE: %.3fs vs %.3fs  speedup=%.2fx\n",
            median(times_cfie), separate_time,
            separate_time / median(times_cfie))

    println("\n=== Summary ===")
    println("Plate EFIE:  $(round(median(times_plate),digits=3))s  (vs 1.02s baseline → $(round(1.02/median(times_plate),digits=2))×)")
    println("Sphere CFIE: $(round(median(times_cfie),digits=3))s  vs separate≈$(round(separate_time,digits=3))s → $(round(separate_time/median(times_cfie),digits=2))× speedup")
else
    println("\nSphere mesh not found, skipping CFIE test: $sphere_mesh_path")
    println("\n=== Summary ===")
    println("Plate EFIE:  $(round(median(times_plate),digits=3))s  (vs 1.02s baseline → $(round(1.02/median(times_plate),digits=2))×)")
end

