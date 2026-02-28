# Phase 8.1 benchmark: Measure Z assembly improvement (per-row lock + batched accumulation)
using EMSuite, LinearAlgebra, Printf, Statistics

const MESH_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles")

# Warmup
set_frequency!(3e8)
mesh_plate = read_nas_mesh(joinpath(@__DIR__, "plate_benchmark.nas"), scale=1.0)
basis_plate = RWGBasis(mesh_plate)
Z_w = assemble_impedance_matrix(EFIE(3e8), basis_plate)
println("Warmup done.  N_plate = ", num_basis(basis_plate))

# ===== Case 1: Plate EFIE Direct =====
println("\n===== Case 1: Plate EFIE Direct (N=$(num_basis(basis_plate))) =====")
times1 = Float64[]
for trial in 1:3
    t = @elapsed begin
        Z = assemble_impedance_matrix(EFIE(3e8), basis_plate)
    end
    push!(times1, t)
    @printf("  Trial %d: Z assembly = %.3fs\n", trial, t)
end
@printf("  → Median Z_assembly = %.3fs  (baseline: 1.02s)\n", median(times1))

# ===== Case 2: Jet EFIE Direct =====
println("\n===== Case 2: Jet EFIE Direct =====")
mesh_jet = read_nas_mesh(joinpath(MESH_DIR, "jet_100MHz.nas"), scale=0.001)
set_frequency!(1e8)
basis_jet = RWGBasis(mesh_jet)
println("  N_jet = ", num_basis(basis_jet))

# warmup
Z_jw = assemble_impedance_matrix(EFIE(1e8), basis_jet)

times2 = Float64[]
for trial in 1:3
    t = @elapsed begin
        Z = assemble_impedance_matrix(EFIE(1e8), basis_jet)
    end
    push!(times2, t)
    @printf("  Trial %d: Z assembly = %.3fs\n", trial, t)
end
@printf("  → Median Z_assembly = %.3fs  (baseline: 20.70s)\n", median(times2))

# ===== Case 2b: Jet CFIE Direct =====
println("\n===== Case 2b: Jet CFIE Direct (N=$(num_basis(basis_jet))) =====")

# warmup CFIE
Z_cw = assemble_impedance_matrix(CFIE(1e8, 0.5), basis_jet)

times2b = Float64[]
for trial in 1:3
    t = @elapsed begin
        Z = assemble_impedance_matrix(CFIE(1e8, 0.5), basis_jet)
    end
    push!(times2b, t)
    @printf("  Trial %d: Z assembly = %.3fs\n", trial, t)
end
@printf("  → Median Z_assembly = %.3fs  (baseline: 168.29s)\n", median(times2b))

# ===== Summary =====
println("\n" * "="^60)
println("Phase 8.1 Z Assembly Benchmark Summary")
println("="^60)
@printf("Case 1  Plate EFIE:  %.3fs → %.3fs  (%.1f%%)\n", 1.02, median(times1), (median(times1)/1.02 - 1)*100)
@printf("Case 2  Jet EFIE:    %.3fs → %.3fs  (%.1f%%)\n", 20.70, median(times2), (median(times2)/20.70 - 1)*100)
@printf("Case 2b Jet CFIE:    %.3fs → %.3fs  (%.1f%%)\n", 168.29, median(times2b), (median(times2b)/168.29 - 1)*100)
println("="^60)
