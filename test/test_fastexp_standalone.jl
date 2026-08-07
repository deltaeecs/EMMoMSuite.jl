#!/usr/bin/env julia
"""
test_fastexp_standalone.jl

Standalone test for FastExp module without loading full EMMoMSuite.
"""

# Activate environment for dependencies
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# Load FastExp module directly
include(joinpath(@__DIR__, "..", "src", "IntegralEquations", "FastExp.jl"))
using .FastExpModule
using Printf
using Test

println("="^70)
println("  FastExp Lookup Table Validation (Standalone)")
println("="^70)

# Test parameters
freq = 1e8  # 100 MHz
k = 2π * freq / 299792458.0
λ = 2π / k

println("\n[1/3] Accuracy Test")

# Create lookup table
table = FastExpTable(k)
println("  Table created: $(table.n_entries) entries, R_max = $(round(table.R_max/λ, digits=2))λ")

# Test R values from 0.01λ to 20λ
test_Rs = λ .* [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]

errors = Float64[]
for R in test_Rs
    G_exact = exp(-im * k * R) / (4π * R)
    G_fast = fast_green_func(table, R)
    
    rel_error = abs(G_fast - G_exact) / abs(G_exact)
    push!(errors, rel_error)
    
    passed = rel_error < 1e-3 ? "✓" : "✗"
    @printf("  %s R = %6.3f λ: rel_error = %.2e\n", passed, R/λ, rel_error)
end

max_rel_error = maximum(errors)

println("\n  Maximum relative error: $(round(max_rel_error*100, digits=4))%")

println("\n[2/3] Singular Case Test")

# Test R → 0 (should return zero)
G_small = fast_green_func(table, 1e-12)
println("  G(R≈0) = $(abs(G_small)) (correctly handled: $(abs(G_small) < 1e-10 ? "✓" : "✗"))")

println("\n[3/3] Performance Benchmark")

# Generate random R values in typical range [0.1λ, 10λ]
N_samples = 100_000
R_vals = λ .* (0.1 .+ 9.9 .* rand(N_samples))

# Benchmark FastExp
GC.gc()
t_fast = @elapsed begin
    G_fast = [fast_green_func(table, R) for R in R_vals]
end

# Benchmark direct exp()
GC.gc()
t_direct = @elapsed begin
    G_direct = [exp(-im * k * R) / (4π * R) for R in R_vals]
end

speedup = t_direct / t_fast

@printf("  Direct exp():  %.3f ms  (%.2f ns/call)\n", t_direct*1000, t_direct/N_samples*1e9)
@printf("  FastExp LUT:   %.3f ms  (%.2f ns/call)\n", t_fast*1000, t_fast/N_samples*1e9)
@printf("  Speedup:       %.2f×\n", speedup)

# Accuracy check
max_error = maximum(abs.(G_fast .- G_direct) ./ abs.(G_direct))
@printf("  Max rel error: %.2e (%.3f%%)\n", max_error, max_error*100)

println("\n" * "="^70)
if speedup > 1.5 && max_error < 1e-3
    println("  ✓✓✓ All Tests Passed! ✓✓✓")
    println("  FastExp is $(round(speedup, digits=1))× faster with $(round(max_error*100, digits=3))% max error")
else
    println("  ⚠ Warning: Performance or accuracy issues detected")
end
println("="^70)
