#!/usr/bin/env julia
"""
test_fastexp.jl

Test FastExp lookup table correctness and performance.
可被 runtests.jl include，也可独立运行。
"""

if !isdefined(Main, :EMSuite)
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using EMSuite
end
using EMSuite.IntegralEquations.VEFIEModule: FastExpModule
using Printf
using Statistics
using Test

# Test parameters
freq = 1e8  # 100 MHz
k = 2π * freq / 299792458.0
λ = 2π / k

@testset "FastExp Correctness" begin

    # Create lookup table
    table = FastExpModule.FastExpTable(k)
    
    # Test R values from 0.01λ to 20λ
    test_Rs = λ .* [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    
    max_rel_error = 0.0
    for R in test_Rs
        G_exact = exp(-im * k * R) / (4π * R)
        G_fast = FastExpModule.fast_green_func(table, R)
        
        rel_error = abs(G_fast - G_exact) / abs(G_exact)
        max_rel_error = max(max_rel_error, rel_error)
        
        @test rel_error < 1e-3  # < 0.1% error
    end
end

@testset "FastExp Singular Case" begin

    table = FastExpModule.FastExpTable(k)
    
    # Test R → 0 (should return zero)
    G_small = FastExpModule.fast_green_func(table, 1e-12)
    @test abs(G_small) < 1e-10
end

@testset "FastExp Performance" begin

    table = FastExpModule.FastExpTable(k)
    
    # Generate random R values in typical range [0.1λ, 10λ]
    N_samples = 100_000
    R_vals = λ .* (0.1 .+ 9.9 .* rand(N_samples))
    
    # Benchmark FastExp
    GC.gc()
    t_fast = @elapsed begin
        G_fast = [FastExpModule.fast_green_func(table, R) for R in R_vals]
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
    
    @test speedup > 1.5  # Should be at least 1.5× faster
    
    # Accuracy check
    max_error = maximum(abs.(G_fast .- G_direct) ./ abs.(G_direct))
    
    @test max_error < 1e-3
end
