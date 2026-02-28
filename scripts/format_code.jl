#!/usr/bin/env julia
# Format all Julia source files in src/ using JuliaFormatter
# Config: .JuliaFormatter.toml in project root (loaded automatically)
# Run with: julia --startup-file=no scripts/format_code.jl
# Note: test/ is excluded due to JuliaFormatter v1.x compatibility issues with
# certain Julia 1.12 syntax (@test a ≈ b atol=c patterns).

using JuliaFormatter   # loads from global environment

project_root = joinpath(@__DIR__, "..")
src_dir = joinpath(project_root, "src")

println("Formatting Julia files in src/ ...")
changed = format(src_dir)

if changed
    println("Done: some files were reformatted.")
else
    println("Done: all files were already well-formatted.")
end
