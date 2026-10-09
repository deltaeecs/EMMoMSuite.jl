# run_case.jl — entry point for the gmsh-driven RCS case family.
#
# Usage (from anywhere; activates the benchmark project so Plots is available):
#   julia --project=benchmark benchmark/cases/run_case.jl                # all cases
#   julia --project=benchmark benchmark/cases/run_case.jl G1_sphere_CFIE_Mie G3_box_EFIE

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))   # benchmark project: EMMoMSuite + Plots

include(joinpath(@__DIR__, "CaseRunner.jl"))
using .CaseRunner

names = String[]
for a in ARGS
    if a == "-h" || a == "--help"
        println("""
        Usage: julia --project=benchmark benchmark/cases/run_case.jl [case names...]
        No arguments runs every case registered in cases.toml.""")
        exit(0)
    end
    push!(names, a)
end

results = CaseRunner.run_cases(names)
failed = count(r -> !isfile(joinpath(r.outdir, "report.md")), results)
exit(failed == 0 ? 0 : 1)
