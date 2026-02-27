using Test

@testset "EMSuite" begin
    include("test_materials.jl")
    include("test_geometry.jl")
    include("test_basis_functions.jl")
    include("test_integral_equations.jl")
    include("test_solvers.jl")
    include("test_mlfma.jl")
    include("test_parallel.jl")
    include("test_mpiarray.jl")
    include("test_lebedev.jl")
    include("test_postprocessing.jl")
    include("test_io.jl")
    include("test_integration.jl")
    include("test_workflow.jl")
    include("test_scfie.jl")
end
