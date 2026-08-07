using Test
using EMMoMSuite
using TOML

@testset "Workflow Integration" begin
    # 1. Create a dummy mesh file
    mesh_file = "test_plate.nas"
    open(mesh_file, "w") do f
        write(f, """
GRID, 1, , 0.0, 0.0, 0.0
GRID, 2, , 1.0, 0.0, 0.0
GRID, 3, , 0.0, 1.0, 0.0
GRID, 4, , 1.0, 1.0, 0.0
CTRIA3, 1, 1, 1, 2, 3
CTRIA3, 2, 1, 2, 4, 3
""")
    end
    
    # 2. Create a config file
    config_file = "test_config.toml"
    config_data = Dict(
        "simulation" => Dict(
            "name" => "TestPlate",
            "frequency" => 300e6,
            "solver_type" => "GMRES"
        ),
        "geometry" => Dict(
            "mesh_file" => mesh_file,
            "unit_scale" => 1.0
        ),
        "basis" => Dict(
            "type" => "RWG"
        ),
        "excitation" => Dict(
            "type" => "PlaneWave",
            "theta" => 0.0,
            "phi" => 0.0,
            "polarization" => [1.0, 0.0, 0.0]
        ),
        "solver" => Dict(
            "tolerance" => 1e-3,
            "max_iter" => 100,
            "restart" => 10
        ),
        "output" => Dict(
            "directory" => "test_results",
            "save_fields" => true,
            "save_rcs" => true,
            "format" => "h5"
        )
    )
    
    open(config_file, "w") do f
        TOML.print(f, config_data)
    end
    
    # 3. Run simulation
    try
        result = run_simulation(config_file)
        
        @test result isa SimulationResult
        @test length(result.currents) > 0
        @test isfile(joinpath("test_results", "TestPlate_result.h5"))
        @test isfile(joinpath("test_results", "TestPlate.log"))
        
    finally
        # Cleanup
        # Close logger to release file lock
        # This is tricky because global_logger is global.
        # In a real app, we might not need to delete logs immediately.
        # For tests, we can try to force GC or just ignore the error if we can't delete.
        
        rm(mesh_file, force=true)
        rm(config_file, force=true)
        try
            rm("test_results", recursive=true, force=true)
        catch e
            @warn "Could not remove test_results directory: $e"
        end
    end
end
