# test_cov_core.jl — Core 模块覆盖率补测
using Test
using EMMoMSuite
using EMMoMSuite.CoreModule
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using LinearAlgebra
using TOML

@testset "源与激励" begin
    pw = PlaneWave(300e6, 0.0, 0.0, [1.0, 0.0, 0.0])
    @test incident_field(pw, [0.0, 0.0, 0.0]) isa AbstractVector
    @test length(incident_field(pw, [0.0, 0.0, 0.0])) == 3
    mesh = generate_sphere_mesh(0.3, 4, 8)
    basis = RWGBasis(mesh)
    V = excitation_vector(pw, basis)
    @test length(V) == num_basis(basis)
    dg = DeltaGapSource(300e6, [1], 1.0 + 0im)
    V2 = excitation_vector(dg, basis)
    @test length(V2) == num_basis(basis)
end

@testset "配置与结果" begin
    sim = SimulationConfig(; name = "t", frequency = 1e9)
    @test sim.frequency == 1e9
    geo = GeometryConfig(; mesh_file = "x.msh")
    @test geo.unit_scale == 1.0
    basis_c = BasisConfig(; type = "RWG")
    @test basis_c.type == "RWG"
    exc = ExcitationConfig(; type = "PlaneWave")
    @test exc.theta == 0.0
    solv = SolverConfig(; tolerance = 1e-6)
    @test solv.max_iter == 1000
    out = OutputConfig()
    @test out isa OutputConfig
    cfg = EMMoMSuiteConfig(
        ; simulation = sim, geometry = geo, basis = basis_c,
        excitation = exc, solver = solv, output = out,
    )
    @test cfg.simulation.frequency == 1e9
    res = SimulationResult(cfg, [1.0 + 0im, 2.0 + 0im])
    @test length(res.currents) == 2
    @test res.metrics == Dict{String,Any}()
    # TOML 配置加载
    path = joinpath(tempdir(), "cfg_$(getpid()).toml")
    open(path, "w") do io
        write(io, """
        [simulation]
        name = "cov"
        frequency = 2e9
        [geometry]
        mesh_file = "m.msh"
        [basis]
        type = "RWG"
        [excitation]
        type = "PlaneWave"
        [solver]
        tolerance = 1e-5
        [output]
        directory = "res"
        """)
    end
    cfg2 = load_config(path)
    @test cfg2 isa EMMoMSuiteConfig
    rm(path; force = true)
end
