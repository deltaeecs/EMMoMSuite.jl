# test_cov_ports.jl — Ports 模块覆盖率补测
using Test
using EMMoMSuite
using EMMoMSuite.PortsModule
using LinearAlgebra, SparseArrays

@testset "LumpedPort" begin
    p = lumped_port(1, 3; impedance = 50.0)
    @test port_id(p) == 1
    I = zeros(ComplexF64, 10)
    I[3] = 0.1 + 0.0im
    @test port_current(p, I) == 0.1
    @test port_voltage(p, I) == 50.0 * 0.1
    @test port_power(p, I) ≈ 0.5 * real(50.0 * 0.1 * conj(0.1))
    Z = zeros(ComplexF64, 10, 10)
    assemble_lumped_port_impedance!(Z, p)
    @test Z[3, 3] == 50.0
    V = zeros(ComplexF64, 10)
    add_port_excitation!(V, p)
    @test V[3] == 0  # :load 类型无激励贡献
    # 源型端口（voltage 激励）
    ps = lumped_port(2, 5; type = :voltage_source, voltage = 1.0 + 0im)
    @test ps.port_type == :voltage_source
    V2 = zeros(ComplexF64, 10)
    add_port_excitation!(V2, ps)
    @test V2[5] == 1.0
    @test port_voltage(ps, zeros(ComplexF64, 10)) == 0.0
end

@testset "S 参数提取与 Touchstone" begin
    ports = [lumped_port(1, 1), lumped_port(2, 2)]
    Zefie = ComplexF64[75.0 5.0; 5.0 25.0]
    S = extract_s_matrix(ports, Zefie)
    @test size(S) == (2, 2)
    data = sparam_data([1e9, 2e9], [S, S], [50.0, 50.0])
    @test data.frequencies == [1e9, 2e9]
    path = joinpath(tempdir(), "sp_$(getpid()).s2p")
    write_touchstone(data, path; version = 1)
    d2 = read_touchstone(path)
    @test length(d2.frequencies) == 2
    @test size(d2.s_matrices[1]) == (2, 2)
    rm(path; force = true)
end

@testset "WavePort 与模式" begin
    wp = WavePort(1; mode = :TE10, a = 1.0, b = 0.5)
    @test wp.mode == :TE10
    pm = compute_port_modes(wp, 3e9, 2)
    @test pm isa PortModes
    @test length(pm.beta) >= 1
    @test length(pm.mode_types) == length(pm.beta)
    @test length(pm.Z_c) == length(pm.beta)
    # 模式解析
    @test EMMoMSuite.PortsModule._parse_mode_indices(:TM21) == (:TM, 2, 1)
    @test_throws ErrorException port_current(wp, zeros(ComplexF64, 2))
end

@testset "差分/同轴/离散端口与混合模式" begin
    dp = DifferentialPairPort(1, 1, 2)
    @test port_id(dp) == 1
    cp = CoaxPort(1, 0.5, 2.0)
    @test coax_impedance(cp) ≈ (376.730313461 / (2π)) * log(4.0)
    disc = DiscretePort(1, 4, 5, :resistor, 100.0)
    @test port_id(disc) == 1
    pairs = [DifferentialPairPort(1, 1, 2)]
    M = mixed_mode_transform_matrix(pairs)
    @test size(M) == (2, 2)
    # 混合模式 S 参数转换
    S = [0.1 0.9; 0.9 0.1]
    @test size(sdd_matrix(S, pairs)) == (1, 1)
    @test size(scc_matrix(S, pairs)) == (1, 1)
    @test size(scd_matrix(S, pairs)) == (1, 1)
end
