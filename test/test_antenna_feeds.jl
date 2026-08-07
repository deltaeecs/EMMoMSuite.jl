using Test
using EMMoMSuite

@testset "Antenna feed edge selection" begin
    freq = 3e8
    L = 0.5
    mesh = generate_cylinder_mesh(0.005, L, 8, 20; closed=false)
    set_frequency!(freq)
    basis = RWGBasis(mesh)

    feed_edges = select_gap_feed_edges(basis; axis=3, center=0.0)

    @test length(feed_edges) == 8
    @test all(!basis.functions[idx].is_boundary for idx in feed_edges)
    @test all(isapprox(basis.functions[idx].center[3], 0.0; atol=1e-12) for idx in feed_edges)

    efie = EFIE(freq)
    source = DeltaGapSource(freq, feed_edges, 1.0 + 0.0im)
    I = assemble_impedance_matrix(efie, basis) \ excitation_vector(efie, source, basis)
    Zin = input_impedance(source, I, basis)

    @test real(Zin) > 20.0
    @test imag(Zin) > 10.0
end