using Test
using EMMoMSuite

@testset "N-Muller DeltaGap 激励 + input_impedance heuristic" begin
    mesh = generate_sphere_mesh(0.1, 4, 6)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    @test N > 0

    freq = 300e6
    eps_r = 4.0
    pmchw = PMCHW(freq, eps_r)
    nmuller = NMuller(freq, eps_r)

    V_feed = 1.0 + 0im
    feed = DeltaGapSource(freq, [1], V_feed)
    edge_len = basis.functions[1].edge_length

    @testset "DeltaGap vector preserves current heuristic layout" begin
        V_pmchw = excitation_vector(pmchw, feed, basis)
        V_nmuller = excitation_vector(nmuller, feed, basis)

        @test length(V_nmuller) == 2N
        @test all(iszero, V_nmuller[1:N])
        @test V_nmuller[N + 1] ≈ EMMoMSuite.Constants.eps0 * V_feed * edge_len
        @test V_nmuller[N + 1:2N] ≈ EMMoMSuite.Constants.eps0 .* V_pmchw[1:N]
    end

    @testset "input_impedance currently uses front-half heuristic only" begin
        I_2N = zeros(ComplexF64, 2N)
        I_2N[1] = 2.0 / edge_len
        I_2N[N + 1] = 999.0 + 999.0im
        Z_in = input_impedance(nmuller, feed, I_2N, basis)
        @test Z_in ≈ ComplexF64(V_feed) / 2.0
    end

    @testset "input_impedance throws when front-half heuristic current is zero" begin
        I_zero = zeros(ComplexF64, 2N)
        I_zero[N + 1] = 1.0 + 0im
        @test_throws ErrorException input_impedance(nmuller, feed, I_zero, basis)
    end

    @testset "direct solve still produces finite heuristic impedance" begin
        Z_nmuller = assemble_impedance_matrix(nmuller, basis)
        V_nmuller = excitation_vector(nmuller, feed, basis)
        I_nmuller = Z_nmuller \ V_nmuller
        Zin_nmuller = input_impedance(nmuller, feed, I_nmuller, basis)

        @test isfinite(real(Zin_nmuller))
        @test isfinite(imag(Zin_nmuller))
        @test abs(Zin_nmuller) > 0
    end
end