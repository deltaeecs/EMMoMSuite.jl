using Test
using EMMoMSuite
using LinearAlgebra

function make_nmuller_fixture(; freq = 120e6, eps_r = 4.0, mu_r = 1.0, radius = 0.1, n_theta = 4, n_phi = 6)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    op = NMuller(freq, eps_r, mu_r)
    return op, basis
end

@testset "N-Muller dense baseline" begin
    op, basis = make_nmuller_fixture()
    pmchw = PMCHW(op.freq, op.eps_r, op.mu_r)
    Z = assemble_impedance_matrix(op, basis)
    Z_pmchw = assemble_impedance_matrix(pmchw, basis)
    N = num_basis(basis)

    @test size(Z) == (2N, 2N)
    @test eltype(Z) == ComplexF64
    @test norm(Z) > 0
    @test !any(isnan, Z)
    @test !any(isinf, Z)

    @testset "block structure is non-trivial and formulation differs from PMCHW" begin
        Z11 = Z[1:N, 1:N]
        Z12 = Z[1:N, N+1:2N]
        Z21 = Z[N+1:2N, 1:N]
        Z22 = Z[N+1:2N, N+1:2N]

        @test norm(Z11) > 0
        @test norm(Z12) > 0
        @test norm(Z21) > 0
        @test norm(Z22) > 0
        @test norm(Z - Z_pmchw) / (norm(Z_pmchw) + 1e-30) > 1e-6
    end

    @testset "plane-wave excitation length matches 2N system" begin
        src = PlaneWave(op.freq, π / 2, 0.0, [0.0, 0.0, 1.0])
        V = excitation_vector(op, src, basis)
        @test length(V) == 2N
        @test norm(V) > 0
    end

    @testset "dense direct solve is executable" begin
        src = PlaneWave(op.freq, π / 2, 0.0, [0.0, 0.0, 1.0])
        V = excitation_vector(op, src, basis)
        I = solve!(LUSolver(), Z, V)
        @test length(I) == 2N
        @test norm(I) > 0
        @test !any(isnan, I)
    end
end