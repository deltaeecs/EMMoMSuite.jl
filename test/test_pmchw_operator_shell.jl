using Test
using EMSuite
using LinearAlgebra

function make_pmchw_shell_fixture(; freq = 150e6, eps_r = 4.0, radius = 0.1, n_theta = 4, n_phi = 6)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r)
    shell = PMCHWBlockOperator(pmchw, basis)
    return pmchw, basis, shell
end

@testset "PMCHW block/operator shell" begin
    pmchw, basis, shell = make_pmchw_shell_fixture()
    N = num_basis(basis)
    Z = assemble_impedance_matrix(pmchw, basis)

    @test size(shell) == (2N, 2N)
    @test eltype(shell) == ComplexF64

    @testset "dense weak form matches direct PMCHW matrix" begin
        @test weak_form(shell) ≈ Z
        @test Matrix(shell) ≈ Z
    end

    @testset "block split reconstructs dense matrix" begin
        blocks = pmchw_blocks(shell)
        Z_reconstructed = [blocks.EJ blocks.EM; blocks.HJ blocks.HM]
        @test Z_reconstructed ≈ Z
        @test norm(blocks.EM + blocks.HJ) / (norm(blocks.EM) + 1e-30) < 1e-10
    end

    @testset "shell matvec matches weak form matvec" begin
        x = ComplexF64.(randn(2N), randn(2N))
        y_shell = shell * x
        y_dense = Z * x
        @test y_shell ≈ y_dense rtol=1e-12 atol=1e-12

        y_mul = zeros(ComplexF64, 2N)
        mul!(y_mul, shell, x)
        @test y_mul ≈ y_dense rtol=1e-12 atol=1e-12
    end

    @testset "strong form lives at shell layer" begin
        pairing = pmchw_block_pairing_matrix(basis)
        expected_default = pairing \ (Z * pairing)
        @test strong_form(shell) ≈ expected_default
        @test norm(strong_form(shell) - weak_form(shell)) / (norm(weak_form(shell)) + 1e-30) > 1e-10

        D = Diagonal(fill(2.0 + 0.0im, 2N))
        paired_shell = PMCHWBlockOperator(pmchw, basis; test_pairing = Matrix(D), trial_pairing = Matrix(D))
        expected = Matrix(D) \ (Z * Matrix(D))
        @test strong_form(paired_shell) ≈ expected
    end

    @testset "GMRES can solve through shell interface" begin
        source = DeltaGapSource(pmchw.freq, [1], 1.0 + 0.0im)
        V = excitation_vector(pmchw, source, basis)
        solver = GMRESSolver(restart = 40, maxiter = 6, tol = 1e-10, verbose = false)

        I_shell, hist_shell = solve!(solver, shell, V)
        I_dense, hist_dense = solve!(solver, Z, V)

        @test norm(I_shell - I_dense) / (norm(I_dense) + 1e-30) < 1e-8
        @test !isempty(hist_shell)
        @test !isempty(hist_dense)
    end
end