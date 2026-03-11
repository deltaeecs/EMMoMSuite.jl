using Test
using EMSuite
using StaticArrays
using LinearAlgebra

@testset "EFIE Matrix Assembly" begin
    # Create a simple mesh with 2 triangles sharing an edge
    # Nodes:
    # 1: (0,0,0)
    # 2: (1,0,0)
    # 3: (0,1,0)
    # 4: (1,1,0)

    nodes = [
        0.0 1.0 0.0 1.0
        0.0 0.0 1.0 1.0
        0.0 0.0 0.0 0.0
    ]

    # Elements (Triangles)
    # T1: 1-2-3
    # T2: 2-4-3 (Note: 2-3 is the shared edge)
    elements = [
        1 2
        2 4
        3 3
    ]

    tags = [1, 1]

    mesh = TriangleMesh(2, nodes, elements, tags)

    # Create RWG basis
    basis = RWGBasis(mesh)

    # 1 internal edge -> 1 basis function
    @test num_basis(basis) == 1

    # Create EFIE operator
    freq = 3e6 # 3 MHz (Even lower frequency)
    efie = EFIE(freq)

    # Assemble matrix
    Z = assemble_impedance_matrix(efie, basis)

    # Check dimensions
    @test size(Z) == (1, 1)

    # Check value is non-zero
    @test abs(Z[1, 1]) > 0.0

    # Check symmetry (trivial for 1x1)

    # Create a larger mesh (2x2 square, 8 triangles) to check symmetry
    # Or just manually add more triangles

    # Let's try a slightly larger mesh
    # 3 triangles in a strip
    # 1-2-3, 2-4-3, 4-5-3?

    # Nodes:
    # 1:(0,0), 2:(1,0), 3:(0,1), 4:(1,1), 5:(2,0), 6:(2,1)
    nodes2 = [
        0.0 1.0 0.0 1.0 2.0 2.0
        0.0 0.0 1.0 1.0 0.0 1.0
        0.0 0.0 0.0 0.0 0.0 0.0
    ]

    # T1: 1-2-3
    # T2: 2-4-3
    # T3: 2-5-4
    # T4: 5-6-4
    elements2 = [
        1 2 2 5
        2 4 5 6
        3 3 4 4
    ]

    mesh2 = TriangleMesh(4, nodes2, elements2, [1, 1, 1, 1])
    basis2 = RWGBasis(mesh2)

    # Edges:
    # T1-T2: 2-3 (internal)
    # T2-T3: 2-4 (internal)
    # T3-T4: 5-4 (internal)
    # Plus boundary edges

    # Should have 3 basis functions?
    # Let's check
    nb = num_basis(basis2)
    # println("Num basis: ", nb)

    Z2 = assemble_impedance_matrix(efie, basis2)

    @test size(Z2) == (nb, nb)

    # Check symmetry
    # Note: Due to semi-analytical integration (analytical source, numerical test),
    # perfect symmetry is not expected with low-order quadrature.
    # Relative error of 1e-3 is acceptable.
    @test isapprox(Z2, transpose(Z2), rtol = 1e-3)
end

@testset "MFIE Matrix Assembly" begin
    # Use the same mesh as EFIE
    nodes = [
        0.0 1.0 0.0 1.0
        0.0 0.0 1.0 1.0
        0.0 0.0 0.0 0.0
    ]
    elements = [
        1 2
        2 4
        3 3
    ]
    tags = [1, 1]
    mesh = TriangleMesh(2, nodes, elements, tags)
    basis = RWGBasis(mesh)

    freq = 300e6
    mfie = MFIE(freq)

    Z = assemble_impedance_matrix(mfie, basis)

    @test size(Z) == (1, 1)

    # MFIE matrix is generally not symmetric
    # But for 1x1 it is symmetric (scalar).

    # Check value
    # For planar triangles, self term is 0.5 * Gram matrix.
    # Interaction term is 0 (coplanar).
    # So Z should be 0.5 * Gram.

    # Let's calculate Gram manually for this case.
    # Edge length l = sqrt(2).
    # Area A = 0.5.
    # f = l/2A rho = sqrt(2) rho.
    # int f.f dS = 2 int rho.rho dS.
    # rho.rho on unit triangle: int (x^2+y^2) dx dy ?
    # Depends on rho definition.
    # Anyway, it should be positive real.

    @test real(Z[1, 1]) > 0.0
    @test abs(imag(Z[1, 1])) < 1e-10 # Should be real for self-term
end

@testset "CFIE Matrix Assembly" begin
    # Use the same mesh
    nodes = [
        0.0 1.0 0.0 1.0
        0.0 0.0 1.0 1.0
        0.0 0.0 0.0 0.0
    ]
    elements = [
        1 2
        2 4
        3 3
    ]
    tags = [1, 1]
    mesh = TriangleMesh(2, nodes, elements, tags)
    basis = RWGBasis(mesh)

    freq = 300e6
    cfie = CFIE(freq, 0.5)

    Z = assemble_impedance_matrix(cfie, basis)

    @test size(Z) == (1, 1)
    @test abs(Z[1, 1]) > 0.0
end

@testset "PlaneWave Polarization" begin
    src = PlaneWave(1.0e9, π / 4, π, [0.0, 0.0, 1.0])
    k_hat = [sin(src.theta) * cos(src.phi), sin(src.theta) * sin(src.phi), cos(src.theta)]

    @test isapprox(norm(src.polarization), 1.0; atol = 1e-12)
    @test isapprox(dot(src.polarization, k_hat), 0.0; atol = 1e-12)

    @test_throws ArgumentError PlaneWave(1.0e9, 0.0, 0.0, [0.0, 0.0, 2.0])
end
