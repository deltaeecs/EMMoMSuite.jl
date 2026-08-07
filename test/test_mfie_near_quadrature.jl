using Test

@testset "MFIE Near Quadrature" begin
    using EMMoMSuite
    using LinearAlgebra

    nodes = [
        0.0  1.0   0.0   0.05
        0.0  0.0   1.0   0.05
        0.0  0.0   0.0   0.01
    ]
    elements = [
        1  2
        2  4
        3  3
    ]
    tags = [1, 1]

    mesh = TriangleMesh(2, nodes, elements, tags)
    basis = RWGBasis(mesh)
    mfie = MFIE(300e6)

    tri1 = get_triangle_info(mesh, basis, 1)
    tri2 = get_triangle_info(mesh, basis, 2)

    z_actual = zeros(ComplexF64, 3, 3)
    EMMoMSuite.IntegralEquations.MFIEModule.mfie_interaction!(z_actual, mfie, tri1, tri2)

    r_test_near = EMMoMSuite.Geometry.get_global_quad_points(tri1, mfie.gq_near)
    r_src_near = EMMoMSuite.Geometry.get_global_quad_points(tri2, mfie.gq_near)
    z_near = zeros(ComplexF64, 3, 3)
    EMMoMSuite.IntegralEquations.MFIEModule.calc_k_term_fast!(
        z_near,
        mfie,
        tri1,
        tri2,
        r_test_near,
        r_src_near,
        mfie.gq_near.weight,
    )

    r_test_far = EMMoMSuite.Geometry.get_global_quad_points(tri1, mfie.gq_far)
    r_src_far = EMMoMSuite.Geometry.get_global_quad_points(tri2, mfie.gq_far)
    z_far = zeros(ComplexF64, 3, 3)
    EMMoMSuite.IntegralEquations.MFIEModule.calc_k_term_fast!(
        z_far,
        mfie,
        tri1,
        tri2,
        r_test_far,
        r_src_far,
        mfie.gq_far.weight,
    )

    @test isapprox(z_actual, z_near; rtol = 1e-12, atol = 1e-12)
    rel_diff = norm(z_near - z_far) / max(norm(z_near), eps())
    @test rel_diff > 1e-3
end

@testset "MFIE Close Nonadjacent Quadrature" begin
    using EMMoMSuite
    using LinearAlgebra

    nodes = [
        0.0  1.0  0.0  1.0  0.0  1.0  0.0
        0.0  0.0  1.0  1.0  0.0  0.0  1.0
        0.0  0.0  0.0  0.0  0.02 0.02 0.02
    ]
    elements = [
        1  2  5
        2  4  6
        3  3  7
    ]
    tags = [1, 1, 1]

    mesh = TriangleMesh(3, nodes, elements, tags)
    basis = RWGBasis(mesh)
    mfie = MFIE(300e6)

    tri1 = get_triangle_info(mesh, basis, 1)
    tri2 = get_triangle_info(mesh, basis, 3)

    @test EMMoMSuite.IntegralEquations.MFIEModule.needs_near_quadrature(tri1, tri2)

    z_actual = zeros(ComplexF64, 3, 3)
    EMMoMSuite.IntegralEquations.MFIEModule.mfie_interaction!(z_actual, mfie, tri1, tri2)

    r_test_near = EMMoMSuite.Geometry.get_global_quad_points(tri1, mfie.gq_near)
    r_src_near = EMMoMSuite.Geometry.get_global_quad_points(tri2, mfie.gq_near)
    z_near = zeros(ComplexF64, 3, 3)
    EMMoMSuite.IntegralEquations.MFIEModule.calc_k_term_fast!(
        z_near,
        mfie,
        tri1,
        tri2,
        r_test_near,
        r_src_near,
        mfie.gq_near.weight,
    )

    @test isapprox(z_actual, z_near; rtol = 1e-12, atol = 1e-12)
end

@testset "MFIE Gap-Based Near Quadrature" begin
    using EMMoMSuite
    using LinearAlgebra

    h = 0.35 * sqrt(3)
    x_tip = 2h / 3
    x_base = -h / 3
    half_side = 0.35

    nodes = [
        x_tip   x_base  x_base  1.3 - x_tip  1.3 - x_base  1.3 - x_base
        0.0     half_side  -half_side  0.0  half_side  -half_side
        0.0     0.0    0.0    0.02 0.02 0.02
    ]
    elements = [
        1  4
        2  5
        3  6
    ]
    tags = [1, 1]

    mesh = TriangleMesh(2, nodes, elements, tags)
    basis = RWGBasis(mesh)
    mfie = MFIE(300e6)

    tri1 = get_triangle_info(mesh, basis, 1)
    tri2 = get_triangle_info(mesh, basis, 2)

    r1 = maximum(norm(tri1.vertices[:, i] - tri1.center) for i = 1:3)
    r2 = maximum(norm(tri2.vertices[:, i] - tri2.center) for i = 1:3)
    center_dist = norm(tri1.center - tri2.center)
    @test center_dist > 1.5 * (r1 + r2)

    min_gap = EMMoMSuite.IntegralEquations.MFIEModule._triangle_triangle_distance(tri1, tri2)
    @test min_gap < 0.75 * (r1 + r2)
    @test EMMoMSuite.IntegralEquations.MFIEModule.needs_near_quadrature(tri1, tri2)

    z_actual = zeros(ComplexF64, 3, 3)
    EMMoMSuite.IntegralEquations.MFIEModule.mfie_interaction!(z_actual, mfie, tri1, tri2)

    r_test_near = EMMoMSuite.Geometry.get_global_quad_points(tri1, mfie.gq_near)
    r_src_near = EMMoMSuite.Geometry.get_global_quad_points(tri2, mfie.gq_near)
    z_near = zeros(ComplexF64, 3, 3)
    EMMoMSuite.IntegralEquations.MFIEModule.calc_k_term_fast!(
        z_near,
        mfie,
        tri1,
        tri2,
        r_test_near,
        r_src_near,
        mfie.gq_near.weight,
    )

    @test isapprox(z_actual, z_near; rtol = 1e-12, atol = 1e-12)
end

