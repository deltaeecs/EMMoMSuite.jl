# P3 性能基准：AIM/IE-FFT vs MLFMA（matvec 时间 + 精度）
using EMMoMSuite
using LinearAlgebra, Random
using EMMoMSuite.Geometry, EMMoMSuite.BasisFunctions, EMMoMSuite.IntegralEquations

function run_case(r, nm, h_ratio)
    freq = 300e6
    λ = 299792458.0 / freq
    mesh = generate_sphere_mesh(r, nm[1], nm[2])
    basis = RWGBasis(mesh)
    efie = EFIE(freq)
    N = num_basis(basis)
    Z = assemble_impedance_matrix(efie, basis)
    op = AIMOperator(efie, basis; h_ratio = h_ratio, near_radius = 0.35)
    Random.seed!(3)
    x = randn(ComplexF64, N)
    y_aim = op * x
    y_d = Z * x
    err = norm(y_aim - y_d) / norm(y_d)
    op * x
    a_aim = @allocated op * x
    t_aim = minimum(@elapsed(op * x) for _ in 1:10)
    opm = MLFMAOperator(efie, basis, 0.25 * λ)
    opm * x
    t_mlfma = minimum(@elapsed(opm * x) for _ in 1:5)
    println("case r=", r, "λ mesh=", nm, " N=", N, " grid=", op.grid.n,
            " | AIM err=", round(err, digits=4),
            " | matvec aim=", round(t_aim*1000, digits=2), " ms (alloc ", round(a_aim/1024, digits=1), " KB), mlfma=", round(t_mlfma*1000, digits=2), " ms",
            " | speedup=", round(t_mlfma / t_aim, digits=1))
end

run_case(0.6, (12, 24), 0.1)
