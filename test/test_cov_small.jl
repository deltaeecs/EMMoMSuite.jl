# test_cov_small.jl — Core/Utilities/Solvers/Materials 覆盖率补测（轻量冒烟）
using Test
using EMMoMSuite
using EMMoMSuite.CoreModule
using EMMoMSuite.MaterialsModule
using EMMoMSuite.Utilities
using EMMoMSuite.Solvers
using LinearAlgebra, SparseArrays, Random
using HDF5

@testset "Core 常量" begin
    @test Constants.c0 ≈ 299792458.0
    @test Constants.eta0 ≈ 120π atol = 1.0
    @test Constants.mu0 > 0
    @test Constants.eps0 > 0
end

@testset "Materials 基础与色散模型" begin
    @test permittivity(PEC()) == Inf
    @test permeability(PEC()) == Constants.mu0
    @test permittivity(PMC()) == Constants.eps0
    @test permeability(PMC()) == Inf
    d = Dielectric(4.0)
    @test permittivity(d) ≈ 4.0 * Constants.eps0
    @test permeability(d) ≈ Constants.mu0
    @test impedance(d) ≈ sqrt(Constants.mu0 / (4.0 * Constants.eps0))
    @test Dielectric(4.0, 2.0) isa Dielectric
    @test ImpedanceSurface(complex(50.0)) isa ImpedanceSurface
    d = DebyeModel(4.0, [DebyePole(0.5, 1e-10)])
    dr = DrudeModel(1.0, 1e13, 1e8)
    lo = LorentzModel(1.0, [LorentzPole(0.1, 1e12, 1e8)])
    ω = 2π * 1e9
    @test eval_permittivity(d, ω) isa ComplexF64
    @test eval_permittivity(dr, ω) isa ComplexF64
    @test eval_permittivity(lo, ω) isa ComplexF64
    @test length(eval_permittivity(d, [ω, 2ω])) == 2
    @test eval_permittivity(lo, [ω, 2ω]) isa Vector
end

@testset "MaterialLibrary" begin
    lib = MaterialLibrary()
    entry = MaterialEntry("vac", Isotropic(1.0), (1e6, 1e12), v"1.0.0", "test")
    add_material!(lib, entry)
    @test get_material(lib, "vac", 1e9) isa MaterialModel
    @test_throws DomainError get_material(lib, "vac", 1e15)
    @test_throws KeyError get_material(lib, "nope", 1e9)
    aniso = Anisotropic([2.0 0 0; 0 2.0 0; 0 0 2.0], [1.0 0 0; 0 1.0 0; 0 0 1.0])
    @test aniso isa Anisotropic
    path = joinpath(tempdir(), "matlib_$(getpid()).jld2")
    save_library(lib, path)
    lib2 = load_library(path)
    @test get_material(lib2, "vac", 1e9) isa MaterialModel
    rm(path; force = true)
end

@testset "Utilities: 参数/Progress/数值工具" begin
    set_frequency!(1e9)
    @test get_k0() > 0
    @test get_eta0() > 0
    @test get_omega() > 0
    p = Progress(5)
    for _ in 1:5
        next!(p)
    end
    @test p.current == 5
    f(x) = x^2 - 2
    @test abs(find_zero_bisection(f, 1.0) - sqrt(2)) < 1e-6
    pts = rand(3, 20)
    idxs, dists = knn_bruteforce(pts, pts, 1)
    @test length(idxs) == 20
    @test all(d -> all(>=(0.0), d), dists)
end

@testset "Utilities: 稀疏矩阵存取" begin
    Random.seed!(7)
    S = sprand(ComplexF64, 6, 6, 0.4)
    path = joinpath(tempdir(), "spmat_$(getpid()).h5")
    h5open(path, "w") do fid
        g = create_group(fid, "sparse")
        save_sparse_matrix(g, "M", S)
    end
    S2 = h5open(path, "r") do fid
        load_sparse_matrix(fid["sparse"], "M")
    end
    @test size(S2) == (6, 6)
    @test nnz(S2) == nnz(S)
    rm(path; force = true)
end

@testset "Utilities: Mie 级数" begin
    th = collect(0.0:0.2:pi)
    r = calculate_mie_rcs_pec_sphere(0.5, 300e6, th)
    @test length(r) == length(th)
    @test all(isfinite, r)
    re, rh, ru = calculate_mie_rcs_dielectric_sphere(0.5, 300e6, th, 4.0, 1.0)
    @test length(re) == length(th)
    @test length(rh) == length(th)
    @test length(ru) == length(th)
end

@testset "Solvers: 直接与迭代求解" begin
    Random.seed!(11)
    A = randn(ComplexF64, 10, 10) + 20I
    b = randn(ComplexF64, 10)
    x = solve!(LUSolver(), A, b)
    @test norm(A * x - b) / norm(b) < 1e-10
    xg = solve!(GMRESSolver(restart = 10, maxiter = 200, tol = 1e-10), A, b)
    @test norm(A * xg - b) / norm(b) < 1e-8
    xb = solve!(BiCGSTABSolver(maxiter = 200, tol = 1e-10), A, b)
    @test norm(A * xb - b) / norm(b) < 1e-8
    # 预条件
    P = DiagonalPreconditioner(A)
    @test norm(P \ b - b ./ diag(A)) < 1e-10
    Pid = IdentityPreconditioner()
    @test Pid \ b == b
end
