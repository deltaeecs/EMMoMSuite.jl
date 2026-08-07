using Test
using EMMoMSuite
using LinearAlgebra
using SparseArrays

@testset "Direct Solvers" begin
    # Test dense system
    A = [2.0 1.0; 1.0 3.0]
    b = [3.0, 4.0]
    # x should be [1.0, 1.0]
    
    solver = LUSolver()
    x = solve!(solver, A, b)
    
    @test isapprox(x, [1.0, 1.0])
    
    # Test complex system
    Ac = [1.0+im 0.0; 0.0 1.0-im]
    bc = [1.0+im, 1.0-im]
    
    xc = solve!(solver, Ac, bc)
    @test isapprox(xc, [1.0, 1.0])
    
    # Test sparse system
    As = sparse([1, 2, 3], [1, 2, 3], [1.0, 2.0, 3.0], 3, 3)
    bs = [1.0, 2.0, 3.0]
    
    xs = solve!(solver, As, bs)
    @test isapprox(xs, [1.0, 1.0, 1.0])
end

@testset "Iterative Solvers" begin
    # Test GMRES
    A = [2.0 1.0; 1.0 3.0]
    b = [3.0, 4.0]
    
    solver = GMRESSolver(verbose=false)
    x = solve!(solver, A, b)
    @test isapprox(x, [1.0, 1.0], atol=1e-5)
    
    # Test BiCGSTAB
    solver_bicg = BiCGSTABSolver(verbose=false)
    x_bicg = solve!(solver_bicg, A, b)
    @test isapprox(x_bicg, [1.0, 1.0], atol=1e-5)
    
    # Test with initial guess
    x0 = [0.0, 0.0]
    x_gmres_x0 = solve!(solver, A, b, x0)
    @test isapprox(x_gmres_x0, [1.0, 1.0], atol=1e-5)

    # Test with Preconditioner
    P = Diagonal([1.0/2.0, 1.0/3.0])
    x_precond = solve!(solver, A, b; Pl=P)
    @test isapprox(x_precond, [1.0, 1.0], atol=1e-5)
end

