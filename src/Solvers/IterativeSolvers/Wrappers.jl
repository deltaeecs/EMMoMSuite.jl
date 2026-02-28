module IterativeWrappers

using ..CoreModule
using LinearAlgebra
using IterativeSolvers

import ..CoreModule: solve!

export GMRESSolver, BiCGSTABSolver, solve!

"""
    GMRESSolver <: AbstractSolver

Iterative solver using the Generalized Minimal Residual (GMRES) method.

Minimizes the residual norm \$\\min_x ||\\mathbf{b} - \\mathbf{A}\\mathbf{x}||_2\$ over the Krylov subspace \$\\mathcal{K}_m(\\mathbf{A}, \\mathbf{b})\$.

# Fields
- `restart`: Number of iterations before restarting (controls memory usage).
- `maxiter`: Maximum number of outer iterations.
- `tol`: Relative tolerance for convergence.
- `verbose`: If true, prints convergence history.
"""
struct GMRESSolver <: AbstractSolver
    restart::Int
    maxiter::Int
    tol::Float64
    verbose::Bool
end

function GMRESSolver(; restart = 30, maxiter = 100, tol = 1e-6, verbose = false)
    GMRESSolver(restart, maxiter, tol, verbose)
end

function solve!(solver::GMRESSolver, A, b, x0 = nothing; Pl = Identity(), Pr = Identity())
    kwargs = (
        restart = solver.restart,
        maxiter = solver.maxiter,
        reltol = solver.tol,
        verbose = solver.verbose,
        Pl = Pl,
        Pr = Pr,
    )

    if x0 !== nothing
        x = copy(x0)
        return gmres!(x, A, b; kwargs...)
    else
        return gmres(A, b; kwargs...)
    end
end

"""
    BiCGSTABSolver <: AbstractSolver

Iterative solver using the Biconjugate Gradient Stabilized (BiCGSTAB) method.

A variant of BiCG that provides smoother convergence and does not require the transpose of the matrix. Suitable for non-symmetric linear systems arising in MoM.

# Fields
- `maxiter`: Maximum number of iterations.
- `tol`: Relative tolerance for convergence.
- `verbose`: If true, prints convergence history.
"""
struct BiCGSTABSolver <: AbstractSolver
    maxiter::Int
    tol::Float64
    verbose::Bool
end

function BiCGSTABSolver(; maxiter = 100, tol = 1e-6, verbose = false)
    BiCGSTABSolver(maxiter, tol, verbose)
end

function solve!(solver::BiCGSTABSolver, A, b, x0 = nothing; Pl = Identity(), Pr = Identity())
    kwargs = (
        max_mv_products = solver.maxiter * 4,
        reltol = solver.tol,
        verbose = solver.verbose,
        Pl = Pl,
    )

    if x0 !== nothing
        x = copy(x0)
        return bicgstabl!(x, A, b; kwargs...)
    else
        return bicgstabl(A, b; kwargs...)
    end
end

end
