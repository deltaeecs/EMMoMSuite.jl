#!/usr/bin/env julia
# benchmark/mlfma_nearrange_tradeoff.jl
#
# issue #22 (problem 3) — (leaf size, adaptive near_range) efficiency trade-off map.
#
# MLFMAOperator derives the tree-uniform near_range adaptively from the leaf
# level (adaptive_near_range: kR_min >= 0.55*L_leaf guarantees real-spectrum
# M2L convergence). Electrically SMALL leaf cubes therefore force a large
# near_range, and the direct near-field matrix grows as (2nr+1)^3 cubes per
# cube. This script sweeps leaf sizes in wavelengths on ONE fixed mesh and
# reports, per configuration:
#   - adaptive near_range and Z_near nnz,
#   - construction (setup) time,
#   - matvec time and heap allocations,
#   - MLFMA-vs-dense matvec relative error.
#
# Usage: julia --project=. benchmark/mlfma_nearrange_tradeoff.jl [n_theta n_phi]

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.MLFMA
using EMMoMSuite.FastAlgorithms.MLFMA.Level: adaptive_near_range
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncationLCal
using LinearAlgebra
using Printf
using SparseArrays

function main(n_theta::Int = 20, n_phi::Int = 40)
    freq = 300e6
    λ = 299792458.0 / freq
    radius = 0.5

    println("="^74)
    println(" MLFMA (leaf size, adaptive near_range) trade-off — issue #22 problem 3")
    println("="^74)
    @printf("freq = %.0f MHz, λ = %.3f m, sphere r = %.2f m\n", freq / 1e6, λ, radius)

    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    efie = EFIE(freq)
    N = num_basis(basis)
    @printf("mesh %d×%d, unknowns N = %d\n\n", n_theta, n_phi, N)

    # Dense reference matvec (small N keeps this cheap)
    x = ones(ComplexF64, N)
    Z = assemble_impedance_matrix(efie, basis)
    y_dense = Z * x

    @printf("%-10s %-5s %-12s %-10s %-10s %-12s %-10s %-12s\n",
            "leaf/λ", "nr", "setup(s)", "Z_near nnz", "matvec(s)", "matvec MiB", "rel err", "L_leaf")
    println("-"^88)

    for leaf_rel in (1.0, 0.5, 0.25, 0.1)
        leaf = leaf_rel * λ
        L_leaf = max(truncationLCal(leaf; λ = λ), 0)
        nr = adaptive_near_range(λ, leaf, L_leaf)

        t_setup = @elapsed op = MLFMAOperator(efie, basis, leaf)
        nnz = nnz(op.Z_near)

        y = op * x                     # warm-up (compile)
        t_mv = @elapsed op * x
        alloc = @allocated op * x
        rel = norm(y - y_dense) / norm(y_dense)

        @printf("%-10.2f %-5d %-12.2f %-10d %-10.4f %-12.1f %-10.3e %-12d\n",
                leaf_rel, nr, t_setup, nnz, t_mv, alloc / 1024 / 1024, rel, L_leaf)
    end

    println("-"^88)
    println("Guidance: electrically larger leaf cubes keep the adaptive near_range")
    println("small (the direct near-field matrix grows as (2nr+1)^3 cubes per cube);")
    println("finer leaves raise near_range for M2L convergence (kR_min >= 0.55*L) and")
    println("should be paired with BlockJacobiPreconditioner. See issue #22.")
end

let args = length(ARGS) >= 2 ? (parse(Int, ARGS[1]), parse(Int, ARGS[2])) : ()
    main(args...)
end
