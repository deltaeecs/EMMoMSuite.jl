# Advanced Guide

This guide focuses on medium/large-scale workflows and fidelity controls used in EMSuite Phase 15.

## 1. PMCHW Block/Operator Strong-Form Solve

Use PMCHW shell with strong form when comparing dense and matrix-free backends under the same physical probe.

```julia
using EMSuite
using IterativeSolvers

freq = 300e6
mesh = generate_sphere_mesh(0.5, 10, 20)
basis = RWGBasis(mesh)
pmchw = PMCHW(freq, 4.0)

# Dense reference path
Z_dense = assemble_impedance_matrix(pmchw, basis)
dense_shell = PMCHWBlockOperator(pmchw, basis, DensePMCHWBackend(Z_dense))

feed = DeltaGapSource(freq, [1], 1.0 + 0im)
rhs = excitation_vector(pmchw, feed, basis)

coeffs_sf, hist = gmres(
    strong_form(dense_shell),
    strong_form_rhs(dense_shell, rhs);
    reltol = 1e-4,
    maxiter = 200,
    log = true,
)
I = recover_trial_coefficients(dense_shell, coeffs_sf)
Zin = input_impedance(pmchw, feed, I, basis)
```

## 2. PMCHW MLFMA Budget Controls

`PMCHWMLFMAErrorBudget` exposes near/far and truncation controls explicitly.

```julia
using EMSuite
using EMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator

budget = PMCHWMLFMAErrorBudget(Float64;
    near_range_scale = 4.0,
    min_near_range = 4,
    max_near_range = 16,
    L_min = 0,
)

op = PMCHWMLFMAOperator(pmchw, basis, 0.10; budget = budget)
mlfma_shell = PMCHWBlockOperator(
    pmchw,
    basis,
    MatrixFreePMCHWBackend(op);
    block_source = Z_dense,
)
```

Recommended workflow:
1. Keep the same `rhs` and the same strong-form GMRES settings.
2. Sweep budget settings only.
3. Track `near_density`, residual, and impedance gap versus dense shell.

## 3. Medium Regression Gates

High-value, medium-scale specialized tests:

- `test/test_pmchw_gate_s_mlfma_medium.jl`
- `test/test_pmchw_mlfma_budget_medium.jl`
- `test/test_pmchw_mlfma_budget_krylov_medium.jl`
- `test/test_pmchw_block_fidelity_medium.jl`
- `test/test_nmuller_comparison_medium.jl`

Run examples:

```bash
julia --project=. test/test_pmchw_gate_s_mlfma_medium.jl
julia --project=. test/test_pmchw_mlfma_budget_krylov_medium.jl
```

## 4. Benchmark Entrypoints

- `benchmark/compare_pmchw_mlfma_budget.jl`
- `benchmark/compare_pmchw_mlfma_budget_krylov.jl`
- `benchmark/compare_pmchw_block_fidelity_medium.jl`
- `benchmark/compare_pmchw_gmres_trajectory_medium.jl`
- `benchmark/compare_pmchw_nmuller_sphere.jl`

Typical run:

```bash
julia --project=. benchmark/compare_pmchw_mlfma_budget.jl medium
julia --project=. benchmark/compare_pmchw_mlfma_budget_krylov.jl medium
```
