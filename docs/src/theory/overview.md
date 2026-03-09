# Theory Overview

This section summarizes the modeling and solver stack used in EMSuite.

## 1. Problem Classes

- PEC surface scattering: EFIE, MFIE, CFIE.
- Dielectric volume/surface modeling: VEFIE, SCFIE, PMCHW, N-Muller.
- Ports and excitation-driven radiation/scattering workflows.

## 2. Discretization

- Surface bases: RWG.
- Volume bases: SWG, PWC, PWCHex, RBF.
- Matrix assembly follows MoM Galerkin testing with geometry-aware quadrature.

## 3. Operator Layer

- Dense operators provide reference behavior and baseline truth probes.
- PMCHW is structured as 2N block systems with explicit EJ/EM/HJ/HM semantics.
- Strong form is defined with pairing transforms to align trial/test spaces.

## 4. Fast Algorithms

- MLFMA accelerates far interactions through octree hierarchy.
- Near interactions are assembled sparsely and merged with matrix-free far passes.
- Error budgets control near-range, effective leaf size, and truncation controls.

## 5. Krylov and Preconditioning

- LU, GMRES, BiCGSTAB are available depending on scale.
- Block-Jacobi preconditioners are supported for near-field sparse blocks.
- Strong-form solves are used as primary comparison route for dense vs matrix-free backends.

## 6. Validation Philosophy

- Legacy implementation is the source of truth for parity checks.
- Dense-vs-fast probes and medium regression gates are used to separate solver error from operator fidelity.
- For dielectric transmission, formulation-level comparison (PMCHW vs N-Muller) is tracked explicitly.
