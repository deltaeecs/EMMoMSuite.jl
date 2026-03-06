# Phase 15 Governance Refresh: Theory -> Implementation -> Test

## Scope

This document defines enforceable rules to keep PMCHW/MLFMA implementation strictly aligned with theory and verified by reproducible tests.

## Applicable Principles

- Principle 1: TDD workflow (RED -> GREEN -> REFACTOR)
- Principle 2: Legacy alignment is mandatory for comparable paths
- Principle 3: Difference diagnosis must follow geometry/constants/integration/assembly order
- Principle 5: Every successful logical change must be committed
- Principle 7: End-of-phase review iteration is mandatory
- Principle 9: Plan and DoD must be explicit and testable

## Why Detailed Docs Still Failed

- Spec text was descriptive, not executable: no machine-checkable invariants were attached.
- Kernel approximation policy was not locked per path (real-k vs complex-k usage drifted by stage).
- Coefficient chain checkpoints existed, but not as mandatory pass-level tests.
- Near/far partition assumptions were changed without a synchronized acceptance matrix.
- Debug scripts and production path occasionally diverged in level flow, causing false conclusions.

## Legacy Alignment Baseline

Use these as single source references for algorithm behavior and coefficient chain:

- `MoM_Kernels/*` (kernel-level constants and singular handling)
- `MoM_Basics/*` (basis, geometry conventions, sign conventions)
- `MoM_AllinOne/*` (end-to-end assembly and solver behavior)

For PMCHW-specific blocks where no legacy PMCHW exists:

- Use EMSuite direct PMCHW assembly as the internal reference baseline.
- Enforce consistency between MLFMA decomposition and direct block decomposition.

## Enforceable Consistency Contract

For each PMCHW block and pass, lock the following:

- Phase convention: source and receive exponent signs must match one declared convention.
- Kernel medium convention: if direct uses lossless-kernel approximation, MLFMA phase kernel must use the same convention.
- Coefficient chain: all external factors are declared once and checked via pass-level ratio tests.
- Near/far partition: one declared neighborhood width; octree neighbors, translation index range, and near extraction must match.
- Dataflow parity: diagnostic pass flow must match production `mul!` flow exactly.

## Test Gates (Mandatory)

### Gate A: Structural Invariants

- `HJ + EM` near-field block infinity norm ratio < 1e-8
- Non-trivial near/far split: `nnz(Z_near) < (2N)^2`
- Two octrees keep identical topology ordering assumptions where required by implementation

### Gate B: Pass-Level Alignment

Use fixed seed and normalized vectors.

- EJ k0 pass vs direct-far EJ-k0: relative error and correlation tracked
- EJ k1 pass vs direct-far EJ-k1: relative error and correlation tracked
- EM/HM k0 and k1 pass checks tracked independently
- Sum of passes must match corresponding direct-far block sum

### Gate C: End-to-End Accuracy

- `test/test_pmchw_mlfma_operator.jl` case `15.11 MLFMA mul! vs Direct`: `rel_err < 0.10`
- PMCHW input impedance gate: `Re` error < 5% and physically valid sign constraints

### Gate D: Regression Safety

Before merge, run:

- PMCHW targeted test file
- EFIE/SCFIE smoke checks that share MLFMA core
- Diagnostic scripts used for decisions must be committed or archived with exact command and seed

## Review Iteration Plan

Perform at least two clean review rounds after implementation appears complete.

Round 1 checklist:

- Architecture: no hidden divergence between debug flow and production flow
- Algorithm: all constants and kernel conventions match declared contract
- Engineering: no temporary toggles, no ambiguous comments, no dead branches

Round 2 checklist:

- Re-run all mandatory gates from clean environment
- Confirm no new issues introduced in shared MLFMA components

Exit rule:

- Two consecutive rounds with no new findings

## Definition of Done (DoD)

A task is done only if all conditions are met:

- RED/GREEN/REFACTOR evidence exists in commit history
- Consistency contract items are explicitly satisfied
- Gate A/B/C/D all pass and results are recorded
- Roadmap and progress files are updated with date and outcomes

## Trace Links

- Roadmap anchor: `REFACTORING_ROADMAP.md` -> Phase 15 -> task `15.G1`
- Progress anchor: `REFACTORING_PROGRESS.md` -> Phase 15 update log -> governance refresh entry
