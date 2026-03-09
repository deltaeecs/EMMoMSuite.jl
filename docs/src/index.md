# EMSuite.jl

EMSuite is a Julia-based computational electromagnetics toolkit centered on Method of Moments workflows for scattering, radiation, and transmission problems.

## What You Get

- Surface and volume formulations: EFIE, MFIE, CFIE, VEFIE, SCFIE, PMCHW, and N-Muller.
- Basis support: RWG, SWG, PWC, PWCHex, RBF.
- Fast backend: MLFMA with medium/large-scale regression gates.
- Solver stack: LU, GMRES, BiCGSTAB, and preconditioners.
- Post-processing: RCS, near/far fields, antenna metrics, and export helpers.

## Documentation Map

- [Installation](guide/installation.md)
- [Quick Start](guide/quick_start.md)
- [Advanced Guide](guide/advanced.md)
- [Examples](guide/examples.md)
- [Theory Overview](theory/overview.md)
- [API Reference](api/public_api.md)

## Local Development Install

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

For package usage details, continue with the [User Guide](guide/quick_start.md).
