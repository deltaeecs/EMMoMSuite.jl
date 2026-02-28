module CFIEModule

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..EFIEModule
using ..MFIEModule
using ..Impedance
using LinearAlgebra

import ..CoreModule: assemble_impedance_matrix

export CFIE, assemble_impedance_matrix

"""
    CFIE{FT, CT} <: AbstractIntegralOperator

Combined Field Integral Equation (CFIE) operator.

A linear combination of EFIE and MFIE used to eliminate internal resonance problems that plague EFIE and MFIE for closed structures at specific frequencies.

# Mathematical Formulation

The CFIE is defined as:
```math
\\text{CFIE} = \\alpha \\cdot \\text{EFIE} + (1-\\alpha) \\eta \\cdot \\text{MFIE}
```

The impedance matrix is constructed as:
```math
\\mathbf{Z}_{CFIE} = \\alpha \\mathbf{Z}_{EFIE} + (1-\\alpha) \\eta \\mathbf{Z}_{MFIE}
```

where:
- \$\\alpha\$ is the weighting factor (typically 0.5).
- \$\\eta\$ is the intrinsic impedance of the medium (used to balance units).

# Fields
- `freq`: Operating frequency.
- `alpha`: Weighting factor \$\\alpha\$ (0 to 1).
- `efie`: Underlying `EFIE` operator.
- `mfie`: Underlying `MFIE` operator.
"""
struct CFIE{FT<:AbstractFloat, CT<:Complex} <: AbstractIntegralOperator
    freq::FT
    alpha::FT
    efie::EFIE{FT, CT}
    mfie::MFIE{FT, CT}
end

function CFIE(freq::FT, alpha::FT = 0.5) where {FT}
    efie = EFIE(freq)
    mfie = MFIE(freq)
    return CFIE{FT, Complex{FT}}(freq, alpha, efie, mfie)
end

"""
    assemble_impedance_matrix(cfie::CFIE, basis::RWGBasis)

Assemble the CFIE impedance matrix as a linear combination of EFIE and MFIE matrices:

    Z = 伪 * Z_EFIE + (1-伪) * Z_MFIE

Each sub-assembly uses the independently optimised parallel kernel.
Note: Z_MFIE is already scaled by 畏 in `assemble_impedance_matrix(MFIE)`.
"""
function assemble_impedance_matrix(cfie::CFIE{FT, CT}, basis::RWGBasis{IT, FT}) where {IT, FT, CT}
    Z_efie = assemble_impedance_matrix(cfie.efie, basis)
    Z_mfie = assemble_impedance_matrix(cfie.mfie, basis)
    # Z = alpha * Z_efie + (1-alpha) * Z_mfie  (in-place: avoids third N×N allocation)
    # Note: Z_mfie is already weighted by eta in assemble_impedance_matrix(MFIE)
    α = cfie.alpha
    β = FT(1.0) - α
    @. Z_efie = α * Z_efie + β * Z_mfie
    return Z_efie
end

end

