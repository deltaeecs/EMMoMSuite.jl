module CFIEModule

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..EFIEModule
using ..MFIEModule
using StaticArrays
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

Assemble the impedance matrix Z for the CFIE.
"""
function assemble_impedance_matrix(cfie::CFIE{FT, CT}, basis::RWGBasis{IT, FT}) where {IT, FT, CT}
    Z_efie = assemble_impedance_matrix(cfie.efie, basis)
    Z_mfie = assemble_impedance_matrix(cfie.mfie, basis)
    
    # Z = alpha * Z_efie + (1-alpha) * Z_mfie
    # Note: Z_mfie is already scaled by eta in assemble_impedance_matrix(MFIE)
    factor_efie = cfie.alpha
    factor_mfie = (1.0 - cfie.alpha)
    
    return factor_efie * Z_efie + factor_mfie * Z_mfie
end

end
