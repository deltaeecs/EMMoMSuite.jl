"""
    NMullerModule — N-Muller 均匀介质体表面积分方程

实现 N-Muller transmission formulation 的 dense 基线装配。
"""
module NMullerModule

using ..CoreModule
using ..BasisFunctions
using ..PMCHWModule: PMCHW, _l_block_operator, _k_real, assemble_K_pmchw_offdiag
using ..PMCHWBlockOperatorsModule: pmchw_surface_gram_matrix
using LinearAlgebra
using ...Utilities.Parameters: set_frequency!

using ..CoreModule: DeltaGapSource

import ..CoreModule: assemble_impedance_matrix
import ..CoreModule: excitation_vector

export NMuller, assemble_impedance_matrix

struct NMuller{FT<:AbstractFloat,CT<:Complex} <: AbstractIntegralOperator
    freq::FT
    k0::FT
    eta0::FT
    k1::CT
    eta1::CT
    eps_r::CT
    mu_r::CT
end

function NMuller(freq::FT, eps_r_in, mu_r_in = 1.0) where {FT<:AbstractFloat}
    CT = Complex{FT}
    k0 = FT(2π * freq / Constants.c0)
    eta0 = FT(Constants.eta0)

    eps_r = CT(eps_r_in)
    mu_r = CT(mu_r_in)

    k1 = CT(k0) * sqrt(eps_r * mu_r)
    eta1 = CT(eta0) * sqrt(mu_r / eps_r)

    set_frequency!(Float64(freq))
    return NMuller{FT,CT}(freq, k0, eta0, k1, eta1, eps_r, mu_r)
end

function _medium_l_blocks(op::NMuller{FT,CT}) where {FT,CT}
    k0_c = CT(op.k0)
    eta0_c = CT(op.eta0)
    k1_c = op.k1
    eta1_c = op.eta1

    ej0 = _l_block_operator(op.k0, op.eta0, k0_c, eta0_c, :EJ)
    ej1 = _l_block_operator(k1_c, eta1_c, k1_c, eta1_c, :EJ)
    hm0 = _l_block_operator(op.k0, op.eta0, k0_c, eta0_c, :HM)
    hm1 = _l_block_operator(k1_c, eta1_c, k1_c, eta1_c, :HM)
    return ej0, ej1, hm0, hm1, k1_c
end

"""
    assemble_impedance_matrix(op::NMuller, basis::RWGBasis) -> Matrix{Complex}

按 scuff-em 文档中的 N-Muller `M^(1) - M^(2)` 结构构造 dense baseline。
当前实现只作为 Phase 15 alternative formulation 的第一版组装基线。
"""
function assemble_impedance_matrix(op::NMuller{FT,CT}, basis::RWGBasis{IT,FT}) where {IT,FT,CT}
    N = num_basis(basis)
    Z = zeros(CT, 2N, 2N)

    ej0_op, ej1_op, hm0_op, hm1_op, k1_c = _medium_l_blocks(op)
    ej0 = assemble_impedance_matrix(ej0_op, basis)
    ej1 = assemble_impedance_matrix(ej1_op, basis)
    hm0 = assemble_impedance_matrix(hm0_op, basis)
    hm1 = assemble_impedance_matrix(hm1_op, basis)
    k0 = assemble_K_pmchw_offdiag(basis, op.k0)
    k1 = assemble_K_pmchw_offdiag(basis, k1_c)

    mass = Matrix{CT}(pmchw_surface_gram_matrix(basis))

    mu0 = CT(Constants.mu0)
    eps0 = CT(Constants.eps0)
    mu1 = mu0 * op.mu_r
    eps1 = eps0 * op.eps_r

    Z11 = CT(0.5) * (mu0 + mu1) .* mass .+ mu0 .* k0 .- mu1 .* k1
    Z12 = .-mu0 .* hm0 .+ mu1 .* hm1
    Z21 = .-eps0 .* ej0 .+ eps1 .* ej1
    Z22 = CT(0.5) * (eps0 + eps1) .* mass .- mu0 .* k0 .+ mu1 .* k1

    Z[1:N, 1:N] .= Z11
    Z[1:N, N+1:2N] .= Z12
    Z[N+1:2N, 1:N] .= Z21
    Z[N+1:2N, N+1:2N] .= Z22
    return Z
end

function excitation_vector(op::NMuller{FT,CT}, source::PlaneWave, basis::RWGBasis{IT,FT}) where {IT,FT,CT}
    pmchw = PMCHW(op.freq, op.eps_r, op.mu_r)
    V_pmchw = excitation_vector(pmchw, source, basis)
    N = num_basis(basis)
    V = zeros(CT, 2N)
    V[1:N] .= CT(Constants.mu0) .* V_pmchw[N+1:2N]
    V[N+1:2N] .= CT(Constants.eps0) .* V_pmchw[1:N]
    return V
end

function excitation_vector(op::NMuller{FT,CT}, source::DeltaGapSource, basis::RWGBasis{IT,FT}) where {IT,FT,CT}
    pmchw = PMCHW(op.freq, op.eps_r, op.mu_r)
    V_pmchw = excitation_vector(pmchw, source, basis)
    N = num_basis(basis)
    V = zeros(CT, 2N)
    # This is only a minimal heuristic bridge for Phase 15 diagnostics.
    # Delta-gap semantics for N-Muller have not yet been validated against
    # a formulation-specific port model or an external reference implementation.
    V[1:N] .= CT(Constants.mu0) .* V_pmchw[N+1:2N]
    V[N+1:2N] .= CT(Constants.eps0) .* V_pmchw[1:N]
    return V
end

end # module NMullerModule