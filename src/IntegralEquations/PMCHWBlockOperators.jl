module PMCHWBlockOperatorsModule

using LinearAlgebra
using StaticArrays
using SparseArrays

using ..CoreModule
using ..BasisFunctions
using ..PMCHWModule: PMCHW, assemble_impedance_matrix
using ..Geometry: GaussQuadratureInfo, get_global_quad_points, TriangleInfo
using ..Impedance: assemble_generic

export AbstractPMCHWBackend,
    DensePMCHWBackend,
    MatrixFreePMCHWBackend,
    PMCHWBlockOperator,
    pmchw_blocks,
    pmchw_surface_gram_matrix,
    pmchw_block_pairing_matrix,
    strong_form_rhs,
    recover_trial_coefficients,
    weak_form,
    strong_form

abstract type AbstractPMCHWBackend end

struct DensePMCHWBackend{CT<:Complex,MT<:AbstractMatrix{CT}} <: AbstractPMCHWBackend
    matrix::MT
end

struct MatrixFreePMCHWBackend{A} <: AbstractPMCHWBackend
    operator::A
end

struct PairingTransformedOperator{A,TP,RP} <: AbstractIntegralOperator
    base::A
    test_pairing::TP
    trial_pairing::RP
end

struct RWGSurfaceGramOperator{FT<:AbstractFloat}
    gq_info
end

function _rwg_surface_gram_interaction!(
    Z_local::AbstractMatrix{CT},
    gram::RWGSurfaceGramOperator{FT},
    tri_test::TriangleInfo{IT,FT},
    tri_source::TriangleInfo{IT,FT},
) where {IT,FT,CT<:Complex}
    tri_test.triID == tri_source.triID || return nothing

    r_quad = get_global_quad_points(tri_test, gram.gq_info)
    w_quad = gram.gq_info.weight
    n_pts = length(w_quad)

    v1 = SVector{3,FT}(tri_test.vertices[:, 1])
    v2 = SVector{3,FT}(tri_test.vertices[:, 2])
    v3 = SVector{3,FT}(tri_test.vertices[:, 3])

    inv_4A2 = inv(4 * tri_test.area^2)

    @inbounds for m = 1:3
        lm = tri_test.edgel[m]
        tri_test.inBfsID[m] == 0 && continue

        for n = 1:3
            ln = tri_test.edgel[n]
            tri_test.inBfsID[n] == 0 && continue

            sum_val = zero(CT)
            for k = 1:n_pts
                rk = r_quad[k]
                rho_m = rk - (m == 1 ? v1 : m == 2 ? v2 : v3)
                rho_n = rk - (n == 1 ? v1 : n == 2 ? v2 : v3)
                sum_val += dot(rho_m, rho_n) * w_quad[k]
            end

            Z_local[m, n] += (lm * ln * inv_4A2) * sum_val
        end
    end

    return nothing
end

function pmchw_surface_gram_matrix(basis::RWGBasis{IT,FT}) where {IT,FT}
    gram = RWGSurfaceGramOperator{FT}(GaussQuadratureInfo(:Triangle, 4, FT))
    wrapper = (Z, op, t1, t2) -> _rwg_surface_gram_interaction!(Z, op, t1, t2)
    return assemble_generic(gram, basis, wrapper, symmetric = false)
end

function pmchw_block_pairing_matrix(basis::RWGBasis{IT,FT}) where {IT,FT}
    surface_mass = pmchw_surface_gram_matrix(basis)
    N = num_basis(basis)
    CT = Complex{FT}
    block_mass = zeros(CT, 2N, 2N)
    block_mass[1:N, 1:N] .= surface_mass
    block_mass[N+1:2N, N+1:2N] .= surface_mass
    return block_mass
end

struct PMCHWBlockOperator{FT<:AbstractFloat,CT<:Complex,B<:AbstractPMCHWBackend,MB<:AbstractMatrix{CT},TP,RP} <:
       AbstractIntegralOperator
    pmchw::PMCHW{FT,CT}
    basis::RWGBasis
    backend::B
    ej::MB
    em::MB
    hj::MB
    hm::MB
    test_pairing::TP
    trial_pairing::RP
end

function _split_pmchw_blocks(Z::AbstractMatrix{CT}, basis::RWGBasis) where {CT<:Complex}
    N = num_basis(basis)
    return (
        EJ = Matrix{CT}(Z[1:N, 1:N]),
        EM = Matrix{CT}(Z[1:N, N+1:2N]),
        HJ = Matrix{CT}(Z[N+1:2N, 1:N]),
        HM = Matrix{CT}(Z[N+1:2N, N+1:2N]),
    )
end

function PMCHWBlockOperator(
    pmchw::PMCHW{FT,CT},
    basis::RWGBasis;
    test_pairing = nothing,
    trial_pairing = nothing,
) where {FT<:AbstractFloat,CT<:Complex}
    Z = assemble_impedance_matrix(pmchw, basis)
    backend = DensePMCHWBackend(Z)
    default_pairing = pmchw_block_pairing_matrix(basis)
    test_pairing === nothing && (test_pairing = default_pairing)
    trial_pairing === nothing && (trial_pairing = default_pairing)
    return PMCHWBlockOperator(pmchw, basis, backend; test_pairing = test_pairing, trial_pairing = trial_pairing)
end

function PMCHWBlockOperator(
    pmchw::PMCHW{FT,CT},
    basis::RWGBasis,
    backend::DensePMCHWBackend{CT};
    test_pairing = nothing,
    trial_pairing = nothing,
) where {FT<:AbstractFloat,CT<:Complex}
    block_data = _split_pmchw_blocks(backend.matrix, basis)
    return PMCHWBlockOperator{FT,CT,typeof(backend),Matrix{CT},typeof(test_pairing),typeof(trial_pairing)}(
        pmchw,
        basis,
        backend,
        block_data.EJ,
        block_data.EM,
        block_data.HJ,
        block_data.HM,
        test_pairing,
        trial_pairing,
    )
end

function PMCHWBlockOperator(
    pmchw::PMCHW{FT,CT},
    basis::RWGBasis,
    backend::MatrixFreePMCHWBackend;
    test_pairing = nothing,
    trial_pairing = nothing,
    block_source = nothing,
) where {FT<:AbstractFloat,CT<:Complex}
    test_pairing === nothing && (test_pairing = pmchw_block_pairing_matrix(basis))
    trial_pairing === nothing && (trial_pairing = test_pairing)
    block_source === nothing && (block_source = assemble_impedance_matrix(pmchw, basis))
    block_data = _split_pmchw_blocks(block_source, basis)
    return PMCHWBlockOperator{FT,CT,typeof(backend),Matrix{CT},typeof(test_pairing),typeof(trial_pairing)}(
        pmchw,
        basis,
        backend,
        block_data.EJ,
        block_data.EM,
        block_data.HJ,
        block_data.HM,
        test_pairing,
        trial_pairing,
    )
end

pmchw_blocks(op::PMCHWBlockOperator) = (EJ = op.ej, EM = op.em, HJ = op.hj, HM = op.hm)

weak_form(backend::DensePMCHWBackend) = backend.matrix
weak_form(backend::MatrixFreePMCHWBackend) = backend.operator
weak_form(op::PMCHWBlockOperator) = weak_form(op.backend)

function strong_form(op_backend::PairingTransformedOperator)
    return op_backend
end

function _strong_form_matrix(base::AbstractMatrix, test_pairing, trial_pairing)
    strong = base
    trial_pairing !== nothing && (strong = strong * trial_pairing)
    test_pairing !== nothing && (strong = test_pairing \ strong)
    return strong
end

function strong_form(op::PMCHWBlockOperator)
    base = weak_form(op)
    if base isa AbstractMatrix
        return _strong_form_matrix(base, op.test_pairing, op.trial_pairing)
    end
    return PairingTransformedOperator(base, op.test_pairing, op.trial_pairing)
end

function strong_form_rhs(op::PMCHWBlockOperator, rhs::AbstractVector)
    op.test_pairing === nothing && return rhs
    return op.test_pairing \ rhs
end

function recover_trial_coefficients(op::PMCHWBlockOperator, coeffs::AbstractVector)
    op.trial_pairing === nothing && return coeffs
    return op.trial_pairing * coeffs
end

Base.eltype(::PMCHWBlockOperator{FT,CT}) where {FT,CT} = CT
Base.size(op::PMCHWBlockOperator) = size(weak_form(op))
Base.size(op::PMCHWBlockOperator, dim::Int) = size(weak_form(op), dim)
Base.Matrix(op::PMCHWBlockOperator) = Matrix(weak_form(op))

Base.eltype(op::PairingTransformedOperator) = eltype(op.base)
Base.size(op::PairingTransformedOperator) = size(op.base)
Base.size(op::PairingTransformedOperator, dim::Int) = size(op.base, dim)

function LinearAlgebra.mul!(y::AbstractVector, op::PairingTransformedOperator, x::AbstractVector)
    x_trial = op.trial_pairing === nothing ? x : op.trial_pairing * x
    y_base = similar(y, eltype(op.base), length(y))
    mul!(y_base, op.base, x_trial)
    if op.test_pairing === nothing
        copyto!(y, y_base)
    else
        y .= op.test_pairing \ y_base
    end
    return y
end

Base.:*(op::PairingTransformedOperator, x::AbstractVector) = (y = similar(x, eltype(op), length(x)); mul!(y, op, x); y)

function LinearAlgebra.mul!(y::AbstractVector, op::PMCHWBlockOperator, x::AbstractVector)
    mul!(y, weak_form(op), x)
    return y
end

Base.:*(op::PMCHWBlockOperator, x::AbstractVector) = weak_form(op) * x

end # module PMCHWBlockOperatorsModule