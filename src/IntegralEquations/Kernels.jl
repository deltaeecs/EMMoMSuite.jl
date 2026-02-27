module Kernels

using StaticArrays
using LinearAlgebra

export green_function_free_space, grad_green_function_free_space

"""
    green_function_free_space(r, r_prime, k)

Calculate the free space Green's function G(r, r') = exp(-j*k*R) / R.
"""
function green_function_free_space(r::AbstractVector{FT}, r_prime::AbstractVector{FT}, k::FT) where {FT<:Real}
    R_vec = r - r_prime
    R = norm(R_vec)
    if R < 1e-12
        return Complex{FT}(0.0) 
    end
    return exp(-im * k * R) / R
end

"""
    grad_green_function_free_space(r, r_prime, k)

Calculate the gradient of the free space Green's function with respect to r.
Returns a vector.
"""
function grad_green_function_free_space(r::AbstractVector{FT}, r_prime::AbstractVector{FT}, k::FT) where {FT<:Real}
    R_vec = r - r_prime
    R = norm(R_vec)
    if R < 1e-12
        return SVector{3, Complex{FT}}(0, 0, 0)
    end
    G = exp(-im * k * R) / R
    # grad G = -(jk + 1/R) * G * R_hat
    #        = -(jk + 1/R) * G * R_vec / R
    factor = -(im * k + 1/R) / R
    return factor * G * R_vec
end

end
