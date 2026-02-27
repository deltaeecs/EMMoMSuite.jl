# ==================== Mesh Interface ====================
"""
    AbstractMesh

Abstract type for all mesh types.
Must implement:
- `vertices(mesh)`: Return vertex coordinates matrix
- `elements(mesh)`: Return element connectivity matrix
- `dimension(mesh)`: Return mesh dimension
- `num_vertices(mesh)`: Return number of vertices
- `num_elements(mesh)`: Return number of elements
"""
abstract type AbstractMesh end

vertices(mesh::AbstractMesh) = error("Not implemented")
elements(mesh::AbstractMesh) = error("Not implemented")
dimension(mesh::AbstractMesh) = error("Not implemented")
num_vertices(mesh::AbstractMesh) = size(vertices(mesh), 2)
num_elements(mesh::AbstractMesh) = size(elements(mesh), 2)

# ==================== Basis Function Interface ====================
"""
    AbstractBasisFunction

Abstract type for all basis functions.
Must implement:
- `num_basis(bf)`: Return number of basis functions
- `support(bf, i)`: Return support elements for i-th basis function
- `evaluate(bf, i, r)`: Evaluate i-th basis function at point r
"""
abstract type AbstractBasisFunction end

num_basis(bf::AbstractBasisFunction) = error("Not implemented")
support(bf::AbstractBasisFunction, i::Int) = error("Not implemented")
evaluate(bf::AbstractBasisFunction, i::Int, r::AbstractVector) = error("Not implemented")

# ==================== Integral Operator Interface ====================
"""
    AbstractIntegralOperator

Abstract type for integral operators (EFIE, MFIE, etc.).
Must implement:
- `kernel(op, r, r_prime)`: Calculate Green's function kernel
- `impedance_element(op, bf_test, i, bf_trial, j)`: Calculate impedance matrix element
"""
abstract type AbstractIntegralOperator end

kernel(op::AbstractIntegralOperator, r, r_prime) = error("Not implemented")
impedance_element(op::AbstractIntegralOperator, bf_test, i, bf_trial, j) = error("Not implemented")

"""
    assemble_impedance_matrix(op::AbstractIntegralOperator, basis::AbstractBasisFunction)

Assemble the impedance matrix for the given operator and basis.
"""
function assemble_impedance_matrix(op::AbstractIntegralOperator, basis::AbstractBasisFunction)
    error("Not implemented")
end

# ==================== Solver Interface ====================
"""
    AbstractSolver

Abstract type for solvers.
Must implement:
- `solve!(solver, A, b, x0)`: Solve linear system Ax = b
"""
abstract type AbstractSolver end

solve!(solver::AbstractSolver, A, b, x0=nothing) = error("Not implemented")

# ==================== Source Interface ====================
"""
    AbstractSource

Abstract type for excitation sources.
Must implement:
- `incident_field(source, r)`: Calculate incident field at point r
- `excitation_vector(source, basis)`: Calculate excitation vector
"""
abstract type AbstractSource end

incident_field(source::AbstractSource, r) = error("Not implemented")
excitation_vector(source::AbstractSource, basis) = error("Not implemented")

# ==================== Fast Algorithm Interface ====================
"""
    AbstractFastAlgorithm

Abstract type for fast algorithms (MLFMA, ACA, etc.).
Must implement:
- `setup!(algo, basis, operator)`: Initialize the algorithm
- `multiply(algo, x)`: Compute A*x efficiently
"""
abstract type AbstractFastAlgorithm end

setup!(algo::AbstractFastAlgorithm, basis, operator) = error("Not implemented")
multiply(algo::AbstractFastAlgorithm, x) = error("Not implemented")

