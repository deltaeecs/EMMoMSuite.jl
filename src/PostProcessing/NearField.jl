module NearField

using StaticArrays
using LinearAlgebra
using ...CoreModule
using ...CoreModule: Constants
using ...Geometry
using ...BasisFunctions
using ...IntegralEquations.Kernels
using ...Utilities.Parameters

export calculate_near_field

"""
    calculate_near_field(points, basis, I_coeffs) -> Vector{SVector{3, Complex}}

Calculate electric field at observation points from surface currents.

# Mathematical Background

The electric field is computed via vector and scalar potentials:

```math
\\mathbf{E}(\\mathbf{r}) = -j\\omega\\mu_0 \\mathbf{A}(\\mathbf{r}) - \\nabla\\Phi(\\mathbf{r})
```

where:
```math
\\mathbf{A}(\\mathbf{r}) = \\sum_n I_n \\int_{S_n} \\mathbf{f}_n(\\mathbf{r}') G(\\mathbf{r}, \\mathbf{r}') dS'
```
```math
\\Phi(\\mathbf{r}) = \\frac{1}{j\\omega\\varepsilon_0} \\sum_n I_n \\int_{S_n} \\nabla' \\cdot \\mathbf{f}_n G(\\mathbf{r}, \\mathbf{r}') dS'
```

# Arguments

- `points::Vector{SVector{3, FT}}`: Observation point coordinates [x, y, z] in meters
- `basis::RWGBasis`: RWG basis functions on surface mesh
- `I_coeffs::Vector{Complex{FT}}`: Current coefficients (from MoM solve)

# Returns

- `Vector{SVector{3, Complex{FT}}}`: Electric field [Ex, Ey, Ez] at each point [V/m]

# Examples

```julia
# After solving EFIE
Z = assemble_impedance_matrix(efie)
V = compute_excitation_vector(efie, plane_wave)
I = Z \\ V

# Define observation grid (10×10 grid at z=1m above antenna)
x_obs = range(-0.5, 0.5, length=10)
y_obs = range(-0.5, 0.5, length=10)
z_obs = 1.0

points = [SVector(x, y, z_obs) for x in x_obs, y in y_obs] |> vec

# Compute near field
E_near = calculate_near_field(points, basis, I)

# Extract magnitude
E_mag = [norm(E) for E in E_near]
E_dB = 20 .* log10.(E_mag ./ maximum(E_mag))

# Visualize (requires plotting package)
using Plots
heatmap(x_obs, y_obs, reshape(E_dB, 10, 10), 
        xlabel="x [m]", ylabel="y [m]", 
        title="Near Field (dB)")
```

# Performance Notes

- Uses 3-point Gauss quadrature per triangle
- Complexity: O(N_points × N_basis × 2 triangles/basis)
- No acceleration (direct integration)
- For large meshes (>1000 unknowns) or many points (>10000), 
  consider MLFMA-accelerated near-field

# Validity Conditions

- **Near-field**: Valid at any distance (r > 0)
- **Observation points inside mesh**: Gives interior field (if penetrable)
- **Observation points on surface**: Singular behavior (use field averaging)

# See Also

- [`farField`](@ref): Far-field radiation pattern
- [`RWGBasis`](@ref): Surface basis functions
- [`assemble_impedance_matrix`](@ref): MoM solver

# Notes

- Returned field includes both incident and scattered components if I_coeffs from total-field solve
- For scattered-field only, use I_coeffs from scattered-field formulation
- Quadrature accuracy: 3-point rule sufficient for smooth integrands
"""
function calculate_near_field(
    points::Vector{SVector{3,FT}},
    basis::RWGBasis{IT,FT},
    I_coeffs::Vector{Complex{FT}},
) where {IT,FT}
    num_points = length(points)
    E_field = zeros(SVector{3,Complex{FT}}, num_points)

    omega = get_omega()
    k = omega / Constants.c0

    # Precompute constants
    const_A = -im * omega * Constants.mu0
    const_Phi = 1.0 / (im * omega * Constants.eps0)

    # Quadrature for triangles
    gq = GaussQuadratureInfo(:Triangle, 3, FT)

    mesh = basis.mesh

    # Loop over basis functions
    for n = 1:num_basis(basis)
        In = I_coeffs[n]
        if abs(In) < 1e-12
            continue
        end

        bf = basis.functions[n]

        # Loop over support triangles
        for k_supp = 1:2
            tri_idx = bf.support[k_supp]
            if tri_idx == 0
                continue
            end

            sign = bf.signs[k_supp]

            # Get triangle info
            tri = get_triangle_info(mesh, basis, tri_idx)

            # Divergence of f_n is constant on triangle
            # div(f_n) = +/- l_n / A_n
            div_f = sign * bf.edge_length / tri.area

            # Quadrature points on triangle
            q_points = get_global_quad_points(tri, gq)
            q_weights = gq.weight .* tri.area

            # Contribution to A and Phi
            for (qp, w) in zip(q_points, q_weights)
                # Value of f_n at qp
                # f_n = sign * l_n / (2*A) * rho
                # rho = r - v_op
                # We need v_op (opposite vertex)
                # In get_triangle_info, we don't explicitly identify v_op for this basis
                # But we can deduce it.
                # bf.local_edge_idx[k] is the local index of the edge.
                # The opposite vertex is the one with that local index.
                local_op_idx = bf.local_edge_idx[k_supp]
                v_op = tri.vertices[:, local_op_idx]

                rho = qp - v_op
                f_val = (sign * bf.edge_length / (2 * tri.area)) * rho

                for i = 1:num_points
                    obs = points[i]
                    G = green_function_free_space(obs, qp, k)
                    grad_G = grad_green_function_free_space(obs, qp, k)

                    # E = -j*w*mu*A - grad(Phi)
                    # Phi = (1/eps) * int(rho * G)
                    # rho = (-1/j*w) * div(J)
                    # Phi = (-1/j*w*eps) * int(div(J) * G)
                    # grad(Phi) = (-1/j*w*eps) * int(div(J) * grad(G))
                    # -grad(Phi) = (1/j*w*eps) * int(div(J) * grad(G))
                    # const_Phi = 1/(j*w*eps)
                    # So -grad(Phi) = const_Phi * div_f * grad_G

                    term_A = const_A * f_val * G
                    term_Phi = const_Phi * div_f * grad_G

                    E_field[i] += In * (term_A + term_Phi) * w
                end
            end
        end
    end

    return E_field
end

# ─── Phase 17.1: Volume basis functions ────────────────────────────────────────

"""
    calculate_near_field(points, basis::SWGBasis, I_coeffs, permittivities)

Calculate the scattered electric field at observation points for a VEFIE SWG
solution using the free-space dyadic Green's function:

    E(r) = -jωμ₀ Σₙ Iₙ κ Σ_{T∈support(n)} ∫_T f_n(r') G(r,r') V dV'
           + (1/jωε₀) Σₙ Iₙ κ Σ_{T∈support(n)} ∇·f_n ∫_T ∇G(r,r') V dV'

`permittivities[t]` is the complex relative permittivity of tetrahedron `t`.
"""
function calculate_near_field(
    points::Vector{SVector{3,FT}},
    basis::SWGBasis{IT,FT},
    I_coeffs::Vector{Complex{FT}},
    permittivities::Vector{Complex{FT}},
) where {IT,FT}
    num_points = length(points)
    E_field = zeros(SVector{3,Complex{FT}}, num_points)

    omega = get_omega()
    k = omega / Constants.c0

    const_A = -im * omega * Constants.mu0         # coefficient of vector potential term
    const_Phi = 1.0 / (im * omega * Constants.eps0) # coefficient of scalar potential term

    # Tetrahedral Gauss quadrature (4 points)
    gq_pts, gq_wts = Geometry.gaussQuadratureTet(4, FT)
    nqp = length(gq_wts)

    mesh = basis.mesh
    nodes = mesh.node
    tetras = mesh.tetras

    for n = 1:length(basis.functions)
        In = I_coeffs[n]
        abs(In) < 1e-12 && continue

        bf = basis.functions[n]

        for k_supp = 1:2
            t_idx = bf.support[k_supp]
            t_idx == 0 && continue

            eps_r = permittivities[t_idx]
            kappa = (eps_r - 1.0) / eps_r
            factor = kappa * In

            # Tet vertices
            vi = tetras[:, t_idx]
            r1 = nodes[:, vi[1]]
            r2 = nodes[:, vi[2]]
            r3 = nodes[:, vi[3]]
            r4 = nodes[:, vi[4]]

            # Tet volume
            vol = abs(det(hcat(r2 - r1, r3 - r1, r4 - r1))) / 6.0

            # Free vertex (vertex opposite the shared face)
            lf = bf.local_face_idx[k_supp]
            v_free = nodes[:, vi[lf]]

            # Constant amplitude factor for f_n = ±(A/(3V)) * (r - v_free)
            sign_k = k_supp == 1 ? 1.0 : -1.0   # + or − tet
            const_bf = bf.area / (3.0 * vol)       # A/(3V)
            # Plus tet: f_n = const_bf*(r-v_free), Minus tet: f_n = -const_bf*(r-v_free)
            # → sign_k already handled by basis definition,
            #   for the minus tet we use v_free - r, so sign = -1

            # div(f_n): +A/V for plus tet, -A/V for minus tet
            div_fn = sign_k * bf.area / vol   # A/V with sign

            # Quadrature
            for gi = 1:nqp
                u, v, w, x = gq_pts[1, gi], gq_pts[2, gi], gq_pts[3, gi], gq_pts[4, gi]
                rgi = u * r1 + v * r2 + w * r3 + x * r4

                # f_n value at rgi (linear vector field)
                ρ = rgi .- v_free
                f_val = sign_k * const_bf .* ρ

                wvol = gq_wts[gi] * vol

                for i = 1:num_points
                    obs = points[i]
                    G = green_function_free_space(obs, SVector{3,FT}(rgi), k)
                    grad_G = grad_green_function_free_space(obs, SVector{3,FT}(rgi), k)

                    term_A = const_A * SVector{3,Complex{FT}}(f_val) * G
                    term_Phi = const_Phi * div_fn * grad_G

                    E_field[i] += factor * wvol * (term_A + term_Phi)
                end
            end
        end
    end

    return E_field
end

"""
    calculate_near_field(points, basis::PWCBasis, I_coeffs, permittivities)

Near-field E calculation for PWC body basis functions.

For PWC, J_eq is piecewise constant (∇·J_eq = 0 within each tet), so only the
vector potential term contributes:

    E(r) ≈ -jωμ₀ Σ_t κ_t J_t ∫_T G(r,r') dV'
"""
function calculate_near_field(
    points::Vector{SVector{3,FT}},
    basis::PWCBasis{IT,FT},
    I_coeffs::Vector{Complex{FT}},
    permittivities::Vector{Complex{FT}},
) where {IT,FT}
    num_points = length(points)
    E_field = zeros(SVector{3,Complex{FT}}, num_points)

    omega = get_omega()
    k = omega / Constants.c0
    const_A = -im * omega * Constants.mu0

    gq_pts, gq_wts = Geometry.gaussQuadratureTet(4, FT)
    nqp = length(gq_wts)

    mesh = basis.mesh
    nodes = mesh.node
    tetras = mesh.tetras

    for t = 1:length(basis.functions)
        pwc = basis.functions[t]
        eps_r = permittivities[t]
        kappa = (eps_r - 1.0) / eps_r
        vol = pwc.volume

        # Current vector for this tet
        Jt = SVector{3,Complex{FT}}(
            I_coeffs[pwc.inBfsID[1]],
            I_coeffs[pwc.inBfsID[2]],
            I_coeffs[pwc.inBfsID[3]],
        )

        vi = tetras[:, t]
        r1 = nodes[:, vi[1]]
        r2 = nodes[:, vi[2]]
        r3 = nodes[:, vi[3]]
        r4 = nodes[:, vi[4]]

        for gi = 1:nqp
            u, v, w, x = gq_pts[1, gi], gq_pts[2, gi], gq_pts[3, gi], gq_pts[4, gi]
            rgi = u * r1 + v * r2 + w * r3 + x * r4
            wvol = gq_wts[gi] * vol

            for i = 1:num_points
                obs = points[i]
                G = green_function_free_space(obs, SVector{3,FT}(rgi), k)

                E_field[i] += kappa * wvol * const_A * Jt * G
            end
        end
    end

    return E_field
end

end
